program post_lammps

use mData_Proc
use mString
implicit none

integer, parameter :: maxnarg = 15
integer            :: narg
character(sl)      :: arg(maxnarg)

integer            :: iarg, n

character(5), parameter :: keyword(1) = ["-skip"]

integer                    :: np, nlines, npoints
character(sl), allocatable :: property(:)
integer,       allocatable :: indx(:)
real(rb),      allocatable :: value(:,:)

integer       :: i, j, first, nbins, window, skip
character(sl) :: infile, action, line
logical       :: exist
character(12) :: C

! Defaults:
skip = 0

! Read options:
iarg = 1
if (iargc() < 1) call Usage_Message
call getarg( iarg, line )
do while (any(keyword == line))
  select case (trim(line))
    case ("-skip")
      iarg = iarg + 1
      call getarg( iarg, line )
      skip = str2int( line )
      if (skip < 0) stop "Unacceptable value of skip parameter"
  end select
  iarg = iarg + 1
  call getarg( iarg, line )
end do

! First argument is the specified action:
call getarg( iarg, action )

! Check specified action:
select case (trim(action))
  case ("acf","eacf","histogram"); n = 1
  case ("block","ineff","extract","list","stats"); n = 0
  case default; call Usage_Message
end select

! Read action arguments:
if (iargc() < iarg + n) call Usage_Message
select case (trim(action))
  case ("acf","eacf")
    iarg = iarg + 1
    call getarg( iarg, line )
    window = str2int( line )
    if (window <= 0) stop "Unacceptable maximum window size"
  case ("histogram")
    iarg = iarg + 1
    call getarg( iarg, line )
    nbins = str2int( line )
    if (nbins <= 0) stop "Unacceptable number of bins"
end select

! Read file name:
iarg = iarg + 1
if (iargc() < iarg) call Usage_Message
call getarg( iarg, infile )
inquire( file = infile, exist = exist )
if (.not.exist) then
  write(6,'("Error: specified file <",A,"> does not exist.")') trim(infile)
  stop
end if

! Perform listing action, if specified:
if (trim(action) == "list") then
  call List_Properties( infile )
  stop
end if

! Read property names:
np = iargc() - iarg
if (np == 0) stop "No properties specified."
allocate( property(np), indx(np) )
do i = 1, np
  call getarg( iarg + i, property(i) )
end do

! Search for the last run containing the specified properties:
call find_last_run( infile, np, property, indx, first, nlines )

! Calculate the actual number of points:
npoints = 0
do j = 1, nlines
  if (mod(j-1,skip+1) == 0) npoints = npoints + 1
end do
allocate( value(np,npoints) )

! Read property values from the log file:
open( unit = 10, file = infile, status = "old" )
write(C,*) first-2
read(10,'('//C//'/)')
npoints = 0
do j = 1, nlines
  read(unit=10,fmt='(A'//csl//')') line
  if (mod(j-1,skip+1) == 0) then
    call split( line, narg, arg )
    npoints = npoints + 1
    do i = 1, np
      read(arg(indx(i)),*) value(i,npoints)
    end do
  end if
end do
close(10)

! Perform specified action:
select case (trim(action))
  case ("acf")
    call Compute_ACF( np, property, npoints, value, window, print = .true., norm = .false. )
  case ("block")
    call Block_Analysis( np, property, npoints, value )
  case ("eacf")
    call Compute_ACF( np, property, npoints, value, window, print = .true., norm = .true. )
  case ("extract")
    call Print_Properties( np, property, npoints, value )
  case ("histogram")
    call Build_Histograms( np, property, npoints, value, nbins )
  case ("stats")
    call Statistics( np, property, npoints, value, print = .true. )
end select

contains

  !=================================================================================================

  subroutine Usage_Message
    write(6,'("Usage: post_lammps [options] action [args] file-name property-1 [property-2] ...")')
    write(6,'("  action = acf or block or extract or histogram or list or stats")')
    write(6,'("    acf args = window")')
    write(6,'("    block args = none")')
    write(6,'("    extract args = none")')
    write(6,'("    histogram args = nbins")')
    write(6,'("    list args = none")')
    write(6,'("    stats args = none")')
    stop
  end subroutine Usage_Message

  !=================================================================================================

  subroutine List_Properties( file )
    character(*), intent(in)  :: file
    integer :: iline, ierr, last
    character(sl) :: line, titles
    logical :: titles_found
    open( unit = 10, file = file, status = "old", iostat = ierr )
    if (ierr /= 0) call error( "could not open file", file )
    iline = 0
    first = 0
    last = 0
    titles_found = .false.
    do while (ierr == 0)
      iline = iline + 1
      read(unit=10,fmt='(A'//csl//')',iostat=ierr) line
      if (titles_found) then
        titles = line
        first = iline + 1
        titles_found = .false.
      end if
      select case (line(1:12))
        case ("Memory usage")
          titles_found = .true.
        case ("Loop time of")
          last = iline - 1
      end select
    end do
    close(10)
    if (last > first) then
      write(6,'("Properties: ",A)') trim(titles)
      write(6,'("Number of points: ",A)') trim(int2str(last - first + 1))
    else
      call error( "could not find a complete last run in file", file )
    end if
  end subroutine List_Properties

  !=================================================================================================
  !> Search for the last run output in a lammps log file and returns the index of the first line
  !! of property svals and the total number of lines of property svals.

  subroutine find_last_run( file, np, property, indx, first, nlines )
    character(*), intent(in)  :: file
    integer,      intent(in)  :: np
    character(*), intent(in)  :: property(np)
    integer,      intent(out) :: indx(np)
    integer,      intent(out) :: first, nlines
    integer :: iline, ierr, narg, i, j, last
    character(sl) :: line, arg(maxnarg)
    logical :: titles_found, run_found
    open( unit = 10, file = file, status = "old", iostat = ierr )
    if (ierr /= 0) call error( "could not open file", file )
    iline = 0
    first = 0
    last = 0
    titles_found = .false.
    run_found = .false.
    do while (ierr == 0)
      iline = iline + 1
      read(unit=10,fmt='(A'//csl//')',iostat=ierr) line
      if (titles_found) then
        call split( line, narg, arg )
        indx = 0
        do i = 1, np
          do j = 1, narg
            if (arg(j) == property(i)) indx(i) = j
          end do
        end do
        run_found = all(indx > 0)
        if (run_found) first = iline + 1
        titles_found = .false.
      end if
      select case (line(1:12))
        case ("Memory usage")
          titles_found = .true.
        case ("Loop time of")
          if (run_found) then
            last = iline - 1
            run_found = .false.
          end if
      end select
    end do
    close(10)
    if (last <= first) &
      call error( "could not find a complete last run with specified properties in file", file )
    nlines = last - first + 1
  end subroutine find_last_run

  !=================================================================================================

  subroutine Block_Analysis( np, property, npoints, value )
    integer,       intent(in) :: np, npoints
    character(sl), intent(in) :: property(np)
    real(rb),      intent(in) :: value(np,npoints)
    integer, parameter   :: nbloc = 5
    integer, parameter   :: b0(nbloc) = [2,3,5,7,11]
    integer :: i, j
    type(Block_Analyzer) :: Block(nbloc)
    ! Initialize block analyzers:
    do i = 1, nbloc
      call Block(i) % Setup( interval = 1, first = b0(i), others = 2 )
      call Block(i) % Props % Add( property )
    end do
    ! Gather properties and execute block analysis:
    do j = 1, npoints
      do i = 1, nbloc
        call Block(i) % Exec( j, value(:,j) )
      end do
    end do
    ! Flush block analysis results:
    call Block(1) % Flush( 6, separator = sep )
    do i = 2, nbloc
      write(6,'("")')
      call Block(i) % Flush( 6, separator = sep, titles = .false. )
    end do
  end subroutine Block_Analysis

  !=================================================================================================

  subroutine Print_Properties( np, property, npoints, value )
    integer,       intent(in) :: np, npoints
    character(sl), intent(in) :: property(np)
    real(rb),      intent(in) :: value(np,npoints)
    integer :: j
    type(Data_Output) :: Output
    call Output % Setup( unit = 6, interval = 1, separator = sep )
    call Output % Props % Add( property )
    do j = 1, npoints
      call Output % Exec( j, value(:,j) )
    end do
  end subroutine Print_Properties

  !=================================================================================================

  subroutine Statistics( np, property, npoints, value, print, max_value, min_value, mean, variance )
    integer,       intent(in) :: np, npoints
    character(sl), intent(in) :: property(np)
    real(rb),      intent(in) :: value(np,npoints)
    logical,       intent(in) :: print
    real(rb),      intent(out), optional :: max_value(np), min_value(np), mean(np), variance(np)
    real(rb) :: vmax(np), vmin(np), avg(np), var(np)
    real(rb) :: acc(np), acc2(np)
    integer :: i
    acc = value(:,1)
    acc2 = acc*acc
    vmin = acc
    vmax = acc
    do i = 2, npoints
      vmin = min(vmin,value(:,i))
      vmax = max(vmax,value(:,i))
      acc = acc + value(:,i)
      acc2 = acc2 + value(:,i)**2
    end do
    avg = acc/real(npoints,rb)
    var = acc2/real(npoints,rb) - avg**2
    if (print) then
      write(6,'("property'//sep//'mininum'//sep//'maximum'//sep//'mean'//sep//'std_dev")')
      do i = 1, np
        write(6,'(A,4("'//sep//'",A))') trim(property(i)), trim(real2str(vmin(i))), &
          trim(real2str(vmax(i))), trim(real2str(avg(i))), trim(real2str(sqrt(var(i))))
      end do
    end if
    if (present(min_value)) min_value = vmin 
    if (present(max_value)) max_value = vmax 
    if (present(mean)) mean = avg 
    if (present(variance)) variance = var
  end subroutine Statistics

  !=================================================================================================

  subroutine Build_Histograms( np, property, npoints, value, nbins )
    integer,       intent(in)  :: np, npoints, nbins
    character(sl), intent(in)  :: property(np)
    real(rb),      intent(in)  :: value(np,npoints)
    integer :: i, j, bin
    real(rb) :: vmin(np), vmax(np), delta(np), val
    integer :: Histo(np,nbins)
    call Statistics( np, property, npoints, value, .false., vmax, vmin )
    vmax = vmax + epsilon(vmax)
    delta = (vmax - vmin)/real(nbins,rb)
    Histo = 0
    do j = 1, npoints
      do i = 1, np
        bin = int((value(i,j) - vmin(i))/delta(i)) + 1
        histo(i,bin) = histo(i,bin) + 1
      end do
    end do
    write(6,'(3A)',advance="no") trim(property(1)), sep, "H("//trim(property(1))//")"
    do i = 2, np
      write(6,'(4A)',advance="no") sep,trim(property(i)), sep, "H("//trim(property(i))//")"
    end do
    write(6,'()')
    do bin = 1, nbins
      val = vmin(1) + (bin - 0.5_rb)*delta(1)
      write(6,'(3A)',advance="no") trim(real2str(val)), sep, trim(int2str(histo(1,bin)))
      do i = 2, np
        val = vmin(i) + (bin - 0.5_rb)*delta(i)
        write(6,'(4A)',advance="no") sep,trim(real2str(val)), sep, trim(int2str(histo(i,bin)))
      end do
      write(6,'()')
    end do
    stop
  end subroutine Build_Histograms

  !=================================================================================================

  subroutine Compute_ACF( np, property, npoints, value, window, print, norm, auto_corr_fun )
    integer,       intent(in)  :: np, npoints, window
    character(sl), intent(in)  :: property(np)
    real(rb),      intent(in)  :: value(np,npoints)
    logical,       intent(in)  :: print, norm
    real(rb),      intent(out), optional :: auto_corr_fun(np,0:window)
    integer  :: i, j, delta
    real(rb) :: acf(np,0:window), avg(np)
    if (window > npoints-1) call error( "Window size cannot be larger than", int2str(npoints-1) )
    if (norm) call Statistics( np, property, npoints, value, .false., mean = avg )
    acf = 0.0_rb
    do delta = 0, window
      do i = 1, npoints-delta
        if (norm) then
          acf(:,delta) = acf(:,delta) + (value(:,i) - avg)*(value(:,i+delta) - avg)
        else
          acf(:,delta) = acf(:,delta) + value(:,i)*value(:,i+delta)
        end if
      end do
      acf(:,delta) = acf(:,delta)/real(npoints-delta,rb)
    end do
    if (norm) forall(delta=0:window) acf(:,delta) = acf(:,delta)/acf(:,0)
    if (print) then
      write(6,'("delta")',advance="no")
      do i = 1, np
        write(6,'(2A)',advance="no") sep, "acf<"//trim(property(i))//">"
      end do
      write(6,'()')
      do j = 0, window
        write(6,'(A)',advance="no") trim(int2str(j))
        do i = 1, np
          write(6,'(2A)',advance="no") sep, trim(real2str(acf(i,j)))
        end do
        write(6,'()')
      end do
    end if
    if (present(auto_corr_fun)) auto_corr_fun = acf
  end subroutine Compute_ACF

  !=================================================================================================

end program post_lammps
