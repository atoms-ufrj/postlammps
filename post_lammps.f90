program post_lammps

use mData_Proc
use mString
implicit none

integer,   parameter :: maxnarg = 15
integer              :: narg
character(sl)        :: arg(maxnarg)

integer                    :: np, nlines, npoints
character(sl), allocatable :: property(:)
integer,       allocatable :: indx(:)
real(rb),      allocatable :: value(:,:)

real(rb),      allocatable :: vmin(:), vmax(:), mean(:), variance(:)

integer       :: i, j, first, nbins, skip
character(sl) :: infile, action, line
character(12) :: C

! Check provided command-line arguments:
if (iargc() < 3) call Usage_Message

! First argument is the specified action:
call getarg( 1, action )

! Read action args:
select case (trim(action))
  case ("histogram")
    narg = 1
    call getarg( 2, line )
    nbins = str2int( line )
    if (nbins <= 0) stop "Unacceptable number of bins"
  case ("block","extract","stats")
    narg = 0
  case default
    call Usage_Message
end select

! Read file name:
call getarg( narg+2, infile )

! Read skip parameter:
call getarg( narg+3, line )
skip = str2int( line )
if (skip < 0) stop "Unacceptable value of skip parameter"

! Read property names:
np = iargc() - (narg+3)
allocate( property(np), indx(np) )
do i = 1, np
  call getarg( i+narg+3, property(i) )
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
      value(i,npoints) = str2real( arg(indx(i)) )
    end do
  end if
end do
close(10)

! Perform specified action:
select case (trim(action))
  case ("block")
    call Block_Analysis( np, property, npoints, value )
  case ("extract")
    call Print_Properties( np, property, npoints, value )
  case ("histogram")
    call Build_Histograms( np, property, npoints, value, nbins )
  case ("stats")
    allocate( vmax(np), vmin(np), mean(np), variance(np) )
    call Statistics( np, property, npoints, value, vmax, vmin, mean, variance )
    write(6,'("property'//sep//'mininum'//sep//'maximum'//sep//'mean'//sep//'std_dev")')
    do i = 1, np
      write(6,'(A,4("'//sep//'",A))') trim(property(i)), trim(real2str(vmin(i))), &
        trim(real2str(vmax(i))), trim(real2str(mean(i))), trim(real2str(sqrt(variance(i))))
    end do
end select

contains

  !=================================================================================================

  subroutine Usage_Message
    write(6,'("Usage: post_lammps action [args] file-name skip property-1 [property-2] ...")')
    write(6,'("  action = block or extract or histogram or stats")')
    write(6,'("    block args = none")')
    write(6,'("    extract args = none")')
    write(6,'("    histogram args = nbins")')
    write(6,'("    stats args = none")')
    stop
  end subroutine Usage_Message

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
    integer :: i, j
    write(6,'(A)',advance="no") trim(property(1))
    do i = 2, np
      write(6,'(A,A)',advance="no") sep, trim(property(i))
    end do
    write(6,'()')
    do j = 1, npoints
      write(6,'(A)',advance="no") trim(real2str(value(1,j)))
      do i = 2, np
        write(6,'(A,A)',advance="no") sep, trim(real2str(value(i,j)))
      end do
      write(6,'()')
    end do
  end subroutine Print_Properties

  !=================================================================================================

  subroutine Statistics( np, property, npoints, value, max_value, min_value, mean, variance )
    integer,       intent(in)  :: np, npoints
    character(sl), intent(in)  :: property(np)
    real(rb),      intent(in)  :: value(np,npoints)
    real(rb),      intent(out) :: max_value(np), min_value(np), mean(np), variance(np)
    real(rb) :: acc(np), acc2(np)
    integer :: i
    acc = value(:,1)
    acc2 = acc*acc
    min_value = acc
    max_value = acc
    do i = 2, npoints
      min_value = min(min_value,value(:,i))
      max_value = max(max_value,value(:,i))
      acc = acc + value(:,i)
      acc2 = acc2 + value(:,i)**2
    end do
    mean = acc/real(npoints,rb)
    variance = acc2/real(npoints,rb) - mean**2
  end subroutine Statistics

  !=================================================================================================

  subroutine Build_Histograms( np, property, npoints, value, nbins )
    integer,       intent(in)  :: np, npoints, nbins
    character(sl), intent(in)  :: property(np)
    real(rb),      intent(in)  :: value(np,npoints)
    integer :: i, j, bin
    real(rb) :: vmin(np), vmax(np), mean(np), variance(np), delta(np), val
    integer :: Histo(np,nbins)
    call Statistics( np, property, npoints, value, vmax, vmin, mean, variance )
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
  end subroutine Build_Histograms

  !=================================================================================================

end program post_lammps
