program postlammps

!   This file is part of Foobar.

!    Foobar is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.

!    Foobar is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

use mData_Proc
use mString
implicit none

integer, parameter :: maxnarg = 15
integer            :: narg
character(sl)      :: arg(maxnarg)

integer            :: iarg, n

character(6), parameter :: keyword(2) = ["-every","-delim"]
integer :: every
character :: delim

integer                    :: np, nlines, npoints
character(sl), allocatable :: property(:)
integer,       allocatable :: indx(:)
real(rb),      allocatable :: value(:,:)

integer       :: i, j, first, nbins, window
character(sl) :: infile, action, line
logical       :: exist
character(12) :: C

! Defaults:
every = 1
delim = " "

! Read options:
iarg = 1
if (iargc() < 1) call Usage_Message
call getarg( iarg, line )
do while (any(keyword == line))
  select case (trim(line))
    case ("-every")
      iarg = iarg + 1
      call getarg( iarg, line )
      every = str2int( line )
      if (every < 1) call error( "Unacceptable parameter 'every'" )
    case ("-delim")
      iarg = iarg + 1
      call getarg( iarg, line )
      select case (trim(line))
        case ("comma"); delim = ","
        case ("space"); delim = " "
        case ("semicolon"); delim = ";"
        case default; call error( "Unacceptable delimiter" )
      end select
  end select
  iarg = iarg + 1
  call getarg( iarg, line )
end do

! First argument is the specified action:
call getarg( iarg, action )

! Check specified action:
select case (trim(action))
  case ("acf","acfn","histo"); n = 1
  case ("block","ineff","print","props","stats"); n = 0
  case default; call error( "Unrecognized action", action )
end select

! Read action arguments:
if (iargc() < iarg + n) call Usage_Message
select case (trim(action))
  case ("acf","acfn")
    iarg = iarg + 1
    call getarg( iarg, line )
    window = str2int( line )
    if (window <= 0) call error( "Unacceptable maximum window size" )
  case ("histo")
    iarg = iarg + 1
    call getarg( iarg, line )
    nbins = str2int( line )
    if (nbins <= 0) call error( "Unacceptable number of bins" )
end select

! Read file name:
iarg = iarg + 1
if (iargc() < iarg) call Usage_Message
call getarg( iarg, infile )
inquire( file = infile, exist = exist )
if (.not.exist) call error( "specified file <", infile, "> does not exist." )

! Perform listing action, if specified:
if (trim(action) == "props") then
  call List_Properties( infile )
  stop
end if

! Read property names:
np = iargc() - iarg
if (np == 0) call error( "No properties specified." )
allocate( property(np), indx(np) )
do i = 1, np
  call getarg( iarg + i, property(i) )
end do

! Search for the last run containing the specified properties:
call find_last_run( infile, np, property, indx, first, nlines )

! Calculate the actual number of points:
npoints = 0
do j = 1, nlines
  if (mod(j-1,every) == 0) npoints = npoints + 1
end do
allocate( value(np,npoints) )

! Read property values from the log file:
open( unit = 10, file = infile, status = "old" )
write(C,*) first-2
read(10,'('//C//'/)')
npoints = 0
do j = 1, nlines
  read(unit=10,fmt='(A'//csl//')') line
  if (mod(j-1,every) == 0) then
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
    call Compute_ACF( np, property, npoints, value, window, norm = .false. )
  case ("acfn")
    call Compute_ACF( np, property, npoints, value, window, norm = .true. )
  case ("block")
    call Block_Analysis( np, property, npoints, value )
  case ("histo")
    call Build_Histograms( np, property, npoints, value, nbins )
  case ("ineff")
    call Correlation_Analisys( np, property, npoints, value, print = .true. )
  case ("print")
    call Print_Properties( np, property, npoints, value )
  case ("stats")
    call Statistics( np, property, npoints, value, print = .true. )
end select

contains

  !=================================================================================================

  subroutine Usage_Message
    write(6,'("Usage: post_lammps [options] action [args] file-name property-1 [property-2 ...]")')
    write(6,'("  action = acf or acfn or block or histo or ineff or print or props or stats")')
    write(6,'("    acf   args = window")')
    write(6,'("    acfn  args = window")')
    write(6,'("    block args = none")')
    write(6,'("    histo args = nbins")')
    write(6,'("    ineff args = none")')
    write(6,'("    print args = none")')
    write(6,'("    props args = none")')
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
    call Block(1) % Flush( 6, separator = delim )
    do i = 2, nbloc
      write(6,'("")')
      call Block(i) % Flush( 6, separator = delim, titles = .false. )
    end do
  end subroutine Block_Analysis

  !=================================================================================================

  subroutine Print_Properties( np, property, npoints, value )
    integer,       intent(in) :: np, npoints
    character(sl), intent(in) :: property(np)
    real(rb),      intent(in) :: value(np,npoints)
    integer :: j
    type(Data_Output) :: Output
    call Output % Setup( unit = 6, interval = 1, separator = delim )
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
      write(6,'("property'//delim//'mininum'//delim//'maximum'//delim//'mean'//delim//'std_dev")')
      do i = 1, np
        write(6,'(A,4("'//delim//'",A))') trim(property(i)), trim(real2str(vmin(i))), &
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
    write(6,'(3A)',advance="no") trim(property(1)), delim, "H("//trim(property(1))//")"
    do i = 2, np
      write(6,'(4A)',advance="no") delim,trim(property(i)), delim, "H("//trim(property(i))//")"
    end do
    write(6,'()')
    do bin = 1, nbins
      val = vmin(1) + (bin - 0.5_rb)*delta(1)
      write(6,'(3A)',advance="no") trim(real2str(val)), delim, trim(int2str(histo(1,bin)))
      do i = 2, np
        val = vmin(i) + (bin - 0.5_rb)*delta(i)
        write(6,'(4A)',advance="no") delim,trim(real2str(val)), delim, trim(int2str(histo(i,bin)))
      end do
      write(6,'()')
    end do
  end subroutine Build_Histograms

  !=================================================================================================

  subroutine Compute_ACF( np, property, npoints, value, window, norm )
    integer,       intent(in)  :: np, npoints, window
    character(sl), intent(in)  :: property(np)
    real(rb),      intent(in)  :: value(np,npoints)
    logical,       intent(in)  :: norm
    integer  :: i, j, delta
    real(rb) :: acf(np,0:window), avg(np)
    logical  :: allpos(np)
    if (window > npoints-1) call error( "Maximum allowable window size is", int2str(npoints) )
    if (norm) call Statistics( np, property, npoints, value, .false., mean = avg )
    acf = 0.0_rb
    allpos = .true.
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
    write(6,'("delta")',advance="no")
    do i = 1, np
      write(6,'(2A)',advance="no") delim, "acf<"//trim(property(i))//">"
    end do
    write(6,'()')
    do j = 0, window
      write(6,'(A)',advance="no") trim(int2str(j))
      do i = 1, np
        write(6,'(2A)',advance="no") delim, trim(real2str(acf(i,j)))
       end do
       write(6,'()')
    end do
  end subroutine Compute_ACF

  !=================================================================================================

  subroutine Correlation_Analisys( np, property, npoints, value, print )
    integer,       intent(in)  :: np, npoints
    character(sl), intent(in)  :: property(np)
    real(rb),      intent(in)  :: value(np,npoints)
    logical,       intent(in)  :: print
    integer  :: i, delta
    real(rb) :: acf(np,0:npoints), avg(np), inv_n, g(np), error(np)
    logical  :: positive(np)
    call Statistics( np, property, npoints, value, .false., mean = avg )
    acf = zero
    positive = .true.
    delta = -1
    do while ((delta < npoints).and.any(positive))
      delta = delta + 1
      do i = 1, npoints-delta
        acf(:,delta) = acf(:,delta) + (value(:,i) - avg)*(value(:,i+delta) - avg)
      end do
      acf(:,delta) = acf(:,delta)/real(npoints-delta,rb)
      where (positive) positive = acf(:,delta) > zero
      where (.not.positive) acf(:,delta) = zero
    end do
    forall (i=1:delta) acf(:,i) = acf(:,i)/acf(:,0)
    g = zero
    inv_n = one/real(npoints,rb)
    do i = 1, delta
      g = g + (one - delta*inv_n)*acf(:,i)
    end do
    g = one + two*g
    error = sqrt(g*acf(:,0)*inv_n)
    if (print) then
      do i = 1, np
        write(6,'(50("-"))')
        write(6,'("Property............: ",A)') trim(property(i))
        write(6,'("Stat. inefficency...: ",A)') trim(real2str(g(i)))
        write(6,'("Average.............: ",A)') trim(real2str(avg(i)))
        write(6,'("Uncertainty.........: ",A)') trim(real2str(error(i)))
      end do
      write(6,'(50("-"))')
    end if
  end subroutine Correlation_Analisys

  !=================================================================================================

end program postlammps
