program post_lammps
use mData_Proc
use mString
implicit none

character, parameter :: sep = ","
integer, parameter   :: nbloc = 5
integer, parameter   :: b0(nbloc) = [2,3,5,7,11]
type(Block_Analyzer) :: Block(nbloc)

integer, parameter :: maxnarg = 15
integer            :: narg
character(sl)      :: arg(maxnarg)

integer                    :: np
character(sl), allocatable :: property(:)
integer,       allocatable :: indx(:)
integer                    :: npoints
real(rb),      allocatable :: value(:)
real(rb),      allocatable :: vmin(:), vmax(:), delta(:)
integer,       allocatable :: histo(:,:), bin

integer       :: i, j, first, nbins
character(sl) :: infile, action, line
character(12) :: C
real(rb)      :: rval

! Check provided command-line arguments:
if (iargc() < 3) then
  write(6,'("Usage: post_lammps action <args> <file name> <property 1> <property 2> ...")')
  write(6,'("  action = block or extract or histogram")')
  write(6,'("    block args = none")')
  write(6,'("    extract args = none")')
  write(6,'("    histogram args = number-of-bins")')
  stop
end if

call getarg( 1, action )

if (trim(action) == "histogram") then
  narg = 1
  call getarg( 2, line )
  nbins = str2int( line )
else if ( (trim(action) == "block").or.(trim(action) == "extract") ) then
  narg = 0
else
  write(*,'("Error: specified action was not recognized")')
  stop
end if

call getarg( narg+2, infile )
np = iargc() - (narg+2)
allocate( property(np), indx(np) )
do i = 1, np
  call getarg( i+narg+2, property(i) )
end do

! Search for the last run containing the specified properties:
call find_last_run( infile, np, property, indx, first, npoints )

select case (trim(action))
  case ("block")
    ! Get ready for block analysis:
    if (trim(action)=="block") then
      do i = 1, nbloc
        call Block(i) % Setup( interval = 1, first = b0(i), others = 2 )
        call Block(i) % Props % Add( property )
      end do
    end if
  case ("extract")
    write(6,'(A)',advance="no") trim(property(1))
    do i = 2, np
      write(6,'(A,A)',advance="no") sep, trim(property(i))
    end do
    write(6,'()')
  case ("histogram")
    allocate( vmax(np), vmin(np), delta(np), histo(np,nbins) )
    vmin = huge(0.0_rb)
    vmax = -huge(0.0_rb)
    histo = 0
end select

! Perform the blocking analysis:
open( unit = 10, file = infile, status = "old" )
write(C,*) first-2
read(10,'('//C//'/)')
allocate( value(np) )
do j = 1, npoints
  read(unit=10,fmt='(A'//csl//')') line
  call split( line, narg, arg )
  do i = 1, np
    value(i) = str2real( arg(indx(i)) )
  end do
  select case (trim(action))
    case ("block")
      do i = 1, nbloc
        call Block(i)%Exec( 1, value )
      end do
    case ("extract")
      write(6,'(A)',advance="no") trim(arg(indx(1)))
      do i = 2, np
        write(6,'(A,A)',advance="no") sep, trim(arg(indx(i)))
      end do
      write(6,'()')
     case ("histogram")
       vmin = min(vmin,value)
       vmax = max(vmax,value)
   end select
end do
close(10)

select case (trim(action))
  case ("block")

    ! Flush block analysis results:
    call Block(1)%Flush( 6, separator = sep )
    do i = 2, nbloc
      write(6,'("")')
      call Block(i)%Flush( 6, separator = sep, titles = .false. )
    end do

  case ("histogram")

    vmax = vmax + epsilon(vmax)
    delta = (vmax - vmin)/real(nbins,rb)
    open( unit = 10, file = infile, status = "old" )
    write(C,*) first-2
    read(10,'('//C//'/)')
    do j = 1, npoints
      read(unit=10,fmt='(A'//csl//')') line
      call split( line, narg, arg )
      do i = 1, np
        rval = str2real( arg(indx(i)) )
        bin = int((rval - vmin(i))/delta(i)) + 1
        histo(i,bin) = histo(i,bin) + 1
      end do
    end do
    close(10)
    write(6,'(3A)',advance="no") trim(property(1)), sep, "H("//trim(property(1))//")"
    do i = 2, np
      write(6,'(4A)',advance="no") sep,trim(property(i)), sep, "H("//trim(property(i))//")"
    end do
    write(6,'()')
    do bin = 1, nbins
      write(6,'(3A)',advance="no") trim(real2str(vmin(1)+(bin-0.5_rb)*delta(1))), sep, trim(int2str(histo(1,bin)))
      do i = 2, np
        write(6,'(4A)',advance="no") sep,trim(real2str(vmin(i)+(bin-0.5_rb)*delta(i))), sep, trim(int2str(histo(i,bin)))
      end do
      write(6,'()')
    end do
end select

contains

  !=================================================================================================
  !> Search for the last run output in a lammps log file and returns the index of the first line
  !! of property values and the total number of lines of property values.

  subroutine find_last_run( file, np, property, indx, first, npoints )
    character(*), intent(in)  :: file
    integer,      intent(in)  :: np
    character(*), intent(in)  :: property(np)
    integer,      intent(out) :: indx(np)
    integer,      intent(out) :: first, npoints
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
    npoints = last - first + 1
  end subroutine find_last_run

  !=================================================================================================

end program post_lammps
