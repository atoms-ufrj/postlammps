program postlammps

!   This file is part of Postlammps.
!
!    Postlammps is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    Postlammps is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with Postlammps.  If not, see <http://www.gnu.org/licenses/>.

use mData_Proc
use mString
implicit none

integer, parameter :: maxnarg = 50
integer            :: narg
character(sl)      :: arg(maxnarg)

integer            :: input, iarg, argcount, n

character(3), parameter :: keyword(8) = [character(3) :: "-e","-d","-in","-c","-p","-nt","-r","-mm"]
integer :: every
character :: delim

integer                    :: np, nlines, npoints
character(sl), allocatable :: property(:)
integer,       allocatable :: indx(:)
real(rb),      allocatable :: value(:,:)

integer       :: i, j, nbins, window, initial, final, bmin, bmax, bsize
character(sl) :: infile, action, line
logical       :: read_from_file, props, plain, print_titles, range = .false., openmm
real(rb)      :: percentage

type tLine
  character(sl) :: line = ""
  type(tLine), pointer :: next => null()
end type tLine
type(tLine), pointer :: StdIn => null()
type(tLine), pointer :: titles => null()
type(tLine), pointer :: current => null()

! Defaults:
plain = .false.
input = 5
every = 1
delim = " "
print_titles = .true.
read_from_file = .false.
percentage = 100.0_rb
openmm = .false.

! Read options:
argcount = command_argument_count()
iarg = 1
if (argcount < 1) call Usage_Message
call get_command_argument( iarg, line )
do while (any(keyword == line))
  select case (trim(line))
    case ("-e")
      iarg = iarg + 1
      call get_command_argument( iarg, line )
      every = str2int( line )
      if (every < 1) call error( "Unacceptable parameter 'every'" )
    case ("-d")
      iarg = iarg + 1
      call get_command_argument( iarg, line )
      select case (trim(line))
        case ("comma"); delim = ","
        case ("space"); delim = " "
        case ("semicolon"); delim = ";"
        case ("tab"); delim = achar(9)
        case default; call error( "Unacceptable delimiter" )
      end select
    case ("-in")
      iarg = iarg + 1
      call get_command_argument( iarg, infile )
      inquire( file = infile, exist = read_from_file )
      if (.not.read_from_file) call error( "Specified input file ", infile, "does not exist" )
    case ("-c")
      iarg = iarg + 1
      call get_command_argument( iarg, line )
      percentage = str2real(line)
      if ((percentage <= 0.0_rb).or.(percentage > 100.0_rb)) call error( "wrong -c option definition" )
    case ("-r")
      iarg = iarg + 1
      call get_command_argument( iarg, line )
      initial = str2int(line)
      iarg = iarg + 1
      call get_command_argument( iarg, line )
      final = str2int(line)
      range = .true.
    case ("-p")
      plain = .true.
    case ("-nt")
      print_titles = .false.
    case ("-mm")
      openmm = .true.
      plain = .true.
  end select
  iarg = iarg + 1
  call get_command_argument( iarg, line )
end do

! First argument is the specified action:
call get_command_argument( iarg, action )

! Check specified action:
select case (trim(action))
  case ("batch","obm"); n = 2
  case ("acfun","fluct","histo","mvavg"); n = 1
  case ("block","ineff","print","props","sampl","stats"); n = 0
  case default; call error( "Unrecognized action", action )
end select
props = trim(action) == "props"

! Read action arguments:
if (argcount < iarg + n) call Usage_Message
select case (trim(action))
  case ("batch","obm")
    call get_command_argument( iarg+1, line )
    bmin = str2int( line )
    call get_command_argument( iarg+2, line )
    bmax = str2int( line )
    if ((bmin <= 0).or.(bmax <= 0).or.(bmin > bmax)) then
      call error( "Unacceptable minimum and maximum block sizes" )
    end if
  case ("acfun","fluct","mvavg")
    call get_command_argument( iarg+1, line )
    window = str2int( line )
    if (window <= 0) call error( "Unacceptable maximum window size" )
  case ("histo")
    call get_command_argument( iarg+1, line )
    nbins = str2int( line )
    if (nbins <= 0) call error( "Unacceptable number of bins" )
end select
iarg = iarg + n

! Read property names:
np = argcount - iarg
if ((np == 0).and.(.not.props)) call error( "no properties have been specified")
allocate( property(np), indx(np) )
do i = 1, np
  call get_command_argument( iarg + i, property(i) )
end do

! Read data from the standard input or the specified input file:
if (read_from_file) open( newunit = input, file = infile, status = "old" )
if (plain) then
  call read_plain_file( input, np, property, indx, nlines, props )
else
  call read_lammps_log( input, np, property, indx, nlines, props )
end if
if (read_from_file) close( input )
if (props) then
  write(6,'("Properties: ",A)') trim(titles % line)
  write(6,'("Number of points: ",A)') trim(int2str(nlines))
  stop
end if

! Calculate the actual number of points:
npoints = 0
if (range) then
  if ((initial < 1).or.(final > nlines)) call error( "invalid range" )
else
  initial = int((1.0_rb - 0.01_rb*percentage)*nlines) + 1
  final = nlines
end if
do j = initial, final
  if (mod(j-1,every) == 0) npoints = npoints + 1
end do
allocate( value(np,npoints) )

! Read property values:
current => titles
npoints = 0
do j = 1, initial-1
  current => current % next
end do
do j = initial, final
  current => current % next
  if (mod(j-1,every) == 0) then
    call split( current % line, narg, arg )
    npoints = npoints + 1
    do i = 1, np
      read(arg(indx(i)),*) value(i,npoints)
    end do
  end if
end do

! Perform specified action:
select case (trim(action))
  case ("acfun")
    call Compute_ACF( np, property, npoints, value, window, norm = .false. )
  case ("block")
    call Block_Analysis( np, property, npoints, value )
  case ("fluct")
    call Compute_ACF( np, property, npoints, value, window, norm = .true. )
  case ("histo")
    call Build_Histograms( np, property, npoints, value, nbins )
  case ("ineff")
    call Correlation_Analisys( np, property, npoints, value, print = .true. )
  case ("print")
    call Print_Properties( np, property, npoints, value )
  case ("sampl")
    call Subsample( np, property, npoints, value )
  case ("stats")
    call Statistics( np, property, npoints, value, print = .true. )
  case ("mvavg")
    call Print_Properties( np, property, npoints, value, moving_average = .true. )
  case ("batch")
    do bsize = bmin, bmax
      call Batch_Means( np, property, npoints, value, bsize )
    end do
  case ("obm")
    do bsize = bmin, bmax
!    bsize = floor(sqrt(real(npoints,rb)))
      call Overlapping_Batch_Means( np, property, npoints, value, bsize )
    end do
end select

contains

  !=================================================================================================

  subroutine Usage_Message
    write(6,'("Usage: postlammps [options] action [args] property-1 [property-2 ...]")')
    write(6,'()')
    write(6,'("  action = acfun / block / fluct / histo / ineff / print / props / sampl / stats")')
    write(6,'("    acfun <maxlag>: Computes autocorrelation functions (ACF) from 0 to maxlag")')
    write(6,'("    block: Performs normalization-group blocking analysis")')
    write(6,'("    fluct <maxlag>: Computes normalized fluctuation ACF from 0 to maxlag")')
    write(6,'("    histo <nbins>: Builds histograms with specified number of bins")')
    write(6,'("    ineff: Computes statistical inefficiencies and uncertainties")')
    write(6,'("    print: Prints the values of the selected properties")')
    write(6,'("    props: Lists all properties available in the log file")')
    write(6,'("    sampl: Samples uncorrelated points from the original data")')
    write(6,'("    stats: Computes basic statistics")')
    write(6,'("    mvavg <window>: Prints moving-window averages of the selected properties")')
    write(6,'("    obm <bmin> <bmax>: Performs overlapping batch mean analysis")')
    write(6,'()')
    write(6,'("  options = -e / -d / -in / -c / -p / -nt / -r")')
    write(6,'("    -in <file>: Specifies the name of the log file to be processed")')
    write(6,'("    -p: Tells postlammps to read a plain data file instead of a lammps log file")')
    write(6,'("    -mm: Tells postlammps to read an OpenMM state-data report")')
    write(6,'("    -e <n>: Skips n lines between property inputs")')
    write(6,'("    -d <delim>: Specifies the item delimiter used for output")')
    write(6,'("       delim = space or comma or semicolon or tab")')
    write(6,'("    -nt: Does not print property titles")')
    write(6,'("    -c <X>: Consider only the last X% of data")')
    write(6,'("    -r <X> <Y>: Consider only data within a specified range from X to Y")')
    stop
  end subroutine Usage_Message

  !=================================================================================================

  subroutine read_plain_file( file, np, property, indx, nlines, props )
    integer,      intent(in)  :: file
    integer,      intent(in)  :: np
    character(*), intent(in)  :: property(np)
    integer,      intent(out) :: indx(np)
    integer,      intent(out) :: nlines
    logical,      intent(in)  :: props
    integer :: ierr, narg, i, j
    character(sl) :: line, arg(maxnarg)
    logical :: titles_found
    allocate( titles )
    read(input,'(A'//csl//')') titles % line
    if (openmm) then
      do while (titles % line(1:2) /= '#"')
        read(input,'(A'//csl//')') titles % line
      end do
      titles % line = translateFromOpenMM(titles % line)
    end if
    call split( titles % line, narg, arg )
    indx = 0
    do i = 1, np
      do j = 1, narg
        if (arg(j) == property(i)) indx(i) = j
      end do
    end do
    titles_found = all(indx > 0).or.props
    if (titles_found) then
      current => titles
      read(input,'(A'//csl//')',iostat=ierr) line
      nlines = 0
      do while (ierr == 0)
        nlines = nlines + 1
        allocate( current % next )
        current => current % next
        current % line = line
        read(input,'(A'//csl//')',iostat=ierr) line
      end do
    else
      call error( "could not find the specified properties" )
    end if
  end subroutine read_plain_file

  !=================================================================================================

  subroutine read_lammps_log( file, np, property, indx, nlines, props )
    integer,      intent(in)  :: file
    integer,      intent(in)  :: np
    character(*), intent(in)  :: property(np)
    integer,      intent(out) :: indx(np)
    integer,      intent(out) :: nlines
    logical,      intent(in)  :: props
    integer :: iline, ierr, narg, i, j, first, last
    character(sl) :: line, arg(maxnarg)
    logical :: titles_found, run_found
    iline = 0
    first = 0
    last = 0
    titles_found = .false.
    run_found = .false.
    ierr = 0
    do while (ierr == 0)
      iline = iline + 1
      read(input,'(A'//csl//')',iostat=ierr) line
      if (ierr == 0) then
        if (associated(current)) then
          allocate( current % next )
          current => current % next
        else
          allocate( current )
          StdIn => current
        end if
        current % line = line
        if (titles_found) then
          titles => current
          call split( line, narg, arg )
          indx = 0
          do i = 1, np
            do j = 1, narg
              if (arg(j) == property(i)) indx(i) = j
            end do
          end do
          run_found = all(indx > 0).or.props
          if (run_found) first = iline + 1
          titles_found = .false.
        end if
        select case (line(1:12))
          case ("Memory usage","Per MPI rank")
            titles_found = .true.
          case ("Loop time of")
            if (run_found) then
              last = iline - 1
              run_found = .false.
            end if
        end select
      end if
    end do
    if (last <= first) &
      call error( "could not find a complete last run with specified properties" )
    nlines = last - first + 1
  end subroutine read_lammps_log

  !=================================================================================================

  function translateFromOpenMM(props) result( new )
    character(sl), intent(in) :: props
    character(sl)             :: new
    new = replace(props, '"', '')
    new = replace(new, 'Progress (%)', 'Progress')
    new = replace(new, '#Step', 'Step')
    new = replace(new, 'Time (ps)', 'Time')
    new = replace(new, 'Potential Energy (kJ/mole)', 'PotEng')
    new = replace(new, 'Kinetic Energy (kJ/mole)', 'KinEng')
    new = replace(new, 'Total Energy (kJ/mole)', 'TotEng')
    new = replace(new, 'Temperature (K)', 'Temp')
    new = replace(new, 'Box Volume (nm^3)', 'Vol')
    new = replace(new, 'Density (g/mL)', 'Density')
    new = replace(new, 'Speed (ns/day)', 'Speed')
    new = replace(new, 'Elapsed Time (s)', 'Elapsed')
    new = replace(new, 'Time Remaining', 'Remaining')
    new = replace(new, 'Atomic Virial (kJ/mole)', 'Virial')
    new = replace(new, 'Nonbonded Virial (kJ/mole)', 'NBVirial')
    new = replace(new, 'Atomic Pressure (atm)', 'Press')
    new = replace(new, 'Molecular Virial (kJ/mole)', 'MolVirial')
    new = replace(new, 'Molecular Pressure (atm)', 'MolPress')
    new = replace(new, 'Molecular Kinetic Energy (kJ/mole)', 'MolKinEng')
    new = replace(new, 'alchemical_vdw_energy', 'AlchemVdwEng')
    new = replace(new, 'alchemical_coulomb_energy', 'AlchemCoulEng')
    new = replace(new, 'Coulomb Energy (kJ/mole)', 'Ecoul')
    new = replace(new, achar(9), ' ')
    new = replace(new, '"', '')
    new = replace(new, ' (', '_')
    new = replace(new, ')', '')
  end function translateFromOpenMM

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

  subroutine Print_Properties( np, property, npoints, value, interval, moving_average)
    integer,       intent(in) :: np, npoints
    character(sl), intent(in) :: property(np)
    real(rb),      intent(in) :: value(np,npoints)
    integer,       intent(in), optional :: interval
    logical,       intent(in), optional :: moving_average
    integer :: j, n, nacc
    logical :: avg
    real(rb) :: acc(np)

    if (present(interval)) then
      n = interval
    else
      n = 1
    endif
    avg = present(moving_average)
    if (avg) avg = moving_average
    if (print_titles) call write_str( 6, property, delim )

    if (avg) then
      acc = 0.0_rb
      nacc = 0
      do j = 1, npoints, n
        if (nacc < window) then
          acc = acc + value(:,j)
          nacc = nacc + 1
        else
          acc = acc + value(:,j) - value(:,j-window)
        endif
        call write_str(6, real2str(acc/nacc), delim )
      end do
    else
      do j = 1, npoints, n
        call write_str(6, real2str(value(:,j)), delim )
      end do
    endif
  end subroutine Print_Properties

  !=================================================================================================

  subroutine Statistics( np, property, npoints, value, print, max_value, min_value, mean, variance )
    integer,       intent(in) :: np, npoints
    character(sl), intent(in) :: property(np)
    real(rb),      intent(in) :: value(np,npoints)
    logical,       intent(in) :: print
    real(rb),      intent(out), optional :: max_value(np), min_value(np), mean(np), variance(np)
    real(rb) :: vmax(np), vmin(np), avg(np), var(np), slope(np)
    real(rb) :: acc(np), acc2(np), accprog(np), N
    integer :: i
    acc = value(:,1)
    acc2 = acc*acc
    accprog = acc
    vmin = acc
    vmax = acc
    do i = 2, npoints
      vmin = min(vmin,value(:,i))
      vmax = max(vmax,value(:,i))
      acc = acc + value(:,i)
      acc2 = acc2 + value(:,i)**2
      accprog = accprog + real(i,rb)*value(:,i)
    end do
    N = real(npoints,rb)
    avg = acc/N
    var = acc2/N - avg**2
    slope = 12.0_rb*(accprog/N - 0.5_rb*(N+1)*avg)/(N*N - 1.0_rb)
    if (print) then
      if (print_titles) call write_str( 6, ["property","mean    ","std_dev ",        &
                                            "minimum ","maximum ","slope   "], delim )
      do i = 1, np
        write(6,'(A,5("'//delim//'",A))') trim(property(i)), trim(real2str(avg(i))), &
          trim(real2str(sqrt(var(i)))), trim(real2str(vmin(i))), trim(real2str(vmax(i))), &
          trim(real2str(slope(i)))
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

  subroutine Correlation_Analisys( np, property, npoints, value, print, stat_ineff )
    integer,       intent(in)  :: np, npoints
    character(sl), intent(in)  :: property(np)
    real(rb),      intent(in)  :: value(np,npoints)
    logical,       intent(in)  :: print
    real(rb),      intent(out), optional :: stat_ineff(np)
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
      g = g + (one - i*inv_n)*acf(:,i)
    end do
    g = one + two*g
    error = sqrt(g*acf(:,0)*inv_n)
    if (print) then
      if (print_titles) call write_str( 6, ["property   ","average    ","uncertainty","stat_ineff "], delim )
      do i = 1, np
        write(6,'(A,3("'//delim//'",A))') trim(property(i)),        &
                                          trim(real2str(avg(i))),   &
                                          trim(real2str(error(i))), &
                                          trim(real2str(g(i)))
      end do
    end if
    if (present(stat_ineff)) stat_ineff = g
  end subroutine Correlation_Analisys

  !=================================================================================================

  subroutine Subsample( np, property, npoints, value)
    integer,       intent(in)  :: np, npoints
    character(sl), intent(in)  :: property(np)
    real(rb),      intent(in)  :: value(np,npoints)
    real(rb) :: g(np)
    call Correlation_Analisys( np, property, npoints, value, .false., g )
    call Print_Properties( np, property, npoints, value, ceiling(maxval(g)) )
  end subroutine Subsample

  !=================================================================================================

  subroutine Batch_Means( np, property, npoints, value, bsize )
    integer,       intent(in)  :: np, npoints
    character(sl), intent(in)  :: property(np)
    real(rb),      intent(in)  :: value(np,npoints)
    integer,       intent(in)  :: bsize
    integer :: i, nblocks
    real(rb) :: avg(np), acc2(np), var(np)
    avg = sum(value,2)/npoints
    nblocks = npoints/bsize
    acc2 = 0.0_rb
    do i = 1, nblocks
      acc2 = acc2 + (sum(value(:,(i-1)*bsize+1:i*bsize),2)/bsize - avg)**2
    end do
    var = acc2/(real(nblocks,rb)*real(nblocks - 1,rb))
    write(6,'(A)') trim(join([int2str(bsize),real2str(var)],delim))
  end subroutine Batch_Means

  !=================================================================================================

  subroutine Overlapping_Batch_Means( np, property, npoints, value, bsize )
    integer,       intent(in)  :: np, npoints
    character(sl), intent(in)  :: property(np)
    real(rb),      intent(in)  :: value(np,npoints)
    integer,       intent(in)  :: bsize
    integer :: i
    real(rb) :: avg(np), blocksum(np), acc2(np), var(np)
    avg = sum(value,2)/npoints
    blocksum = sum(value(:,1:bsize),2)
    acc2 = (blocksum/bsize - avg)**2
    do i = 1, npoints-bsize
      blocksum = blocksum - value(:,i) + value(:,i+bsize)
      acc2 = acc2 + (blocksum/bsize - avg)**2
    end do
    var = bsize*acc2/(real(npoints - bsize,rb)*real(npoints - bsize + 1,rb))
    write(6,'(A)') trim(join([int2str(bsize),real2str(var)],delim))
  end subroutine Overlapping_Batch_Means

  !=================================================================================================

end program postlammps
