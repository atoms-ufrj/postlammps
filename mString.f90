module mString

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

use mConstants

implicit none

character(3), parameter, private :: delimiters = " ,;"//achar(9)
character,    parameter, private :: comment_mark = "#"

contains

  !=================================================================================================

  elemental function has_macros( str ) result( yes )
    character(sl), intent(in) :: str
    logical                   :: yes
    yes = scan(trim(str),"*?") > 0
  end function has_macros

  !=================================================================================================

  elemental function match_str( a, b ) result( match )
    character(sl), intent(in) :: a, b
    logical                   :: match
    integer :: n, i
    i = 1
    match = .true.
    n = min(len_trim(a),len_trim(b))
    do while (match.and.(i <= n))
      if ((a(i:i) == "*").or.(b(i:i) == "*")) return
      match = (a(i:i) == b(i:i)).or.(a(i:i) == "?").or.(b(i:i) == "?")
      i = i + 1
    end do
    if (match) match = len_trim(a) == len_trim(b)
  end function match_str

  !=================================================================================================

  subroutine clean( str )
    character(*), intent(inout) :: str
    integer :: m, n
    m = scan(trim(str),comment_mark)
    if (m > 0) str = str(1:m-1)
    m = verify(str,delimiters)
    if (m == 0) then
      str = ""
    else
      n = verify(str,delimiters,back=.true.)
      str = str(m:n)
    end if
  end subroutine clean

  !=================================================================================================

  subroutine split( str, narg, arg )
    character(*), intent(in)  :: str
    integer,      intent(out) :: narg
    character(*), intent(out) :: arg(:)
    logical :: letter, word
    integer :: i, wlen
    narg = 0
    wlen = 0
    word = .false.
    do i = 1, len_trim(str)
      letter = scan(str(i:i),delimiters) == 0
      if (word) then
        if (letter) then
          wlen = wlen + 1
          arg(narg)(wlen:wlen) = str(i:i)
        else
          arg(narg) = arg(narg)(1:wlen)
          if (narg == size(arg)) return
          word = .false.
        end if
      else
        if (letter) then
          narg = narg + 1
          wlen = 1
          arg(narg)(wlen:wlen) = str(i:i)
          word = .true.
        end if
      end if
    end do
    if (word) arg(narg) = arg(narg)(1:wlen)
  end subroutine split

  !=================================================================================================

  function join( narg, arg ) result( str )
    integer,      intent(in)  :: narg
    character(*), intent(in)  :: arg(:)
    character(sl)             :: str
    integer :: i
    if (narg == 0) then
      str = ""
    else
      str = arg(1)
      do i = 2, narg
        str = trim(str)//" "//arg(i)
      end do
    end if
  end function join

  !=================================================================================================

  subroutine next_command( unit, narg, arg )
    integer,      intent(in)  :: unit
    integer,      intent(out) :: narg
    character(*), intent(out) :: arg(:)
    integer       :: ioerr
    character(sl) :: line
    read(unit,'(A'//csl//')',iostat=ioerr) line
    call clean( line )
    do while ((ioerr == 0).and.(line == ""))
      read(unit,'(A'//csl//')',iostat=ioerr) line
      call clean( line )
    end do
    if (ioerr == 0) then
      call split( line, narg, arg )
    else
      narg = 0
    end if
  end subroutine next_command

  !=================================================================================================

  function str2int( str ) result( i )
    character(*), intent(in) :: str
    integer                  :: i
    integer :: ioerr
    read(str,*,iostat=ioerr) i
    if (ioerr /= 0) call error( str, "is not an integer number" )
  end function

  !=================================================================================================

  function str2real( str ) result( r )
    character(*), intent(in) :: str
    real(rb)                 :: r
    integer :: ioerr
    read(str,*,iostat=ioerr) r
    if (ioerr /= 0) call error( str, "is not a real number" )
  end function

  !=================================================================================================

  function int2str( i ) result( str )
    integer, intent(in) :: i
    character(sl)       :: str
    integer :: ioerr
    write(str,*,iostat=ioerr) i
    str = adjustl(str)
  end function int2str

  !=================================================================================================

  function real2str( r ) result( str )
    real(rb), intent(in) :: r
    character(sl)        :: str
    integer :: ioerr
    write(str,*,iostat=ioerr) r
    str = adjustl(str)
  end function real2str

  !=================================================================================================

  subroutine warning( msg, msg1, msg2, msg3, msg4, msg5, msg6, msg7, msg8, msg9 )
    character(*), intent(in)           :: msg
    character(*), intent(in), optional :: msg1, msg2, msg3, msg4, msg5, msg6, msg7, msg8, msg9
    write(*,'("WARNING: ",A)',advance='no') trim(msg)
    if (present(msg1)) write(*,'(" ",A)',advance='no') trim(msg1)
    if (present(msg2)) write(*,'(" ",A)',advance='no') trim(msg2)
    if (present(msg3)) write(*,'(" ",A)',advance='no') trim(msg3)
    if (present(msg4)) write(*,'(" ",A)',advance='no') trim(msg4)
    if (present(msg5)) write(*,'(" ",A)',advance='no') trim(msg5)
    if (present(msg6)) write(*,'(" ",A)',advance='no') trim(msg6)
    if (present(msg7)) write(*,'(" ",A)',advance='no') trim(msg7)
    if (present(msg8)) write(*,'(" ",A)',advance='no') trim(msg8)
    if (present(msg9)) write(*,'(" ",A)',advance='no') trim(msg9)
    write(*,'()')
  end subroutine warning

  !=================================================================================================

  subroutine error( msg, msg1, msg2, msg3, msg4, msg5, msg6, msg7, msg8, msg9 )
    character(*), intent(in)           :: msg
    character(*), intent(in), optional :: msg1, msg2, msg3, msg4, msg5, msg6, msg7, msg8, msg9
    write(*,'("ERROR: ",A)',advance='no') trim(msg)
    if (present(msg1)) write(*,'(" ",A)',advance='no') trim(msg1)
    if (present(msg2)) write(*,'(" ",A)',advance='no') trim(msg2)
    if (present(msg3)) write(*,'(" ",A)',advance='no') trim(msg3)
    if (present(msg4)) write(*,'(" ",A)',advance='no') trim(msg4)
    if (present(msg5)) write(*,'(" ",A)',advance='no') trim(msg5)
    if (present(msg6)) write(*,'(" ",A)',advance='no') trim(msg6)
    if (present(msg7)) write(*,'(" ",A)',advance='no') trim(msg7)
    if (present(msg8)) write(*,'(" ",A)',advance='no') trim(msg8)
    if (present(msg9)) write(*,'(" ",A)',advance='no') trim(msg9)
    write(*,'()')
    stop
  end subroutine error

  !=================================================================================================

end module mString
