!> This module define a classes for storage and management of property names and dimensions.
!!
!! @author Charlles R. A. Abreu (abreu@eq.ufrj.br)
!! @date Sept 21, 2013
module mProp_List
use mConstants
implicit none

!> A structured type for storing and managing property names:
type Prop_List

  !> The number of properties in the list.
  integer :: Number = 0

  !> The names of the properties in the list.
  character(sl), allocatable :: Name(:)

  !> The dimension of each property in the list.
  integer, allocatable :: Dim(:)

  !> The sum of dimensions of the properties in the list.
  integer :: Total = 0

  !> The initial position of each property in a contiguous array of data.
  integer, allocatable :: First(:)

  !> The final position of each property in a contiguous array of data.
  integer, allocatable :: Last(:)

  !> A character to be used as separator of property names. Default = " " (space).
  character :: separator = " "

  contains

    procedure, private :: Prop_List_Add_Scalar
    procedure, private :: Prop_List_Add_Array
    !> Adds a new property or an array of new properties to the list.
    !! @param[in] Name (string scalar/array) name of each new property to be added.
    !! @param[in] Dim (integer scalar/array, optional) dimension of each new property to be added
    !!            (Default = 1 for each property).
    generic :: Add => Prop_List_Add_Scalar, Prop_List_Add_Array

    !> Returns a string with the names of the properties in the list.
    !! @param[in] begin (string, optional) a heading string. Default = "" (none).
    !! @param[in] middle (string, optional) a separation string. Default = Prop_List::separator.
    !! @param[in] end (string, optional) a trailing string. Default = "" (none).
    procedure :: Titles => Prop_List_Titles

    !! Sums the values of all dimensions of each property.
    !! @param[in] data (real array) an array with dimension Prop_List::Total.
    !! @returns an array with dimension Prop_List::Number with the sums over dimensions of each
    !!          property.
    procedure :: Sum => Prop_List_Sum

end type Prop_List

private :: Prop_List_Add_Scalar, Prop_List_Add_Array, Prop_List_Titles, Prop_List_Sum

contains
  !-------------------------------------------------------------------------------------------------
  subroutine Prop_List_Add_Scalar( a, Name, Dim )
    class(Prop_List), intent(inout)       :: a
    character(*),    intent(in)           :: Name
    integer,         intent(in), optional :: Dim
    if (present(Dim)) then
      call Prop_List_Add_Array( a, [Name], [Dim] )
    else
      call Prop_List_Add_Array( a, [Name] )
    end if
  end subroutine Prop_List_Add_Scalar
  !-------------------------------------------------------------------------------------------------
  subroutine Prop_List_Add_Array( a, Name, Dim )
    class(Prop_List), intent(inout)        :: a
    character(*),     intent(in)           :: Name(:)
    integer,          intent(in), optional :: Dim(:)
    integer :: i
    character(sl) :: caux(a%Number + size(Name))
    integer       :: iaux(a%Number + size(Name))
    caux(1:a%Number) = a%Name
    iaux(1:a%Number) = a%Dim
    forall(i=1:size(Name)) caux(a%Number+i) = trim(adjustl(Name(i)))
    if (present(Dim)) then
      iaux(a%Number+1:size(iaux)) = Dim
    else
      iaux(a%Number+1:size(iaux)) = 1
    end if
    a%Number = a%Number + size(Name)
    if (allocated(a%Name)) deallocate(a%Name,a%Dim,a%First,a%Last)
    allocate( a%Name(a%Number), a%Dim(a%Number), a%First(a%Number), a%Last(a%Number) )
    a%Name = caux
    a%Dim = iaux
    forall (i=1:a%Number)
      a%First(i) = sum(a%Dim(1:i-1)) + 1
      a%Last(i)  = sum(a%Dim(1:i))
    end forall
    a%Total = a%Last(a%Number)
  end subroutine Prop_List_Add_Array
  !-------------------------------------------------------------------------------------------------
  function Prop_List_Titles( Props, begin, middle, end ) result( Titles )
    class(Prop_List), intent(inout)        :: Props
    character(*),     intent(in), optional :: begin, middle, end
    character(sl)                          :: Titles
    character(10) :: C
    integer :: i, j, k, n
    if (present(begin)) then
      Titles = begin
      n = len(begin)
    else
      Titles = ""
      n = 0
    end if
    k = 0
    do i = 1, props%Number
      do j = 1, props%Dim(i)
        Titles = Titles(1:n)//trim(props%Name(i))
        n = n + len_trim(props%Name(i))
        if (props%Dim(i) > 1) then
          write(C,'(I10)') j
          C = adjustl(C)
          Titles = Titles(1:n)//"["//trim(C)//"]"
          n = n + len_trim(C) + 2
        end if
        k = k + 1
        if (k < props%Total) then
          if (present(middle)) then
            Titles = Titles(1:n)//middle
            n = n + len(middle)
          else
            Titles = Titles(1:n)//props%separator
            n = n + 1
          end if
        end if
      end do
    end do
    if (present(end)) Titles = Titles(1:n)//end
  end function Prop_List_Titles
  !-------------------------------------------------------------------------------------------------
  function Prop_List_Sum( Props, Data ) result( S )
    class(Prop_List), intent(in) :: Props
    real(rb),         intent(in) :: Data(Props%Total)
    real(rb)                     :: S(Props%Number)
    integer :: i
    forall (i=1:Props%Number) S(i) = sum(Data(Props%First(i):Props%Last(i)))
  end function Prop_List_Sum
  !-------------------------------------------------------------------------------------------------
end module mProp_List
