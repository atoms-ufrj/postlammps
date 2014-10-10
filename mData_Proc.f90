!> This module defines classes for managing simulation data processing tasks.
!!
!! @author Charlles R. A. Abreu (abreu@eq.ufrj.br)
!! @date Sept 17, 2013
!!
!! @todo Use polymorphic arrays (feature not implemented in gfortran 4.6) to manage several data
!!       processors simultaneously.
module mData_Proc
use mConstants
use mProp_List
implicit none

type, private :: Container
  real(rb), pointer :: Data(:)
end type Container

!===================================================================================================
!                                            Data_Proc
!===================================================================================================

!> An abstract class for general simulation data processors.
type, abstract :: Data_Proc

  !> True if data processing has not been executed for the first time.
  logical, private :: not_executed = .true.

  !> The list of properties associated with the data processor.
  type(Prop_List) :: Props

  !> Interval of action of the data processor.
  integer, private :: Interval = 1

  contains

    !> Executes the defined data processing task with a provided set of system property values.
    !! @param[in] data (real vector) set of system properties to be handled by the data processor.
    !! @remark The size of vector "data" must match the sum of dimensions of all properties
    !!         associated with the data processor.
    procedure :: Exec => Data_Proc_Exec

    procedure(Data_Proc_Init), deferred, private :: Init
    procedure(Data_Proc_Act),  deferred, private :: Act

end type Data_Proc

private :: Data_Proc_Init, Data_Proc_Act, Data_Proc_Exec

!> Abstract interface for deferred procedures:
abstract interface

  subroutine Data_Proc_Init( a )
    import :: Data_Proc
    class(Data_Proc), intent(inout) :: a
  end subroutine Data_Proc_Init

  subroutine Data_Proc_Act( a, Step, Data )
    import :: Data_Proc, rb
    class(Data_Proc), intent(inout) :: a
    integer,          intent(in)    :: Step
    real(rb),         intent(in)    :: Data(a%Props%Total)
  end subroutine Data_Proc_Act

end interface

!===================================================================================================
!                                          Data_Output
!===================================================================================================
!> A class for managing output of data to external files or other output units.
type, extends(Data_Proc) :: Data_Output

  integer, private :: unit = 6

  contains

    !> Deferred binding procedure for performing initial step.
    procedure, private :: Init => Data_Output_Init

    !> Deferred binding procedure for the data output processor action.
    procedure, private :: Act => Data_Output_Act

    !> Performs initialization of the data output processor.
    !! @param[in] unit (integer) an output unit.
    !! @param[in] interval (integer, optional) interval between actual executions of data output
    !!            calls. Default = 1.
    !! @param[in] separator (character, optional) a character to separate output values. Default =
    !!            Props%separator.
    procedure :: Setup => Data_Output_Setup

end type Data_Output

private :: Data_Output_Init, Data_Output_Act, Data_Output_Setup

!===================================================================================================
!                                         Block_Analyzer
!===================================================================================================
!> A class for performing data analysis based on the Block Average algorithm of Flyvbjerg and
!! Petersen (1989). The implementation is similar to the on-the-fly version described by
!! Kent et al. (2007).
!!
!! @cite Flyvbjerg_1989a Flyvbjerg and Petersen, J. Chem. Phys. 91, 461, 1989. <br>
!! @cite Kent_et_al_2007a Kent et al., J. Comput. Chem. 8, 2309–2316, 2007.
type, extends(Data_Proc) :: Block_Analyzer

  integer, private :: entries                   !< Number of sampled values of each property
  integer, private :: levels                    !< Number of levels of blocking operations
  integer, private :: B0 = 2                    !< Block size at level 0
  integer, private :: BS = 2                    !< Block size at other levels
  type(container), allocatable, private :: Acc(:,:)  !< Accumulators at each level

  contains

    !> Deferred binding procedure for performing initial step.
    procedure, private :: Init => Block_Analyzer_Init

    !> Deferred binding procedure for the block analyzer action.
    procedure, private :: Act => Block_Analyzer_Act

    !> Performs initialization of the block analyzer.
    !! @param[in] interval (integer) interval between actual executions of the block analyzer
    !!            calls (Default = 1).
    !! @param[in] first (integer) the size of blocks at level zero of blocking operations
    !!            (Default = 2).
    !! @param[in] others (integer) size of blocks at other levels of blocking operations
    !!            (Default = 2).
    procedure :: Setup => Block_Analyzer_Setup

    !> Performs block analysis computations and saves the results in a provided output unit.
    !! @param[in] unit (integer) an output unit for saving the block analysis results.
    !! @param[in] separator (character, optional) a character to separate output values.
    !!            Default = Props%separator.
    procedure :: Flush => Block_Analyzer_Flush

    !> Adds a new level to accummulators.
    procedure, private :: Add_Level => Block_Analyzer_Add_Level

end type Block_Analyzer

private :: Block_Analyzer_Init, Block_Analyzer_Act, Block_Analyzer_Setup, &
           Block_Analyzer_Add_Level, Block_Analyzer_Flush

!===================================================================================================
!                                            Correlator
!===================================================================================================

!> A class for performing mean square displacement (MSD) and autocorrelation function (ACF)
!! calculations based on the Multiple Window Algorithm of Dubbeldam et al. (2009).
!!
!! @cite Dubbeldam_2009a Dubbeldam et al., Molecular Simulation 35(1), 1084-1097, 2009.
type, abstract, extends(Data_Proc) :: Correlator

  integer, private :: entries                   !< Number of sampled values of each property
  integer, private :: levels                    !< Number of levels of MSD computations
  integer, private :: M = 5                     !< Number of MSD computations per level
  integer, private :: Mm1 = 4                   !< M - 1 (saved to avoid repeated computations)
  character(3), private :: Type                 !< Type of operation carried out by the correlator

  type(Container), allocatable, private :: Record(:,:)   !< Saved data at each level
  type(Container), allocatable, private :: Acc(:,:)      !< Accumulators at all levels/displacements
  type(Container), private :: Acc0                       !< Accumulator at displacement=0

  contains

    !> Deferred binding procedure for performing initial step.
    procedure, private :: Init => Correlator_Init

    !> Deferred binding procedure for the data processor action.
    procedure, private :: Act => Correlator_Act

    !> Performs initialization of the mean square displacement computer.
    !! @param[in] interval (integer) interval between actual executions of the MSD computer
    !!            calls (Default = 1).
    !! @param[in] M (integer) number of entries per level of correlation computations.
    procedure :: Setup => Correlator_Setup

    !> Performs MSD or ACF computations and saves the results in a provided output unit.
    !! @param[in] unit (integer) an output unit for saving the MSD or ACF results.
    !! @param[in] separator (character, optional) a character to separate output values.
    !!            Default = Props%separator.
    procedure :: Flush => Correlator_Flush

    !> Adds a new level to the accummulators:
    procedure, private :: Add_Level => Correlator_Add_Level

    !> Performs the required operation for MSD or ACF calculation:
    procedure(Correlator_Operation), nopass, deferred, private :: Operation

    !> Defines the type of operator performed by the data correlator:
    procedure(Correlation_Define_Type), deferred, private :: Define_Type

end type Correlator

private :: Correlator_Init, Correlator_Act, Correlator_Setup, Correlator_Add_Level, &
           Correlator_Operation, Correlation_Define_Type

abstract interface

  function Correlator_Operation( x, y ) result( z )
    import :: rb
    real(rb), intent(in) :: x(:), y(:)
    real(rb)             :: z(size(x))
  end function Correlator_Operation

  subroutine Correlation_Define_Type( a )
    import :: Correlator
    class(Correlator), intent(inout) :: a
  end subroutine Correlation_Define_Type

end interface

!===================================================================================================
!                                            MSD_Comp
!===================================================================================================

type, extends(Correlator) :: MSD_Comp
  contains
    procedure, nopass, private :: Operation => MSD_Comp_Operation
    procedure, private :: Define_Type => MSD_Comp_Define_Type
end type MSD_Comp

private :: MSD_Comp_Operation, MSD_Comp_Define_Type

!===================================================================================================
!                                            ACF_Comp
!===================================================================================================

type, extends(Correlator) :: ACF_Comp
  contains
    procedure, nopass, private :: Operation => ACF_Comp_Operation
    procedure, private :: Define_Type => ACF_Comp_Define_Type
end type ACF_Comp

private :: ACF_Comp_Operation, ACF_Comp_Define_Type

!===================================================================================================

contains

!===================================================================================================
!                                           Data_Proc
!===================================================================================================
  subroutine Data_Proc_Exec( a, Step, Data )
    class(Data_Proc), intent(inout) :: a
    integer,          intent(in)    :: Step
    real(rb),         intent(in)    :: Data(a%Props%Total)
    if (mod(Step,a%Interval) == 0) then
      if (a%not_executed) then
        if (a%Props%Number > 0) then
          call a%Init()
        else
          stop "Error: trying to execute a data processor with no associated properties."
        end if
        a%not_executed = .false.
      end if
      call a%Act( Step, Data )
    end if
  end subroutine Data_Proc_Exec
!===================================================================================================
!                                           Data_Output
!===================================================================================================
  subroutine Data_Output_Init( a )
    class(Data_Output), intent(inout) :: a
    write(a%unit,'(A)') "Step"//a%props%separator//trim( a%Props%Titles() )
  end subroutine Data_Output_Init
  !-------------------------------------------------------------------------------------------------
  subroutine Data_Output_Act( a, Step, Data )
    class(Data_Output), intent(inout) :: a
    integer,          intent(in)    :: Step
    real(rb),         intent(in)    :: Data(a%Props%Total)
    integer :: i
    character(sl) :: C
    write(C,*) Data(i)
    do i = 2, a%Props%Total
      write(a%unit,'(A,A)',advance="no") trim(C), a%Props%separator
      write(C,*) Data(i)
    end do
    write(a%unit,'(A)') trim(C)
  end subroutine Data_Output_Act
  !-------------------------------------------------------------------------------------------------
  subroutine Data_Output_Setup( a, unit, interval, separator )
    class(Data_Output), intent(inout)        :: a
    integer,            intent(in), optional :: unit, interval
    character,          intent(in), optional :: separator
    if (present(unit)) a%unit = unit
    if (present(interval)) a%interval = interval
    if (present(separator)) a%props%separator = separator
  end subroutine Data_Output_Setup
!===================================================================================================
!                                        Block_Analyzer
!===================================================================================================
  subroutine Block_Analyzer_Init( a )
    class(Block_Analyzer), intent(inout) :: a
    integer :: j
    a%entries = 0
    a%levels = 0
    allocate( a%Acc(0:2,0:a%levels) )
    do j = 0, 2
      allocate( a%Acc(j,0)%Data(a%Props%Total) )
      a%Acc(j,0)%Data = zero
    end do
  end subroutine Block_Analyzer_Init
  !-------------------------------------------------------------------------------------------------
  subroutine Block_Analyzer_Act( a, Step, Data )
    class(Block_Analyzer), intent(inout) :: a
    integer,               intent(in)    :: Step
    real(rb),              intent(in)    :: Data(a%Props%Total)
    integer  :: level, milestone
    real(rb) :: Average(size(Data))
    a%entries = a%entries + 1                          ! Update sample size
    call Accumulate( Data, 0 )                         ! Accumulate data at level 0
    if (mod(a%entries,a%B0) == 0) then                 ! Check necessity to act at level 1
      Average = a%Acc(0,0)%Data/a%B0                   ! Average data of level 0
      call Accumulate( Average, 1 )                    ! Accumulate averages at level 1
      a%Acc(0,0)%Data = zero                           ! Reset volatile accumulator of level 0
    end if
    level = 2                                          ! Set level = 2
    milestone = a%B0*a%BS                              ! Set milestone for level 2
    do while (mod(a%entries,milestone)  == 0)          ! Check necessity to act at current level
      Average = a%Acc(0,level-1)%Data/a%BS             ! Average data of level below current level
      call Accumulate( Average, level )                ! Accumulate averages at current level
      a%Acc(0,level-1)%Data = zero                     ! Reset volatile accumulator of level below
      level = level + 1                                ! Raise to next level
      milestone = milestone*a%BS                       ! Set milestone for next level
    end do
    contains
      subroutine Accumulate( Data, level )
        real(rb), intent(in) :: Data(:)
        integer,  intent(in) :: level
        if (level > a%levels) call a%Add_Level
        a%Acc(0,level)%Data = a%Acc(0,level)%Data + Data
        a%Acc(1,level)%Data = a%Acc(1,level)%Data + Data
        a%Acc(2,level)%Data = a%Acc(2,level)%Data + Data*Data
      end subroutine Accumulate
  end subroutine Block_Analyzer_Act
  !-------------------------------------------------------------------------------------------------
  subroutine Block_Analyzer_Add_Level( a )
    class(Block_Analyzer), intent(inout) :: a
    integer :: NL, j
    type(Container), allocatable :: caux(:,:)
    NL = a%levels + 1
    allocate( caux(0:2,0:NL) )
    caux(:,0:a%levels) = a%Acc
    call move_alloc( caux, a%Acc )
    do j = 0, 2
      allocate( a%Acc(j,NL)%Data(a%Props%Total) )
      a%Acc(j,NL)%Data = zero
    end do
    a%levels = NL
  end subroutine Block_Analyzer_Add_Level
  !-------------------------------------------------------------------------------------------------
  subroutine Block_Analyzer_Setup( a, interval, first, others )
    class(Block_Analyzer), intent(inout)          :: a
    integer,               intent(in),   optional :: interval, first, others
    if (a%not_executed) then
      if (present(interval)) a%interval = interval
      if (present(first)) a%B0 = first
      if(present(others)) a%BS = others
    else
      stop "Error: trying to setup a block analyzer which has already been executed."
    end if
  end subroutine Block_Analyzer_Setup
  !-------------------------------------------------------------------------------------------------
  subroutine Block_Analyzer_Flush( a, unit, separator, titles )
    class(Block_Analyzer), intent(inout)        :: a
    integer,               intent(in), optional :: unit
    character,             intent(in), optional :: separator
    logical,               intent(in), optional :: titles

    integer   :: level, i, out, nBlocks(0:a%levels), Block_Size(0:a%levels)
    real(rb)  :: Mean(a%Props%Total,0:a%levels), Var(a%Props%Total,0:a%levels)
    character :: sep
    logical   :: print_titles
    character(1000) :: C
    character(sl)   :: D

    print_titles = .true.; if (present(titles)) print_titles = titles
    out = 6; if (present(unit)) out = unit
    sep = a%Props%separator; if (present(separator)) sep = separator

    if (print_titles) write(out,'(A)') "size"//sep//"number"//sep//"log(size)"// &
      trim( a%Props%Titles( sep//"mean(",  ")"//sep//"mean(",  ")" ) ) // &
      trim( a%Props%Titles( sep//"stdev(", ")"//sep//"stdev(", ")" ) ) // &
      trim( a%Props%Titles( sep//"error(", ")"//sep//"error(", ")" ) ) // &
      trim( a%Props%Titles( sep//"ineff(", ")"//sep//"ineff(", ")" ) )

    Block_Size(0) = 1
    forall (level=1:a%levels) Block_Size(level) = a%B0*a%BS**(level-1)
    do level = 0, a%levels
      nBlocks(level) = a%entries/Block_Size(level)
      if (nBlocks(level) > 1) then
        Mean(:,level) = a%Acc(1,level)%Data/nBlocks(level)
        Var(:,level) = a%Acc(2,level)%Data/nBlocks(level) - Mean(:,level)**2
      end if
    end do

    do level = 0, a%levels
      if (nBlocks(level) > 1) then
        write(C,*) Block_Size(level)
        write(D,*) nBlocks(level)
        C = trim(adjustl(C))//sep//trim(adjustl(D))
        C = join( C, log10(real(Block_Size(level),rb)) )
        do i = 1, size(Mean,1)
          C = join( C, Mean(i,level) )
        end do
        do i = 1, size(Var,1)
          C = join( C, sqrt(Var(i,level)/(nBlocks(level) - one)) )
        end do
        do i = 1, size(Var,1)
          C = join( C, sqrt(half*Var(i,level))/(nBlocks(level) - one)**3)
        end do
        do i = 1, size(Var,1)
          C = join( C, Block_Size(level)*Var(i,level)/Var(i,0) )
        end do
        if (print_titles.or.(level > 0)) write(unit,'(A)') trim(C)
      end if
    end do

    contains
      function join( a, b ) result( c )
        character(*), intent(in) :: a
        real(rb),     intent(in) :: b
        character(1000)          :: c
        write(c,*) b
        c = trim(adjustl(a))//sep//trim(adjustl(c))
      end function join
  end subroutine Block_Analyzer_Flush
!===================================================================================================
!                                            Correlator
!===================================================================================================
  subroutine Correlator_Init( a )
    class(Correlator), intent(inout) :: a
    integer :: i
    a%entries = 0
    a%levels = 0
    allocate( a%Record(0:a%Mm1,0:0), a%Acc(1:a%Mm1,0:0)  )
    allocate( a%Record(0,0)%Data(a%Props%Total) )
    do i = 1, a%Mm1
      allocate( a%Record(i,0)%Data(a%Props%Total) )
      allocate( a%Acc(i,0)%Data(a%Props%Number) )
      a%Acc(i,0)%Data = zero
    end do
    allocate( a%Acc0%Data(a%Props%Number) )
    a%Acc0%Data = zero
    call a%Define_Type
  end subroutine Correlator_Init
  !-------------------------------------------------------------------------------------------------
  subroutine Correlator_Act( a, Step, Data )
    class(Correlator), intent(inout) :: a
    integer,         intent(in)    :: Step
    real(rb),        intent(in)    :: Data(a%Props%Total)
    integer  :: pos, bin, level, unity, delta, origin
    real(rb) :: SD(a%Props%Total)
    a%entries = a%entries + 1
    bin = 0
    level = 0
    unity = 1
    SD = a%Operation( Data, Data )
    a%Acc0%Data = a%Acc0%Data + a%Props%Sum( SD )
    do while (bin == 0)
      if (level > a%levels) call a%Add_Level
      pos = a%entries/unity
      do delta = 1, min(a%Mm1,pos-1)
        origin = mod(pos-delta,a%M)
        SD = a%Operation( Data, a%Record(origin,level)%Data )
        a%Acc(delta,level)%Data = a%Acc(delta,level)%Data + a%Props%Sum( SD )
      end do
      bin = mod(pos,a%M)
      a%Record(bin,level)%Data = Data
      level = level + 1
      unity = unity*a%M
    end do
  end subroutine Correlator_Act
  !-------------------------------------------------------------------------------------------------
  subroutine Correlator_Add_Level( a )
    class(Correlator), intent(inout) :: a
    integer :: j, NL
    type(container), allocatable :: aux(:,:)
    NL = a%levels+1
    allocate( aux(0:a%Mm1,0:NL) )
    aux(:,0:a%levels) = a%Record
    call move_alloc( aux, a%Record )
    allocate( aux(1:a%Mm1,0:NL) )
    aux(:,0:a%levels) = a%Acc
    call move_alloc( aux, a%Acc )
    allocate( a%Record(0,NL)%Data(a%Props%Total) )
    do j = 1, a%Mm1
      allocate( a%Record(j,NL)%Data(a%Props%Total) )
      allocate( a%Acc(j,NL)%Data(a%Props%Number) )
      a%Acc(j,NL)%Data = zero
    end do
    a%levels = NL
  end subroutine Correlator_Add_Level
  !-------------------------------------------------------------------------------------------------
  subroutine Correlator_Setup( a, interval, M )
    class(Correlator), intent(inout)        :: a
    integer,           intent(in), optional :: interval, M
    if (a%not_executed) then
      if (present(interval)) a%interval = interval
      if (present(M)) then
        a%M = M
        a%Mm1 = M - 1
      end if
    else
      stop "Error: trying to setup a MSD computer which has already been executed."
    end if
  end subroutine Correlator_Setup
  !-------------------------------------------------------------------------------------------------
  subroutine Correlator_Flush( a, unit, separator )
    class(Correlator), intent(inout)        :: a
    integer,           intent(in), optional :: unit
    character,         intent(in), optional :: separator
    integer :: out, i, level, prop, ndata
    character :: sep
    character(1000) :: C
    if (present(unit)) then
      out = unit
    else
      out = 6
    end if
    if (present(separator)) then
      sep = separator
    else
      sep = a%Props%separator
    end if
    write(out,'(A)',advance='no') "steps"//sep//a%Type//"["//trim(a%Props%Name(1))//"]"
    do i = 2, a%Props%Number
      write(out,'(A)',advance='no') sep//a%Type//"["//trim(a%Props%Name(i))//"]"
    end do
    write(out,*)

    write(C,*) 0
    write(out,'(A)',advance="no") trim(adjustl(C))
    do prop = 1, a%Props%Number
      write(C,*) a%Acc0%Data(prop)/(a%entries * a%Props%Dim(prop))
      write(out,'(A)',advance="no") sep//trim(adjustl(C))
    end do
    write(out,'()')

    do level = 0, a%levels
      ndata = a%entries/a%M**level
      do i = 1, a%Mm1
        ndata = ndata - 1
        if (ndata > 0) then
          write(C,*) i*a%M**level
          write(out,'(A)',advance="no") trim(adjustl(C))
          do prop = 1, a%Props%Number
            write(C,*) a%Acc(i,level)%Data(prop)/(ndata * a%Props%Dim(prop))
            write(out,'(A)',advance="no") sep//trim(adjustl(C))
          end do
          write(out,'()')
        end if
      end do
    end do

  end subroutine Correlator_Flush
!===================================================================================================
!                                       MSD_Comp
!===================================================================================================
  function MSD_Comp_Operation( x, y ) result( z )
    real(rb), intent(in) :: x(:), y(:)
    real(rb)             :: z(size(x))
    z = (x - y)**2
  end function MSD_Comp_Operation
  !-------------------------------------------------------------------------------------------------
  subroutine MSD_Comp_Define_Type( a )
    class(MSD_Comp), intent(inout) :: a
    a%Type = "MSD"
  end subroutine MSD_Comp_Define_Type
!===================================================================================================
!                                       ACF_Comp
!===================================================================================================
  function ACF_Comp_Operation( x, y ) result( z )
    real(rb), intent(in) :: x(:), y(:)
    real(rb)             :: z(size(x))
    z = x*y
  end function ACF_Comp_Operation
  !-------------------------------------------------------------------------------------------------
  subroutine ACF_Comp_Define_Type( a )
    class(ACF_Comp), intent(inout) :: a
    a%Type = "ACF"
  end subroutine ACF_Comp_Define_Type
!===================================================================================================
end module mData_Proc
