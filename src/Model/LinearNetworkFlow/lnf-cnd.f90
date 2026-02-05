!> @brief Linear Network Flow Static Conductance Package
!!
!! This module contains the static conductance package for the LNF model.
!! It provides conductance values for DISU connections and contributes
!! to the matrix during assembly.
!!
!! Phase 1: Static only - no PERIOD blocks, no stress-period updates.
!<
module LnfCndModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DONE, LENMEMPATH, LENVARNAME, LINELENGTH
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error, count_errors
  use NumericalPackageModule, only: NumericalPackageType
  use BaseDisModule, only: DisBaseType
  use MemoryManagerModule, only: mem_allocate, mem_deallocate, mem_setptr
  use MemoryHelperModule, only: create_mem_path
  use MatrixBaseModule

  implicit none

  private
  public :: LnfCndType
  public :: cnd_cr

  type, extends(NumericalPackageType) :: LnfCndType

    ! -- conductance data
    real(DP), pointer :: cond_scalar => null() !< scalar conductance value
    real(DP), dimension(:), pointer, contiguous :: cond => null() !< conductance array (NJA)
    logical(LGP), pointer :: use_scalar => null() !< flag to use scalar conductance

    ! -- pointers to model data
    integer(I4B), dimension(:), pointer, contiguous :: ia => null() !< CSR row pointer
    integer(I4B), dimension(:), pointer, contiguous :: ja => null() !< CSR column indices
    integer(I4B), pointer :: nja => null() !< number of connections

  contains

    procedure :: cnd_df
    procedure :: cnd_ar
    procedure :: cnd_fc
    procedure :: cnd_cq
    procedure :: cnd_save_model_flows
    procedure :: cnd_da
    procedure :: allocate_scalars
    procedure :: allocate_arrays
    procedure, private :: source_options
    procedure, private :: source_griddata
    procedure, private :: log_options

  end type LnfCndType

contains

  !> @brief Create a new CND package
  !<
  subroutine cnd_cr(cndobj, name_model, input_mempath, inunit, iout, dis)
    ! -- dummy
    type(LnfCndType), pointer :: cndobj
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: input_mempath
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    class(DisBaseType), pointer, intent(in) :: dis
    !
    ! -- Create the object
    allocate (cndobj)
    !
    ! -- create name and memory path
    call cndobj%set_names(1, name_model, 'CND', 'CND', input_mempath)
    !
    ! -- Allocate scalars
    call cndobj%allocate_scalars()
    !
    ! -- Set variables
    cndobj%inunit = inunit
    cndobj%iout = iout
    cndobj%dis => dis
    !
  end subroutine cnd_cr

  !> @brief Define the CND package
  !<
  subroutine cnd_df(this, dis)
    ! -- dummy
    class(LnfCndType) :: this
    class(DisBaseType), pointer :: dis
    !
    ! -- Set pointers to discretization data
    this%dis => dis
    this%ia => dis%con%ia
    this%ja => dis%con%ja
    this%nja => dis%con%nja
    !
    ! -- Source options
    call this%source_options()
    !
    ! -- Allocate arrays
    call this%allocate_arrays()
    !
    ! -- Source griddata
    call this%source_griddata()
    !
  end subroutine cnd_df

  !> @brief Allocate and read
  !<
  subroutine cnd_ar(this)
    ! -- dummy
    class(LnfCndType) :: this
    !
    ! -- Nothing to do for static conductance
    !
  end subroutine cnd_ar

  !> @brief Fill coefficients - add conductance to matrix
  !!
  !! For each connection (n, m), add:
  !!   A(n,n) -= C  (subtract from diagonal)
  !!   A(n,m) += C  (add to off-diagonal)
  !!
  !<
  subroutine cnd_fc(this, kiter, matrix_sln, idxglo, rhs, hnew)
    ! -- dummy
    class(LnfCndType) :: this
    integer(I4B), intent(in) :: kiter
    class(MatrixBaseType), pointer :: matrix_sln
    integer(I4B), dimension(:), intent(in) :: idxglo
    real(DP), dimension(:), intent(inout) :: rhs
    real(DP), dimension(:), intent(inout) :: hnew
    ! -- local
    integer(I4B) :: n, m, ipos, idiag, isympos
    real(DP) :: cond
    !
    ! -- Loop over all nodes
    do n = 1, this%dis%nodes
      !
      ! -- Loop over connections for this node
      do ipos = this%ia(n) + 1, this%ia(n + 1) - 1
        m = this%ja(ipos)
        !
        ! -- Skip lower triangle (handle each connection once)
        if (m < n) cycle
        !
        ! -- Get conductance for this connection
        if (this%use_scalar) then
          cond = this%cond_scalar
        else
          cond = this%cond(ipos)
        end if
        !
        ! -- Add to row n: A(n,m) += cond, A(n,n) -= cond
        idiag = this%ia(n)
        call matrix_sln%add_value_pos(idxglo(ipos), cond)
        call matrix_sln%add_value_pos(idxglo(idiag), -cond)
        !
        ! -- Add to row m: A(m,n) += cond, A(m,m) -= cond
        isympos = this%dis%con%isym(ipos)
        call matrix_sln%add_value_pos(idxglo(isympos), cond)
        call matrix_sln%add_value_pos(idxglo(this%ia(m)), -cond)
        !
      end do
    end do
    !
  end subroutine cnd_fc

  !> @brief Calculate flows for CND package
  !!
  !! Calculate flow between connected cells and store in flowja.
  !! Flow is positive from n to m when h(n) > h(m).
  !!
  !<
  subroutine cnd_cq(this, hnew, flowja)
    ! -- dummy
    class(LnfCndType) :: this
    real(DP), dimension(:), intent(in) :: hnew
    real(DP), dimension(:), intent(inout) :: flowja
    ! -- local
    integer(I4B) :: n, m, ipos, isympos
    real(DP) :: cond, qnm
    !
    ! -- Loop over all nodes
    do n = 1, this%dis%nodes
      !
      ! -- Loop over connections for this node
      do ipos = this%ia(n) + 1, this%ia(n + 1) - 1
        m = this%ja(ipos)
        !
        ! -- Skip lower triangle (handle each connection once)
        if (m < n) cycle
        !
        ! -- Get conductance for this connection
        if (this%use_scalar) then
          cond = this%cond_scalar
        else
          cond = this%cond(ipos)
        end if
        !
        ! -- Calculate flow from n to m: qnm = cond * (hn - hm)
        qnm = cond * (hnew(n) - hnew(m))
        !
        ! -- Add to flowja: positive qnm means flow out of n into m
        flowja(ipos) = flowja(ipos) + qnm
        !
        ! -- Add symmetric flow (flow from m to n is negative of qnm)
        isympos = this%dis%con%isym(ipos)
        flowja(isympos) = flowja(isympos) - qnm
        !
      end do
    end do
    !
  end subroutine cnd_cq

  !> @brief Record flowja to binary file
  !<
  subroutine cnd_save_model_flows(this, flowja, icbcfl, icbcun)
    ! -- dummy
    class(LnfCndType) :: this
    real(DP), dimension(:), intent(in) :: flowja
    integer(I4B), intent(in) :: icbcfl
    integer(I4B), intent(in) :: icbcun
    ! -- local
    integer(I4B) :: ibinun
    !
    ! -- Set unit number for binary output
    if (this%ipakcb < 0) then
      ibinun = icbcun
    elseif (this%ipakcb == 0) then
      ibinun = 0
    else
      ibinun = this%ipakcb
    end if
    if (icbcfl == 0) ibinun = 0
    !
    ! -- Write the face flows if requested
    if (ibinun /= 0) then
      call this%dis%record_connection_array(flowja, ibinun, this%iout)
    end if
    !
  end subroutine cnd_save_model_flows

  !> @brief Deallocate
  !<
  subroutine cnd_da(this)
    ! -- dummy
    class(LnfCndType) :: this
    !
    ! -- Deallocate arrays
    call mem_deallocate(this%cond)
    !
    ! -- Deallocate scalars
    call mem_deallocate(this%cond_scalar)
    call mem_deallocate(this%use_scalar)
    !
    ! -- Deallocate parent
    call this%NumericalPackageType%da()
    !
  end subroutine cnd_da

  !> @brief Allocate scalars
  !<
  subroutine allocate_scalars(this)
    ! -- dummy
    class(LnfCndType) :: this
    !
    ! -- Allocate scalars in parent
    call this%NumericalPackageType%allocate_scalars()
    !
    ! -- Allocate scalars for this package
    call mem_allocate(this%cond_scalar, 'COND_SCALAR', this%memoryPath)
    call mem_allocate(this%use_scalar, 'USE_SCALAR', this%memoryPath)
    !
    ! -- Initialize
    this%cond_scalar = DZERO
    this%use_scalar = .true.
    !
  end subroutine allocate_scalars

  !> @brief Allocate arrays
  !<
  subroutine allocate_arrays(this)
    ! -- dummy
    class(LnfCndType) :: this
    ! -- local
    integer(I4B) :: ipos
    !
    ! -- Allocate conductance array sized by NJA
    call mem_allocate(this%cond, this%nja, 'COND', this%memoryPath)
    !
    ! -- Initialize to zero
    do ipos = 1, this%nja
      this%cond(ipos) = DZERO
    end do
    !
  end subroutine allocate_arrays

  !> @brief Source options from input context
  !<
  subroutine source_options(this)
    use MemoryManagerExtModule, only: mem_set_value
    use LnfCndInputModule, only: LnfCndParamFoundType
    ! -- dummy
    class(LnfCndType) :: this
    ! -- local
    type(LnfCndParamFoundType) :: found
    !
    ! -- Source options (currently none required)
    !
    ! -- Log options
    if (this%iout > 0) then
      call this%log_options(found)
    end if
    !
  end subroutine source_options

  !> @brief Log options
  !<
  subroutine log_options(this, found)
    use LnfCndInputModule, only: LnfCndParamFoundType
    ! -- dummy
    class(LnfCndType) :: this
    type(LnfCndParamFoundType), intent(in) :: found
    !
    write (this%iout, '(1x,a)') 'PROCESSING CND OPTIONS'
    write (this%iout, '(1x,a)') 'END OF CND OPTIONS'
    !
  end subroutine log_options

  !> @brief Source griddata from input context
  !<
  subroutine source_griddata(this)
    use MemoryManagerExtModule, only: mem_set_value
    ! -- dummy
    class(LnfCndType) :: this
    ! -- local
    real(DP), dimension(:), pointer, contiguous :: cond_input => null()
    integer(I4B) :: ipos
    !
    ! -- Try to get CONDUCTANCE from input context
    call mem_setptr(cond_input, 'CONDUCTANCE', this%input_mempath)
    !
    if (associated(cond_input)) then
      ! -- Check if it's a single value (scalar) or array
      if (size(cond_input) == 1) then
        this%use_scalar = .true.
        this%cond_scalar = cond_input(1)
        write (this%iout, '(1x,a,g15.7)') &
          'USING SCALAR CONDUCTANCE = ', this%cond_scalar
      else
        ! -- Array form
        this%use_scalar = .false.
        if (size(cond_input) /= this%nja) then
          write (errmsg, '(a,i0,a,i0)') &
            'CONDUCTANCE array size (', size(cond_input), &
            ') does not match NJA (', this%nja, ')'
          call store_error(errmsg)
        else
          do ipos = 1, this%nja
            this%cond(ipos) = cond_input(ipos)
          end do
          write (this%iout, '(1x,a)') 'CONDUCTANCE READ FROM ARRAY'
        end if
      end if
    else
      write (errmsg, '(a)') 'CONDUCTANCE not found in input'
      call store_error(errmsg)
    end if
    !
    ! -- Check for errors
    if (count_errors() > 0) then
      call store_error('Error in CND GRIDDATA', terminate=.true.)
    end if
    !
  end subroutine source_griddata

end module LnfCndModule
