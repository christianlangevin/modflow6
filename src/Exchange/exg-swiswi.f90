!> @brief This module contains the SwiExchangeModule Module
!!
!! This module contains the code for connecting freshwater and
!! saltwater GWF models.
!!
!<
module SwiSwiExchangeModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, LINELENGTH
  use SimModule, only: count_errors, store_error, store_error_filename, &
                       store_error_unit
  use SimVariablesModule, only: errmsg, model_loc_idx
  use MemoryHelperModule, only: create_mem_path
  use BaseModelModule, only: BaseModelType, GetBaseModelFromList
  use BaseExchangeModule, only: BaseExchangeType, AddBaseExchangeToList
  use ListsModule, only: basemodellist, baseexchangelist
  use NumericalExchangeModule, only: NumericalExchangeType
  use GwfModule, only: GwfModelType
  use MatrixBaseModule

  implicit none

  private
  public :: swiswi_cr
  public :: SwiSwiExchangeType

  !> @brief Derived type for SwiSwiExchangeType
  !!
  !! This derived type contains information and methods for
  !! connecting freshwater and saltwater GWF models.
  !<
  type, extends(NumericalExchangeType) :: SwiSwiExchangeType

    ! model pointers
    class(GwfModelType), pointer :: gwf_fresh => null() !< pointer to GWF Model 1
    class(GwfModelType), pointer :: gwf_salt => null() !< pointer to GWF Model 2

    ! options
    character(len=LINELENGTH), pointer :: filename => null() !< name of the input file
    integer(I4B), pointer :: ipr_input => null() !< flag to print input
    integer(I4B), pointer :: ipr_flow => null() !< print flag for cell by cell flows

    ! number of connections
    integer(I4B), pointer :: nexg => null() !< number of connections (number of cells in fresh model or salt model)

    ! matrix position index arrays
    integer(I4B), dimension(:), pointer, contiguous :: idxglo => null() !< mapping to global (solution) amat
    integer(I4B), dimension(:), pointer, contiguous :: idxsymglo => null() !< mapping to global (solution) symmetric amat

  contains

    procedure :: exg_df => swi_swi_df
    procedure :: exg_ac => swi_swi_ac
    procedure :: exg_mc => swi_swi_mc
    procedure :: exg_ar => swi_swi_ar
    procedure :: exg_rp => swi_swi_rp
    procedure :: exg_ad => swi_swi_ad
    procedure :: exg_cf => swi_swi_cf
    procedure :: exg_fc => swi_swi_fc
    ! procedure :: exg_fn => swi_swi_fn
    ! ! procedure :: exg_cq => swi_swi_cq
    ! ! procedure :: exg_bd => swi_swi_bd
    ! procedure :: exg_ot => swi_swi_ot
    procedure :: exg_da => swi_swi_da
    ! procedure :: exg_fp => swi_swi_fp
    ! procedure :: get_iasym => swi_swi_get_iasym
    procedure :: allocate_scalars
    procedure :: allocate_arrays
    procedure :: source_options
    ! procedure :: source_dimensions
    ! procedure :: source_data
    ! procedure :: noder
    ! procedure :: cellstr
    procedure :: connects_model => swi_swi_connects_model

  end type SwiSwiExchangeType

contains

  subroutine swiswi_cr(filename, name, id, m1_id, m2_id, input_mempath)
    ! dummy
    character(len=*), intent(in) :: filename !< filename for reading
    character(len=*) :: name !< exchange name
    integer(I4B), intent(in) :: id !< id for the exchange
    integer(I4B), intent(in) :: m1_id !< id for model 1
    integer(I4B), intent(in) :: m2_id !< id for model 2
    character(len=*), intent(in) :: input_mempath
    ! local
    type(SwiSwiExchangeType), pointer :: exchange
    class(BaseModelType), pointer :: mb
    class(BaseExchangeType), pointer :: baseexchange
    integer(I4B) :: m1_index, m2_index

    ! Create a new exchange and add it to the baseexchangelist container
    allocate (exchange)
    baseexchange => exchange
    call AddBaseExchangeToList(baseexchangelist, baseexchange)

    ! Assign id and name
    exchange%id = id
    exchange%name = name
    exchange%memoryPath = create_mem_path(exchange%name)
    exchange%input_mempath = input_mempath

    ! allocate scalars and set defaults
    call exchange%allocate_scalars()
    exchange%filename = filename
    exchange%typename = 'SWI-SWI'

    ! -- set gwf_fresh
    m1_index = model_loc_idx(m1_id)
    if (m1_index > 0) then
      mb => GetBaseModelFromList(basemodellist, m1_index)
      select type (mb)
      type is (GwfModelType)
        ! exchange%model1 => mb
        exchange%gwf_fresh => mb
      end select
    end if
    ! exchange%v_model1 => get_virtual_model(m1_id)
    ! exchange%is_datacopy = .not. exchange%v_model1%is_local
    !
    ! -- set gwf_salt
    m2_index = model_loc_idx(m2_id)
    if (m2_index > 0) then
      mb => GetBaseModelFromList(basemodellist, m2_index)
      select type (mb)
      type is (GwfModelType)
        ! exchange%model2 => mb
        exchange%gwf_salt => mb
      end select
    end if
    ! exchange%v_model2 => get_virtual_model(m2_id)
    !
    ! -- Verify that gwf model1 is of the correct type
    if (.not. associated(exchange%gwf_fresh) .and. m1_index > 0) then
      write (errmsg, '(3a)') 'Problem with SWI-SWI exchange ', &
        trim(exchange%name), &
        '.  First specified GWF Model does not appear to be of the &
        &correct type.'
      call store_error(errmsg, terminate=.true.)
    end if
    !
    ! -- Verify that gwf model2 is of the correct type
    if (.not. associated(exchange%gwf_salt) .and. m2_index > 0) then
      write (errmsg, '(3a)') 'Problem with SWI-SWI exchange ', &
        trim(exchange%name), &
        '.  Second specified GWF Model does not appear to be of the &
        &correct type.'
      call store_error(errmsg, terminate=.true.)
    end if

  end subroutine swiswi_cr

  !> @ brief Define SWI SWI exchange
  !!
  !! Define SWI to SWI exchange object.
  !<
  subroutine swi_swi_df(this)
    ! modules
    use SimVariablesModule, only: iout
    ! dummy
    class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType
    ! local

    ! log the exchange
    write (iout, '(/a,a)') ' Creating exchange: ', this%name

    ! -- Ensure models are in same solution
    if (associated(this%gwf_fresh) .and. associated(this%gwf_salt)) then
      if (this%gwf_fresh%idsoln /= this%gwf_salt%idsoln) then
        call store_error('Two models are connected in a GWF '// &
                         'exchange but they are in different solutions. '// &
                         'GWF models must be in same solution: '// &
                         trim(this%gwf_fresh%name)//' '// &
                         trim(this%gwf_salt%name))
        call store_error_filename(this%filename)
      end if
    end if

    ! Ensure fresh model has active SWI Package
    if (this%gwf_fresh%inswi == 0) then
      call store_error( &
        'A SWI-SWI exchange is active, but the freshwater GWF model &
        &does not have an active Seawater Intrusion (Package).  Activate &
        &the SWI Package for the freshwater GWF model.')
      call store_error_filename(this%filename)
    end if

    ! Ensure salt model has active SWI Package
    if (this%gwf_salt%inswi == 0) then
      call store_error( &
        'A SWI-SWI exchange is active, but the saltwater GWF model &
        &does not have an active Seawater Intrusion (Package).  Activate &
        &the SWI Package for the saltwater GWF model.')
      call store_error_filename(this%filename)
    end if

    ! Set the configuration for each SWI package
    call this%gwf_fresh%swi%set_configuration("FRESHWATER_TWO_FLUID")
    call this%gwf_salt%swi%set_configuration("SALTWATER_TWO_FLUID")

    ! Check to make sure fresh and salt models have same number of nodes
    if (this%gwf_fresh%dis%nodes /= this%gwf_salt%dis%nodes) then
      call store_error( &
        'A SWI-SWI exchange is active, but the fresh GWF model &
        &does not have the same number of nodes as the salt &
        &GWF model.  Both models must have the same discretization &
        & properties.')
      call store_error_filename(this%filename)
    end if

    ! Set the number of exchanges equal to the number of equations in
    ! the fresh or saltwater models
    this%nexg = this%gwf_fresh%dis%nodes

    ! source options
    call this%source_options(iout)

    ! ! -- source dimensions
    ! call this%source_dimensions(iout)

    ! allocate arrays
    call this%allocate_arrays()

    ! ! -- source exchange data
    ! call this%source_data(iout)

    return
  end subroutine swi_swi_df

  !> @ brief Add connections
  !!
  !! Override parent exg_ac so that connections can be added here.
  !<
  subroutine swi_swi_ac(this, sparse)
    ! modules
    use SparseModule, only: sparsematrix
    ! dummy
    class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType
    type(sparsematrix), intent(inout) :: sparse
    ! local
    integer(I4B) :: n, iglo, jglo

    ! add exchange connections
    do n = 1, this%nexg
      iglo = n + this%gwf_fresh%moffset
      jglo = n + this%gwf_salt%moffset
      call sparse%addconnection(iglo, jglo, 1)
      call sparse%addconnection(jglo, iglo, 1)
    end do

  end subroutine swi_swi_ac

  !> @ brief Map connections
  !!
  !! Map the connections in the global matrix
  !<
  subroutine swi_swi_mc(this, matrix_sln)
    ! modules
    use SparseModule, only: sparsematrix
    ! dummy
    class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType
    class(MatrixBaseType), pointer :: matrix_sln !< the system matrix
    ! local
    integer(I4B) :: n, iglo, jglo

    ! map exchange connections
    do n = 1, this%nexg
      iglo = n + this%gwf_fresh%moffset
      jglo = n + this%gwf_salt%moffset
      this%idxglo(n) = matrix_sln%get_position(iglo, jglo)
      this%idxsymglo(n) = matrix_sln%get_position(jglo, iglo)
    end do

  end subroutine swi_swi_mc

  !> @ brief Allocate and read
  !<
  subroutine swi_swi_ar(this)
    ! dummy
    class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType

    ! Set head pointers for fresh and salt model swi packages.
    call this%gwf_fresh%swi%set_head_pointers(this%gwf_fresh%x, &
                                              this%gwf_salt%x, &
                                              this%gwf_fresh%xold, &
                                              this%gwf_salt%xold)
    call this%gwf_salt%swi%set_head_pointers(this%gwf_fresh%x, &
                                             this%gwf_salt%x, &
                                             this%gwf_fresh%xold, &
                                             this%gwf_salt%xold)
  end subroutine swi_swi_ar

  !> @ brief Read and prepare
  !<
  subroutine swi_swi_rp(this)
    ! modules
    use TdisModule, only: readnewdata
    ! dummy
    class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType

    ! -- Check with TDIS on whether or not it is time to RP
    if (.not. readnewdata) return

  end subroutine swi_swi_rp

  !> @ brief Advance
  !<
  subroutine swi_swi_ad(this)
    ! dummy
    class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType
  end subroutine swi_swi_ad

  !> @ brief Calculate coefficients
  !<
  subroutine swi_swi_cf(this, kiter)
    ! dummy
    class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType
    integer(I4B), intent(in) :: kiter
  end subroutine swi_swi_cf

  !> @ brief Fill coefficients
  !!
  !! Fill coefficient matrix with cross exchange terms
  !<
  subroutine swi_swi_fc(this, kiter, matrix_sln, rhs_sln, inwtflag)
    ! modules
    ! dummy
    class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType
    integer(I4B), intent(in) :: kiter
    class(MatrixBaseType), pointer :: matrix_sln
    real(DP), dimension(:), intent(inout) :: rhs_sln
    integer(I4B), optional, intent(in) :: inwtflag
    ! local
    integer(I4B) :: n, iglo, jglo
    real(DP) :: term

    ! -- TODO: Since these terms are all Newton, do they go on the
    !    here in _fc or should we put in _fn
    !  ** this is not true.  they are not all newton terms.  The
    !     storage term for the freshwater model has a saltwater model
    !     dependency and vice versa.

    ! -- Put terms into amatsln and rhs
    do n = 1, this%nexg

      iglo = n + this%gwf_fresh%moffset
      jglo = n + this%gwf_salt%moffset
      term = DZERO ! todo: need to calculate this term

      ! fill off-diagonal matrix coefficient
      call matrix_sln%set_value_pos(this%idxglo(n), term)
      call matrix_sln%set_value_pos(this%idxsymglo(n), term)

      ! fill diagonal position for fresh model and salt model nodes
      call matrix_sln%add_diag_value(iglo, -term)
      call matrix_sln%add_diag_value(jglo, -term)
    end do
    ! !
    ! ! -- Set inwt to exchange newton, but shut off if requested by caller
    ! inwt = this%inewton
    ! if (present(inwtflag)) then
    !   if (inwtflag == 0) inwt = 0
    ! end if
    ! if (inwt /= 0) then
    !   call this%exg_fn(kiter, matrix_sln)
    ! end if
    !
    ! -- Return
    return
  end subroutine swi_swi_fc

  ! !> @ brief Fill Newton
  ! !!
  ! !! Fill amatsln with Newton terms
  ! !<
  ! subroutine swi_swi_fn(this, kiter, matrix_sln)
  !   ! -- modules
  !   use SmoothingModule, only: sQuadraticSaturationDerivative
  !   ! -- dummy
  !   class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType
  !   integer(I4B), intent(in) :: kiter
  !   class(MatrixBaseType), pointer :: matrix_sln
  !   ! -- local
  !   !
  !   ! -- Return
  !   return
  ! end subroutine swi_swi_fn

  ! !> @ brief Output
  ! !!
  ! !! Write output
  ! !<
  ! subroutine swi_swi_ot(this)
  !   ! -- modules
  !   ! -- dummy
  !   class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType
  !   ! -- local
  !   !
  !   ! -- Return
  !   return
  ! end subroutine swi_swi_ot

  !> @ brief Source options
  !!
  !! Source the options block
  !<
  subroutine source_options(this, iout)
    ! -- modules
    use MemoryManagerExtModule, only: mem_set_value
    use ExgSwiswiInputModule, only: ExgSwiswiParamFoundType
    use SourceCommonModule, only: filein_fname
    ! -- dummy
    class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType
    integer(I4B), intent(in) :: iout
    ! -- local
    type(ExgSwiswiParamFoundType) :: found
    !
    ! -- update defaults with idm sourced values
    call mem_set_value(this%ipr_input, 'IPR_INPUT', this%input_mempath, &
                       found%ipr_input)
    call mem_set_value(this%ipr_flow, 'IPR_FLOW', this%input_mempath, &
                       found%ipr_flow)
    !
    write (iout, '(1x,a)') 'Processing SWI-SWI exchange options'
    !
    if (found%ipr_input) then
      write (iout, '(4x,a)') &
        'Exchange information will be printed to the listing file.'
    end if
    !
    if (found%ipr_flow) then
      write (iout, '(4x,a)') &
        'Exchange flows will be printed to list file.'
    end if
    !
    write (iout, '(1x,a)') 'End processing SWI-SWI exchange options'
    !
    ! -- Return
    return
  end subroutine source_options

  ! !> @ brief Source dimensions
  ! !!
  ! !! Source the dimensions block
  ! !<
  ! subroutine source_dimensions(this, iout)
  !   ! -- modules
  !   use ConstantsModule, only: LENVARNAME, DEM6
  !   use InputOutputModule, only: getunit, openfile
  !   use MemoryManagerExtModule, only: mem_set_value
  !   use CharacterStringModule, only: CharacterStringType
  !   use ExgSwiswiInputModule, only: ExgSwiswiParamFoundType
  !   use SourceCommonModule, only: filein_fname
  !   ! -- dummy
  !   class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType
  !   integer(I4B), intent(in) :: iout
  !   ! -- local
  !   type(ExgSwiswiParamFoundType) :: found
  !   !
  !   ! -- update defaults with idm sourced values
  !   call mem_set_value(this%nexg, 'NEXG', this%input_mempath, found%nexg)
  !   !
  !   write (iout, '(1x,a)') 'PROCESSING SWI-SWI EXCHANGE DIMENSIONS'
  !   !
  !   if (found%nexg) then
  !     write (iout, '(4x,a,i0)') &
  !       'Number of exchanges (NEXG) set to', this%nexg
  !   else
  !     write (errmsg, '(a)') 'Error in DIMENSIONS block: NEXG not found.'
  !     call store_error(errmsg, terminate=.false.)
  !     call store_error_filename(this%filename)
  !   end if
  !   !
  !   write (iout, '(1x,a)') 'END OF SWI-SWI EXCHANGE DIMENSIONS'
  !   !
  !   ! -- Return
  !   return
  ! end subroutine source_dimensions

  ! !> @brief Source exchange data from input context
  ! !<
  ! subroutine source_data(this, iout)
  !   ! -- modules
  !   use MemoryManagerModule, only: mem_setptr
  !   ! -- dummy
  !   class(SwiSwiExchangeType) :: this !< instance of exchange object
  !   integer(I4B), intent(in) :: iout !< the output file unit
  !   ! -- local
  !   integer(I4B), dimension(:, :), contiguous, pointer :: cellidm1
  !   integer(I4B), dimension(:, :), contiguous, pointer :: cellidm2
  !   real(DP), dimension(:), contiguous, pointer :: amat_fresh
  !   real(DP), dimension(:), contiguous, pointer :: rhs_fresh
  !   real(DP), dimension(:), contiguous, pointer :: amat_salt
  !   real(DP), dimension(:), contiguous, pointer :: rhs_salt
  !   character(len=20) :: cellstr1, cellstr2
  !   integer(I4B) :: nerr
  !   integer(I4B) :: iexg, nodem1, nodem2
  !   ! -- format
  !   character(len=*), parameter :: fmtexglabel = "(1x, 3a10, 50(a16))"
  !   character(len=*), parameter :: fmtexgdata = &
  !                                  "(5x, a, 1x, a ,I10, 50(1pg16.6))"
  !   !
  !   call mem_setptr(cellidm1, 'CELLIDM1', this%input_mempath)
  !   call mem_setptr(cellidm2, 'CELLIDM2', this%input_mempath)
  !   call mem_setptr(amat_fresh, 'AMAT_FRESH', this%input_mempath)
  !   call mem_setptr(rhs_fresh, 'RHS_FRESH', this%input_mempath)
  !   call mem_setptr(amat_salt, 'AMAT_SALT', this%input_mempath)
  !   call mem_setptr(rhs_salt, 'RHS_SALT', this%input_mempath)
  !   !
  !   write (iout, '(1x,a)') 'PROCESSING EXCHANGEDATA'
  !   !
  !   if (this%ipr_input /= 0) then
  !     write (iout, fmtexglabel) 'NODEM1', 'NODEM2', 'AMAT_FRESH', &
  !       'RHS_FRESH', 'AMAT_SALT', 'RHS_SALT'
  !   end if
  !   !
  !   do iexg = 1, this%nexg
  !     !
  !     if (associated(this%gwf_fresh)) then
  !       !
  !       ! -- Determine user node number
  !       nodem1 = this%noder(this%gwf_fresh, cellidm1(:, iexg), iout)
  !       this%nodem1(iexg) = nodem1
  !       !
  !     else
  !       this%nodem1(iexg) = -1
  !     end if
  !     !
  !     if (associated(this%gwf_salt)) then
  !       !
  !       ! -- Determine user node number
  !       nodem2 = this%noder(this%gwf_salt, cellidm2(:, iexg), iout)
  !       this%nodem2(iexg) = nodem2
  !       !
  !     else
  !       this%nodem2(iexg) = -1
  !     end if
  !     !
  !     ! -- Read rest of input line
  !     this%amat_fresh(iexg) = amat_fresh(iexg)
  !     this%rhs_fresh(iexg) = rhs_fresh(iexg)
  !     this%amat_salt(iexg) = amat_salt(iexg)
  !     this%rhs_salt(iexg) = rhs_salt(iexg)
  !     !
  !     ! -- Write the data to listing file if requested
  !     if (this%ipr_input /= 0) then
  !       cellstr1 = this%cellstr(this%gwf_fresh, cellidm1(:, iexg), iout)
  !       cellstr2 = this%cellstr(this%gwf_salt, cellidm2(:, iexg), iout)
  !       write (iout, fmtexgdata) trim(cellstr1), trim(cellstr2), &
  !         this%amat_fresh(iexg), this%rhs_fresh(iexg), &
  !         this%amat_salt(iexg), this%rhs_salt(iexg)
  !     end if
  !     !
  !     ! -- Check to see if nodem1 is outside of active domain
  !     if (associated(this%gwf_fresh)) then
  !       if (nodem1 <= 0) then
  !         cellstr1 = this%cellstr(this%gwf_fresh, cellidm1(:, iexg), iout)
  !         write (errmsg, *) &
  !           trim(adjustl(this%gwf_fresh%name))// &
  !           ' Cell is outside active grid domain ('// &
  !           trim(adjustl(cellstr1))//').'
  !         call store_error(errmsg)
  !       end if
  !     end if
  !     !
  !     ! -- Check to see if nodem2 is outside of active domain
  !     if (associated(this%gwf_salt)) then
  !       if (nodem2 <= 0) then
  !         cellstr2 = this%cellstr(this%gwf_salt, cellidm2(:, iexg), iout)
  !         write (errmsg, *) &
  !           trim(adjustl(this%gwf_salt%name))// &
  !           ' Cell is outside active grid domain ('// &
  !           trim(adjustl(cellstr2))//').'
  !         call store_error(errmsg)
  !       end if
  !     end if
  !   end do
  !   !
  !   write (iout, '(1x,a)') 'END OF EXCHANGEDATA'
  !   !
  !   ! -- Stop if errors
  !   nerr = count_errors()
  !   if (nerr > 0) then
  !     call store_error('Errors encountered in exchange input file.')
  !     call store_error_filename(this%filename)
  !   end if
  !   !
  !   ! -- Return
  !   return
  ! end subroutine source_data

  !> @ brief Allocate scalars
  !!
  !! Allocate scalar variables
  !<
  subroutine allocate_scalars(this)
    ! modules
    use MemoryManagerModule, only: mem_allocate
    use ConstantsModule, only: DZERO
    ! dummy
    class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType

    allocate (this%filename)
    this%filename = ''

    call mem_allocate(this%ipr_input, 'IPR_INPUT', this%memoryPath)
    call mem_allocate(this%ipr_flow, 'IPR_FLOW', this%memoryPath)
    call mem_allocate(this%nexg, 'NEXG', this%memoryPath)

    this%ipr_input = 0
    this%ipr_flow = 0
    this%nexg = 0

  end subroutine allocate_scalars

  !> @ brief Deallocate
  !!
  !! Deallocate memory associated with this object
  !<
  subroutine swi_swi_da(this)
    ! modules
    use MemoryManagerModule, only: mem_deallocate
    use MemoryManagerExtModule, only: memorystore_remove
    use SimVariablesModule, only: idm_context
    ! dummy
    class(SwiSwiExchangeType) :: this !< SwiSwiExchangeType

    ! deallocate IDM memory
    call memorystore_remove(this%name, '', idm_context)

    ! arrays
    call mem_deallocate(this%idxglo)
    call mem_deallocate(this%idxsymglo)

    ! scalars
    deallocate (this%filename)
    call mem_deallocate(this%ipr_input)
    call mem_deallocate(this%ipr_flow)
    call mem_deallocate(this%nexg)

  end subroutine swi_swi_da

  !> @ brief Allocate arrays
  !!
  !! Allocate arrays
  !<
  subroutine allocate_arrays(this)
    ! modules
    use MemoryManagerModule, only: mem_allocate
    ! dummy
    class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType
    ! local
    integer(I4B) :: i

    ! call mem_allocate(this%nodem1, this%nexg, 'NODEM1', this%memoryPath)
    ! call mem_allocate(this%nodem2, this%nexg, 'NODEM2', this%memoryPath)
    ! call mem_allocate(this%amat_fresh, this%nexg, 'AMAT_FRESH', this%memoryPath)
    ! call mem_allocate(this%rhs_fresh, this%nexg, 'RHS_FRESH', this%memoryPath)
    ! call mem_allocate(this%amat_salt, this%nexg, 'AMAT_SALT', this%memoryPath)
    ! call mem_allocate(this%rhs_salt, this%nexg, 'RHS_SALT', this%memoryPath)
    call mem_allocate(this%idxglo, this%nexg, 'IDXGLO', this%memoryPath)
    call mem_allocate(this%idxsymglo, this%nexg, 'IDXSYMGLO', this%memoryPath)
    ! call mem_allocate(this%simvals, this%nexg, 'SIMVALS', this%memoryPath)

    ! Initialize
    do i = 1, this%nexg
      this%idxglo(i) = 0
      this%idxsymglo(i) = 0
      ! this%amat_fresh(i) = DNODATA
      ! this%rhs_fresh(i) = DNODATA
      ! this%amat_salt(i) = DNODATA
      ! this%rhs_salt(i) = DNODATA
    end do
    !
    ! -- Return
    return
  end subroutine allocate_arrays

  ! !> @ brief Final processing
  ! !!
  ! !! Conduct any final processing
  ! !<
  ! subroutine swi_swi_fp(this)
  !   ! -- dummy
  !   class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType
  !   !
  !   ! -- Return
  !   return
  ! end subroutine swi_swi_fp

  ! !> @ brief Set symmetric flag
  ! !!
  ! !! Return flag indicating whether or not this exchange will cause the
  ! !! coefficient matrix to be asymmetric.
  ! !<
  ! function swi_swi_get_iasym(this) result(iasym)
  !   ! -- dummy
  !   class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType
  !   ! -- local
  !   integer(I4B) :: iasym
  !   !
  !   ! -- Newton formulation results in asymmetric matrix
  !   iasym = 1
  !   !
  !   ! -- Return
  !   return
  ! end function swi_swi_get_iasym

  ! !> @brief
  ! !<
  ! function noder(this, model, cellid, iout)
  !   ! -- modules
  !   use GeomUtilModule, only: get_node
  !   ! -- dummy
  !   class(SwiSwiExchangeType) :: this !< instance of exchange object
  !   class(GwfModelType), pointer, intent(in) :: model
  !   integer(I4B), dimension(:), pointer, intent(in) :: cellid
  !   integer(I4B), intent(in) :: iout !< the output file unit
  !   integer(I4B) :: noder, node
  !   !
  !   if (model%dis%ndim == 1) then
  !     node = cellid(1)
  !   elseif (model%dis%ndim == 2) then
  !     node = get_node(cellid(1), 1, cellid(2), &
  !                     model%dis%mshape(1), 1, &
  !                     model%dis%mshape(2))
  !   else
  !     node = get_node(cellid(1), cellid(2), cellid(3), &
  !                     model%dis%mshape(1), &
  !                     model%dis%mshape(2), &
  !                     model%dis%mshape(3))
  !   end if
  !   noder = model%dis%get_nodenumber(node, 0)
  !   !
  !   ! -- return
  !   return
  ! end function noder

  ! !> @brief
  ! !<
  ! function cellstr(this, model, cellid, iout)
  !   ! -- modules
  !   ! -- dummy
  !   class(SwiSwiExchangeType) :: this !< instance of exchange object
  !   class(GwfModelType), pointer, intent(in) :: model
  !   integer(I4B), dimension(:), pointer, intent(in) :: cellid
  !   integer(I4B), intent(in) :: iout !< the output file unit
  !   character(len=20) :: cellstr
  !   character(len=*), parameter :: fmtndim1 = &
  !                                  "('(',i0,')')"
  !   character(len=*), parameter :: fmtndim2 = &
  !                                  "('(',i0,',',i0,')')"
  !   character(len=*), parameter :: fmtndim3 = &
  !                                  "('(',i0,',',i0,',',i0,')')"
  !   !
  !   cellstr = ''
  !   !
  !   select case (model%dis%ndim)
  !   case (1)
  !     write (cellstr, fmtndim1) cellid(1)
  !   case (2)
  !     write (cellstr, fmtndim2) cellid(1), cellid(2)
  !   case (3)
  !     write (cellstr, fmtndim3) cellid(1), cellid(2), cellid(3)
  !   case default
  !   end select
  !   !
  !   ! -- return
  !   return
  ! end function cellstr

  ! !> @ brief Cast polymorphic object as exchange
  ! !!
  ! !! Cast polymorphic object as exchange
  ! !<
  ! function CastAsSwiExchange(obj) result(res)
  !   implicit none
  !   ! -- dummy
  !   class(*), pointer, intent(inout) :: obj
  !   ! -- return
  !   class(SwiSwiExchangeType), pointer :: res
  !   !
  !   res => null()
  !   if (.not. associated(obj)) return
  !   !
  !   select type (obj)
  !   class is (SwiSwiExchangeType)
  !     res => obj
  !   end select
  !   !
  !   ! -- Return
  !   return
  ! end function CastAsSwiExchange

  ! !> @ brief Get exchange from list
  ! !!
  ! !! Return an exchange from the list for specified index
  ! !<
  ! function GetSwiExchangeFromList(list, idx) result(res)
  !   implicit none
  !   ! -- dummy
  !   type(ListType), intent(inout) :: list
  !   integer(I4B), intent(in) :: idx
  !   ! -- return
  !   class(SwiSwiExchangeType), pointer :: res
  !   ! -- local
  !   class(*), pointer :: obj
  !   !
  !   obj => list%GetItem(idx)
  !   res => CastAsSwiExchange(obj)
  !   !
  !   ! -- Return
  !   return
  ! end function GetSwiExchangeFromList

  !> @brief Should return true when the exchange should be added to the
  !! solution where the model resides
  !<
  function swi_swi_connects_model(this, model) result(is_connected)
    ! dummy
    class(SwiSwiExchangeType) :: this !< the instance of the exchange
    class(BaseModelType), pointer, intent(in) :: model !< the model to which the exchange might hold a connection
    ! return
    logical(LGP) :: is_connected !< true, when connected

    is_connected = .false.
    select type (model)
    class is (GwfModelType)
      if (associated(this%gwf_fresh, model)) then
        is_connected = .true.
      else if (associated(this%gwf_salt, model)) then
        is_connected = .true.
      end if
    end select

  end function

end module SwiSwiExchangeModule
