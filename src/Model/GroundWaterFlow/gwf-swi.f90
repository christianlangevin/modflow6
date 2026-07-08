module GwfSwiModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LINELENGTH, DONE, DZERO, LENBUDTXT, &
                             MNORMAL, DHNOFLO, C3D_VERTICAL, LENMODELNAME, &
                             DEM12
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error, store_error_filename
  use NumericalPackageModule, only: NumericalPackageType
  use BlockParserModule, only: BlockParserType
  use BaseDisModule, only: DisBaseType
  use GwfNpfModule, only: GwfNpfType
  use GwfStoModule, only: GwfStoType
  use SmoothingModule, only: sQuadraticSaturation, &
                             sQuadraticSaturationDerivative
  use MemoryManagerModule, only: mem_setptr, mem_allocate
  use MemoryHelperModule, only: create_mem_path
  use MatrixBaseModule
  use GwfConductanceUtilsModule, only: hcond, vcond
  use UtlTvaModule, only: UtlTvaType
  use GwfNpfFormulationModule, only: GwfNpfFormulationType, SWI_FLOW
  use GwfStoExtModule, only: GwfStoFormulationType, SWI_STORAGE

  implicit none
  private
  public :: GwfSwiType
  public :: swi_cr
  public :: swi_thksat

  character(len=LENBUDTXT), dimension(1) :: budtxt = & !< text labels for budget terms
    &['         STORAGE']
  integer(I4B), parameter :: FRESHWATER_ONE_FLUID = 0 !< freshwater model without swi-swi exchange
  integer(I4B), parameter :: FRESHWATER_TWO_FLUID = 1 !< freshwater model with swi-swi exchange
  integer(I4B), parameter :: SALTWATER_ONE_FLUID = 2 !< saltwater model without swi-swi exchange, not supported
  integer(I4B), parameter :: SALTWATER_TWO_FLUID = 3 !< saltwater model with swi-swi exchange
  real(DP), parameter :: ZONE_MINIMUM_THICKNESS = 1.d-3 !< minimum thickness of saltwater or freshwater zone to avoid numerical issues

  type, extends(NumericalPackageType) :: GwfSwiType

    character(len=LENMODELNAME) :: freshwater_model_name = ''
    character(len=LENMODELNAME) :: saltwater_model_name = ''

    integer(I4B), pointer :: iconfiguration => null() !< 0 is freshwater, 1 is freshwater two fluid, 2 is saltwater, 3 is saltwater two fluid
    integer(I4B), pointer :: isaltwater => null() !< 0 is freshwater, 1 is saltwater
    integer(I4B), pointer :: izetaout => null() !< unit number for binary zeta output file
    integer(I4B), pointer :: intva => null() !< unit number for reading hsalt
    real(DP), dimension(:), pointer, contiguous :: zeta => null() !< starting zeta
    real(DP), dimension(:), pointer, contiguous :: hcof => null() !< hcof contribution to amat
    real(DP), dimension(:), pointer, contiguous :: rhs => null() !< rhs contribution
    real(DP), dimension(:), pointer, contiguous :: storage => null() !< calculated swi storage

    ! information needed for full implementation
    real(DP), pointer :: alphaf => null() !< rhof / (rhos - rhof), default is 40
    real(DP), pointer :: alphas => null() !< rhos / (rhos - rhof), default is 41
    real(DP), pointer :: hsalt_user => null() !< user specified head of saltwater (zero by default)
    real(DP), dimension(:), pointer, contiguous :: hfresh => null() !< head of freshwater
    real(DP), dimension(:), pointer, contiguous :: hsalt => null() !< head of saltwater
    real(DP), dimension(:), pointer, contiguous :: hfreshold => null() !< old head of freshwater
    real(DP), dimension(:), pointer, contiguous :: hsaltold => null() !< old head of saltwater
    integer(I4B), dimension(:), pointer, contiguous :: ibound => null() !< pointer to model ibound

    ! pointers for transient simulations
    integer(I4B), pointer :: insto => null() !< pointer to check of storage package is on
    integer(I4B), pointer :: iss => null() !< pointer to gwf steady state flag
    real(DP), dimension(:), pointer, contiguous :: sy => null() !< pointer to storage package specific yield
    real(DP), dimension(:), pointer, contiguous :: oldsy => null() !< pointer to storage package specific yield
    real(DP), dimension(:), pointer, contiguous :: oldss => null() !< pointer to storage package specific yield

    ! objects
    type(UtlTvaType) :: tva_reader

    ! flow-extension formulation (spike: single-fluid freshwater horizontal flow)
    class(GwfNpfFormulationType), pointer :: npf_form => null()
    ! storage-extension formulation (spike: SWI storage budget/cbc via bd/save_flows)
    class(GwfStoFormulationType), pointer :: sto_form => null()

  contains

    procedure :: swi_df
    procedure :: swi_ar
    procedure :: swi_rp
    procedure :: swi_ad
    procedure :: swi_fc
    procedure :: swi_fn
    procedure :: swi_cc
    procedure :: swi_cq
    procedure :: swi_bd
    procedure :: swi_save_model_flows
    procedure :: swi_ot_dv
    procedure :: swi_da
    procedure, private :: swi_load
    procedure, private :: source_options
    procedure, private :: log_options
    procedure :: allocate_scalars
    procedure, private :: allocate_arrays
    procedure, private :: source_griddata
    procedure, private :: update_zeta
    procedure, private :: swi_fc_storage
    procedure, private :: swi_fn_storage
    procedure, public :: swi_saturated_thickness

    procedure :: set_configuration
    procedure :: set_head_pointers
    procedure :: get_zetanew
    procedure :: get_zetaold

  end type GwfSwiType

  !> @brief NPF flow formulation for SWI (spike: single-fluid freshwater,
  !! horizontal flow). Overrides the default NPF connection conductance with the
  !! freshwater-column conductance S^f = S^w - S^s.
  !<
  type, extends(GwfNpfFormulationType) :: SwiNpfFormulationType
    class(GwfSwiType), pointer :: swi => null() !< owning SWI package
    type(GwfNpfType), pointer :: npf => null() !< the model NPF package
  contains
    procedure :: is_active => swinpf_is_active
    procedure :: cf => swinpf_cf
    procedure :: fc => swinpf_fc
    procedure :: fn => swinpf_fn
    procedure :: cq => swinpf_cq
  end type SwiNpfFormulationType

  !> @brief STO storage formulation for SWI (spike: provides the SWI storage
  !! budget term and cell-by-cell flows via bd/save_flows. The storage matrix
  !! fill (fc/fn/cq) is still done by the existing swi_fc_storage path for now,
  !! so no cells are flagged SWI_STORAGE and fc/fn/cq are not dispatched here.)
  !<
  type, extends(GwfStoFormulationType) :: SwiStoFormulationType
    class(GwfSwiType), pointer :: swi => null() !< owning SWI package
    type(GwfStoType), pointer :: sto => null() !< the model STO package
  contains
    procedure :: is_active => swisto_is_active
    procedure :: fc => swisto_fc
    procedure :: fn => swisto_fn
    procedure :: cq => swisto_cq
    procedure :: bd => swisto_bd
    procedure :: save_flows => swisto_save_flows
  end type SwiStoFormulationType

contains

  !> @brief Create a new swi package object
  !<
  subroutine swi_cr(swi, name_model, input_mempath, inunit, iout, dis)
    ! modules
    use MemoryManagerExtModule, only: mem_set_value
    ! dummy
    type(GwfSwiType), pointer :: swi
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: input_mempath
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    class(DisBaseType), pointer, intent(in) :: dis
    ! formats
    character(len=*), parameter :: fmtswi = &
      "(1x, /1x, 'SWI -- Seawater Intrusion Package, Version 8, 1/8/2024', &
      &' input read from mempath: ', A, //)"
    !
    ! create SWI object
    allocate (swi)
    !
    ! create name and memory path
    call swi%set_names(1, name_model, 'SWI', 'SWI', input_mempath)
    !
    ! allocate scalars
    call swi%allocate_scalars()
    !
    ! set variables
    swi%inunit = inunit
    swi%iout = iout
    !
    ! set pointers
    swi%dis => dis
    !
    ! check if pkg is enabled,
    if (inunit > 0) then
      ! print message identifying pkg
      write (swi%iout, fmtswi) input_mempath
    end if
    !
    ! return
    return
  end subroutine swi_cr

  !> @brief Allocate arrays, load from IDM, and assign head
  !<
  subroutine swi_df(this)
    ! dummy
    class(GwfSwiType) :: this
    ! local

    ! allocate arrays
    call this%allocate_arrays(this%dis%nodes)

    ! load from IDM
    call this%swi_load()

    ! setup the tva reader
    if (this%intva > 0) then
      call this%tva_reader%initialize(this%dis, this%dis%nodesuser, &
                                      this%intva, this%iout, &
                                      this%name_model, 'SWI')
      if (.not. this%tva_reader%has('SALTWATER_HEAD')) then
        write (errmsg, '(a)') 'TVA Input file does not have SALTWATER_HEAD &
          &variable.  Specify SALTWATER_HEAD as a time variable input array in &
          &the TVA input file.'
        call store_error(errmsg, terminate=.false.)
        call store_error_filename(this%input_fname)
      end if
    end if

  end subroutine swi_df

  !> @brief Setup pointers
  !<
  subroutine swi_ar(this, ibound, npf, sto)
    ! dummy
    class(GwfSwiType), target :: this
    integer(I4B), dimension(:), pointer, contiguous :: ibound !< model ibound
    type(GwfNpfType), pointer, intent(inout) :: npf !< model NPF package
    type(GwfStoType), pointer, intent(inout) :: sto !< model STO package
    ! local
    integer(I4B) :: n
    integer(I4B) :: ipos, ihc
    type(SwiNpfFormulationType), pointer :: npf_form
    type(SwiStoFormulationType), pointer :: sto_form

    ! set pointer to ibound
    this%ibound => ibound

    ! set pointer to gwf steady state flag
    call mem_setptr(this%insto, 'INSTO', &
                    create_mem_path(this%name_model))
    call mem_setptr(this%iss, 'ISS', &
                    create_mem_path(this%name_model))
    if (this%insto > 0) then
      call mem_setptr(this%sy, 'SY', &
                      create_mem_path(this%name_model, 'STO'))
    end if

    ! If hfresh is not associated, then this swi package must not
    ! be part of a model that has a swi-swi exchange with another
    ! model.  In this case, isaltwater must be 0 because we assume
    ! that if a model is not connected with a swi-swi exchange then
    ! it must be a freshwater model.
    if (.not. associated(this%hfresh)) then
      call mem_setptr(this%hfresh, 'X', &
                      create_mem_path(this%name_model))
    end if
    if (.not. associated(this%hfreshold)) then
      call mem_setptr(this%hfreshold, 'XOLD', &
                      create_mem_path(this%name_model))
    end if

    ! If this is a single-fluid freshwater model
    ! then memory needs to be allocated for hsalt
    if (this%iconfiguration == FRESHWATER_ONE_FLUID) then

      ! allocate saltwater head array
      call mem_allocate(this%hsalt, size(this%hfresh), 'HSALT', this%memoryPath)
      do n = 1, size(this%hsalt)
        this%hsalt(n) = this%hsalt_user
      end do

      ! allocate old saltwater head array
      call mem_allocate(this%hsaltold, size(this%hfresh), 'HSALTOLD', &
                        this%memoryPath)
      do n = 1, size(this%hsaltold)
        this%hsaltold(n) = this%hsalt_user
      end do
    end if

    ! register the SWI NPF flow formulation and flag its (horizontal)
    ! connections so NPF dispatches them to the SWI override
    allocate (npf_form)
    npf_form%swi => this
    npf_form%npf => npf
    this%npf_form => npf_form
    call npf%add_flow_formulation(this%npf_form, SWI_FLOW)
    do n = 1, npf%dis%nodes
      do ipos = npf%dis%con%ia(n) + 1, npf%dis%con%ia(n + 1) - 1
        ihc = npf%dis%con%ihc(npf%dis%con%jas(ipos))
        if (ihc /= C3D_VERTICAL) npf%iformulation(ipos) = SWI_FLOW
      end do
    end do

    ! register the SWI STO storage formulation and flag its cells so STO
    ! dispatches storage fc/fn/cq (and bd/save_flows) to the SWI override. The
    ! override assembles freshwater storage as the water column (S^w, via the STO
    ! default) minus the salt column (S^s); the SWI storage budget term is the
    ! interface-movement (fresh<->salt conversion) volume.
    if (this%insto > 0) then
      allocate (sto_form)
      sto_form%swi => this
      sto_form%sto => sto
      this%sto_form => sto_form
      call sto%add_sto_formulation(this%sto_form, SWI_STORAGE)
      do n = 1, sto%dis%nodes
        sto%iformulation(n) = SWI_STORAGE
      end do
    end if

  end subroutine swi_ar

  !> @brief Read and prepare
  !<
  subroutine swi_rp(this)
    ! dummy
    class(GwfSwiType) :: this
    ! local

    if (this%intva > 0) then
      call this%tva_reader%tva_rp()
    end if

  end subroutine swi_rp

  !> @brief Advance
  !<
  subroutine swi_ad(this, irestore)
    ! dummy
    class(GwfSwiType) :: this
    integer(I4B), intent(in) :: irestore
    ! local
    integer(I4B) :: n

    if (this%iconfiguration == FRESHWATER_ONE_FLUID) then
      if (irestore == 0) then
        ! advance hsaltold
        do n = 1, size(this%hsalt)
          this%hsaltold(n) = this%hsalt(n)
        end do
      else
        ! redoing time step so reset hsalt to hsaltold
        do n = 1, size(this%hsalt)
          this%hsalt(n) = this%hsaltold(n)
        end do
      end if
    end if

    ! if time array series used, then update hsalt with new values
    if (this%intva > 0) then
      call this%tva_reader%tva_ad()
      if (this%iconfiguration == FRESHWATER_ONE_FLUID) then
        call this%tva_reader%fill(this%hsalt, 'SALTWATER_HEAD')
      end if
    end if

  end subroutine swi_ad

  !> @ brief Fill A and right-hand side for the package
  !!
  !!  Fill the coefficient matrix and right-hand side
  !!
  !<
  subroutine swi_fc(this, kiter, hold, hnew, matrix_sln, idxglo, rhs, npf, &
                    sto, inwt)
    ! modules
    ! dummy variables
    class(GwfSwiType) :: this !< GwfSwiType object
    integer(I4B), intent(in) :: kiter !< outer iteration number
    real(DP), intent(in), dimension(:) :: hold !< previous heads
    real(DP), intent(in), dimension(:) :: hnew !< current heads
    class(MatrixBaseType), pointer :: matrix_sln !< A matrix
    integer(I4B), intent(in), dimension(:) :: idxglo !< global index model to solution
    real(DP), intent(inout), dimension(:) :: rhs !< right-hand side
    type(GwfNpfType), intent(in) :: npf
    type(GwfStoType), intent(in) :: sto
    ! local variables
    integer(I4B) :: inwt
    ! formats
    !
    ! fill zeta
    call this%update_zeta()
    !
    ! SPIKE: horizontal freshwater flow is now filled by the SwiNpfFormulationType
    ! flow formulation (dispatched from npf_fc), not the npf_fc_swi correction.
    ! call npf_fc_swi(npf, kiter, matrix_sln, idxglo, &
    !                 rhs, hnew, this%zeta, this%isaltwater)
    if (this%isaltwater == 1) then
      call npf_fc_vflow(npf, kiter, matrix_sln, idxglo, &
                        rhs, hnew)
    end if
    !
    ! Storage is now assembled by the SWI STO formulation (swisto_fc), dispatched
    ! from sto_fc for the flagged SWI_STORAGE cells -- no bolt-on correction here.
    !
    ! return
    return
  end subroutine swi_fc

  !> @ brief Calculate fresh storage terms
  !<
  subroutine swi_fc_storage(this, kiter, hold, hnew, matrix_sln, idxglo, &
                            rhs, sto, inwt, isaltwater)
    ! modules
    use TdisModule, only: delt
    use GwfStorageUtilsModule, only: SsCapacity, SyCapacity, SsTerms, SyTerms
    !
    class(GwfStoType) :: sto
    ! dummy variables
    class(GwfSwiType) :: this !< GwfSwiType object
    integer(I4B), intent(in) :: kiter !< outer iteration number
    real(DP), intent(in), dimension(:) :: hold !< previous heads
    real(DP), intent(in), dimension(:) :: hnew !< current heads
    class(MatrixBaseType), pointer :: matrix_sln !< A matrix
    integer(I4B), intent(in), dimension(:) :: idxglo !< global index model to solution
    real(DP), intent(inout), dimension(:) :: rhs !< right-hand side
    integer(I4B) :: inwt, isaltwater
    ! local variables
    integer(I4B) :: n
    integer(I4B) :: idiag
    real(DP) :: tled
    real(DP) :: sc1
    real(DP) :: sc2
    real(DP) :: rho1
    real(DP) :: rho2
    real(DP) :: sc1old
    real(DP) :: sc2old
    real(DP) :: rho1old
    real(DP) :: rho2old
    real(DP) :: tp
    real(DP) :: bt
    real(DP) :: snold
    real(DP) :: snnew
    real(DP) :: aterm
    real(DP) :: rhsterm
    real(DP) :: zetanew
    real(DP) :: zetaold

    ! set variables
    tled = DONE / delt

    ! Calculate storage terms for freshwater model by subtracting out saltwater side
    do n = 1, this%dis%nodes

      ! initialize
      this%hcof(n) = DZERO
      this%rhs(n) = DZERO

      if (this%ibound(n) < 1) cycle

      ! calculate zetanew and zetaold
      zetanew = this%get_zetanew(n)
      zetaold = this%get_zetaold(n)

      ! aquifer elevations and thickness
      tp = this%dis%top(n)
      bt = this%dis%bot(n)

      ! aquifer saturation
      if (isaltwater == 0) then
        ! Freshwater model, freshwater saturation is considered to be 1
        ! if iconvert == 0, otherwise old and new saturation are calculated
        ! using the old and new head values.  In this case snnew and snold
        ! are the total saturated cell fractions.
        if (sto%iconvert(n) == 0) then
          snold = DONE
          snnew = DONE
        else
          snold = sQuadraticSaturation(tp, bt, hold(n), sto%satomega)
          snnew = sQuadraticSaturation(tp, bt, hnew(n), sto%satomega)
        end if
      else
        ! Saltwater model, calculate old and new saturations using the
        ! old and new zeta values. In this case snnew and snold
        ! are the saltwater saturated cell fractions.
        snold = sQuadraticSaturation(tp, bt, zetaold, sto%satomega)
        snnew = sQuadraticSaturation(tp, bt, zetanew, sto%satomega)
      end if

      ! storage coefficients
      sc1 = SsCapacity(sto%istor_coef, tp, bt, this%dis%area(n), sto%ss(n))
      rho1 = sc1 * tled

      ! Handle time-varying storage case
      if (sto%integratechanges /= 0) then
        ! Integration of storage changes (e.g. when using TVS):
        !    separate the old (start of time step) and new (end of time step)
        !    primary storage capacities
        sc1old = SsCapacity(sto%istor_coef, tp, bt, this%dis%area(n), &
                            this%oldss(n))
        rho1old = sc1old * tled
      else
        ! No integration of storage changes: old and new values are
        !    identical => normal MF6 storage formulation
        rho1old = rho1
      end if

      ! calculate specific storage terms
      call SsTerms(sto%iconvert(n), sto%iorig_ss, sto%iconf_ss, tp, bt, &
                   rho1, rho1old, snnew, snold, hnew(n), hold(n), &
                   aterm, rhsterm)

      ! scale down the saltwater part
      ! THIS APPEARS TO RECALCULATE THE SALTWATER FRACTION OF THE CELL
      ! IN SNNEW.  THEN THE HCOF AND RHS TERMS ARE MULTIPLIED BY THIS
      ! SALTWATER FRACTION.  FOR THE FRESHWATER MODEL, THE COMPRESSIBLE
      ! STORAGE CHANGE IS SUBTRACTED AS IT IS A CORRECTION.  FOR THE
      ! SALTWATER MODEL, WHICH IS NOT TREATED AS A CORRECTION, THE ENTIRE
      ! COMPRESSIBLE STORAGE CHANGE IS ADDED TO HCOF AND RHS.
      ! IS IT OKAY TO USE A FULLY FORWARD-WEIGHTED SNNEW HERE FOR THE
      ! SALTWATER FRACTION OR SHOULD IT BE TIME AVERAGED?
      snnew = sQuadraticSaturation(tp, bt, zetanew, sto%satomega)
      aterm = aterm * snnew
      rhsterm = rhsterm * snnew

      ! add specific storage terms to amat and rhs -
      ! subtract out aterm and rhsterm from the saltwater zone for freshwater
      if (isaltwater == 1) then
        ! If this is saltwater rmodel, then flip sign on aterm and rhsterm
        aterm = -aterm
        rhsterm = -rhsterm
      end if
      idiag = this%dis%con%ia(n)
      call matrix_sln%add_value_pos(idxglo(idiag), -aterm)
      rhs(n) = rhs(n) - rhsterm

      ! THIS IS WHERE THERE WAS LIKELY A BUG IN THE IMPLEMENTATION.  HERE WE SEE THAT THE SY
      ! CONTRIBUTION IS ONLY APPLIED FOR ICONVERT /= 0.  HOWEVER, WE STILL NEED
      ! TO APPLY AN SY PART EVEN IF ICONVERT == 0 TO HANDLE FOR STORAGE CHANGES
      ! RESULTING FROM CHANGE IN THE ZETA ELEVATION.  COMMENTING OUT THIS IF STATEMENT
      ! SO THAT THE SY CONTRIBUTION IS ALWAYS APPLIED.
      ! specific yield
      ! cdl -- if (sto%iconvert(n) /= 0) then
      rhsterm = DZERO

      ! calculate the saltwater new and old saturations
      snold = sQuadraticSaturation(tp, bt, zetaold, sto%satomega)
      snnew = sQuadraticSaturation(tp, bt, zetanew, sto%satomega)

      ! secondary storage coefficient
      sc2 = SyCapacity(this%dis%area(n), this%sy(n))
      rho2 = sc2 * tled

      if (sto%integratechanges /= 0) then
        ! Integration of storage changes (e.g. when using TVS):
        ! separate the old (start of time step) and new (end of time step)
        ! secondary storage capacities
        sc2old = SyCapacity(this%dis%area(n), this%oldsy(n))
        rho2old = sc2old * tled
      else
        ! No integration of storage changes: old and new values are
        ! identical => normal MF6 storage formulation
        rho2old = rho2
      end if

      ! Calculate specific yield terms from bot to zeta. This will
      ! represent the change in storage due to change in zeta elevation.
      ! This storage term needs to occur for both iconvert == 0 and iconvert /= 0
      if (inwt /= 0) then
        call SyTerms(tp, bt, rho2, rho2old, snnew, snold, &
                     aterm, rhsterm)

        ! add specific yield terms to amat and rhs -
        ! subtract out aterm and rhsterm from the saltwater zone for freshwater

        ! THEREFORE, for saltwater, flip sign on aterm and rhsterm
        if (isaltwater == 1) then
          aterm = -aterm
          rhsterm = -rhsterm
        end if

        idiag = this%dis%con%ia(n)
        call matrix_sln%add_value_pos(idxglo(idiag), -aterm)
        rhs(n) = rhs(n) - rhsterm
      else
        !
        rho2 = -rho2 * this%alphaf
        rho2old = -rho2old * this%alphaf
        ! add specific yield terms to hcof and rhs -
        ! subtract out terms from the saltwater zone for freshwater
        !
        ! THEREFORE, for saltwater, flip sign on terms
        if (isaltwater == 1) then
          rho2 = -rho2
          rho2old = -rho2old
        end if
        !-------------------------------------------------------------
        if (zetanew > this%dis%bot(n) .and. zetaold > this%dis%bot(n)) then
          ! new and old zeta above bottom
          this%hcof(n) = rho2
          this%rhs(n) = rho2 * hold(n)
        else if (zetanew > this%dis%bot(n) .and. zetaold < this%dis%bot(n)) then
          ! zetanew above bottom but zetaold is not
          this%hcof(n) = rho2
          this%rhs(n) = -rho2 * this%dis%bot(n) / this%alphaf
        else if (zetanew < this%dis%bot(n) .and. zetaold > this%dis%bot(n)) then
          ! zetanew is below bottom, zetaold above bottom
          this%hcof(n) = DZERO
          this%rhs(n) = rho2 * (this%dis%bot(n) / this%alphaf + hold(n))
        end if
      end if
      !
      ! cdl -- end if
      !
    end do
    !
    ! return
    return
  end subroutine swi_fc_storage

  subroutine swi_fn(this, kiter, matrix_sln, idxglo, rhs, hnew, hold, npf, sto)
    ! dummy
    class(GwfSwiType) :: this
    integer(I4B) :: kiter
    class(MatrixBaseType), pointer :: matrix_sln
    integer(I4B), intent(in), dimension(:) :: idxglo
    real(DP), intent(inout), dimension(:) :: rhs
    real(DP), intent(inout), dimension(:) :: hnew
    real(DP), intent(inout), dimension(:) :: hold
    type(GwfNpfType) :: npf
    type(GwfStoType) :: sto
    ! local variables
    real(DP) :: dssdh
    !
    dssdh = -this%alphaf ! derivative of saltwater saturation with respect to freshwater head
    if (this%isaltwater == 1) dssdh = this%alphas ! derivative of saltwater saturation with respect to saltwater head
    !
    call npf_fn_swi(npf, kiter, matrix_sln, idxglo, rhs, hnew, &
                    this%zeta, dssdh, this%isaltwater)
    !
    ! Storage Newton terms are now assembled by the SWI STO formulation
    ! (swisto_fn), dispatched from sto_fn for the flagged SWI_STORAGE cells.

  end subroutine swi_fn

  !> @ brief Calculate fresh storage terms
  !<
  subroutine swi_fn_storage(this, kiter, hold, hnew, matrix_sln, idxglo, &
                            rhs, sto, dssdh, isaltwater)
    ! modules
    use TdisModule, only: delt
    use GwfStorageUtilsModule, only: SsCapacity, SyCapacity, SsTerms, SyTerms
    !
    class(GwfStoType) :: sto
    ! dummy variables
    class(GwfSwiType) :: this !< GwfSwiType object
    integer(I4B), intent(in) :: kiter !< outer iteration number
    real(DP), intent(in), dimension(:) :: hold !< previous heads
    real(DP), intent(in), dimension(:) :: hnew !< current heads
    class(MatrixBaseType), pointer :: matrix_sln !< A matrix
    integer(I4B), intent(in), dimension(:) :: idxglo !< global index model to solution
    real(DP), intent(inout), dimension(:) :: rhs !< right-hand side
    integer(I4B), intent(in) :: isaltwater !< index for saltwater equation
    real(DP) :: dssdh
    ! local variables
    integer(I4B) :: n
    integer(I4B) :: idiag
    real(DP) :: tled
    real(DP) :: sc1
    real(DP) :: sc2
    real(DP) :: rho1
    real(DP) :: rho2
    real(DP) :: tp
    real(DP) :: bt
    real(DP) :: tthk
    real(DP) :: h
    real(DP) :: snnew
    real(DP) :: derv
    real(DP) :: rterm
    real(DP) :: drterm
    real(DP) :: zetanew
    real(DP) :: zetaold
    ! real(DP) :: sew,dereps
    !
    ! set variables
    tled = DONE / delt
    !
    ! loop through and calculate storage contribution to hcof and rhs
    do n = 1, this%dis%nodes
      idiag = this%dis%con%ia(n)
      if (this%ibound(n) <= 0) cycle
      !
      ! calculate zetanew and zetaold
      zetanew = this%get_zetanew(n)
      zetaold = this%get_zetaold(n)
      !
      ! aquifer elevations and thickness
      tp = this%dis%top(n)
      bt = this%dis%bot(n)
      tthk = tp - bt
      h = hnew(n)
      !
      ! aquifer saturation
      snnew = sQuadraticSaturation(tp, bt, zetanew)
      !
      ! storage coefficients
      sc1 = SsCapacity(sto%istor_coef, tp, bt, this%dis%area(n), sto%ss(n))
      sc2 = SyCapacity(this%dis%area(n), this%sy(n))
      rho1 = sc1 * tled
      rho2 = sc2 * tled
      !
      ! calculate newton terms for specific storage
      !    and specific yield
      ! cdl -- if (sto%iconvert(n) /= 0) then
      ! ----------------------------------------------------
      ! calculate saturation derivative as dS/dzeta * dzeta/dh_fresh
      !    derv = sQuadraticSaturationDerivative(tp, bt, zetanew)
      !    derv = derv * dssdh
      ! -----------------------------------------------------
      ! calculate saturation derivative directly as dS / dh
      !    sew = sQuadraticSaturation(tp, bt, zetanew, sto%satomega)
      !    dereps = 1e-6
      !    if (isaltwater == 0) then
      !      zetanew = this%get_zetanew(n, eps_fresh=dereps)
      !    else
      !      zetanew = this%get_zetanew(n, eps_salt=dereps)
      !    endif
      !    derv = sQuadraticSaturation(tp, bt, zetanew, sto%satomega)
      !    derv = (derv-sew)/dereps
      ! ----------------------------------------------------
      ! calculate saturation derivative for saltwater as rhof/(rhos-rhof)/TOTTHICK and freshwater as rhos/(rhos-rhof)/TOTTHICK
      derv = dssdh / tthk
      !      if(isaltwater.eq.0) derv = -0.5
      !      if(isaltwater.eq.1) derv = 0.5125
      ! ----------------------------------------------------
      !
      ! newton terms for specific storage
      if (sto%iconf_ss == 0) then
        if (sto%iorig_ss == 0) then
          drterm = -rho1 * derv * (h - bt) + rho1 * tthk * snnew * derv
        else
          drterm = -(rho1 * derv * h)
        end if
        !sp**         call matrix_sln%add_value_pos(idxglo(idiag), drterm)
        !sp**         rhs(n) = rhs(n) + drterm * h
      end if
      !
      ! newton terms for specific yield
      !    only calculated if the current saturation
      !    is less than one
      if (snnew < DONE) then
        ! calculate newton terms for specific yield
        if (snnew > DZERO) then
          rterm = -rho2 * tthk * snnew
          drterm = -rho2 * tthk * derv
          ! swi correction on freshwater flips sign so for saltwater equation flip sign back
          if (isaltwater == 1) then
            rterm = -rterm
            drterm = -drterm
            rho2 = -rho2
          end if
          ! subtract saltwater part from total flow terms
          call matrix_sln%add_value_pos(idxglo(idiag), -drterm - rho2)
          !     rhs(n) = rhs(n) -(- rterm + drterm * hnew(n) + rho2 * bt)  !csp**** check with single equation
          rhs(n) = rhs(n) - (-rterm + drterm * hnew(n) + rho2 * bt)
        end if
      end if
      ! cdl -- end if
    end do
    !
    ! return
    return
  end subroutine swi_fn_storage
  !
  !

  !-----------------------------------------------
  !
  !> @brief convergence check
  !<
  subroutine swi_cc(this)
    ! dummy
    class(GwfSwiType) :: this
    ! local

    ! recalculate zeta
    call this%update_zeta()

  end subroutine swi_cc

! CSP***--------------------------------------------------------------------------------
  !> @ brief Calculate flows for package
  !!
  !!  Flow calculation for the STO package components. Components include
  !!  specific storage and specific yield storage.
  !!
  !<
! subroutine swi_cq(this, flowja, hnew, hold)
!   ! modules
!   ! dummy variables
!   class(GwfSwiType) :: this !< GwfStoType object
!   real(DP), dimension(:), contiguous, intent(inout) :: flowja !< connection flows
!   real(DP), dimension(:), contiguous, intent(in) :: hnew !< current head
!   real(DP), dimension(:), contiguous, intent(in) :: hold !< previous head
!   ! local variables
!   integer(I4B) :: n
!   integer(I4B) :: idiag
!   real(DP) :: rate
!   !
!   ! initialize strg arrays
!   do n = 1, this%dis%nodes
!     this%storage(n) = DZERO
!   end do
!   !
!   ! Loop through cells
!   do n = 1, this%dis%nodes
!     !
!     ! Calculate change in freshwater storage
!     rate = this%hcof(n) * hnew(n) - this%rhs(n)
!     this%storage(n) = rate
!     !
!     ! Add storage term to flowja
!     idiag = this%dis%con%ia(n)
!     flowja(idiag) = flowja(idiag) + rate
!   end do
!   !
!   ! return
!   return
! end subroutine swi_cq
! CSP***--------------------------------------------------------------------------------

  !> @ brief Calculate flows and storages for SWI package and adjust flowja
  !!
  !<
  subroutine swi_cq(this, hnew, hold, flowja, npf, sto)
    ! modules
    !
    ! dummy variables
    class(GwfSwiType) :: this !< GwfSwiType object
    real(DP), intent(in), dimension(:) :: hnew
    real(DP), dimension(:), contiguous, intent(in) :: hold !< previous head
    real(DP), dimension(:), contiguous, intent(inout) :: flowja !< connection flows
    type(GwfNpfType), intent(in) :: npf
    type(GwfStoType), intent(in) :: sto
    ! local variables
    integer(I4B) :: n, ii, m, ictn, ictm
    real(DP) :: qnm
    !
    ! todo: need to issue error if xt3d is active
    ! if (npf%ixt3d /= 0) then
    !   call npf%xt3d%xt3d_fc(kiter, matrix_sln, idxglo, rhs, hnew)
    ! else

    ! loop over nodes and connections and call swi_qcalc to subtract qnm from flowja
    ictn = 1
    ictm = 1
    do n = 1, npf%dis%nodes
      do ii = npf%dis%con%ia(n) + 1, npf%dis%con%ia(n + 1) - 1
        if (npf%dis%con%mask(ii) == 0) cycle

        ! Calculate terms only for upper triangle but insert into
        ! upper and lower parts of amat.
        m = npf%dis%con%ja(ii)
        if (m < n) cycle
        !
        call swi_qcalc(npf, n, m, ii, ictn, ictm, hnew(n), hnew(m), qnm, &
                       this%zeta)
        !
        ! For saltwater equation, flip sign on rate to just add it
        if (this%isaltwater == 1) qnm = -qnm
        !
        !       change sign to subtract out qnm from salt side
        flowja(ii) = flowja(ii) - qnm
        flowja(this%dis%con%isym(ii)) = flowja(this%dis%con%isym(ii)) + qnm
      end do
    end do
    ! endif ! xt3d if-check
    !
    ! Storage flows are now assembled by the SWI STO formulation (swisto_cq),
    ! dispatched from sto_cq for the flagged SWI_STORAGE cells; the SWI storage
    ! budget/cbc term is reported via swisto_bd / swisto_save_flows.
    !
    ! return
    return
  end subroutine swi_cq

  !> @ brief Calculate flows for SWI package and adjust flowja
  !!
  !<
  subroutine swi_qcalc(npf, n, m, ii, ictn, ictm, hn, hm, qnm, zeta)
    ! modules
    use ConstantsModule, only: DONE
    ! dummy
    class(GwfNpfType) :: npf
    integer(I4B) :: n, m, ii, ihc, ictn, ictm
    real(DP), intent(in) :: hn
    real(DP), intent(in) :: hm
    real(DP), intent(out) :: qnm
    real(DP), intent(in), dimension(:) :: zeta
    ! local
    real(DP) :: hyn, hym
    real(DP) :: cond
    real(DP) :: satn, satm

    ! --Calculate freshwater flow between nodes n and m on the saltwater side
    ! For SWI correction is only in horizontal direction
    ihc = npf%dis%con%ihc(npf%dis%con%jas(ii))
    if (ihc == C3D_VERTICAL) return

    ! Use NPF to get the effective hydraulic conductivity
    hyn = npf%hy_eff(n, m, ihc, ipos=ii)
    hym = npf%hy_eff(m, n, ihc, ipos=ii)

    ! Horizontal conductance
    ! calculate saturation based on zeta, so that hcond is for the
    ! region from zeta down to bottom; hnew is passed in so that
    ! upstream is based on head and not zeta
    call swi_thksat(n, npf%dis%top(n), npf%dis%bot(n), &
                    zeta(n), satn, npf%inewton)
    call swi_thksat(m, npf%dis%top(m), npf%dis%bot(m), &
                    zeta(m), satm, npf%inewton)
    cond = hcond(npf%ibound(n), npf%ibound(m), &
                 ictn, ictm, &
                 npf%inewton, &
                 npf%dis%con%ihc(npf%dis%con%jas(ii)), &
                 npf%icellavg, &
                 npf%condsat(npf%dis%con%jas(ii)), &
                 hn, hm, &
                 satn, satm, &
                 hyn, hym, &
                 npf%dis%top(n), npf%dis%top(m), &
                 npf%dis%bot(n), npf%dis%bot(m), &
                 npf%dis%con%cl1(npf%dis%con%jas(ii)), &
                 npf%dis%con%cl2(npf%dis%con%jas(ii)), &
                 npf%dis%con%hwva(npf%dis%con%jas(ii)))
    !
    ! Calculate flow positive into cell n
    qnm = cond * (hm - hn)
    !
    ! return
    return
  end subroutine swi_qcalc

  !> @ brief Model budget calculation for package
  !!
  !!  Budget calculation for the STO package components. Components include
  !!  specific storage and specific yield storage.
  !!
  !<
  subroutine swi_bd(this, isuppress_output, model_budget)
    ! modules
    use TdisModule, only: delt
    use BudgetModule, only: BudgetType, rate_accumulator
    ! dummy variables
    class(GwfSwiType) :: this !< GwfSwiType object
    integer(I4B), intent(in) :: isuppress_output !< flag to suppress model output
    type(BudgetType), intent(inout) :: model_budget !< model budget object
    ! local variables
    real(DP) :: rin
    real(DP) :: rout
    !
    ! Add swi storage rates to model budget
    call rate_accumulator(this%storage, rin, rout)
    call model_budget%addentry(rin, rout, delt, budtxt(1), &
                               isuppress_output, '     SWI')
    !
    ! return
    return
  end subroutine swi_bd

  !> @ brief Save model flows for package
  !!
  !!  Save cell-by-cell budget terms for the STO package.
  !!
  !<
  subroutine swi_save_model_flows(this, icbcfl, icbcun)
    ! dummy variables
    class(GwfSwiType) :: this !< GwfSwiType object
    integer(I4B), intent(in) :: icbcfl !< flag to output budget data
    integer(I4B), intent(in) :: icbcun !< cell-by-cell file unit number
    ! local variables
    integer(I4B) :: ibinun
    integer(I4B) :: iprint, nvaluesp, nwidthp
    character(len=1) :: cdatafmp = ' ', editdesc = ' '
    real(DP) :: dinact
    !
    ! Set unit number for binary output
    if (this%ipakcb < 0) then
      ibinun = icbcun
    elseif (this%ipakcb == 0) then
      ibinun = 0
    else
      ibinun = this%ipakcb
    end if
    if (icbcfl == 0) ibinun = 0
    !
    ! Record the storage rates if requested
    if (ibinun /= 0) then
      iprint = 0
      dinact = DZERO
      !
      ! swi storage
      call this%dis%record_array(this%storage, this%iout, iprint, -ibinun, &
                                 budtxt(1), cdatafmp, nvaluesp, &
                                 nwidthp, editdesc, dinact)
    end if
    !
    ! return
    return
  end subroutine swi_save_model_flows

  !> @brief Save density array to binary file
  !<
  subroutine swi_ot_dv(this, idvfl)
    ! dummy
    class(GwfSwiType) :: this
    integer(I4B), intent(in) :: idvfl
    ! local
    character(len=1) :: cdatafmp = ' ', editdesc = ' '
    integer(I4B) :: ibinun
    integer(I4B) :: iprint
    integer(I4B) :: nvaluesp
    integer(I4B) :: nwidthp
    real(DP) :: dinact
    !
    ! Set unit number for density output
    if (this%izetaout /= 0) then
      ibinun = 1
    else
      ibinun = 0
    end if
    if (idvfl == 0) ibinun = 0
    !
    ! save density array
    if (ibinun /= 0) then
      iprint = 0
      dinact = DHNOFLO
      !
      ! write density to binary file
      if (this%izetaout /= 0) then
        ibinun = this%izetaout
        call this%dis%record_array(this%zeta, this%iout, iprint, ibinun, &
                                   '            ZETA', cdatafmp, nvaluesp, &
                                   nwidthp, editdesc, dinact)
      end if
    end if
    !
    ! Return
    return
  end subroutine swi_ot_dv

  !> @brief Deallocate
  !<
  subroutine swi_da(this)
    ! modules
    use MemoryManagerModule, only: mem_deallocate
    use MemoryManagerExtModule, only: memorystore_remove
    use SimVariablesModule, only: idm_context
    ! dummy
    class(GwfSwiType) :: this

    ! deallocate IDM memory
    call memorystore_remove(this%name_model, 'SWI', idm_context)

    if (this%intva > 0) then
      call this%tva_reader%tva_da()
    end if

    ! deallocate arrays
    call mem_deallocate(this%zeta)
    call mem_deallocate(this%hcof)
    call mem_deallocate(this%rhs)
    call mem_deallocate(this%storage)

    if (this%iconfiguration == FRESHWATER_ONE_FLUID) then
      call mem_deallocate(this%hsalt)
      call mem_deallocate(this%hsaltold)
    end if

    ! deallocate scalars
    call mem_deallocate(this%iconfiguration)
    call mem_deallocate(this%isaltwater)
    call mem_deallocate(this%izetaout)
    call mem_deallocate(this%intva)
    call mem_deallocate(this%alphaf)
    call mem_deallocate(this%alphas)
    call mem_deallocate(this%hsalt_user)

    ! deallocate parent
    call this%NumericalPackageType%da()

    ! nullify pointers
    nullify (this%hfresh)
    nullify (this%hsalt)
    nullify (this%hfreshold)
    nullify (this%hsaltold)

  end subroutine swi_da

  !> @brief Load data from IDM into package
  !<
  subroutine swi_load(this)
    ! modules
    ! dummy
    class(GwfSwiType) :: this

    call this%source_options()
    call this%source_griddata()

  end subroutine swi_load

  !> @ brief Allocate scalars
  !!
  !! Allocate and initialize scalars for the VSC package. The base model
  !! allocate scalars method is also called.
  !<
  subroutine allocate_scalars(this)
    ! modules
    ! dummy
    class(GwfSwiType) :: this

    ! allocate scalars in NumericalPackageType
    call this%NumericalPackageType%allocate_scalars()

    ! Allocate scalars
    call mem_allocate(this%iconfiguration, 'ICONFIGURATION', this%memoryPath)
    call mem_allocate(this%isaltwater, 'ISALTWATER', this%memoryPath)
    call mem_allocate(this%izetaout, 'IZETAOUT', this%memoryPath)
    call mem_allocate(this%intva, 'INTVA', this%memoryPath)
    call mem_allocate(this%alphaf, 'ALPHAF', this%memoryPath)
    call mem_allocate(this%alphas, 'ALPHAS', this%memoryPath)
    call mem_allocate(this%hsalt_user, 'HSALT_USER', this%memoryPath)

    ! Initialize value
    this%izetaout = 0
    this%intva = 0
    this%alphaf = 40.D0
    this%alphas = 41.D0
    this%hsalt_user = DZERO

    ! initialize iconfiguration and isaltwater
    call this%set_configuration("FRESHWATER_ONE_FLUID")

  end subroutine allocate_scalars

  !> @brief Allocate arrays
  !<
  subroutine allocate_arrays(this, nodes)
    ! modules
    use MemoryManagerModule, only: mem_allocate
    ! dummy
    class(GwfSwiType) :: this
    integer(I4B), intent(in) :: nodes
    ! local
    integer(I4B) :: n

    ! Allocate
    call mem_allocate(this%zeta, nodes, 'ZETA', this%memoryPath)
    call mem_allocate(this%hcof, nodes, 'HCOF', this%memoryPath)
    call mem_allocate(this%rhs, nodes, 'RHS', this%memoryPath)
    call mem_allocate(this%storage, nodes, 'STORAGE', this%memoryPath)

    ! initialize
    do n = 1, nodes
      this%zeta(n) = DZERO
      this%hcof(n) = DZERO
      this%rhs(n) = DZERO
      this%storage(n) = DZERO
    end do

  end subroutine allocate_arrays

  !> @brief Update simulation options from input mempath
  !<
  subroutine source_options(this)
    ! modules
    use SourceCommonModule, only: filein_fname
    use OpenSpecModule, only: access, form
    use InputOutputModule, only: getunit, openfile
    use MemoryManagerModule, only: mem_setptr, get_isize
    use MemoryManagerExtModule, only: mem_set_value
    use CharacterStringModule, only: CharacterStringType
    use GwfSwiInputModule, only: GwfSwiParamFoundType
    ! dummy
    class(GwfSwiType) :: this
    character(len=LINELENGTH) :: zeta_fname
    character(len=LINELENGTH) :: tva_fname
    ! locals
    type(GwfSwiParamFoundType) :: found

    ! update defaults with idm sourced values
    zeta_fname = ''
    call mem_set_value(zeta_fname, 'ZETAFILE', this%input_mempath, &
                       found%zetafile)
    call mem_set_value(this%hsalt_user, 'HSALT_USER', this%input_mempath, &
                       found%hsalt_user)

    ! open zeta file
    if (zeta_fname /= '') then
      this%izetaout = getunit()
      call openfile(this%izetaout, this%iout, zeta_fname, 'DATA(BINARY)', &
                    form, access, 'REPLACE', MNORMAL)
      write (this%iout, '(4x,a)') &
        'ZETA information will be written to ', trim(zeta_fname)
    end if

    ! open hsalt file
    tva_fname = ''
    if (filein_fname(tva_fname, 'TVA6_FILENAME', this%input_mempath, &
                     this%input_fname)) then
      this%intva = getunit()
      call openfile(this%intva, this%iout, tva_fname, 'TVA6')
      write (this%iout, '(4x,a)') &
        'Time variable array information will be read from ', trim(tva_fname)
    end if

    ! log options
    if (this%iout > 0) then
      call this%log_options(found)
    end if

  end subroutine source_options

  !> @brief Log options sourced from the input mempath
  !<
  subroutine log_options(this, found)
    ! modules
    use KindModule, only: LGP
    use GwfSwiInputModule, only: GwfSwiParamFoundType
    ! dummy
    class(GwfSwiType) :: this
    ! locals
    type(GwfSwiParamFoundType), intent(in) :: found

    write (this%iout, '(1x,a)') 'Setting SWI Options'
    if (found%zetafile) then
      write (this%iout, '(4x,a)') &
        'Zeta will be written to a binary output file.'
    end if
    if (found%hsalt_user) then
      write (this%iout, '(4x,a,G0,a)') &
        'Saltwater head was set to ', this%hsalt_user, '.'
    end if
    if (found%tva6_filename) then
      write (this%iout, '(4x,a)') &
        'Saltwater head will be read from file.'
    end if
    write (this%iout, '(1x,a,/)') 'End Setting SWI Options'

  end subroutine log_options

  !> @brief Copy grid data from IDM into package
  !<
  subroutine source_griddata(this)
    ! modules
    use MemoryManagerExtModule, only: mem_set_value
    use GwfSwiInputModule, only: GwfSwiParamFoundType
    ! dummy
    class(GwfSwiType) :: this
    ! local
    ! type(GwfSwiParamFoundType) :: found
    ! integer(I4B), dimension(:), pointer, contiguous :: map

  end subroutine source_griddata

  !> @brief Set this model configuration.
  !! This can be called from an exchange, such as the
  !! SwiSwi exchange, which knows which model/swi package
  !! is the saltwater model.
  !<
  subroutine set_configuration(this, config)
    ! dummy
    class(GwfSwiType) :: this
    character(len=*), intent(in) :: config

    select case (config)
    case ("FRESHWATER_ONE_FLUID")
      this%iconfiguration = FRESHWATER_ONE_FLUID
      this%isaltwater = 0
    case ("SALTWATER_ONE_FLUID")
      this%iconfiguration = SALTWATER_ONE_FLUID
      this%isaltwater = 1
    case ("FRESHWATER_TWO_FLUID")
      this%iconfiguration = FRESHWATER_TWO_FLUID
      this%isaltwater = 0
    case ("SALTWATER_TWO_FLUID")
      this%iconfiguration = SALTWATER_TWO_FLUID
      this%isaltwater = 1
    end select
  end subroutine set_configuration

  !> @brief Set head pointers
  !! This will likely be called from the SWI-SWI exchange
  !! which has access to both fresh and salt models.
  !<
  subroutine set_head_pointers(this, hfresh, hsalt, hfreshold, hsaltold)
    class(GwfSwiType) :: this
    real(DP), dimension(:), target, intent(in) :: hfresh
    real(DP), dimension(:), target, intent(in) :: hsalt
    real(DP), dimension(:), target, intent(in) :: hfreshold
    real(DP), dimension(:), target, intent(in) :: hsaltold
    this%hfresh => hfresh
    this%hsalt => hsalt
    this%hfreshold => hfreshold
    this%hsaltold => hsaltold
  end subroutine set_head_pointers

  !> @brief Get new value for zeta
  !!
  !! Include hsalt contribution if available and use perturbations if
  !! provided.
  !!
  !<
  function get_zetanew(this, n, eps_fresh, eps_salt) result(zetanew)
    ! modules
    ! dummy
    class(GwfSwiType) :: this !< this instance
    integer(I4B) :: n !< cell number
    real(DP), optional, intent(in) :: eps_fresh !< perturbation for fresh head
    real(DP), optional, intent(in) :: eps_salt !< perturbation for salt head
    ! return
    real(DP) :: zetanew
    ! locals
    real(DP) :: eps_f
    real(DP) :: eps_s
    if (present(eps_fresh)) then
      eps_f = eps_fresh
    else
      eps_f = DZERO
    end if
    if (present(eps_salt)) then
      eps_s = eps_salt
    else
      eps_s = DZERO
    end if
    if (this%iconfiguration == FRESHWATER_TWO_FLUID .or. &
        this%iconfiguration == SALTWATER_TWO_FLUID) then
      ! two-model swi with fresh and salt
      zetanew = calc_zeta(this%alphaf, &
                          this%hfresh(n) + eps_f, &
                          this%alphas, &
                          this%hsalt(n) + eps_s)

      ! cdl HACK -- keep minimum saltwater thickness
      ! if (zetanew <= this%dis%bot(n)) then
      !   zetanew = this%dis%bot(n) + 1.d-3 ! cdl DEM12
      ! end if

    else
      ! freshwater only simulation
      zetanew = calc_zeta(this%alphaf, this%hfresh(n) + eps_f, &
                          this%alphas, this%hsalt(n))
    end if
  end function get_zetanew

  !> @brief Get old value for zeta
  !!
  !! Include hsalt contribution if available and use perturbations if
  !! provided, though may be able to eliminate perturbations as they
  !! are probably not needed.
  !!
  !<
  function get_zetaold(this, n, eps_fresh, eps_salt) result(zetaold)
    ! modules
    ! dummy
    class(GwfSwiType) :: this !< this instance
    integer(I4B) :: n !< cell number
    real(DP), optional, intent(in) :: eps_fresh !< perturbation for fresh head
    real(DP), optional, intent(in) :: eps_salt !< perturbation for salt head
    ! return
    real(DP) :: zetaold
    ! locals
    real(DP) :: eps_f
    real(DP) :: eps_s
    if (present(eps_fresh)) then
      eps_f = eps_fresh
    else
      eps_f = DZERO
    end if
    if (present(eps_salt)) then
      eps_s = eps_salt
    else
      eps_s = DZERO
    end if
    if (this%iconfiguration == FRESHWATER_TWO_FLUID .or. &
        this%iconfiguration == SALTWATER_TWO_FLUID) then
      ! two-model swi with fresh and salt
      zetaold = calc_zeta(this%alphaf, &
                          this%hfreshold(n) + eps_f, &
                          this%alphas, &
                          this%hsaltold(n) + eps_s)
    else
      ! freshwater only simulation
      zetaold = calc_zeta(this%alphaf, this%hfreshold(n) + eps_f, &
                          this%alphas, this%hsaltold(n))
    end if
  end function get_zetaold

  !> @brief Calculate zeta surface
  !<
  subroutine update_zeta(this)
    ! modules
    ! dummy
    class(GwfSwiType) :: this
    ! locals
    integer(I4B) :: n

    ! Loop through each node and calculate zeta
    do n = 1, this%dis%nodes
      ! skip if inactive
      if (this%ibound(n) == 0) cycle

      ! Calculate zeta
      this%zeta(n) = this%get_zetanew(n)

      ! todo: Do we want to constrain zeta to top and
      !    bot of cell?
      if (this%zeta(n) > this%dis%top(n)) then
        this%zeta(n) = this%dis%top(n)
      end if
      if (this%zeta(n) < this%dis%bot(n)) then
        this%zeta(n) = this%dis%bot(n)
      end if

    end do

  end subroutine update_zeta

  function calc_zeta(alphaf, hf, alphas, hs) result(zeta)
    real(DP), intent(in) :: alphaf
    real(DP), intent(in) :: hf
    real(DP), intent(in) :: alphas
    real(DP), intent(in) :: hs
    real(DP) :: zeta
    zeta = -alphaf * hf + alphas * hs
  end function calc_zeta

  !> @brief Add swi correction term
  !!
  !! This is a retooling of the npf_fc to subtract the saltwater flow
  !> that would occur below zeta.
  subroutine npf_fc_swi(npf, kiter, matrix_sln, idxglo, rhs, hnew, zeta, &
                        isaltwater)
    ! modules
    use ConstantsModule, only: DONE
    ! dummy
    class(GwfNpfType) :: npf
    integer(I4B) :: kiter
    class(MatrixBaseType), pointer :: matrix_sln
    integer(I4B), intent(in), dimension(:) :: idxglo
    real(DP), intent(inout), dimension(:) :: rhs
    real(DP), intent(in), dimension(:) :: hnew
    real(DP), intent(in), dimension(:) :: zeta
    ! local
    integer(I4B) :: n, m, ii, idiag, ihc, isaltwater
    integer(I4B) :: isymcon, idiagm
    real(DP) :: hyn, hym
    real(DP) :: cond
    real(DP) :: satn, satm
    integer(I4B) :: ictn, ictm
    real(DP) :: bs, bf
    !
    ! Calculate conductance and put into amat
    !
    ! todo: need to issue error if xt3d is active
    ! if (npf%ixt3d /= 0) then
    !   call npf%xt3d%xt3d_fc(kiter, matrix_sln, idxglo, rhs, hnew)
    ! else

    ! Set the celltype to be 1 so that zeta is used as the top
    ! surface for evaluation of the saturated thickness in the
    ! saltwater zone
    ictn = 1
    ictm = 1
!
    do n = 1, npf%dis%nodes
      do ii = npf%dis%con%ia(n) + 1, npf%dis%con%ia(n + 1) - 1
        if (npf%dis%con%mask(ii) == 0) cycle

        ! Calculate terms only for upper triangle but insert into
        ! upper and lower parts of amat.
        m = npf%dis%con%ja(ii)
        if (m < n) cycle

        ! For SWI correction is only in horizontal direction
        ihc = npf%dis%con%ihc(npf%dis%con%jas(ii))
        if (ihc == C3D_VERTICAL) cycle

        ! Use NPF to get the effective hydraulic conductivity
        hyn = npf%hy_eff(n, m, ihc, ipos=ii)
        hym = npf%hy_eff(m, n, ihc, ipos=ii)

        ! Horizontal conductance
        ! calculate saturation based on zeta, so that hcond is for the
        ! region from zeta down to bottom; hnew is passed in so that
        ! upstream is based on head and not zeta
        call swi_thksat(n, npf%dis%top(n), npf%dis%bot(n), &
                        zeta(n), satn, npf%inewton)

        bs = satn * (npf%dis%top(n) - npf%dis%bot(n))
        bs = max(bs, 1.d-3)
        if (bs <= 1.d-3) then
          bs = 1.d-3
          satn = bs / (npf%dis%top(n) - npf%dis%bot(n))
        end if
        bf = (DONE - satn) * (npf%dis%top(n) - npf%dis%bot(n))
        if (bf <= 1.d-3) then
          bf = 1.d-3
          satn = 1.d0 - bf / (npf%dis%top(n) - npf%dis%bot(n))
        end if

        call swi_thksat(m, npf%dis%top(m), npf%dis%bot(m), &
                        zeta(m), satm, npf%inewton)

        bs = satm * (npf%dis%top(m) - npf%dis%bot(m))
        if (bs <= 1.d-3) then
          bs = 1.d-3
          satm = bs / (npf%dis%top(m) - npf%dis%bot(m))
        end if
        bf = (DONE - satm) * (npf%dis%top(m) - npf%dis%bot(m))
        if (bf <= 1.d-3) then
          bf = 1.d-3
          satm = 1.d0 - bf / (npf%dis%top(m) - npf%dis%bot(m))
        end if

        cond = hcond(npf%ibound(n), npf%ibound(m), &
                     ictn, ictm, &
                     npf%inewton, &
                     npf%dis%con%ihc(npf%dis%con%jas(ii)), &
                     npf%icellavg, &
                     npf%condsat(npf%dis%con%jas(ii)), &
                     hnew(n), hnew(m), &
                     satn, satm, &
                     hyn, hym, &
                     npf%dis%top(n), npf%dis%top(m), &
                     npf%dis%bot(n), npf%dis%bot(m), &
                     npf%dis%con%cl1(npf%dis%con%jas(ii)), &
                     npf%dis%con%cl2(npf%dis%con%jas(ii)), &
                     npf%dis%con%hwva(npf%dis%con%jas(ii)))

        ! Fill row n, Note signs are flipped in order to
        ! subtract the flow from the saltwater zone for freshwater
        !
        ! THEREFORE, for saltwater, flip sign on cond
        if (isaltwater == 1) cond = -cond
        !
        idiag = npf%dis%con%ia(n)
        call matrix_sln%add_value_pos(idxglo(ii), -cond)
        call matrix_sln%add_value_pos(idxglo(idiag), cond)

        ! Fill row m, Note signs are flipped in order to
        ! subtract the flow from the saltwater zone
        isymcon = npf%dis%con%isym(ii)
        idiagm = npf%dis%con%ia(m)
        call matrix_sln%add_value_pos(idxglo(isymcon), -cond)
        call matrix_sln%add_value_pos(idxglo(idiagm), cond)
      end do
    end do
  end subroutine npf_fc_swi

  ! add vertical connections for saltwater model
  subroutine npf_fc_vflow(npf, kiter, matrix_sln, idxglo, rhs, hnew)
    ! -- modules
    use ConstantsModule, only: DONE
    ! -- dummy
    class(GwfNpfType) :: npf
    integer(I4B) :: kiter
    class(MatrixBaseType), pointer :: matrix_sln
    integer(I4B), intent(in), dimension(:) :: idxglo
    real(DP), intent(inout), dimension(:) :: rhs
    real(DP), intent(in), dimension(:) :: hnew
    ! -- local
    integer(I4B) :: n, m, ii, idiag, ihc
    integer(I4B) :: isymcon, idiagm
    real(DP) :: hyn, hym
    real(DP) :: cond
    ! real(DP) :: satn
    ! real(DP) :: satm
    !
    ! -- Calculate conductance and put into amat
    !
    if (npf%ixt3d /= 0) then
      ! call npf%xt3d%xt3d_fc(kiter, matrix_sln, idxglo, rhs, hnew)
    else
      do n = 1, npf%dis%nodes
        do ii = npf%dis%con%ia(n) + 1, npf%dis%con%ia(n + 1) - 1
          if (npf%dis%con%mask(ii) == 0) cycle

          m = npf%dis%con%ja(ii)
          !
          ! -- Calculate conductance only for upper triangle but insert into
          !    upper and lower parts of amat.
          if (m < n) cycle
          ihc = npf%dis%con%ihc(npf%dis%con%jas(ii))
          hyn = npf%hy_eff(n, m, ihc, ipos=ii)
          hym = npf%hy_eff(m, n, ihc, ipos=ii)
          !
          ! -- Vertical connection
          if (ihc == C3D_VERTICAL) then
            !
            ! -- Calculate vertical conductance
            cond = vcond(npf%ibound(n), npf%ibound(m), &
                         npf%icelltype(n), npf%icelltype(m), npf%inewton, &
                         npf%ivarcv, npf%idewatcv, &
                         npf%condsat(npf%dis%con%jas(ii)), hnew(n), hnew(m), &
                         hyn, hym, &
                         npf%sat(n), npf%sat(m), &
                         npf%dis%top(n), npf%dis%top(m), &
                         npf%dis%bot(n), npf%dis%bot(m), &
                         npf%dis%con%hwva(npf%dis%con%jas(ii)))
            !
            ! -- Vertical flow for perched conditions
            if (npf%iperched /= 0) then
              if (npf%icelltype(m) /= 0) then
                if (hnew(m) < npf%dis%top(m)) then
                  !
                  ! -- Fill row n
                  idiag = npf%dis%con%ia(n)
                  rhs(n) = rhs(n) - cond * npf%dis%bot(n)
                  call matrix_sln%add_value_pos(idxglo(idiag), -cond)
                  !
                  ! -- Fill row m
                  isymcon = npf%dis%con%isym(ii)
                  call matrix_sln%add_value_pos(idxglo(isymcon), cond)
                  rhs(m) = rhs(m) + cond * npf%dis%bot(n)
                  !
                  ! -- cycle the connection loop
                  cycle
                end if
              end if
            end if
            !
          else
            cond = DZERO
            ! satn = this%sat(n)
            ! satm = this%sat(m)
            ! if (this%ihighcellsat /= 0) then
            !   call this%highest_cell_saturation(n, m, &
            !                                     hnew(n), hnew(m), &
            !                                     satn, satm)
            ! end if
            ! !
            ! ! -- Horizontal conductance
            ! cond = hcond(this%ibound(n), this%ibound(m), &
            !              this%icelltype(n), this%icelltype(m), &
            !              this%inewton, &
            !              this%dis%con%ihc(this%dis%con%jas(ii)), &
            !              this%icellavg, &
            !              this%condsat(this%dis%con%jas(ii)), &
            !              hnew(n), hnew(m), satn, satm, hyn, hym, &
            !              this%dis%top(n), this%dis%top(m), &
            !              this%dis%bot(n), this%dis%bot(m), &
            !              this%dis%con%cl1(this%dis%con%jas(ii)), &
            !              this%dis%con%cl2(this%dis%con%jas(ii)), &
            !              this%dis%con%hwva(this%dis%con%jas(ii)))
          end if
          !
          ! -- Fill row n
          idiag = npf%dis%con%ia(n)
          call matrix_sln%add_value_pos(idxglo(ii), cond)
          call matrix_sln%add_value_pos(idxglo(idiag), -cond)
          !
          ! -- Fill row m
          isymcon = npf%dis%con%isym(ii)
          idiagm = npf%dis%con%ia(m)
          call matrix_sln%add_value_pos(idxglo(isymcon), cond)
          call matrix_sln%add_value_pos(idxglo(idiagm), -cond)
        end do
      end do
      !
    end if
  end subroutine npf_fc_vflow

  !> @brief Fill newton terms
  !<
  subroutine npf_fn_swi(npf, kiter, matrix_sln, idxglo, rhs, hnew, zeta, &
                        dzetadh, isaltwater)
    ! dummy
    type(GwfNpfType) :: npf
    integer(I4B) :: kiter
    class(MatrixBaseType), pointer :: matrix_sln
    integer(I4B), intent(in), dimension(:) :: idxglo
    real(DP), intent(inout), dimension(:) :: rhs
    real(DP), intent(inout), dimension(:) :: hnew
    real(DP), intent(inout), dimension(:) :: zeta
    real(DP), intent(in) :: dzetadh
    integer(I4B), intent(in) :: isaltwater !< index for saltwater equation
    ! local
    integer(I4B) :: nodes, nja
    integer(I4B) :: n, m, ii, idiag
    integer(I4B) :: isymcon, idiagm
    integer(I4B) :: iups
    integer(I4B) :: idn
    real(DP) :: cond
    real(DP) :: consterm
    real(DP) :: filledterm
    real(DP) :: derv
    real(DP) :: hds
    real(DP) :: term
    real(DP) :: topup
    real(DP) :: botup
    !
    ! add newton terms to solution matrix
    nodes = npf%dis%nodes
    nja = npf%dis%con%nja
    !
    ! todo: need to issue error if xt3d is active
    !if (npf%ixt3d /= 0) then
    !  call npf%xt3d%xt3d_fn(kiter, nodes, nja, matrix_sln, idxglo, rhs, hnew)
    !else
    !
    do n = 1, nodes
      idiag = npf%dis%con%ia(n)
      do ii = npf%dis%con%ia(n) + 1, npf%dis%con%ia(n + 1) - 1
        if (npf%dis%con%mask(ii) == 0) cycle

        m = npf%dis%con%ja(ii)
        isymcon = npf%dis%con%isym(ii)
        ! work on upper triangle
        if (m < n) cycle
        if (npf%dis%con%ihc(npf%dis%con%jas(ii)) == 0 .and. &
            npf%ivarcv == 0) then
          !call npf%vcond(n,m,hnew(n),hnew(m),ii,cond)
          ! do nothing
        else
          ! determine upstream node
          iups = m
          if (hnew(m) < hnew(n)) iups = n
          idn = n
          if (iups == n) idn = m
          !
          ! no newton terms if upstream cell is confined
          ! for swi, always do newton
          !if (npf%icelltype(iups) == 0) cycle
          !
          ! Set the upstream top and bot, and then recalculate for a
          !    vertically staggered horizontal connection
          topup = npf%dis%top(iups)
          botup = npf%dis%bot(iups)
          if (npf%dis%con%ihc(npf%dis%con%jas(ii)) == 2) then
            topup = min(npf%dis%top(n), npf%dis%top(m))
            botup = max(npf%dis%bot(n), npf%dis%bot(m))
          end if
          !
          ! get saturated conductivity for derivative
          cond = npf%condsat(npf%dis%con%jas(ii))
          !
          ! compute additional term
          consterm = -cond * (hnew(iups) - hnew(idn)) !needs to use hwadi instead of hnew(idn)
          !filledterm = cond
          filledterm = matrix_sln%get_value_pos(idxglo(ii))
          ! use zeta in derivative
          derv = sQuadraticSaturationDerivative(topup, botup, zeta(iups), &
                                                npf%satomega)
          !   derv = 1.25e-2 ! will always be this on linear part where saturation is not near 0 or 1 (away from smoothing)
          derv = derv * dzetadh
          idiagm = npf%dis%con%ia(m)
          ! fill jacobian for n being the upstream node
          if (iups == n) then
            hds = hnew(m)
            !isymcon =  npf%dis%con%isym(ii)
            term = consterm * derv
            ! swi correction on freshwater flips sign so for saltwater equation flip sign back
            if (isaltwater == 1) term = -term
            ! flip signs for swi correction
            rhs(n) = rhs(n) - term * hnew(n) !+ amat(idxglo(isymcon)) * (dwadi * hds - hds) !need to add dwadi
            rhs(m) = rhs(m) + term * hnew(n) !- amat(idxglo(isymcon)) * (dwadi * hds - hds) !need to add dwadi
            ! fill in row of n
            ! flip sign for swi correction
            call matrix_sln%add_value_pos(idxglo(idiag), -term)
            ! fill newton term in off diagonal if active cell
            if (npf%ibound(n) > 0) then
              filledterm = matrix_sln%get_value_pos(idxglo(ii))
              call matrix_sln%set_value_pos(idxglo(ii), filledterm) !* dwadi !need to add dwadi
            end if
            !fill row of m
            filledterm = matrix_sln%get_value_pos(idxglo(idiagm))
            call matrix_sln%set_value_pos(idxglo(idiagm), filledterm) !- filledterm * (dwadi - DONE) !need to add dwadi
            ! fill newton term in off diagonal if active cell
            if (npf%ibound(m) > 0) then
              ! flip sign for swi correction
              call matrix_sln%add_value_pos(idxglo(isymcon), term)
            end if
            ! fill jacobian for m being the upstream node
          else
            hds = hnew(n)
            term = -consterm * derv
            ! swi correction on freshwater flips sign so for saltwater equation flip sign back
            if (isaltwater == 1) term = -term
            ! flip sign for swi correction
            rhs(n) = rhs(n) - term * hnew(m) !+ amat(idxglo(ii)) * (dwadi * hds - hds) !need to add dwadi
            rhs(m) = rhs(m) + term * hnew(m) !- amat(idxglo(ii)) * (dwadi * hds - hds) !need to add dwadi
            ! fill in row of n
            filledterm = matrix_sln%get_value_pos(idxglo(idiag))
            call matrix_sln%set_value_pos(idxglo(idiag), filledterm) !- filledterm * (dwadi - DONE) !need to add dwadi
            ! fill newton term in off diagonal if active cell
            if (npf%ibound(n) > 0) then
              ! flip sign for swi correction
              call matrix_sln%add_value_pos(idxglo(ii), -term)
            end if
            !fill row of m
            ! flip sign for swi correction
            call matrix_sln%add_value_pos(idxglo(idiagm), term)
            ! fill newton term in off diagonal if active cell
            if (npf%ibound(m) > 0) then
              filledterm = matrix_sln%get_value_pos(idxglo(isymcon))
              call matrix_sln%set_value_pos(idxglo(isymcon), filledterm) !* dwadi  !need to add dwadi
            end if
          end if
        end if

      end do
    end do
    !
!    end if
    !
    ! Return
    return
  end subroutine npf_fn_swi

  !> @brief Fractional cell freshwater saturation
  !<
  subroutine swi_thksat(n, top, bot, zeta, thksat, inewton)
    ! dummy
    integer(I4B), intent(in) :: n
    real(DP), intent(in) :: top
    real(DP), intent(in) :: bot
    real(DP), intent(in) :: zeta
    real(DP), intent(inout) :: thksat
    integer(I4B), intent(in) :: inewton
    !
    ! Standard Formulation
    if (zeta >= top) then
      thksat = DONE
    else
      thksat = (zeta - bot) / (top - bot)
    end if
    !
    ! Smoothed thickness
    if (inewton /= 0) then
      thksat = sQuadraticSaturation(top, bot, zeta)
    end if
    !
    ! Return
    return
  end subroutine swi_thksat

  !> @brief Calculated the saturated thickness of either the freshwater
  !< or saltwater zone, depending on which model this package is assigned
  !< to, and return the value.
  !<
  function swi_saturated_thickness(this, n, h, ict) result(zone_thickness)
    ! modules
    ! dummy
    class(GwfSwiType) :: this
    ! locals
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: ict
    real(DP), intent(in) :: h
    real(DP) :: b
    real(DP) :: tp
    real(DP) :: zeta
    real(DP) :: zone_thickness

    zeta = this%get_zetanew(n)
    if (this%isaltwater == 1) then
      b = zeta - this%dis%bot(n)
    else
      tp = this%dis%top(n)
      if (ict == 1) then
        tp = min(tp, h)
      end if
      b = tp - max(zeta, this%dis%bot(n))
    end if
    zone_thickness = max(b, ZONE_MINIMUM_THICKNESS)
  end function swi_saturated_thickness

  !> @brief SWI NPF formulation: is this connection handled by SWI?
  !<
  function swinpf_is_active(this, n, m) result(is_active)
    class(SwiNpfFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: m
    logical(LGP) :: is_active
    is_active = .true.
  end function swinpf_is_active

  !> @brief SWI NPF formulation: per-cell precompute (zeta is updated elsewhere)
  !<
  subroutine swinpf_cf(this, kiter, n)
    class(SwiNpfFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: kiter
    integer(I4B), intent(in) :: n
  end subroutine swinpf_cf

  !> @brief SWI NPF formulation: fill the freshwater conductance for connection
  !! n-m, overriding the default NPF term. condf = cond(S^w) - cond(S^s), i.e. the
  !! conductance based on the freshwater column between the interface and the head.
  !<
  subroutine swinpf_fc(this, n, m, ipos, matrix_sln, rhs, idxglo, hnew)
    class(SwiNpfFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: m
    integer(I4B), intent(in) :: ipos
    class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: idxglo
    real(DP), dimension(:), intent(in) :: hnew
    ! local
    type(GwfNpfType), pointer :: npf
    integer(I4B) :: idiag, ihc, isymcon, idiagm, jas, ictn, ictm
    real(DP) :: hyn, hym, cond, condsw, condf
    real(DP) :: satn, satm, satn_s, satm_s, tn, tm, bn, bm, tthkn, tthkm
    !
    npf => this%npf
    jas = npf%dis%con%jas(ipos)
    ihc = npf%dis%con%ihc(jas)
    ! SWI (single-fluid) overrides horizontal flow only
    if (ihc == C3D_VERTICAL) return
    tn = npf%dis%top(n)
    tm = npf%dis%top(m)
    bn = npf%dis%bot(n)
    bm = npf%dis%bot(m)
    tthkn = tn - bn
    tthkm = tm - bm
    hyn = npf%hy_eff(n, m, ihc, ipos=ipos)
    hym = npf%hy_eff(m, n, ihc, ipos=ipos)
    !
    ! -- full water-column (S^w) conductance, as in fc_default_flow
    cond = hcond(npf%ibound(n), npf%ibound(m), &
                 npf%icelltype(n), npf%icelltype(m), npf%inewton, ihc, &
                 npf%icellavg, npf%condsat(jas), hnew(n), hnew(m), &
                 npf%sat(n), npf%sat(m), hyn, hym, tn, tm, bn, bm, &
                 npf%dis%con%cl1(jas), npf%dis%con%cl2(jas), &
                 npf%dis%con%hwva(jas))
    !
    ! -- saltwater-column (S^s) conductance below the interface (zeta), with the
    !    minimum-thickness clamp used by npf_fc_swi; celltype forced to 1 so zeta
    !    acts as the top surface
    ictn = 1
    ictm = 1
    call swi_thksat(n, tn, bn, this%swi%zeta(n), satn_s, npf%inewton)
    if (satn_s * tthkn <= ZONE_MINIMUM_THICKNESS) &
      satn_s = ZONE_MINIMUM_THICKNESS / tthkn
    call swi_thksat(m, tm, bm, this%swi%zeta(m), satm_s, npf%inewton)
    if (satm_s * tthkm <= ZONE_MINIMUM_THICKNESS) &
      satm_s = ZONE_MINIMUM_THICKNESS / tthkm
    condsw = hcond(npf%ibound(n), npf%ibound(m), ictn, ictm, npf%inewton, ihc, &
                   npf%icellavg, npf%condsat(jas), hnew(n), hnew(m), &
                   satn_s, satm_s, hyn, hym, tn, tm, bn, bm, &
                   npf%dis%con%cl1(jas), npf%dis%con%cl2(jas), &
                   npf%dis%con%hwva(jas))
    !
    ! -- freshwater conductance and fill (same sign pattern as fc_default_flow)
    condf = cond - condsw
    idiag = npf%dis%con%ia(n)
    isymcon = npf%dis%con%isym(ipos)
    idiagm = npf%dis%con%ia(m)
    call matrix_sln%add_value_pos(idxglo(ipos), condf)
    call matrix_sln%add_value_pos(idxglo(idiag), -condf)
    call matrix_sln%add_value_pos(idxglo(isymcon), condf)
    call matrix_sln%add_value_pos(idxglo(idiagm), -condf)
  end subroutine swinpf_fc

  !> @brief SWI NPF formulation: Newton terms for the freshwater flow override.
  !! Mirrors fn_default_flow but uses the freshwater saturation derivative
  !! dS^f/dh = S'(h) + alphaf*S'(zeta) of the upstream cell. Horizontal only.
  !<
  subroutine swinpf_fn(this, n, m, ipos, matrix_sln, rhs, idxglo, hnew)
    use SmoothingModule, only: sQuadraticSaturationDerivative
    class(SwiNpfFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: m
    integer(I4B), intent(in) :: ipos
    class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: idxglo
    real(DP), dimension(:), intent(in) :: hnew
    ! local
    type(GwfNpfType), pointer :: npf
    integer(I4B) :: isymcon, idiag, idiagm, iups, idn, jas
    real(DP) :: cond, consterm, filledterm, derv, term, topup, botup, zetaup
    !
    npf => this%npf
    jas = npf%dis%con%jas(ipos)
    ! horizontal connections only (SWI single-fluid)
    if (npf%dis%con%ihc(jas) == C3D_VERTICAL) return
    idiag = npf%dis%con%ia(n)
    isymcon = npf%dis%con%isym(ipos)
    ! determine upstream node (by freshwater head)
    iups = m
    if (hnew(m) < hnew(n)) iups = n
    idn = n
    if (iups == n) idn = m
    topup = npf%dis%top(iups)
    botup = npf%dis%bot(iups)
    if (npf%dis%con%ihc(jas) == 2) then
      topup = min(npf%dis%top(n), npf%dis%top(m))
      botup = max(npf%dis%bot(n), npf%dis%bot(m))
    end if
    ! saturated conductance and the freshwater saturation derivative
    cond = npf%condsat(jas)
    consterm = -cond * (hnew(iups) - hnew(idn))
    zetaup = this%swi%get_zetanew(iups)
    derv = sQuadraticSaturationDerivative(topup, botup, hnew(iups), npf%satomega) + &
           this%swi%alphaf * &
           sQuadraticSaturationDerivative(topup, botup, zetaup, npf%satomega)
    idiagm = npf%dis%con%ia(m)
    if (iups == n) then
      term = consterm * derv
      rhs(n) = rhs(n) + term * hnew(n)
      rhs(m) = rhs(m) - term * hnew(n)
      call matrix_sln%add_value_pos(idxglo(idiag), term)
      if (npf%ibound(m) > 0) then
        call matrix_sln%add_value_pos(idxglo(isymcon), -term)
      end if
    else
      term = -consterm * derv
      rhs(n) = rhs(n) + term * hnew(m)
      rhs(m) = rhs(m) - term * hnew(m)
      if (npf%ibound(n) > 0) then
        call matrix_sln%add_value_pos(idxglo(ipos), term)
      end if
      call matrix_sln%add_value_pos(idxglo(idiagm), -term)
    end if
  end subroutine swinpf_fn

  !> @brief SWI NPF formulation: intercell flow (spike stub)
  !<
  subroutine swinpf_cq(this, n, m, ipos, flowja, h_new)
    class(SwiNpfFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: m
    integer(I4B), intent(in) :: ipos
    real(DP), dimension(:), intent(inout) :: flowja
    real(DP), dimension(:), intent(in) :: h_new
  end subroutine swinpf_cq

  !> @brief SWI STO formulation: is this cell handled by SWI storage?
  !<
  function swisto_is_active(this, n) result(is_active)
    class(SwiStoFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    logical(LGP) :: is_active
    is_active = .true.
  end function swisto_is_active

  !> @brief SWI STO formulation: freshwater storage fill for cell n, in
  !! override form. Assembled as the difference of two bottom-referenced columns:
  !! the water column (S^w, via the STO default) minus the salt column (S^s).
  !! COMMIT 1 (behavior-preserving relocation): the salt-column terms reproduce
  !! the former swi_fc_storage correction for the single-fluid freshwater model,
  !! branching on sto%inewton (Newton uses SyTerms(S^s) plus the swisto_fn
  !! tangent; Picard uses the -alphaf*rho2 chord). Cleanup into a proper
  !! chord-slope/smoothed-tangent split is Commit 2.
  !<
  subroutine swisto_fc(this, n, matrix_sln, rhs, idxglo, h_old, h_new)
    use TdisModule, only: delt
    use GwfStorageUtilsModule, only: SsCapacity, SyCapacity, SsTerms, SyTerms
    class(SwiStoFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: idxglo
    real(DP), dimension(:), intent(in) :: h_old
    real(DP), dimension(:), intent(in) :: h_new
    ! local
    type(GwfStoType), pointer :: sto
    integer(I4B) :: idiag
    real(DP) :: tled, tp, bt
    real(DP) :: sc1, rho1, rho1old, sc2, rho2, rho2old
    real(DP) :: snold, snnew, sfrac, aterm, rhsterm
    real(DP) :: zetanew, zetaold
    !
    sto => this%sto
    ! reset the Picard-branch storage arrays (read back in swisto_cq)
    this%swi%hcof(n) = DZERO
    this%swi%rhs(n) = DZERO
    if (this%swi%ibound(n) < 1) return
    tled = DONE / delt
    idiag = sto%dis%con%ia(n)
    tp = sto%dis%top(n)
    bt = sto%dis%bot(n)
    !
    ! -- water column (S^w): standard STO default fill (bottom-referenced)
    call sto%fc_default_sto(n, matrix_sln, rhs, idxglo, h_old, h_new)
    !
    ! -- salt column (S^s) correction: freshwater = S^w - S^s (single-fluid)
    zetanew = this%swi%get_zetanew(n)
    zetaold = this%swi%get_zetaold(n)
    !
    ! compressible (Ss): subtract the saltwater fraction of the S^w term
    if (sto%iconvert(n) == 0) then
      snold = DONE
      snnew = DONE
    else
      snold = sQuadraticSaturation(tp, bt, h_old(n), sto%satomega)
      snnew = sQuadraticSaturation(tp, bt, h_new(n), sto%satomega)
    end if
    sc1 = SsCapacity(sto%istor_coef, tp, bt, sto%dis%area(n), sto%ss(n))
    rho1 = sc1 * tled
    if (sto%integratechanges /= 0) then
      rho1old = SsCapacity(sto%istor_coef, tp, bt, sto%dis%area(n), &
                           this%swi%oldss(n)) * tled
    else
      rho1old = rho1
    end if
    call SsTerms(sto%iconvert(n), sto%iorig_ss, sto%iconf_ss, tp, bt, &
                 rho1, rho1old, snnew, snold, h_new(n), h_old(n), aterm, rhsterm)
    sfrac = sQuadraticSaturation(tp, bt, zetanew, sto%satomega)
    aterm = aterm * sfrac
    rhsterm = rhsterm * sfrac
    call matrix_sln%add_value_pos(idxglo(idiag), -aterm)
    rhs(n) = rhs(n) - rhsterm
    !
    ! drainable (Sy): interface-movement storage
    sc2 = SyCapacity(sto%dis%area(n), this%swi%sy(n))
    rho2 = sc2 * tled
    if (sto%integratechanges /= 0) then
      rho2old = SyCapacity(sto%dis%area(n), this%swi%oldsy(n)) * tled
    else
      rho2old = rho2
    end if
    snold = sQuadraticSaturation(tp, bt, zetaold, sto%satomega)
    snnew = sQuadraticSaturation(tp, bt, zetanew, sto%satomega)
    if (sto%inewton /= 0) then
      ! Newton: subtract SyTerms of the salt column (residual-consistent);
      ! the interface-movement tangent is added in swisto_fn
      call SyTerms(tp, bt, rho2, rho2old, snnew, snold, aterm, rhsterm)
      call matrix_sln%add_value_pos(idxglo(idiag), -aterm)
      rhs(n) = rhs(n) - rhsterm
    else
      ! Picard: -alphaf*rho2 stabilizing chord diagonal (stored in hcof/rhs so
      ! swisto_cq can recover the interface-storage rate)
      rho2 = -rho2 * this%swi%alphaf
      rho2old = -rho2old * this%swi%alphaf
      if (zetanew > bt .and. zetaold > bt) then
        this%swi%hcof(n) = rho2
        this%swi%rhs(n) = rho2 * h_old(n)
      else if (zetanew > bt .and. zetaold < bt) then
        this%swi%hcof(n) = rho2
        this%swi%rhs(n) = -rho2 * bt / this%swi%alphaf
      else if (zetanew < bt .and. zetaold > bt) then
        this%swi%hcof(n) = DZERO
        this%swi%rhs(n) = rho2 * (bt / this%swi%alphaf + h_old(n))
      end if
      call matrix_sln%add_value_pos(idxglo(idiag), this%swi%hcof(n))
      rhs(n) = rhs(n) + this%swi%rhs(n)
    end if
  end subroutine swisto_fc

  !> @brief Add the drainable-storage head-derivative (interface movement) as a
  !! Jacobian term (does not change the converged solution).
  !<
  subroutine swisto_sy_deriv(this, n, matrix_sln, rhs, idxglo, h_new, rho2, tthk)
    use SmoothingModule, only: sQuadraticSaturationDerivative
    class(SwiStoFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: idxglo
    real(DP), dimension(:), intent(in) :: h_new
    real(DP), intent(in) :: rho2, tthk
    real(DP) :: tp, bt, omega, zeta, dsfdh, drterm
    integer(I4B) :: idiag
    tp = this%sto%dis%top(n)
    bt = this%sto%dis%bot(n)
    omega = this%sto%satomega
    zeta = this%swi%get_zetanew(n)
    dsfdh = sQuadraticSaturationDerivative(tp, bt, h_new(n), omega) + &
            this%swi%alphaf * sQuadraticSaturationDerivative(tp, bt, zeta, omega)
    drterm = rho2 * tthk * dsfdh
    idiag = this%sto%dis%con%ia(n)
    call matrix_sln%add_value_pos(idxglo(idiag), -drterm)
    rhs(n) = rhs(n) - drterm * h_new(n)
  end subroutine swisto_sy_deriv

  !> @brief Freshwater saturation S^f = S^w(head) - S^s(zeta) for cell n
  !<
  function swisto_satf(this, n, head, old) result(satf)
    class(SwiStoFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    real(DP), intent(in) :: head
    logical(LGP), intent(in), optional :: old
    real(DP) :: satf
    real(DP) :: tp, bt, omega, zeta
    logical(LGP) :: isold
    isold = .false.
    if (present(old)) isold = old
    tp = this%sto%dis%top(n)
    bt = this%sto%dis%bot(n)
    omega = this%sto%satomega
    if (isold) then
      zeta = this%swi%get_zetaold(n)
    else
      zeta = this%swi%get_zetanew(n)
    end if
    satf = sQuadraticSaturation(tp, bt, head, omega) - &
           sQuadraticSaturation(tp, bt, zeta, omega)
  end function swisto_satf

  !> @brief SWI STO formulation: storage Newton terms for cell n. The
  !! interface-movement (drainable) head derivative is a pure Jacobian term and
  !! lives here (added on top of swisto_fc's residual-consistent fill). COMMIT 1:
  !! reproduces the former swi_fn_storage terms for the single-fluid freshwater
  !! model. Called only when Newton is active (sto_fn dispatch).
  !<
  subroutine swisto_fn(this, n, matrix_sln, rhs, idxglo, h_old, h_new)
    use TdisModule, only: delt
    use GwfStorageUtilsModule, only: SyCapacity
    class(SwiStoFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: idxglo
    real(DP), dimension(:), intent(in) :: h_old
    real(DP), dimension(:), intent(in) :: h_new
    ! local
    type(GwfStoType), pointer :: sto
    integer(I4B) :: idiag
    real(DP) :: tled, tp, bt, tthk, h
    real(DP) :: sc2, rho2, snnew, derv, rterm, drterm
    real(DP) :: zetanew
    !
    sto => this%sto
    if (this%swi%ibound(n) < 1) return
    tled = DONE / delt
    idiag = sto%dis%con%ia(n)
    tp = sto%dis%top(n)
    bt = sto%dis%bot(n)
    tthk = tp - bt
    h = h_new(n)
    zetanew = this%swi%get_zetanew(n)
    snnew = sQuadraticSaturation(tp, bt, zetanew)
    sc2 = SyCapacity(sto%dis%area(n), this%swi%sy(n))
    rho2 = sc2 * tled
    ! dS^s/dh_fresh * (1/tthk); dssdh = -alphaf for the freshwater model
    derv = -this%swi%alphaf / tthk
    ! interface-movement (Sy) Newton terms; applied where 0 < S^s < 1
    if (snnew < DONE .and. snnew > DZERO) then
      rterm = -rho2 * tthk * snnew
      drterm = -rho2 * tthk * derv
      call matrix_sln%add_value_pos(idxglo(idiag), -drterm - rho2)
      rhs(n) = rhs(n) - (-rterm + drterm * h + rho2 * bt)
    end if
  end subroutine swisto_fn

  !> @brief SWI STO formulation: freshwater storage rate for cell n, in override
  !! form. The water column (S^w) rate comes from the STO default (fills
  !! strgss/strgsy and flowja); the salt-column (S^s) correction is subtracted
  !! here. The interface-movement (drainable) rate is reported separately as the
  !! SWI storage term (swi%storage) -- the fresh<->salt conversion volume -- and
  !! also added to flowja. COMMIT 1: reproduces the former swi_cq storage block
  !! for the single-fluid freshwater model.
  !<
  subroutine swisto_cq(this, n, flowja, h_new, h_old)
    use TdisModule, only: delt
    use GwfStorageUtilsModule, only: SsCapacity, SyCapacity, SsTerms, SyTerms
    class(SwiStoFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    real(DP), dimension(:), intent(inout) :: flowja
    real(DP), dimension(:), intent(in) :: h_new
    real(DP), dimension(:), intent(in) :: h_old
    ! local
    type(GwfStoType), pointer :: sto
    integer(I4B) :: idiag
    real(DP) :: tled, tp, bt
    real(DP) :: sc1, rho1, rho1old, sc2, rho2, rho2old
    real(DP) :: snold, snnew, sfrac, aterm, rhsterm, rate
    real(DP) :: zetanew, zetaold
    !
    sto => this%sto
    this%swi%storage(n) = DZERO
    ! -- water column (S^w): standard STO default rate (fills strgss/strgsy)
    call sto%cq_default_sto(n, flowja, h_new, h_old)
    if (sto%iss == 1) return
    if (this%swi%ibound(n) < 1) return
    tled = DONE / delt
    idiag = sto%dis%con%ia(n)
    tp = sto%dis%top(n)
    bt = sto%dis%bot(n)
    zetanew = this%swi%get_zetanew(n)
    zetaold = this%swi%get_zetaold(n)
    !
    ! -- compressible (Ss): subtract the saltwater fraction from strgss and flowja
    if (sto%iconvert(n) == 0) then
      snold = DONE
      snnew = DONE
    else
      snold = sQuadraticSaturation(tp, bt, h_old(n), sto%satomega)
      snnew = sQuadraticSaturation(tp, bt, h_new(n), sto%satomega)
    end if
    sc1 = SsCapacity(sto%istor_coef, tp, bt, sto%dis%area(n), sto%ss(n))
    rho1 = sc1 * tled
    if (sto%integratechanges /= 0) then
      rho1old = SsCapacity(sto%istor_coef, tp, bt, sto%dis%area(n), &
                           this%swi%oldss(n)) * tled
    else
      rho1old = rho1
    end if
    call SsTerms(sto%iconvert(n), sto%iorig_ss, sto%iconf_ss, tp, bt, &
                 rho1, rho1old, snnew, snold, h_new(n), h_old(n), &
                 aterm, rhsterm, rate)
    sfrac = sQuadraticSaturation(tp, bt, zetanew, sto%satomega)
    rate = sfrac * rate
    sto%strgss(n) = sto%strgss(n) - rate
    flowja(idiag) = flowja(idiag) - rate
    !
    ! -- drainable (Sy): interface-movement storage -> SWI storage term, flowja
    rate = DZERO
    snold = sQuadraticSaturation(tp, bt, zetaold, sto%satomega)
    snnew = sQuadraticSaturation(tp, bt, zetanew, sto%satomega)
    if (sto%inewton /= 0) then
      sc2 = SyCapacity(sto%dis%area(n), this%swi%sy(n))
      rho2 = sc2 * tled
      if (sto%integratechanges /= 0) then
        rho2old = SyCapacity(sto%dis%area(n), this%swi%oldsy(n)) * tled
      else
        rho2old = rho2
      end if
      call SyTerms(tp, bt, rho2, rho2old, snnew, snold, aterm, rhsterm, rate)
      rate = -rate
    else
      ! Picard: recover the interface-storage rate from the fc chord terms
      rate = this%swi%hcof(n) * h_new(n) - this%swi%rhs(n)
    end if
    this%swi%storage(n) = rate
    flowja(idiag) = flowja(idiag) + rate
  end subroutine swisto_cq

  !> @brief SWI STO formulation: add SWI storage to the model budget
  !<
  subroutine swisto_bd(this, isuppress_output, model_budget)
    use TdisModule, only: delt
    use BudgetModule, only: BudgetType, rate_accumulator
    class(SwiStoFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: isuppress_output
    type(BudgetType), intent(inout) :: model_budget
    real(DP) :: rin, rout
    call rate_accumulator(this%swi%storage, rin, rout)
    call model_budget%addentry(rin, rout, delt, budtxt(1), &
                               isuppress_output, '     SWI')
  end subroutine swisto_bd

  !> @brief SWI STO formulation: save SWI storage cell-by-cell flows
  !<
  subroutine swisto_save_flows(this, iprint, ibinun)
    class(SwiStoFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: iprint
    integer(I4B), intent(in) :: ibinun
    integer(I4B) :: nvaluesp, nwidthp
    character(len=1) :: cdatafmp = ' ', editdesc = ' '
    real(DP) :: dinact
    if (ibinun /= 0) then
      dinact = DZERO
      call this%swi%dis%record_array(this%swi%storage, this%swi%iout, iprint, &
                                     -ibinun, budtxt(1), cdatafmp, nvaluesp, &
                                     nwidthp, editdesc, dinact)
    end if
  end subroutine swisto_save_flows

end module GwfSwiModule
