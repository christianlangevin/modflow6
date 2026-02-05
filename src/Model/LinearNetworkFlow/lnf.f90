!> @brief Linear Network Flow (LNF) Module
!!
!! This module contains the Linear Network Flow model, which is a
!! simplified flow model for network-like systems using DISU
!! discretization with static conductance.
!<
module LnfModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: LENFTYPE, LENMEMPATH, LENPACKAGETYPE, &
                             LINELENGTH, DZERO
  use SimModule, only: store_error
  use SimVariablesModule, only: errmsg, idm_context
  use BaseModelModule, only: BaseModelType
  use ListsModule, only: basemodellist
  use BaseModelModule, only: AddBaseModelToList
  use NumericalModelModule, only: NumericalModelType
  use MemoryManagerModule, only: mem_allocate
  use MemoryHelperModule, only: create_mem_path
  use BudgetModule, only: BudgetType, budget_cr
  use DisuModule, only: DisuType, disu_cr
  use LnfCndModule, only: LnfCndType, cnd_cr
  use LnfOcModule, only: LnfOcType, oc_cr
  use MatrixBaseModule, only: MatrixBaseType
  use BndModule, only: BndType, AddBndToList, GetBndFromList

  implicit none

  private
  public :: lnf_cr
  public :: LnfModelType
  public :: LNF_NBASEPKG, LNF_NMULTIPKG
  public :: LNF_BASEPKG, LNF_MULTIPKG

  type, extends(NumericalModelType) :: LnfModelType
    type(LnfCndType), pointer :: cnd => null() !< conductance package
    type(LnfOcType), pointer :: oc => null() !< output control package
    type(BudgetType), pointer :: budget => null() !< budget object
    integer(I4B), pointer :: incnd => null() !< CND enabled flag
    integer(I4B), pointer :: inoc => null() !< OC enabled flag
    integer(I4B), pointer :: iss => null() !< steady state flag
    integer(I4B), pointer :: inewtonur => null() !< newton under relaxation flag
  contains
    procedure :: model_df => lnf_df
    procedure :: model_ac => lnf_ac
    procedure :: model_mc => lnf_mc
    procedure :: model_ar => lnf_ar
    procedure :: model_rp => lnf_rp
    procedure :: model_ad => lnf_ad
    procedure :: model_fc => lnf_fc
    procedure :: model_cq => lnf_cq
    procedure :: model_bd => lnf_bd
    procedure :: model_ot => lnf_ot
    procedure :: model_da => lnf_da
    procedure :: allocate_scalars
    procedure :: allocate_arrays
    procedure, private :: create_packages
    procedure, private :: create_bndpkgs
    procedure, private :: package_create
    procedure, private :: set_namfile_options
    procedure, private :: log_namfile_options
  end type LnfModelType

  !> @brief LNF base package array descriptors
  !!
  !! LNF model base package types. Only listed packages are candidates
  !! for input and these will be loaded in the order specified.
  !<
  integer(I4B), parameter :: LNF_NBASEPKG = 50
  character(len=LENPACKAGETYPE), dimension(LNF_NBASEPKG) :: LNF_BASEPKG
  data LNF_BASEPKG/'DISU6', 'CND6 ', 'OC6  ', '     ', '     ', & !  5
                  &'     ', '     ', '     ', '     ', '     ', & ! 10
                  &40*'     '/ ! 50

  !> @brief LNF multi package array descriptors
  !!
  !! LNF model multi-instance package types. Only listed packages are
  !! candidates for input and these will be loaded in the order specified.
  !<
  integer(I4B), parameter :: LNF_NMULTIPKG = 50
  character(len=LENPACKAGETYPE), dimension(LNF_NMULTIPKG) :: LNF_MULTIPKG
  data LNF_MULTIPKG/'CHD6 ', '     ', '     ', '     ', '     ', & !  5
                   &'     ', '     ', '     ', '     ', '     ', & ! 10
                   &40*'     '/ ! 50

  ! -- size of supported model package arrays
  integer(I4B), parameter :: NIUNIT_LNF = LNF_NBASEPKG + LNF_NMULTIPKG

contains

  !> @brief Create a new Linear Network Flow model object
  !!
  !! (1) creates model object and add to modellist
  !! (2) assign values
  !!
  !<
  subroutine lnf_cr(filename, id, modelname)
    ! -- dummy
    character(len=*), intent(in) :: filename !< input file
    integer(I4B), intent(in) :: id !< consecutive model number listed in mfsim.nam
    character(len=*), intent(in) :: modelname !< name of the model
    ! -- local
    type(LnfModelType), pointer :: this
    class(BaseModelType), pointer :: model
    !
    ! -- Allocate a new LNF Model (this) and add it to basemodellist
    allocate (this)
    !
    ! -- Set memory path before allocation in memory manager can be done
    this%memoryPath = create_mem_path(modelname)
    !
    call this%allocate_scalars(modelname)
    model => this
    call AddBaseModelToList(basemodellist, model)
    !
    ! -- Assign values
    this%filename = filename
    this%name = modelname
    this%macronym = 'LNF'
    this%id = id
    !
    ! -- Set namfile options
    call this%set_namfile_options()
    !
    ! -- Create utility objects
    call budget_cr(this%budget, this%name)
    !
    ! -- Create model packages
    call this%create_packages()
    !
  end subroutine lnf_cr

  !> @brief Handle namefile options
  !<
  subroutine set_namfile_options(this)
    use MemoryManagerExtModule, only: mem_set_value
    use LnfNamInputModule, only: LnfNamParamFoundType
    class(LnfModelType) :: this
    type(LnfNamParamFoundType) :: found
    character(len=LENMEMPATH) :: input_mempath
    character(len=LINELENGTH) :: lst_fname

    ! -- set input model namfile memory path
    input_mempath = create_mem_path(this%name, 'NAM', idm_context)

    ! -- copy option params from input context
    call mem_set_value(lst_fname, 'LIST', input_mempath, found%list)
    call mem_set_value(this%iprpak, 'PRINT_INPUT', input_mempath, &
                       found%print_input)
    call mem_set_value(this%iprflow, 'PRINT_FLOWS', input_mempath, &
                       found%print_flows)
    call mem_set_value(this%ipakcb, 'SAVE_FLOWS', input_mempath, found%save_flows)

    ! -- create the list file
    call this%create_lstfile(lst_fname, this%filename, found%list, &
                             'LINEAR NETWORK FLOW MODEL (LNF)')

    ! -- activate save_flows if found
    if (found%save_flows) then
      this%ipakcb = -1
    end if

    ! -- log set options
    if (this%iout > 0) then
      call this%log_namfile_options(found)
    end if

  end subroutine set_namfile_options

  !> @brief Write model namfile options to list file
  !<
  subroutine log_namfile_options(this, found)
    use LnfNamInputModule, only: LnfNamParamFoundType
    class(LnfModelType) :: this
    type(LnfNamParamFoundType), intent(in) :: found

    write (this%iout, '(1x,a)') 'BEGIN NAMEFILE OPTIONS'

    if (found%print_input) then
      write (this%iout, '(4x,a)') 'STRESS PACKAGE INPUT WILL BE PRINTED '// &
        'FOR ALL MODEL STRESS PACKAGES'
    end if

    if (found%print_flows) then
      write (this%iout, '(4x,a)') 'PACKAGE FLOWS WILL BE PRINTED '// &
        'FOR ALL MODEL PACKAGES'
    end if

    if (found%save_flows) then
      write (this%iout, '(4x,a)') &
        'FLOWS WILL BE SAVED TO BUDGET FILE SPECIFIED IN OUTPUT CONTROL'
    end if

    write (this%iout, '(1x,a)') 'END NAMEFILE OPTIONS'

  end subroutine log_namfile_options

  !> @brief Define packages of the model
  !<
  subroutine lnf_df(this)
    ! -- dummy
    class(LnfModelType) :: this
    ! -- local
    integer(I4B) :: ip
    class(BndType), pointer :: packobj
    !
    ! -- Define discretization package
    call this%dis%dis_df()
    !
    ! -- Define CND package
    if (this%incnd > 0) then
      call this%cnd%cnd_df(this%dis)
    end if
    !
    ! -- Define OC package
    call this%oc%oc_df()
    !
    ! -- Define budget
    call this%budget%budget_df(NIUNIT_LNF, 'VOLUME', 'L**3')
    !
    ! -- Assign or point model members to dis members
    this%neq = this%dis%nodes
    this%nja = this%dis%nja
    this%ia => this%dis%con%ia
    this%ja => this%dis%con%ja
    !
    ! -- Allocate model arrays, now that neq and nja are known
    call this%allocate_arrays()
    !
    ! -- Define packages
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_df(this%neq, this%dis)
    end do
    !
  end subroutine lnf_df

  !> @brief Add the internal connections of this model to the sparse matrix
  !<
  subroutine lnf_ac(this, sparse)
    ! -- modules
    use SparseModule, only: sparsematrix
    ! -- dummy
    class(LnfModelType) :: this
    type(sparsematrix), intent(inout) :: sparse
    ! -- local
    class(BndType), pointer :: packobj
    integer(I4B) :: ip
    !
    ! -- Add the primary grid connections of this model to sparse
    call this%dis%dis_ac(this%moffset, sparse)
    !
    ! -- Add any package connections
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ac(this%moffset, sparse)
    end do
    !
  end subroutine lnf_ac

  !> @brief Map the positions of this models connections in the
  !! numerical solution coefficient matrix.
  !<
  subroutine lnf_mc(this, matrix_sln)
    ! -- dummy
    class(LnfModelType) :: this
    class(MatrixBaseType), pointer :: matrix_sln
    ! -- local
    class(BndType), pointer :: packobj
    integer(I4B) :: ip
    !
    ! -- Find the position of each connection in the global ia, ja structure
    !    and store them in idxglo.
    call this%dis%dis_mc(this%moffset, this%idxglo, matrix_sln)
    !
    ! -- Map any package connections
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_mc(this%moffset, matrix_sln)
    end do
    !
  end subroutine lnf_mc

  !> @brief Allocate and read
  !<
  subroutine lnf_ar(this)
    ! -- dummy
    class(LnfModelType) :: this
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: ip
    class(BndType), pointer :: packobj
    !
    ! -- Set ibound to 1 for all nodes (no idomain yet)
    do i = 1, this%dis%nodes
      this%ibound(i) = 1
    end do
    !
    ! -- Allocate and read CND package
    if (this%incnd > 0) then
      call this%cnd%cnd_ar()
    end if
    !
    ! -- Set up output control
    call this%oc%oc_ar(this%x, this%dis, DZERO)
    call this%budget%set_ibudcsv(this%oc%ibudcsv)
    !
    ! -- Package input files now open, so allocate and read
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%set_pointers(this%dis%nodes, this%ibound, this%x, &
                                this%xold, this%flowja)
      ! -- Read and allocate package
      call packobj%bnd_ar()
    end do
    !
  end subroutine lnf_ar

  !> @brief Read and prepare (calls package RP routines)
  !<
  subroutine lnf_rp(this)
    ! -- dummy
    class(LnfModelType) :: this
    ! -- local
    class(BndType), pointer :: packobj
    integer(I4B) :: ip
    !
    ! -- Read and prepare output control (OC)
    if (this%inoc > 0) call this%oc%oc_rp()
    !
    ! -- Call boundary package RP routines
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_rp()
      call packobj%bnd_rp_obs()
    end do
    !
  end subroutine lnf_rp

  !> @brief Advance (advance time)
  !<
  subroutine lnf_ad(this)
    ! -- dummy
    class(LnfModelType) :: this
    ! -- local
    class(BndType), pointer :: packobj
    integer(I4B) :: ip
    !
    ! -- Call boundary package AD routines
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ad()
    end do
    !
  end subroutine lnf_ad

  !> @brief Fill coefficients
  !<
  subroutine lnf_fc(this, kiter, matrix_sln, inwtflag)
    ! -- dummy
    class(LnfModelType) :: this
    integer(I4B), intent(in) :: kiter
    class(MatrixBaseType), pointer :: matrix_sln
    integer(I4B), intent(in) :: inwtflag
    ! -- local
    class(BndType), pointer :: packobj
    integer(I4B) :: ip
    !
    ! -- Fill conductance terms
    if (this%incnd > 0) then
      call this%cnd%cnd_fc(kiter, matrix_sln, this%idxglo, this%rhs, this%x)
    end if
    !
    ! -- Call boundary package FC routines
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_fc(this%rhs, this%ia, this%idxglo, matrix_sln)
    end do
    !
  end subroutine lnf_fc

  !> @brief Calculate flows
  !<
  subroutine lnf_cq(this, icnvg, isuppress_output)
    use SparseModule, only: csr_diagsum
    ! -- dummy
    class(LnfModelType) :: this
    integer(I4B), intent(in) :: icnvg
    integer(I4B), intent(in) :: isuppress_output
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: ip
    class(BndType), pointer :: packobj
    !
    ! -- Construct the flowja array.  Flowja is calculated each time, even if
    !    output is suppressed.  (flowja is positive into a cell.)
    do i = 1, this%nja
      this%flowja(i) = DZERO
    end do
    !
    ! -- Calculate CND flows
    if (this%incnd > 0) then
      call this%cnd%cnd_cq(this%x, this%flowja)
    end if
    !
    ! -- Go through packages and call cq routines
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_cf()
      call packobj%bnd_cq(this%x, this%flowja)
    end do
    !
  end subroutine lnf_cq

  !> @brief Model budget
  !<
  subroutine lnf_bd(this, icnvg, isuppress_output)
    use SparseModule, only: csr_diagsum
    ! -- dummy
    class(LnfModelType) :: this
    integer(I4B), intent(in) :: icnvg
    integer(I4B), intent(in) :: isuppress_output
    ! -- local
    integer(I4B) :: ip
    class(BndType), pointer :: packobj
    !
    ! -- Finalize calculation of flowja by adding face flows to the diagonal.
    call csr_diagsum(this%dis%con%ia, this%flowja)
    !
    ! -- Save the solution convergence flag
    this%icnvg = icnvg
    !
    ! -- Budget routines (start by resetting)
    call this%budget%reset()
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_bd(this%budget)
    end do
    !
  end subroutine lnf_bd

  !> @brief Model output
  !<
  subroutine lnf_ot(this)
    use TdisModule, only: kstp, kper, tdis_ot, endofperiod
    use TdisModule, only: totim
    ! -- dummy
    class(LnfModelType) :: this
    ! -- local
    integer(I4B) :: idvsave
    integer(I4B) :: idvprint
    integer(I4B) :: icbcfl
    integer(I4B) :: icbcun
    integer(I4B) :: ibudfl
    integer(I4B) :: ipflag
    integer(I4B) :: ip
    class(BndType), pointer :: packobj
    ! -- formats
    character(len=*), parameter :: fmtnocnvg = &
      "(1X,/9X,'****FAILED TO MEET SOLVER CONVERGENCE CRITERIA IN TIME STEP ', &
      &I0,' OF STRESS PERIOD ',I0,'****')"
    !
    ! -- Set write and print flags
    idvsave = 0
    idvprint = 0
    icbcfl = 0
    ibudfl = 0
    if (this%oc%oc_save('HEAD')) idvsave = 1
    if (this%oc%oc_print('HEAD')) idvprint = 1
    if (this%oc%oc_save('BUDGET')) icbcfl = 1
    if (this%oc%oc_print('BUDGET')) ibudfl = 1
    icbcun = this%oc%oc_save_unit('BUDGET')
    !
    ! -- Override ibudfl and idvprint flags for nonconvergence
    !    and end of period
    ibudfl = this%oc%set_print_flag('BUDGET', this%icnvg, endofperiod)
    idvprint = this%oc%set_print_flag('HEAD', this%icnvg, endofperiod)
    !
    ! -- Calculate and save package observations
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_bd_obs()
      call packobj%bnd_ot_obs()
    end do
    !
    ! -- Save CND flows (FLOW-JA-FACE)
    if (this%incnd > 0) then
      call this%cnd%cnd_save_model_flows(this%flowja, icbcfl, icbcun)
    end if
    !
    ! -- Save package flows
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ot_model_flows(icbcfl=icbcfl, ibudfl=0, icbcun=icbcun)
    end do
    !
    ! -- Print package flows
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ot_model_flows(icbcfl=icbcfl, ibudfl=ibudfl, icbcun=0)
    end do
    !
    ! -- Save and print dependent variables
    ipflag = 0
    call this%oc%oc_ot(ipflag)
    !
    ! -- Print budget summary
    if (ibudfl /= 0) then
      ipflag = 1
      call this%budget%budget_ot(kstp, kper, this%iout)
    end if
    !
    ! -- Write to budget csv every time step
    call this%budget%writecsv(totim)
    !
    ! -- Timing Output
    if (ipflag == 1) call tdis_ot(this%iout)
    !
    ! -- Write non-convergence message
    if (this%icnvg == 0) then
      write (this%iout, fmtnocnvg) kstp, kper
    end if
    !
  end subroutine lnf_ot

  !> @brief Deallocate model memory
  !<
  subroutine lnf_da(this)
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(LnfModelType) :: this
    ! -- local
    class(BndType), pointer :: packobj
    integer(I4B) :: ip
    !
    ! -- Deallocate boundary packages
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_da()
      deallocate (packobj)
    end do
    !
    ! -- Deallocate CND package
    if (this%incnd > 0) then
      call this%cnd%cnd_da()
      deallocate (this%cnd)
    end if
    !
    ! -- Deallocate OC package
    call this%oc%oc_da()
    deallocate (this%oc)
    !
    ! -- Deallocate scalars
    call mem_deallocate(this%incnd)
    call mem_deallocate(this%inoc)
    call mem_deallocate(this%iss)
    call mem_deallocate(this%inewtonur)
    !
    ! -- Deallocate parent
    call this%NumericalModelType%model_da()
    !
    ! -- Deallocate budget
    call this%budget%budget_da()
    deallocate (this%budget)
    !
    ! -- Deallocate discretization object
    call this%dis%dis_da()
    deallocate (this%dis)
    !
  end subroutine lnf_da

  !> @brief Allocate memory for scalar members
  !<
  subroutine allocate_scalars(this, modelname)
    ! -- dummy
    class(LnfModelType) :: this
    character(len=*), intent(in) :: modelname
    !
    ! -- allocate members from parent class
    call this%NumericalModelType%allocate_scalars(modelname)
    !
    ! -- allocate members that are part of model class
    call mem_allocate(this%incnd, 'INCND', this%memoryPath)
    call mem_allocate(this%inoc, 'INOC', this%memoryPath)
    call mem_allocate(this%iss, 'ISS', this%memoryPath)
    call mem_allocate(this%inewtonur, 'INEWTONUR', this%memoryPath)
    !
    ! -- initialize
    this%incnd = 0
    this%inoc = 0
    this%iss = 1 ! default is steady-state
    this%inewtonur = 0
    !
  end subroutine allocate_scalars

  !> @brief Allocate memory for arrays
  !<
  subroutine allocate_arrays(this)
    ! -- dummy
    class(LnfModelType) :: this
    !
    ! -- allocate members from parent class
    call this%NumericalModelType%allocate_arrays()
    !
  end subroutine allocate_arrays

  !> @brief Create model packages
  !<
  subroutine create_packages(this)
    use ConstantsModule, only: LINELENGTH
    use MemoryManagerModule, only: mem_setptr
    use MemoryHelperModule, only: create_mem_path
    use CharacterStringModule, only: CharacterStringType
    use ArrayHandlersModule, only: expandarray
    ! -- dummy
    class(LnfModelType) :: this
    ! -- local
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: pkgtypes => null()
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: pkgnames => null()
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: mempaths => null()
    integer(I4B), dimension(:), contiguous, pointer :: inunits => null()
    character(len=LENMEMPATH) :: model_mempath
    character(len=LINELENGTH) :: pkgtype, pkgname, mempath
    integer(I4B), dimension(:), allocatable :: bndpkgs
    integer(I4B) :: n
    integer(I4B) :: indis = 0 ! DIS enabled flag
    !
    ! -- Set model memory path
    model_mempath = create_mem_path(component=this%name, context=idm_context)
    !
    ! -- Set pointers to model path package attribute arrays
    call mem_setptr(pkgtypes, 'PKGTYPES', model_mempath)
    call mem_setptr(pkgnames, 'PKGNAMES', model_mempath)
    call mem_setptr(mempaths, 'MEMPATHS', model_mempath)
    call mem_setptr(inunits, 'INUNITS', model_mempath)
    !
    ! -- Create packages
    do n = 1, size(pkgtypes)
      pkgtype = pkgtypes(n)
      pkgname = pkgnames(n)
      mempath = mempaths(n)
      !
      select case (pkgtype)
      case ('DISU6')
        indis = 1
        call disu_cr(this%dis, this%name, mempath, indis, this%iout)
      case ('CND6')
        this%incnd = 1
        call cnd_cr(this%cnd, this%name, mempath, this%incnd, this%iout, &
                    this%dis)
      case ('OC6')
        this%inoc = 1
        call oc_cr(this%oc, this%name, mempath, this%inoc, this%iout)
      case ('CHD6')
        call expandarray(bndpkgs)
        bndpkgs(size(bndpkgs)) = n
      case default
        write (errmsg, '(a,a)') &
          'Unknown package type for LNF model: ', trim(pkgtype)
        call store_error(errmsg)
      end select
    end do
    !
    ! -- Create boundary packages
    call this%create_bndpkgs(bndpkgs, pkgtypes, pkgnames, mempaths, inunits)
    !
  end subroutine create_packages

  !> @brief Create boundary condition packages for this model
  !<
  subroutine package_create(this, filtyp, ipakid, ipaknum, pakname, mempath, &
                            inunit, iout)
    ! -- modules
    use ChdModule, only: chd_create
    ! -- dummy
    class(LnfModelType) :: this
    character(len=*), intent(in) :: filtyp
    integer(I4B), intent(in) :: ipakid
    integer(I4B), intent(in) :: ipaknum
    character(len=*), intent(in) :: pakname
    character(len=*), intent(in) :: mempath
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    ! -- local
    class(BndType), pointer :: packobj
    class(BndType), pointer :: packobj2
    integer(I4B) :: ip
    !
    ! -- This part creates the package object
    select case (filtyp)
    case ('CHD6')
      call chd_create(packobj, ipakid, ipaknum, inunit, iout, this%name, &
                      pakname, mempath)
      ! -- Set ictMemPath to empty string since LNF has no NPF
      packobj%ictMemPath = ''
    case default
      write (errmsg, *) 'Invalid package type: ', filtyp
      call store_error(errmsg, terminate=.TRUE.)
    end select
    !
    ! -- Check to make sure that the package name is unique, then store a
    !    pointer to the package in the model bndlist
    do ip = 1, this%bndlist%Count()
      packobj2 => GetBndFromList(this%bndlist, ip)
      if (packobj2%packName == pakname) then
        write (errmsg, '(a,a)') 'Cannot create package.  Package name  '// &
          'already exists: ', trim(pakname)
        call store_error(errmsg, terminate=.TRUE.)
      end if
    end do
    call AddBndToList(this%bndlist, packobj)
    !
  end subroutine package_create

  !> @brief Source package info and begin to process
  !<
  subroutine create_bndpkgs(this, bndpkgs, pkgtypes, pkgnames, &
                            mempaths, inunits)
    ! -- modules
    use ConstantsModule, only: LINELENGTH, LENPACKAGENAME
    use CharacterStringModule, only: CharacterStringType
    ! -- dummy
    class(LnfModelType) :: this
    integer(I4B), dimension(:), allocatable, intent(inout) :: bndpkgs
    type(CharacterStringType), dimension(:), contiguous, &
      pointer, intent(inout) :: pkgtypes
    type(CharacterStringType), dimension(:), contiguous, &
      pointer, intent(inout) :: pkgnames
    type(CharacterStringType), dimension(:), contiguous, &
      pointer, intent(inout) :: mempaths
    integer(I4B), dimension(:), contiguous, &
      pointer, intent(inout) :: inunits
    ! -- local
    integer(I4B) :: ipakid, ipaknum
    character(len=LENFTYPE) :: pkgtype, bndptype
    character(len=LENPACKAGENAME) :: pkgname
    character(len=LENMEMPATH) :: mempath
    integer(I4B), pointer :: inunit
    integer(I4B) :: n

    if (allocated(bndpkgs)) then
      !
      ! -- create stress packages
      ipakid = 1
      bndptype = ''
      do n = 1, size(bndpkgs)
        !
        pkgtype = pkgtypes(bndpkgs(n))
        pkgname = pkgnames(bndpkgs(n))
        mempath = mempaths(bndpkgs(n))
        inunit => inunits(bndpkgs(n))
        !
        if (bndptype /= pkgtype) then
          ipaknum = 1
          bndptype = pkgtype
        end if
        !
        call this%package_create(pkgtype, ipakid, ipaknum, pkgname, mempath, &
                                 inunit, this%iout)
        ipakid = ipakid + 1
        ipaknum = ipaknum + 1
      end do
      !
      ! -- cleanup
      deallocate (bndpkgs)
    end if
  end subroutine create_bndpkgs

end module LnfModule
