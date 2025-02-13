!> @brief Surface Water Flow (SWF) Junction (JNC) Module
!!
!! This module manages junctions where cells connect.
!<
module SwfJncModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENMEMPATH, LINELENGTH, &
                             DZERO, DHALF, DONE, DTWO, &
                             DTWOTHIRDS, DP9, DONETHIRD, &
                             DPREC, DEM10
  use MemoryHelperModule, only: create_mem_path
  use MemoryManagerModule, only: mem_allocate, mem_setptr, get_isize, &
                                 mem_reallocate
  use SimVariablesModule, only: errmsg
  use SimModule, only: count_errors, store_error, store_error_unit, &
                       store_error_filename
  use NumericalPackageModule, only: NumericalPackageType
  use BaseDisModule, only: DisBaseType
  use Disv1dModule, only: Disv1dType
  use SwfCxsModule, only: SwfCxsType
  use SwfDfwModule, only: SwfDfwType
  use ObsModule, only: ObsType, obs_cr
  use ObserveModule, only: ObserveType
  use MatrixBaseModule

  implicit none
  private
  public :: SwfJncType, jnc_cr

  type, extends(NumericalPackageType) :: SwfJncType

    integer(I4B), pointer :: moffset => null() !< equation offset
    integer(I4B), pointer :: njunction => null() !< number of junctions
    integer(I4B), pointer :: njunc_red => null() !< number of reduced junctions (only those without 2 connected cells)
    integer(I4B), pointer :: nvert => null() !< number of vertices
    integer(I4B), pointer :: nq => null() !< number of unique q equations
    integer(I4B), dimension(:), pointer, contiguous :: junc_ivert => null() !< the vertex number corresponding to each junction, size (njunction)
    integer(I4B), dimension(:), pointer, contiguous :: ivert_junc => null() !< the junction number that corresponds to the vertex, size (nvert); 0 if vertex does not correspond to a junction
    integer(I4B), dimension(:), pointer, contiguous :: iajunction_cell => null() !< the index array for jajunction_cell, size (njunction + 1)
    integer(I4B), dimension(:), pointer, contiguous :: jajunction_cell => null() !< csr array to map junction number to cells ()
    integer(I4B), dimension(:), pointer, contiguous :: jstart => null() !< the starting junction number for each cell, size (ncells)
    integer(I4B), dimension(:), pointer, contiguous :: jend => null() !< the ending junction number for each cell, size (ncells)
    integer(I4B), dimension(:), pointer, contiguous :: junc_juncred => null() !< given a junction number, return the reduced number; 0 if not part of reduced set, size (njunctions)
    integer(I4B), dimension(:), pointer, contiguous :: juncred_junc => null() !< given a reduced junction number, return the junction number, size (njunc_red)
    integer(I4B), dimension(:), pointer, contiguous :: irowqup => null() !< given a cell number, return row for the upside q (ncells)
    integer(I4B), dimension(:), pointer, contiguous :: irowqdn => null() !< given a cell number, return row for the downside q (ncells)
    integer(I4B), dimension(:), pointer, contiguous :: icellup => null() !< given a cell number, return number of up side cell if no active junction, otherwise a zero if there is an active junction (ncells)
    integer(I4B), dimension(:), pointer, contiguous :: icelldn => null() !< given a cell number, return number of dn side cell if no active junction, otherwise a zero if there is an active junction (ncells)


    ! user-provided input
    ! integer(I4B), pointer :: is2d => null() !< flag to indicate this model is 2D overland flow and not 1d channel flow
    ! integer(I4B), pointer :: icentral => null() !< flag to use central in space weighting (default is upstream weighting)
    ! integer(I4B), pointer :: iswrcond => null() !< flag to activate the dev SWR conductance formulation
    ! real(DP), pointer :: unitconv !< conversion factor used in mannings equation; calculated from timeconv and lengthconv
    ! real(DP), pointer :: timeconv !< conversion factor from model length units to meters (1.0 if model uses meters for length)
    ! real(DP), pointer :: lengthconv !< conversion factor from model time units to seconds (1.0 if model uses seconds for time)
    ! real(DP), dimension(:), pointer, contiguous :: hnew => null() !< pointer to model xnew
    ! real(DP), dimension(:), pointer, contiguous :: manningsn => null() !< mannings roughness for each reach
    ! integer(I4B), dimension(:), pointer, contiguous :: idcxs => null() !< cross section id for each reach
    ! integer(I4B), dimension(:), pointer, contiguous :: ibound => null() !< pointer to model ibound
    ! integer(I4B), dimension(:), pointer, contiguous :: icelltype => null() !< set to 1 and is accessed by chd for checking

    ! velocity
    ! integer(I4B), pointer :: icalcvelocity => null() !< flag to indicate velocity will be calculated (always on)
    ! integer(I4B), pointer :: isavvelocity => null() !< flag to indicate velocity will be saved
    ! real(DP), dimension(:, :), pointer, contiguous :: vcomp => null() !< velocity components: vx, vy, vz (nodes, 3)
    ! real(DP), dimension(:), pointer, contiguous :: vmag => null() !< velocity magnitude (of size nodes)
    ! integer(I4B), pointer :: nedges => null() !< number of cell edges
    ! integer(I4B), pointer :: lastedge => null() !< last edge number
    ! integer(I4B), dimension(:), pointer, contiguous :: nodedge => null() !< array of node numbers that have edges
    ! integer(I4B), dimension(:), pointer, contiguous :: ihcedge => null() !< edge type (horizontal or vertical)
    ! real(DP), dimension(:, :), pointer, contiguous :: propsedge => null() !< edge properties (Q, area, nx, ny, distance)
    ! real(DP), dimension(:), pointer, contiguous :: grad_dhds_mag => null() !< magnitude of the gradient (of size nodes)
    ! real(DP), dimension(:), pointer, contiguous :: dhdsja => null() !< gradient for each connection (of size njas)

    ! observation data
    integer(I4B), pointer :: inobspkg => null() !< unit number for obs package
    type(ObsType), pointer :: obs => null() !< observation package

    ! pointers
    type(SwfCxsType), pointer :: cxs
    type(SwfDfwType), pointer :: dfw

  contains

    procedure :: jnc_df
    procedure :: allocate_scalars
    procedure :: allocate_arrays
    procedure :: jnc_ac
    procedure :: jnc_mc
    ! procedure :: dfw_load
    ! procedure :: source_options
    ! procedure :: log_options
    ! procedure :: source_griddata
    ! procedure :: log_griddata
    ! procedure :: dfw_ar
    ! procedure :: dfw_rp
    ! procedure :: dfw_ad
    procedure :: jnc_fc
    ! procedure :: dfw_qnm_fc_nr
    ! !procedure :: dfw_qnm_fc
    ! procedure :: dfw_fn
    ! procedure :: dfw_nur
    ! procedure :: dfw_cq
    ! procedure :: dfw_bd
    ! procedure :: dfw_save_model_flows
    ! procedure :: dfw_print_model_flows
    procedure :: jnc_da
    ! procedure :: dfw_df_obs
    ! procedure :: dfw_rp_obs
    ! procedure :: dfw_bd_obs
    ! procedure :: qcalc
    ! procedure :: get_cond
    ! procedure :: get_cond_swr
    ! procedure :: get_cond_n
    ! procedure :: get_flow_area_nm
    ! procedure :: calc_velocity
    ! procedure :: sav_velocity
    ! procedure, public :: increase_edge_count
    ! procedure, public :: set_edge_properties
    ! procedure :: calc_dhds
    ! procedure :: write_cxs_tables
    procedure, private :: jglo_qup
    procedure, private :: jglo_qdn
    procedure, private :: jglo_hup
    procedure, private :: jglo_hdn
    procedure, private :: get_icell_up
    procedure, private :: get_icell_dn

  end type SwfJncType

contains

  !> @brief create package
  !<
  subroutine jnc_cr(jncobj, name_model, input_mempath, inunit, iout, &
                    dfw, cxs)
    ! modules
    use MemoryManagerExtModule, only: mem_set_value
    ! dummy
    type(SwfJncType), pointer :: jncobj !< object to create
    character(len=*), intent(in) :: name_model !< name of the SWF model
    character(len=*), intent(in) :: input_mempath !< memory path
    integer(I4B), intent(in) :: inunit !< flag to indicate if package is active
    integer(I4B), intent(in) :: iout !< unit number for output
    type(SwfDfwType), pointer, intent(in) :: dfw !< the pointer to the dfw package
    type(SwfCxsType), pointer, intent(in) :: cxs !< the pointer to the cxs package
    ! locals
    logical(LGP) :: found_fname
    ! formats
    character(len=*), parameter :: fmtheader = &
      "(1x, /1x, 'DFW --  JUNCTION (JNC) PACKAGE, VERSION 1, 1/31/2025', &
       &' INPUT READ FROM MEMPATH: ', A, /)"
    !
    ! Create the object
    allocate (jncobj)

    ! create name and memory path
    call jncobj%set_names(1, name_model, 'JNC', 'JNC')

    ! Allocate scalars
    call jncobj%allocate_scalars()

    ! Set variables
    jncobj%input_mempath = input_mempath
    jncobj%inunit = inunit
    jncobj%iout = iout

    ! set name of input file
    call mem_set_value(jncobj%input_fname, 'INPUT_FNAME', jncobj%input_mempath, &
                       found_fname)

    ! Set a pointers to passed in objects
    jncobj%cxs => cxs
    jncobj%dfw => dfw

    ! create obs package
    call obs_cr(jncobj%obs, jncobj%inobspkg)

    ! check if jnc is enabled
    if (inunit > 0) then

      ! Print a message identifying the package.
      write (iout, fmtheader) input_mempath

    end if

  end subroutine jnc_cr

  !> @brief Define
  !<
  subroutine jnc_df(this, neq, dis)
    ! dummy
    class(SwfJncType) :: this !< this instance
    integer(I4B), intent(inout) :: neq !< number of equations
    class(DisBaseType), pointer :: dis !< the pointer to the discretization
    ! local
    character(len=*), parameter :: fmtneq = &
      &"(1x, 'The SWF Model is configured for dynamic wave equations',&
      &/3x, 'Number of cell continuity equations: ', I0, &
      &/3x, 'Number of unique motion equations: ', I0, &
      &/3x, 'Number of junction continuity equations: ', I0, &
      &/3x, 'Total number of equations: ', I0&
      &)"

    ! Set a pointers to passed in objects
    this%dis => dis

    ! update number of equations
    select type (dis)
    class is (Disv1dType)
      call count_junctions(dis%nvert, dis%iavert, dis%javert, this%njunction, &
                           this%njunc_red)
      this%nvert = dis%nvert
    end select

    ! allocate arrays
    call this%allocate_arrays()

    ! calculate junction arrays
    select type (dis)
    class is (Disv1dType)

      ! fill arrays to map between vertex number and junction number
      call fill_junc_ivert(dis%nvert, dis%iavert, dis%javert, this%junc_ivert, &
                           this%ivert_junc)

      ! fill csr arrays relating junction to cells
      call fill_junction_cell(this%njunction, dis%iavert, dis%javert, &
                              this%ivert_junc, this%iajunction_cell, &
                              this%jajunction_cell, this%memoryPath)

      ! fill arrays to map between vertex number and junction number
      call fill_jstart_jend(dis%nvert, dis%iavert, dis%javert, &
                           this%ivert_junc, this%jstart, this%jend)

      ! fill arrays to map between junction and reduced junction
      call fill_junc_juncred(this%iajunction_cell, this%jajunction_cell, &
                             this%junc_juncred, this%juncred_junc)

      ! fill icellup with up side cell number or zero if active junction on up side
      call fill_icellupdn(this%jend, this%iajunction_cell, this%jajunction_cell, &
                          this%icellup, this%icelldn)

      ! fill irowqup and irowqdn
      call fill_irow_qupdn(this%jstart, this%jend, this%junc_juncred, &
                           this%icellup, this%nq, this%irowqup, this%irowqdn)

    end select
    print *, "njunction", this%njunction
    print *, "njunc_red", this%njunc_red
    print *, "nvert", this%nvert
    print *, "nq", this%nq
    print *, "nodesuser", dis%nodesuser
    print *, "junc_ivert", this%junc_ivert
    print *, "ivert_junc", this%ivert_junc

    ! for some reason this is not right!
    print *, "ia", this%iajunction_cell
    print *, "ja", this%jajunction_cell
    print *, "jstart", this%jstart
    print *, "jend", this%jend
    print *, "junc_juncred", this%junc_juncred
    print *, "juncred_junc", this%juncred_junc
    print *, "irowqup", this%irowqup
    print *, "irowqdn", this%irowqdn
    print *, "icellup", this%icellup
    print *, "icelldn", this%icelldn


    ! Number of equations include:
    !   reach continuity: dis%nodes
    !   motion equations (flow at reach ends): number of unique qs (nq)
    !   junction continuity: njunctions
    ! Number of unknowns include:
    !   reach stage (nreach)
    !   unique flows (nq)
    !   junction stages (njunctions)
    neq = dis%nodes
    neq = neq + this%nq
    neq = neq + this%njunc_red
    write (this%iout, fmtneq) dis%nodes, this%nq, this%njunc_red, neq

    ! ! Set the distype (either DISV1D or DIS2D)
    ! if (this%dis%is_2d()) then
    !   this%is2d = 1
    ! end if

    ! ! check if dfw is enabled
    ! ! this will need to become if (.not. present(dfw_options)) then
    ! !if (inunit > 0) then


    ! ! load dfw
    ! call this%dfw_load()

    !end if

  end subroutine jnc_df

  !> @brief Add connections to sparse cell connectivity matrix
  !!
  !!                                                    Columns
  !!              h1  h2  h3  ...  hn  Qup1  Qdn1  Qup2  Qdn2  ...  Qupn  Qdnn  hj1  hj2  hj3  ...  hjn 
  !! rows
  !!
  !! Cont cell 1
  !!      cell 2
  !!      ...
  !!      cell n
  !!
  !! Flow Qup_c1
  !!      Qdn_c1
  !!      Qup_c2
  !!      Qdn_c2
  !!      ...
  !!      Qup_cn
  !!      Qdn_cn
  !!
  !! Junc Qjnc1
  !!      Qjnc2
  !!      ...
  !!      Qjncn
  !<
  subroutine jnc_ac(this, moffset, sparse)
    ! modules
    use SparseModule, only: sparsematrix
    ! dummy
    class(SwfJncType) :: this
    integer(I4B), intent(in) :: moffset
    type(sparsematrix), intent(inout) :: sparse
    ! local
    integer(I4B) :: i, j, iglo, jglo, ired, jred

    ! Rows are as follows
    !   nodes: cell continuity (node1, node2, node3, ...)
    !   2 * nodes: motion equations (Qup1, Qdn1, Qup2, Qdn2, ...)
    !   njunctions: junction continuity (j1, j2, ...)

    ! Columns are as follows:
    !   nodes: cell heads (h1, h2, h3, ...)
    !   Q: flows at cell edges (Qup1, Qdn1, Qup2, Qdn2, ...)
    !   njunctions: junction heads (hj1, hj2, hj3, ...)

    ! cell continuity equation -- each row has a diagonal component (hn)
    ! and hup and hdn components
    do i = 1, this%dis%nodes

      ! diagonal component, head in cell
      iglo = i + moffset
      jglo = iglo
      call sparse%addconnection(iglo, jglo, 1)

      ! head on up side (either up side cell or up side junction)
      jglo = this%jglo_hup(i, moffset)
      call sparse%addconnection(iglo, jglo, 1)

      ! head on dn side (either dn side head or dn side junction)
      jglo = this%jglo_hdn(i, moffset)
      call sparse%addconnection(iglo, jglo, 1)

    end do

    ! motion equations -- each row has a coefficient entry for the Q itself,
    ! the cell head, the junction head, and the Q on the other side of the reach
    do i = 1, this%dis%nodes

      ! First, process the down side flow for cell i; there is a dn side Q
      ! equation for every cell.  Sign will be positive into cell i.

      ! diagonal
      iglo = this%jglo_qdn(i, moffset)
      jglo = iglo
      call sparse%addconnection(iglo, jglo, 1)

      ! head for cell i
      jglo = i + moffset
      call sparse%addconnection(iglo, jglo, 1)

      ! head on dn side (either junction or cell)
      jglo = this%jglo_hdn(i, moffset)
      call sparse%addconnection(iglo, jglo, 1)

      ! flow on other side (up)
      jglo = this%jglo_qup(i, moffset)
      call sparse%addconnection(iglo, jglo, 1)

    end do

    do i = 1, this%dis%nodes

      ! Second, process the up side flow for cell i.  Not every cell will have
      ! its own unique Q equation for the up side.  If there is no active
      ! junction on the up side, then Qup for cell i is the same (but reverse 
      ! in sign) for Qdn of the up side cell.

      j = this%jstart(i)
      jred = this%junc_juncred(j)
      if (jred == 0) then
        ! this cell is connected to an upstream cell, so there is no unique
        ! flow term here.  Only formulate this flow equation if there is an
        ! active junction on the up side.
        cycle
      end if

      ! diagonal
      iglo = this%jglo_qup(i, moffset)
      jglo = iglo
      call sparse%addconnection(iglo, jglo, 1)

      ! head for cell i
      jglo = i + moffset
      call sparse%addconnection(iglo, jglo, 1)

      ! head on up side (either junction or cell)
      jglo = this%jglo_hup(i, moffset)
      call sparse%addconnection(iglo, jglo, 1)

      ! flow on other side (dn)
      jglo = this%jglo_qdn(i, moffset)
      call sparse%addconnection(iglo, jglo, 1)

    end do

    ! add diagonal for each junction
    do i = 1, this%njunc_red

      ! global row number for this junction
      iglo = moffset + this%dis%nodes + this%nq + i
      jglo = iglo
      call sparse%addconnection(iglo, jglo, 1)

    end do

    ! construct structure of junction equations
    do j = 1, this%dis%nodes

      ! up side junction for this cell
      i = this%jstart(j)
      ired = this%junc_juncred(i)

      if (ired /= 0) then ! active junction

        ! global row number for this junction
        iglo = moffset + this%dis%nodes + this%nq + ired

        ! global column number for cell head
        jglo = moffset + j
        call sparse%addconnection(iglo, jglo, 1)
      
      end if

      ! ! up side Q
      ! jglo = this%jglo_qup(j, moffset)
      ! call sparse%addconnection(iglo, jglo, 1)

      ! dn side junction for this cell
      i = this%jend(j)
      ired = this%junc_juncred(i)

      if (ired /= 0) then ! active junction

        ! global row number for this junction
        iglo = moffset + this%dis%nodes + this%nq + ired

        ! global column number for cell head
        jglo = moffset + j
        call sparse%addconnection(iglo, jglo, 1)
      
      end if

      ! ! dn side Q
      ! jglo = this%jglo_qdn(j, moffset)
      ! call sparse%addconnection(iglo, jglo, 1)

    end do

  end subroutine jnc_ac

  !> @brief Map cell connections in the numerical solution coefficient matrix.
  subroutine jnc_mc(this, moffset, idxglo, matrix_sln)
    ! dummy
    class(SwfJncType) :: this
    integer(I4B), intent(in) :: moffset
    integer(I4B), dimension(:), intent(inout) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! local
    integer(I4B) :: i, j, ipos, iglo, jglo

    ! todo: need to figure out what global matrix positions to store in
    ! idxglo.  Idxglo is allocated to size swf%nja.  Not sure that is of any
    ! use.
    do i = 1, this%dis%nodes
      iglo = i + moffset
      do ipos = this%dis%con%ia(i), this%dis%con%ia(i + 1) - 1
        j = this%dis%con%ja(ipos)
        jglo = j + moffset
        idxglo(ipos) = matrix_sln%get_position(iglo, jglo)
      end do
    end do

  end subroutine jnc_mc

  !> @brief fill coefficients for jnc
  subroutine jnc_fc(this, kiter, moffset, matrix_sln, idxglo, rhs, stage, &
                    stage_old)
    ! modules
    ! dummy
    class(SwfJncType) :: this !< this instance
    integer(I4B) :: kiter
    integer(I4B), intent(in) :: moffset
    class(MatrixBaseType), pointer :: matrix_sln
    integer(I4B), intent(in), dimension(:) :: idxglo
    real(DP), intent(inout), dimension(:) :: rhs
    real(DP), intent(inout), dimension(:) :: stage
    real(DP), intent(inout), dimension(:) :: stage_old
    ! local
    integer(I4B) :: ipos_sln
    integer(I4B) :: i, j, iglo, jglo, jred, ired
    real(DP) :: val
    real(DP) :: cond
    real(DP) :: depth, dx, width, dhds, cln, clm

    ! half cell conductance for cell i (what head should this be a function of?)
    cond = 1.d0 ! this is half-cell conductance for cell i

    ! add coefficients for cell continuity equations
    ! cell continuity equation -- each row has a diagonal component (hn)
    ! and Qup and Qdn components
    do i = 1, this%dis%nodes

      ! cln = 500.d0
      ! clm = 500.d0

      ! ! need conductance term for cell i in the up direction
      ! ! if iup is a junction, then use half-cell conductance of cell i
      ! ! if iup is a connected cell, then use averaged conductance between cell i and iup
      ! j = this%jstart(i)
      ! jred = this%junc_juncred(j)
      ! if (jred == 0) then
      !   ! average conductance between cell n and cell m
      !   cond = this%dfw%get_cond(i, this%get_icell_up(i), 1, stage(i), stage(this%get_icell_up(i)), cln, clm)
      ! else
      !   ! half cell conductance between cell i and junction jred
      !   depth = stage(i) - this%dis%bot(i)
      !   dx = cln
      !   width = 1.d0
      !   dhds = 1.D0
      !   cond = this%dfw%get_cond_n(i, depth, dx, width, dhds)
      ! end if

      ! diagonal component
      iglo = i + moffset
      jglo = iglo
      ipos_sln = matrix_sln%get_position(iglo, jglo)
      val = -cond ! sum of conductance terms; TODO: and then there will be a storage term
      call matrix_sln%add_value_pos(ipos_sln, val)

      ! head in junction or cell on up side
      jglo = this%jglo_hup(i, moffset)
      ipos_sln = matrix_sln%get_position(iglo, jglo)
      ! todo: will need to know whether this is half cell junction cond or averaged cell-cell cond
      val = cond ! conductance term and then there will be a storage term
      call matrix_sln%add_value_pos(ipos_sln, val)

      ! diagonal component
      iglo = i + moffset
      jglo = iglo
      ipos_sln = matrix_sln%get_position(iglo, jglo)
      val = -cond ! sum of conductance terms; TODO: and then there will be a storage term
      call matrix_sln%add_value_pos(ipos_sln, val)

      ! head in junction or cell on dn side
      jglo = this%jglo_hdn(i, moffset)
      ipos_sln = matrix_sln%get_position(iglo, jglo)
      ! todo: will need to know whether this is half cell junction cond or averaged cell-cell cond
      val = cond ! conductance term and then there will be a storage term
      call matrix_sln%add_value_pos(ipos_sln, val)

      ! ! Qup
      ! jglo = this%jglo_qup(i, moffset)
      ! ipos_sln = matrix_sln%get_position(iglo, jglo)
      ! val = 1.d0  ! this should be reach length / 2
      ! call matrix_sln%add_value_pos(ipos_sln, val)

      ! ! Qdn
      ! jglo = this%jglo_qdn(i, moffset)
      ! ipos_sln = matrix_sln%get_position(iglo, jglo)
      ! val = 1.d0  ! this should be reach length / 2
      ! call matrix_sln%add_value_pos(ipos_sln, val)

    end do

    ! motion equations -- each row has a coefficient entry for the Q itself,
    ! the cell head, the junction head, and the Q on the other side of the reach
    do i = 1, this%dis%nodes
      ! First, process the down side flow for cell i

      ! diagonal
      iglo = this%jglo_qdn(i, moffset)
      jglo = iglo
      ipos_sln = matrix_sln%get_position(iglo, jglo)
      val = -1.d0 ! will also need to add temporal intertial term
      call matrix_sln%add_value_pos(ipos_sln, val)

      ! cell head
      jglo = i + moffset
      ipos_sln = matrix_sln%get_position(iglo, jglo)
      val = -cond
      call matrix_sln%add_value_pos(ipos_sln, val)

      ! junction or cell head on dn side
      jglo = this%jglo_hdn(i, moffset)
      ipos_sln = matrix_sln%get_position(iglo, jglo)
      val = cond
      call matrix_sln%add_value_pos(ipos_sln, val)

      ! flow on other side (up)
      jglo = this%jglo_qup(i, moffset)
      ipos_sln = matrix_sln%get_position(iglo, jglo)
      val = DZERO ! this will be a spatial inertial term
      call matrix_sln%add_value_pos(ipos_sln, val)

    end do


    do i = 1, this%dis%nodes

      ! Second, process the up side flow for cell i
      j = this%jstart(i)
      jred = this%junc_juncred(j)
      if (jred == 0) then
        ! this cell is connected to an upstream cell, so there is no unique
        ! flow term here.  Only formulate this flow equation if there is an
        ! active junction on the up side.
        cycle
      end if

      ! diagonal (Q term itself)
      iglo = this%jglo_qup(i, moffset)
      jglo = iglo
      ipos_sln = matrix_sln%get_position(iglo, jglo)
      val = -1.d0 ! will also need to add temporal intertial term
      call matrix_sln%add_value_pos(ipos_sln, val)

      ! coefficient for cell head
      jglo = moffset + i
      ipos_sln = matrix_sln%get_position(iglo, jglo)
      val = -cond
      call matrix_sln%add_value_pos(ipos_sln, val)

      ! junction head on up side
      jglo = this%jglo_hup(i, moffset)
      ipos_sln = matrix_sln%get_position(iglo, jglo)
      val = cond
      call matrix_sln%add_value_pos(ipos_sln, val)

      ! flow on other side (dn)
      jglo = this%jglo_qdn(i, moffset)
      ipos_sln = matrix_sln%get_position(iglo, jglo)
      val = DZERO ! this will be a spatial inertial term
      call matrix_sln%add_value_pos(ipos_sln, val)
    
    end do

    ! ! add diagonal for each junction
    ! do i = 1, this%njunction

    !   ! global row number for this junction
    !   iglo = moffset + 3 * this%dis%nodes + i
    !   jglo = iglo
    !   ipos_sln = matrix_sln%get_position(iglo, jglo)
    !   val = DZERO ! no diagonal term for junctions
    !   call matrix_sln%add_value_pos(ipos_sln, val)

    ! end do

    do j = 1, this%dis%nodes

      ! up side junction for this cell
      i = this%jstart(j)
      ired = this%junc_juncred(i)

      if (ired /= 0) then

        ! diagonal position
        iglo = moffset + this%dis%nodes + this%nq + ired
        jglo = iglo
        ipos_sln = matrix_sln%get_position(iglo, jglo)
        val = -cond ! coefficient for hj
        call matrix_sln%add_value_pos(ipos_sln, val)

        ! global column number for cell head
        jglo = moffset + j
        ipos_sln = matrix_sln%get_position(iglo, jglo)
        val = cond ! coefficient for hi
        call matrix_sln%add_value_pos(ipos_sln, val)

      end if
      
      ! ! up side Q
      ! jglo = this%jglo_qup(j, moffset)
      ! ipos_sln = matrix_sln%get_position(iglo, jglo)
      ! val = -1.d0 ! 100 percent of flow is for junction; negative to indicate flow into junction
      ! call matrix_sln%add_value_pos(ipos_sln, val)

      ! dn side junction for this cell
      i = this%jend(j)
      ired = this%junc_juncred(i)

      if (ired /= 0) then

        ! diagonal position
        iglo = moffset + this%dis%nodes + this%nq + ired
        jglo = iglo
        ipos_sln = matrix_sln%get_position(iglo, jglo)
        val = -cond ! coefficient for hj
        call matrix_sln%add_value_pos(ipos_sln, val)

        ! global column number for cell head
        jglo = moffset + j
        ipos_sln = matrix_sln%get_position(iglo, jglo)
        val = cond ! coefficient for hi
        call matrix_sln%add_value_pos(ipos_sln, val)

      end if

      ! ! dn side Q
      ! jglo = this%jglo_qdn(j, moffset)
      ! ipos_sln = matrix_sln%get_position(iglo, jglo)
      ! val = -1.d0 ! 100 percent of flow is for junction; negative to indicate flow into junction
      ! call matrix_sln%add_value_pos(ipos_sln, val)

    end do

  end subroutine jnc_fc

  !> @ brief Return position in x array for flow on up side
  function jglo_qup(this, node, moffset) result(jglo)
    ! dummy
    class(SwfJncType) :: this !< this instance
    integer(I4B), intent(in) :: node !< node number in model indices
    integer(I4B), intent(in) :: moffset !< model offset (0 for none)
    integer(I4B) :: jglo
    jglo = moffset + this%dis%nodes + this%irowqup(node)
  end function jglo_qup

  !> @ brief Return position in x array for flow on down side
  function jglo_qdn(this, node, moffset) result(jglo)
    ! dummy
    class(SwfJncType) :: this !< this instance
    integer(I4B), intent(in) :: node !< node number in model indices
    integer(I4B), intent(in) :: moffset !< model offset (0 for none)
    integer(I4B) :: jglo
    jglo = moffset + this%dis%nodes + this%irowqdn(node)
  end function jglo_qdn

  !> @ brief Return position in x array for junction head or head in upside cell
  function jglo_hup(this, node, moffset) result(jglo)
    ! dummy
    class(SwfJncType) :: this !< this instance
    integer(I4B), intent(in) :: node !< node number in model indices
    integer(I4B), intent(in) :: moffset !< model offset (0 for none)
    ! dummy
    integer(I4B) :: j
    integer(I4B) :: jred
    ! return
    integer(I4B) :: jglo
    j = this%jstart(node)
    jred = this%junc_juncred(j)
    if (jred == 0) then ! no junction on upside
      ! calculate position for next up side cell
      jglo = moffset + this%get_icell_up(node)
    else
      ! use location of junction head
      jglo = moffset + this%dis%nodes + this%nq + jred
    end if
  end function jglo_hup

  !> @ brief Return number of up side cell
  function get_icell_up(this, node) result(icell_up)
    class(SwfJncType) :: this !< this instance
    integer(I4B), intent(in) :: node
    integer(I4B) :: icell_up
    icell_up = this%icellup(node)
  end function get_icell_up

  !> @ brief Return number of up side cell
  function get_icell_dn(this, node) result(icell_dn)
    class(SwfJncType) :: this !< this instance
    integer(I4B), intent(in) :: node
    integer(I4B) :: icell_dn
    icell_dn = this%icelldn(node)
  end function get_icell_dn

  !> @ brief Return position in x array for junction head on down side of cell
  function jglo_hdn(this, node, moffset) result(jglo)
    ! dummy
    class(SwfJncType) :: this !< this instance
    integer(I4B), intent(in) :: node !< node number in model indices
    integer(I4B), intent(in) :: moffset !< model offset (0 for none)
    ! dummy
    integer(I4B) :: j
    integer(I4B) :: jred
    ! return
    integer(I4B) :: jglo
    j = this%jend(node)
    jred = this%junc_juncred(j)
    if (jred == 0) then ! no junction on down side
      ! calculate position for next down side cell
      jglo = moffset + node + 1
      if (node + 1 > this%dis%nodes) then
        print *, "i am brokin"
        stop
      end if
    else
      ! use location of junction head
      jglo = moffset + this%dis%nodes + this%nq + jred
    end if
  end function jglo_hdn

  !> @ brief Allocate scalars
  !!
  !! Allocate and initialize scalars for the package. The base model
  !! allocate scalars method is also called.
  !!
  !<
  subroutine allocate_scalars(this)
    ! modules
    ! dummy
    class(SwfJncType) :: this !< this instance
    !
    ! allocate scalars in NumericalPackageType
    call this%NumericalPackageType%allocate_scalars()

    ! Allocate scalars
    call mem_allocate(this%moffset, 'MOFFSET', this%memoryPath)
    call mem_allocate(this%njunction, 'NJUNCTION', this%memoryPath)
    call mem_allocate(this%njunc_red, 'NJUNC_RED', this%memoryPath)
    call mem_allocate(this%nvert, 'NVERT', this%memoryPath)
    call mem_allocate(this%nq, 'NQ', this%memoryPath)
  !   call mem_allocate(this%is2d, 'IS2D', this%memoryPath)
  !   call mem_allocate(this%icentral, 'ICENTRAL', this%memoryPath)
  !   call mem_allocate(this%iswrcond, 'ISWRCOND', this%memoryPath)
  !   call mem_allocate(this%unitconv, 'UNITCONV', this%memoryPath)
  !   call mem_allocate(this%lengthconv, 'LENGTHCONV', this%memoryPath)
  !   call mem_allocate(this%timeconv, 'TIMECONV', this%memoryPath)
  !   call mem_allocate(this%inobspkg, 'INOBSPKG', this%memoryPath)
  !   call mem_allocate(this%icalcvelocity, 'ICALCVELOCITY', this%memoryPath)
  !   call mem_allocate(this%isavvelocity, 'ISAVVELOCITY', this%memoryPath)
  !   call mem_allocate(this%nedges, 'NEDGES', this%memoryPath)
  !   call mem_allocate(this%lastedge, 'LASTEDGE', this%memoryPath)

    this%moffset = 0
    this%njunction = 0
    this%njunc_red = 0
    this%nvert = 0
    this%nq = 0
  !   this%is2d = 0
  !   this%icentral = 0
  !   this%iswrcond = 0
  !   this%unitconv = DONE
  !   this%lengthconv = DONE
  !   this%timeconv = DONE
  !   this%inobspkg = 0
  !   this%icalcvelocity = 0
  !   this%isavvelocity = 0
  !   this%nedges = 0
  !   this%lastedge = 0

  end subroutine allocate_scalars

  !> @brief allocate memory for arrays
  !<
  subroutine allocate_arrays(this)
    ! dummy
    class(SwfJncType) :: this !< this instance
    ! locals
    integer(I4B) :: n

    call mem_allocate(this%junc_ivert, this%njunction, 'JUNC_IVERT', &
                      this%memoryPath)
    call mem_allocate(this%ivert_junc, this%nvert, 'IVERT_JUNC', &
                      this%memoryPath)
    call mem_allocate(this%jstart, this%dis%nodesuser, 'JSTART', &
                      this%memoryPath)
    call mem_allocate(this%jend, this%dis%nodesuser, 'JEND', &
                      this%memoryPath)
    call mem_allocate(this%junc_juncred, this%njunction, 'JUNC_JUNCRED', &
                      this%memoryPath)
    call mem_allocate(this%juncred_junc, this%njunction, 'JUNCRED_JUNC', &
                      this%memoryPath)
    call mem_allocate(this%irowqdn, this%dis%nodesuser, 'IROWQDN', &
                      this%memoryPath)
    call mem_allocate(this%irowqup, this%dis%nodesuser, 'IROWQUP', &
                      this%memoryPath)
    call mem_allocate(this%icellup, this%dis%nodesuser, 'ICELLUP', &
                      this%memoryPath)
    call mem_allocate(this%icelldn, this%dis%nodesuser, 'ICELLDN', &
                      this%memoryPath)

    do n = 1, size(this%junc_ivert)
      this%junc_ivert(n) = 0
    end do
    do n = 1, size(this%ivert_junc)
      this%ivert_junc(n) = 0
    end do
    do n = 1, size(this%jstart)
      this%jstart(n) = 0
    end do
    do n = 1, size(this%jend)
      this%jend(n) = 0
    end do
    do n = 1, size(this%junc_juncred)
      this%junc_juncred(n) = 0
    end do
    do n = 1, size(this%juncred_junc)
      this%juncred_junc(n) = 0
    end do
    do n = 1, size(this%irowqup)
      this%irowqup(n) = 0
    end do
    do n = 1, size(this%irowqdn)
      this%irowqdn(n) = 0
    end do
    do n = 1, size(this%icellup)
      this%icellup(n) = 0
    end do
    do n = 1, size(this%icelldn)
      this%icelldn(n) = 0
    end do

  end subroutine allocate_arrays

  ! !> @brief load data from IDM to package
  ! !<
  ! subroutine dfw_load(this)
  !   ! dummy
  !   class(SwfJncType) :: this !< this instance

  !   ! source input data
  !   call this%source_options()
  !   call this%source_griddata()

  ! end subroutine dfw_load

  ! !> @brief Copy options from IDM into package
  ! !<
  ! subroutine source_options(this)
  !   ! modules
  !   use KindModule, only: LGP
  !   use InputOutputModule, only: getunit, openfile
  !   use MemoryManagerExtModule, only: mem_set_value
  !   use CharacterStringModule, only: CharacterStringType
  !   use SwfDfwInputModule, only: SwfDfwParamFoundType
  !   ! dummy
  !   class(SwfJncType) :: this !< this instance
  !   ! locals
  !   integer(I4B) :: isize
  !   type(SwfDfwParamFoundType) :: found
  !   type(CharacterStringType), dimension(:), pointer, &
  !     contiguous :: obs6_fnames

  !   ! update defaults with idm sourced values
  !   call mem_set_value(this%icentral, 'ICENTRAL', &
  !                      this%input_mempath, found%icentral)
  !   call mem_set_value(this%iswrcond, 'ISWRCOND', &
  !                      this%input_mempath, found%iswrcond)
  !   call mem_set_value(this%lengthconv, 'LENGTHCONV', &
  !                      this%input_mempath, found%lengthconv)
  !   call mem_set_value(this%timeconv, 'TIMECONV', &
  !                      this%input_mempath, found%timeconv)
  !   call mem_set_value(this%iprflow, 'IPRFLOW', &
  !                      this%input_mempath, found%iprflow)
  !   call mem_set_value(this%ipakcb, 'IPAKCB', &
  !                      this%input_mempath, found%ipakcb)
  !   call mem_set_value(this%isavvelocity, 'ISAVVELOCITY', &
  !                      this%input_mempath, found%isavvelocity)

  !   ! save flows option active
  !   if (found%icentral) this%icentral = 1
  !   if (found%ipakcb) this%ipakcb = -1

  !   ! calculate unit conversion
  !   this%unitconv = this%lengthconv**DONETHIRD
  !   this%unitconv = this%unitconv * this%timeconv

  !   ! save velocity active
  !   if (found%isavvelocity) this%icalcvelocity = this%isavvelocity

  !   ! check for obs6_filename
  !   call get_isize('OBS6_FILENAME', this%input_mempath, isize)
  !   if (isize > 0) then
  !     !
  !     if (isize /= 1) then
  !       errmsg = 'Multiple OBS6 keywords detected in OPTIONS block.'// &
  !                ' Only one OBS6 entry allowed.'
  !       call store_error(errmsg)
  !       call store_error_filename(this%input_fname)
  !     end if

  !     call mem_setptr(obs6_fnames, 'OBS6_FILENAME', this%input_mempath)

  !     found%obs6_filename = .true.
  !     this%obs%inputFilename = obs6_fnames(1)
  !     this%obs%active = .true.
  !     this%inobspkg = GetUnit()
  !     this%obs%inUnitObs = this%inobspkg
  !     call openfile(this%inobspkg, this%iout, this%obs%inputFilename, 'OBS')
  !     call this%obs%obs_df(this%iout, this%packName, this%filtyp, this%dis)
  !     call this%dfw_df_obs()
  !   end if

  !   ! log values to list file
  !   if (this%iout > 0) then
  !     call this%log_options(found)
  !   end if

  ! end subroutine source_options

  ! !> @brief Write user options to list file
  ! !<
  ! subroutine log_options(this, found)
  !   use SwfDfwInputModule, only: SwfDfwParamFoundType
  !   class(SwfJncType) :: this !< this instance
  !   type(SwfDfwParamFoundType), intent(in) :: found

  !   write (this%iout, '(1x,a)') 'Setting DFW Options'

  !   if (found%lengthconv) then
  !     write (this%iout, '(4x,a, G0)') 'Mannings length conversion value &
  !                                 &specified as ', this%lengthconv
  !   end if

  !   if (found%timeconv) then
  !     write (this%iout, '(4x,a, G0)') 'Mannings time conversion value &
  !                                 &specified as ', this%timeconv
  !   end if

  !   if (found%lengthconv .or. found%timeconv) then
  !     write (this%iout, '(4x,a, G0)') 'Mannings conversion value calculated &
  !                                 &from user-provided length_conversion and &
  !                                 &time_conversion is ', this%unitconv
  !   end if

  !   if (found%iprflow) then
  !     write (this%iout, '(4x,a)') 'Cell-by-cell flow information will be printed &
  !                                 &to listing file whenever ICBCFL is not zero.'
  !   end if

  !   if (found%ipakcb) then
  !     write (this%iout, '(4x,a)') 'Cell-by-cell flow information will be printed &
  !                                 &to listing file whenever ICBCFL is not zero.'
  !   end if

  !   if (found%obs6_filename) then
  !     write (this%iout, '(4x,a)') 'Observation package is active.'
  !   end if

  !   if (found%isavvelocity) &
  !     write (this%iout, '(4x,a)') 'Velocity will be calculated at cell &
  !                                 &centers and written to DATA-VCOMP in budget &
  !                                 &file when requested.'

  !   if (found%iswrcond) then
  !     write (this%iout, '(4x,a, G0)') 'Conductance will be calculated using &
  !                                      &the SWR development option.'
  !   end if

  !   write (this%iout, '(1x,a,/)') 'End Setting DFW Options'

  ! end subroutine log_options

  ! !> @brief copy griddata from IDM to package
  ! !<
  ! subroutine source_griddata(this)
  !   ! modules
  !   use SimModule, only: count_errors, store_error
  !   use MemoryHelperModule, only: create_mem_path
  !   use MemoryManagerExtModule, only: mem_set_value
  !   use SimVariablesModule, only: idm_context
  !   use SwfDfwInputModule, only: SwfDfwParamFoundType
  !   ! dummy
  !   class(SwfJncType) :: this !< this instance
  !   ! locals
  !   character(len=LENMEMPATH) :: idmMemoryPath
  !   type(SwfDfwParamFoundType) :: found
  !   integer(I4B), dimension(:), pointer, contiguous :: map

  !   ! set memory path
  !   idmMemoryPath = create_mem_path(this%name_model, 'DFW', idm_context)

  !   ! set map to convert user input data into reduced data
  !   map => null()
  !   if (this%dis%nodes < this%dis%nodesuser) map => this%dis%nodeuser

  !   ! update defaults with idm sourced values
  !   call mem_set_value(this%manningsn, 'MANNINGSN', &
  !                      idmMemoryPath, map, found%manningsn)
  !   call mem_set_value(this%idcxs, 'IDCXS', idmMemoryPath, map, found%idcxs)

  !   ! ensure MANNINGSN was found
  !   if (.not. found%manningsn) then
  !     write (errmsg, '(a)') 'Error in GRIDDATA block: MANNINGSN not found.'
  !     call store_error(errmsg)
  !   end if

  !   if (count_errors() > 0) then
  !     call store_error_filename(this%input_fname)
  !   end if

  !   ! log griddata
  !   if (this%iout > 0) then
  !     call this%log_griddata(found)
  !   end if

  ! end subroutine source_griddata

  ! !> @brief log griddata to list file
  ! !<
  ! subroutine log_griddata(this, found)
  !   use SwfDfwInputModule, only: SwfDfwParamFoundType
  !   class(SwfJncType) :: this !< this instance
  !   type(SwfDfwParamFoundType), intent(in) :: found

  !   write (this%iout, '(1x,a)') 'Setting DFW Griddata'

  !   if (found%manningsn) then
  !     write (this%iout, '(4x,a)') 'MANNINGSN set from input file'
  !   end if

  !   if (found%idcxs) then
  !     write (this%iout, '(4x,a)') 'IDCXS set from input file'
  !   end if

  !   call this%write_cxs_tables()

  !   write (this%iout, '(1x,a,/)') 'End Setting DFW Griddata'

  ! end subroutine log_griddata

  ! subroutine write_cxs_tables(this)
  !   ! modules
  !   ! dummy
  !   class(SwfJncType) :: this !< this instance
  !   ! local
  !   ! integer(I4B) :: idcxs
  !   ! integer(I4B) :: n

  !   !-- TODO: write cross section tables
  !   ! do n = 1, this%dis%nodes
  !   !   idcxs = this%idcxs(n)
  !   !   if (idcxs > 0) then
  !   !     call this%cxs%write_cxs_table(idcxs, this%width(n), this%slope(n), &
  !   !                                   this%manningsn(n), this%unitconv)
  !   !   end if
  !   ! end do
  ! end subroutine write_cxs_tables

  ! !> @brief allocate memory
  ! !<
  ! subroutine dfw_ar(this, ibound, hnew)
  !   ! modules
  !   ! dummy
  !   class(SwfJncType) :: this !< this instance
  !   integer(I4B), dimension(:), pointer, contiguous :: ibound !< model ibound array
  !   real(DP), dimension(:), pointer, contiguous, intent(inout) :: hnew !< pointer to model head array
  !   ! local
  !   integer(I4B) :: n

  !   ! store pointer to ibound
  !   this%ibound => ibound
  !   this%hnew => hnew

  !   if (this%icalcvelocity == 1) then
  !     call mem_reallocate(this%vcomp, 3, this%dis%nodes, 'VCOMP', this%memoryPath)
  !     call mem_reallocate(this%vmag, this%dis%nodes, 'VMAG', this%memoryPath)
  !     call mem_reallocate(this%nodedge, this%nedges, 'NODEDGE', this%memoryPath)
  !     call mem_reallocate(this%ihcedge, this%nedges, 'IHCEDGE', this%memoryPath)
  !     call mem_reallocate(this%propsedge, 5, this%nedges, 'PROPSEDGE', &
  !                         this%memoryPath)
  !     do n = 1, this%dis%nodes
  !       this%vcomp(:, n) = DZERO
  !       this%vmag(n) = DZERO
  !     end do
  !   end if

  !   ! observation data
  !   call this%obs%obs_ar()

  ! end subroutine dfw_ar

  ! !> @brief allocate memory
  ! !<
  ! subroutine dfw_rp(this)
  !   ! modules
  !   ! dummy
  !   class(SwfJncType) :: this !< this instance

  !   ! read observations
  !   call this%dfw_rp_obs()

  ! end subroutine dfw_rp

  ! !> @brief advance
  ! !<
  ! subroutine dfw_ad(this, irestore)
  !   class(SwfJncType) :: this !< this instance
  !   integer(I4B), intent(in) :: irestore !< ATS flag for retrying time step (1) or advancing (0)

  !   ! Push simulated values to preceding time/subtime step
  !   call this%obs%obs_ad()

  ! end subroutine dfw_ad

  ! !> @brief fill coefficients
  ! !!
  ! !! The DFW Package is entirely Newton based.  All matrix and rhs terms
  ! !! are added from thish routine.
  ! !!
  ! !<
  ! subroutine dfw_fc(this, kiter, matrix_sln, idxglo, rhs, stage, stage_old)
  !   ! modules
  !   ! dummy
  !   class(SwfJncType) :: this !< this instance
  !   integer(I4B) :: kiter
  !   class(MatrixBaseType), pointer :: matrix_sln
  !   integer(I4B), intent(in), dimension(:) :: idxglo
  !   real(DP), intent(inout), dimension(:) :: rhs
  !   real(DP), intent(inout), dimension(:) :: stage
  !   real(DP), intent(inout), dimension(:) :: stage_old
  !   ! local

  !   ! calculate dhds at cell center for 2d case
  !   if (this%is2d == 1) then
  !     call this%calc_dhds()
  !   end if

  !   ! add qnm contributions to matrix equations
  !   call this%dfw_qnm_fc_nr(kiter, matrix_sln, idxglo, rhs, stage, stage_old)

  ! end subroutine dfw_fc

  ! !> @brief fill coefficients
  ! !!
  ! !< Add qnm contributions to matrix equations
  ! subroutine dfw_qnm_fc_nr(this, kiter, matrix_sln, idxglo, rhs, stage, stage_old)
  !   ! modules
  !   use MathUtilModule, only: get_perturbation
  !   ! dummy
  !   class(SwfJncType) :: this !< this instance
  !   integer(I4B) :: kiter
  !   class(MatrixBaseType), pointer :: matrix_sln
  !   integer(I4B), intent(in), dimension(:) :: idxglo
  !   real(DP), intent(inout), dimension(:) :: rhs
  !   real(DP), intent(inout), dimension(:) :: stage
  !   real(DP), intent(inout), dimension(:) :: stage_old
  !   ! local
  !   integer(I4B) :: n, m, ii, idiag
  !   real(DP) :: qnm
  !   real(DP) :: qeps
  !   real(DP) :: eps
  !   real(DP) :: derv

  !   ! Calculate conductance and put into amat
  !   do n = 1, this%dis%nodes

  !     ! Find diagonal position for row n
  !     idiag = this%dis%con%ia(n)

  !     ! Loop through connections adding matrix terms
  !     do ii = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1

  !       ! skip for masked cells
  !       if (this%dis%con%mask(ii) == 0) cycle

  !       ! connection variables
  !       m = this%dis%con%ja(ii)

  !       ! Fill the qnm term on the right-hand side
  !       qnm = this%qcalc(n, m, stage(n), stage(m), ii)
  !       rhs(n) = rhs(n) - qnm

  !       ! Derivative calculation and fill of n terms
  !       eps = get_perturbation(stage(n))
  !       qeps = this%qcalc(n, m, stage(n) + eps, stage(m), ii)
  !       derv = (qeps - qnm) / eps
  !       call matrix_sln%add_value_pos(idxglo(idiag), derv)
  !       rhs(n) = rhs(n) + derv * stage(n)

  !       ! Derivative calculation and fill of m terms
  !       eps = get_perturbation(stage(m))
  !       qeps = this%qcalc(n, m, stage(n), stage(m) + eps, ii)
  !       derv = (qeps - qnm) / eps
  !       call matrix_sln%add_value_pos(idxglo(ii), derv)
  !       rhs(n) = rhs(n) + derv * stage(m)

  !     end do
  !   end do

  ! end subroutine dfw_qnm_fc_nr

  ! !> @brief fill newton
  ! !<
  ! subroutine dfw_fn(this, kiter, matrix_sln, idxglo, rhs, stage)
  !   ! dummy
  !   class(SwfJncType) :: this !< this instance
  !   integer(I4B) :: kiter
  !   class(MatrixBaseType), pointer :: matrix_sln
  !   integer(I4B), intent(in), dimension(:) :: idxglo
  !   real(DP), intent(inout), dimension(:) :: rhs
  !   real(DP), intent(inout), dimension(:) :: stage
  !   ! local

  !   ! add newton terms to solution matrix
  !   ! todo: add newton terms here instead?
  !   ! this routine is probably not necessary as method is fully newton

  ! end subroutine dfw_fn

  ! !> @brief calculate flow between cells n and m
  ! !<
  ! function qcalc(this, n, m, stage_n, stage_m, ipos) result(qnm)
  !   ! dummy
  !   class(SwfJncType) :: this !< this instance
  !   integer(I4B), intent(in) :: n !< number for cell n
  !   integer(I4B), intent(in) :: m !< number for cell m
  !   real(DP), intent(in) :: stage_n !< stage in reach n
  !   real(DP), intent(in) :: stage_m !< stage in reach m
  !   integer(I4B), intent(in) :: ipos !< connection number
  !   ! local
  !   integer(I4B) :: isympos
  !   real(DP) :: qnm
  !   real(DP) :: cond
  !   real(DP) :: cl1
  !   real(DP) :: cl2

  !   ! Set connection lengths
  !   isympos = this%dis%con%jas(ipos)
  !   if (n < m) then
  !     cl1 = this%dis%con%cl1(isympos)
  !     cl2 = this%dis%con%cl2(isympos)
  !   else
  !     cl1 = this%dis%con%cl2(isympos)
  !     cl2 = this%dis%con%cl1(isympos)
  !   end if

  !   ! Calculate conductance
  !   if (this%iswrcond == 0) then
  !     cond = this%get_cond(n, m, ipos, stage_n, stage_m, cl1, cl2)
  !   else if (this%iswrcond == 1) then
  !     cond = this%get_cond_swr(n, m, ipos, stage_n, stage_m, cl1, cl2)
  !   end if

  !   ! calculate flow between n and m
  !   qnm = cond * (stage_m - stage_n)

  ! end function qcalc

  ! !> @brief calculate effective conductance between cells n and m
  ! !!
  ! !! Calculate half-cell conductances for cell n and cell m and then use
  ! !! harmonic averaging to calculate the effective conductance between the
  ! !< two cells.
  ! function get_cond(this, n, m, ipos, stage_n, stage_m, cln, clm) result(cond)
  !   ! modules
  !   use SmoothingModule, only: sQuadratic
  !   ! dummy
  !   class(SwfJncType) :: this !< this instance
  !   integer(I4B), intent(in) :: n !< number for cell n
  !   integer(I4B), intent(in) :: m !< number for cell m
  !   integer(I4B), intent(in) :: ipos !< connection number
  !   real(DP), intent(in) :: stage_n !< stage in reach n
  !   real(DP), intent(in) :: stage_m !< stage in reach m
  !   real(DP), intent(in) :: cln !< distance from cell n to shared face with m
  !   real(DP), intent(in) :: clm !< distance from cell m to shared face with n
  !   ! local
  !   real(DP) :: depth_n
  !   real(DP) :: depth_m
  !   real(DP) :: dhds_n
  !   real(DP) :: dhds_m
  !   real(DP) :: width_n
  !   real(DP) :: width_m
  !   real(DP) :: range = 1.d-6
  !   real(DP) :: dydx
  !   real(DP) :: smooth_factor
  !   real(DP) :: length_nm
  !   real(DP) :: cond
  !   real(DP) :: cn
  !   real(DP) :: cm

  !   ! we are using a harmonic conductance approach here; however
  !   ! the SWR Process for MODFLOW-2005/NWT uses length-weighted
  !   ! average areas and hydraulic radius instead.
  !   length_nm = cln + clm
  !   cond = DZERO
  !   if (length_nm > DPREC) then

  !     ! Calculate depth in each reach
  !     depth_n = stage_n - this%dis%bot(n)
  !     depth_m = stage_m - this%dis%bot(m)

  !     ! assign gradients
  !     if (this%is2d == 0) then
  !       dhds_n = abs(stage_m - stage_n) / (cln + clm)
  !       dhds_m = dhds_n
  !     else
  !       dhds_n = this%grad_dhds_mag(n)
  !       dhds_m = this%grad_dhds_mag(m)
  !     end if

  !     ! Assign upstream depth, if not central
  !     if (this%icentral == 0) then
  !       ! use upstream weighting
  !       if (stage_n > stage_m) then
  !         depth_m = depth_n
  !       else
  !         depth_n = depth_m
  !       end if
  !     end if

  !     ! Calculate a smoothed depth that goes to zero over
  !     ! the specified range
  !     call sQuadratic(depth_n, range, dydx, smooth_factor)
  !     depth_n = depth_n * smooth_factor
  !     call sQuadratic(depth_m, range, dydx, smooth_factor)
  !     depth_m = depth_m * smooth_factor

  !     ! Get the flow widths for n and m from dis package
  !     call this%dis%get_flow_width(n, m, ipos, width_n, width_m)

  !     ! Calculate half-cell conductance for reach
  !     ! n and m
  !     cn = this%get_cond_n(n, depth_n, cln, width_n, dhds_n)
  !     cm = this%get_cond_n(m, depth_m, clm, width_m, dhds_m)

  !     ! Use harmonic mean to calculate weighted
  !     ! conductance between the centers of reaches
  !     ! n and m
  !     if ((cn + cm) > DPREC) then
  !       cond = cn * cm / (cn + cm)
  !     else
  !       cond = DZERO
  !     end if

  !   end if

  ! end function get_cond

  ! !> @brief Calculate half cell conductance
  ! !!
  ! !! Calculate half-cell conductance for cell n
  ! !< using conveyance and Manning's equation
  ! function get_cond_n(this, n, depth, dx, width, dhds) result(c)
  !   ! modules
  !   ! dummy
  !   class(SwfJncType) :: this !< this instance
  !   integer(I4B), intent(in) :: n !< reach number
  !   real(DP), intent(in) :: depth !< simulated depth (stage - elevation) in reach n for this iteration
  !   real(DP), intent(in) :: dx !< half-cell distance
  !   real(DP), intent(in) :: width !< width of the reach perpendicular to flow
  !   real(DP), intent(in) :: dhds !< gradient
  !   ! return
  !   real(DP) :: c
  !   ! local
  !   real(DP) :: rough
  !   real(DP) :: dhds_sqr
  !   real(DP) :: conveyance

  !   ! Calculate conveyance, which is a * r**DTWOTHIRDS / roughc
  !   rough = this%manningsn(n)
  !   conveyance = this%cxs%get_conveyance(this%idcxs(n), width, depth, rough)
  !   dhds_sqr = dhds**DHALF
  !   if (dhds_sqr < DEM10) then
  !     dhds_sqr = DEM10
  !   end if

  !   ! Multiply by unitconv and divide conveyance by sqrt of friction slope and dx
  !   c = this%unitconv * conveyance / dx / dhds_sqr

  ! end function get_cond_n

  ! !> @brief Calculate effective conductance for cells n and m using SWR method
  ! !!
  ! !! The SWR Process for MODFLOW uses average cell parameters from cell n and
  ! !! m to calculate an effective conductance.  This is different from the
  ! !! default approach used in SWF, which uses harmonic averaging on two half-
  ! !< cell conductances.
  ! function get_cond_swr(this, n, m, ipos, stage_n, stage_m, cln, clm) result(cond)
  !   ! modules
  !   use SmoothingModule, only: sQuadratic
  !   ! dummy
  !   class(SwfJncType) :: this !< this instance
  !   integer(I4B), intent(in) :: n !< number for cell n
  !   integer(I4B), intent(in) :: m !< number for cell m
  !   integer(I4B), intent(in) :: ipos !< connection number
  !   real(DP), intent(in) :: stage_n !< stage in reach n
  !   real(DP), intent(in) :: stage_m !< stage in reach m
  !   real(DP), intent(in) :: cln !< distance from cell n to shared face with m
  !   real(DP), intent(in) :: clm !< distance from cell m to shared face with n
  !   ! local
  !   real(DP) :: depth_n
  !   real(DP) :: depth_m
  !   real(DP) :: dhds_n
  !   real(DP) :: dhds_m
  !   real(DP) :: dhds_nm
  !   real(DP) :: dhds_sqr
  !   real(DP) :: width_n
  !   real(DP) :: width_m
  !   real(DP) :: range = 1.d-6
  !   real(DP) :: dydx
  !   real(DP) :: smooth_factor
  !   real(DP) :: length_nm
  !   real(DP) :: cond
  !   real(DP) :: ravg
  !   real(DP) :: rinv_avg
  !   real(DP) :: area_n, area_m, area_avg
  !   real(DP) :: rhn, rhm, rhavg
  !   real(DP) :: weight_n
  !   real(DP) :: weight_m
  !   real(DP) :: rough_n
  !   real(DP) :: rough_m

  !   ! Use harmonic weighting for 1/manningsn, but using length-weighted
  !   ! averaging for other terms
  !   length_nm = cln + clm
  !   cond = DZERO
  !   if (length_nm > DPREC) then

  !     ! Calculate depth in each reach
  !     depth_n = stage_n - this%dis%bot(n)
  !     depth_m = stage_m - this%dis%bot(m)

  !     ! Assign upstream depth, if not central
  !     if (this%icentral == 0) then
  !       ! use upstream weighting
  !       if (stage_n > stage_m) then
  !         depth_m = depth_n
  !       else
  !         depth_n = depth_m
  !       end if
  !     end if

  !     ! Calculate a smoothed depth that goes to zero over
  !     !    the specified range
  !     call sQuadratic(depth_n, range, dydx, smooth_factor)
  !     depth_n = depth_n * smooth_factor
  !     call sQuadratic(depth_m, range, dydx, smooth_factor)
  !     depth_m = depth_m * smooth_factor

  !     ! Get the flow widths for n and m from dis package
  !     call this%dis%get_flow_width(n, m, ipos, width_n, width_m)

  !     ! linear weight toward node closer to shared face
  !     weight_n = clm / length_nm
  !     weight_m = DONE - weight_n

  !     ! average cross sectional flow area
  !     area_n = this%cxs%get_area(this%idcxs(n), width_n, depth_n)
  !     area_m = this%cxs%get_area(this%idcxs(m), width_m, depth_m)
  !     area_avg = weight_n * area_n + weight_m * area_m

  !     ! average hydraulic radius
  !     if (this%is2d == 0) then
  !       rhn = this%cxs%get_hydraulic_radius(this%idcxs(n), width_n, &
  !                                           depth_n, area_n)
  !       rhm = this%cxs%get_hydraulic_radius(this%idcxs(m), width_m, &
  !                                           depth_m, area_m)
  !       rhavg = weight_n * rhn + weight_m * rhm
  !     else
  !       rhavg = area_avg / width_n
  !     end if
  !     rhavg = rhavg**DTWOTHIRDS

  !     ! average gradient
  !     if (this%is2d == 0) then
  !       dhds_nm = abs(stage_m - stage_n) / (length_nm)
  !     else
  !       dhds_n = this%grad_dhds_mag(n)
  !       dhds_m = this%grad_dhds_mag(m)
  !       dhds_nm = weight_n * dhds_n + weight_m * dhds_m
  !     end if
  !     dhds_sqr = dhds_nm**DHALF
  !     if (dhds_sqr < DEM10) then
  !       dhds_sqr = DEM10
  !     end if

  !     ! weighted harmonic mean for inverse mannings value
  !     weight_n = cln / length_nm
  !     weight_m = DONE - weight_n
  !     rough_n = this%cxs%get_roughness(this%idcxs(n), width_n, depth_n, &
  !                                      this%manningsn(n), dhds_nm)
  !     rough_m = this%cxs%get_roughness(this%idcxs(m), width_m, depth_m, &
  !                                      this%manningsn(m), dhds_nm)
  !     ravg = (weight_n + weight_m) / &
  !            (weight_n / rough_n + weight_m / rough_m)
  !     rinv_avg = DONE / ravg

  !     ! calculate conductance using averaged values
  !     cond = this%unitconv * rinv_avg * area_avg * rhavg / dhds_sqr / length_nm

  !   end if

  ! end function get_cond_swr

  ! !> @brief Calculate flow area between cell n and m
  ! !!
  ! !! Calculate an average flow area between cell n and m.
  ! !! First calculate a flow area for cell n and then for
  ! !! cell m and linearly weight the areas using the connection
  ! !< distances.
  ! function get_flow_area_nm(this, n, m, stage_n, stage_m, cln, clm, &
  !                           ipos) result(area_avg)
  !   ! module
  !   use SmoothingModule, only: sQuadratic
  !   ! dummy
  !   class(SwfJncType) :: this !< this instance
  !   integer(I4B), intent(in) :: n
  !   integer(I4B), intent(in) :: m
  !   real(DP), intent(in) :: stage_n
  !   real(DP), intent(in) :: stage_m
  !   real(DP), intent(in) :: cln
  !   real(DP), intent(in) :: clm
  !   integer(I4B), intent(in) :: ipos
  !   ! local
  !   real(DP) :: depth_n
  !   real(DP) :: depth_m
  !   real(DP) :: width_n
  !   real(DP) :: width_m
  !   real(DP) :: area_n
  !   real(DP) :: area_m
  !   real(DP) :: weight_n
  !   real(DP) :: weight_m
  !   real(DP) :: length_nm
  !   real(DP) :: range = 1.d-6
  !   real(DP) :: dydx
  !   real(DP) :: smooth_factor
  !   ! return
  !   real(DP) :: area_avg

  !   ! depths
  !   depth_n = stage_n - this%dis%bot(n)
  !   depth_m = stage_m - this%dis%bot(m)

  !   ! Assign upstream depth, if not central
  !   if (this%icentral == 0) then
  !     ! use upstream weighting
  !     if (stage_n > stage_m) then
  !       depth_m = depth_n
  !     else
  !       depth_n = depth_m
  !     end if
  !   end if

  !   ! Calculate a smoothed depth that goes to zero over
  !   ! the specified range
  !   call sQuadratic(depth_n, range, dydx, smooth_factor)
  !   depth_n = depth_n * smooth_factor
  !   call sQuadratic(depth_m, range, dydx, smooth_factor)
  !   depth_m = depth_m * smooth_factor

  !   ! Get the flow widths for n and m from dis package
  !   call this%dis%get_flow_width(n, m, ipos, width_n, width_m)

  !   ! linear weight toward node closer to shared face
  !   length_nm = cln + clm
  !   weight_n = clm / length_nm
  !   weight_m = DONE - weight_n

  !   ! average cross sectional flow area
  !   area_n = this%cxs%get_area(this%idcxs(n), width_n, depth_n)
  !   area_m = this%cxs%get_area(this%idcxs(m), width_m, depth_m)
  !   area_avg = weight_n * area_n + weight_m * area_m

  ! end function get_flow_area_nm

  ! !> @brief Calculate average hydraulic gradient magnitude for each cell
  ! !!
  ! !! Go through each cell and calculate the average hydraulic gradient using
  ! !! an xt3d-style gradient interpolation.  This is used for 2D grids in the
  ! !! calculation of an effective conductance.  For 1D grids, gradients are
  ! !< calculated between cell centers.
  ! subroutine calc_dhds(this)
  !   ! modules
  !   use VectorInterpolationModule, only: vector_interpolation_2d
  !   ! dummy
  !   class(SwfJncType) :: this !< this instance
  !   ! local
  !   integer(I4B) :: n
  !   integer(I4B) :: m
  !   integer(I4B) :: ipos
  !   integer(I4B) :: isympos
  !   real(DP) :: cl1
  !   real(DP) :: cl2

  !   do n = 1, this%dis%nodes
  !     this%grad_dhds_mag(n) = DZERO
  !     do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
  !       m = this%dis%con%ja(ipos)
  !       isympos = this%dis%con%jas(ipos)

  !       ! determine cl1 and cl2
  !       if (n < m) then
  !         cl1 = this%dis%con%cl1(isympos)
  !         cl2 = this%dis%con%cl2(isympos)
  !       else
  !         cl1 = this%dis%con%cl2(isympos)
  !         cl2 = this%dis%con%cl1(isympos)
  !       end if

  !       ! store for n < m in upper right triangular part of symmetric dhdsja array
  !       if (n < m) then
  !         if (cl1 + cl2 > DPREC) then
  !           this%dhdsja(isympos) = (this%hnew(m) - this%hnew(n)) / (cl1 + cl2)
  !         else
  !           this%dhdsja(isympos) = DZERO
  !         end if
  !       end if
  !     end do
  !   end do

  !   ! pass dhdsja into the vector interpolation to get the components
  !   ! of the gradient at the cell center
  !   call vector_interpolation_2d(this%dis, this%dhdsja, vmag=this%grad_dhds_mag)

  ! end subroutine calc_dhds

  ! !> @ brief Newton under relaxation
  ! !!
  ! !<  If stage is below the bottom, then pull it up a bit
  ! subroutine dfw_nur(this, neqmod, x, xtemp, dx, inewtonur, dxmax, locmax)
  !   ! dummy
  !   class(SwfJncType) :: this !< this instance
  !   integer(I4B), intent(in) :: neqmod !< number of equations
  !   real(DP), dimension(neqmod), intent(inout) :: x !< dependent variable
  !   real(DP), dimension(neqmod), intent(in) :: xtemp !< temporary dependent variable
  !   real(DP), dimension(neqmod), intent(inout) :: dx !< change in dependent variable
  !   integer(I4B), intent(inout) :: inewtonur !< flag to indication relaxation was applied
  !   real(DP), intent(inout) :: dxmax !< max change in x
  !   integer(I4B), intent(inout) :: locmax !< location of max change
  !   ! local
  !   integer(I4B) :: n
  !   real(DP) :: botm
  !   real(DP) :: xx
  !   real(DP) :: dxx

  !   ! Newton-Raphson under-relaxation
  !   do n = 1, this%dis%nodes
  !     if (this%ibound(n) < 1) cycle
  !     if (this%icelltype(n) > 0) then
  !       botm = this%dis%bot(n)
  !       ! only apply Newton-Raphson under-relaxation if
  !       ! solution head is below the bottom of the model
  !       if (x(n) < botm) then
  !         inewtonur = 1
  !         xx = xtemp(n) * (DONE - DP9) + botm * DP9
  !         dxx = x(n) - xx
  !         if (abs(dxx) > abs(dxmax)) then
  !           locmax = n
  !           dxmax = dxx
  !         end if
  !         x(n) = xx
  !         dx(n) = DZERO
  !       end if
  !     end if
  !   end do

  ! end subroutine dfw_nur

  ! !> @ brief Calculate flow for each connection and store in flowja
  ! !<
  ! subroutine dfw_cq(this, stage, stage_old, flowja)
  !   ! dummy
  !   class(SwfJncType) :: this !< this instance
  !   real(DP), intent(inout), dimension(:) :: stage !< calculated head
  !   real(DP), intent(inout), dimension(:) :: stage_old !< calculated head from previous time step
  !   real(DP), intent(inout), dimension(:) :: flowja !< vector of flows in CSR format
  !   ! local
  !   integer(I4B) :: n, ipos, m
  !   real(DP) :: qnm

  !   do n = 1, this%dis%nodes
  !     do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
  !       m = this%dis%con%ja(ipos)
  !       if (m < n) cycle
  !       qnm = this%qcalc(n, m, stage(n), stage(m), ipos)
  !       flowja(ipos) = qnm
  !       flowja(this%dis%con%isym(ipos)) = -qnm
  !     end do
  !   end do

  ! end subroutine dfw_cq

  ! !> @ brief Model budget calculation for package
  ! !<
  ! subroutine dfw_bd(this, isuppress_output, model_budget)
  !   ! modules
  !   use BudgetModule, only: BudgetType
  !   ! dummy variables
  !   class(SwfJncType) :: this !< this instance
  !   integer(I4B), intent(in) :: isuppress_output !< flag to suppress model output
  !   type(BudgetType), intent(inout) :: model_budget !< model budget object

  !   ! Add any DFW budget terms
  !   ! none

  ! end subroutine dfw_bd

  ! !> @ brief save flows for package
  ! !<
  ! subroutine dfw_save_model_flows(this, flowja, icbcfl, icbcun)
  !   ! dummy
  !   class(SwfJncType) :: this !< this instance
  !   real(DP), dimension(:), intent(in) :: flowja !< vector of flows in CSR format
  !   integer(I4B), intent(in) :: icbcfl !< flag to indicate if flows should be saved
  !   integer(I4B), intent(in) :: icbcun !< unit number for flow output
  !   ! local
  !   integer(I4B) :: ibinun

  !   ! Set unit number for binary output
  !   if (this%ipakcb < 0) then
  !     ibinun = icbcun
  !   elseif (this%ipakcb == 0) then
  !     ibinun = 0
  !   else
  !     ibinun = this%ipakcb
  !   end if
  !   if (icbcfl == 0) ibinun = 0

  !   ! Write the face flows if requested
  !   if (ibinun /= 0) then
  !     ! flowja
  !     call this%dis%record_connection_array(flowja, ibinun, this%iout)
  !   end if

  !   ! Calculate velocities at cell centers and write, if requested
  !   if (this%isavvelocity /= 0) then
  !     if (ibinun /= 0) call this%sav_velocity(ibinun)
  !   end if

  ! end subroutine dfw_save_model_flows

  ! !> @ brief print flows for package
  ! !<
  ! subroutine dfw_print_model_flows(this, ibudfl, flowja)
  !   ! modules
  !   use TdisModule, only: kper, kstp
  !   use ConstantsModule, only: LENBIGLINE
  !   ! dummy
  !   class(SwfJncType) :: this !< this instance
  !   integer(I4B), intent(in) :: ibudfl !< print flag
  !   real(DP), intent(inout), dimension(:) :: flowja !< vector of flows in CSR format
  !   ! local
  !   character(len=LENBIGLINE) :: line
  !   character(len=30) :: tempstr
  !   integer(I4B) :: n, ipos, m
  !   real(DP) :: qnm
  !   ! formats
  !   character(len=*), parameter :: fmtiprflow = &
  !     &"(/,4x,'CALCULATED INTERCELL FLOW FOR PERIOD ', i0, ' STEP ', i0)"

  !   ! Write flowja to list file if requested
  !   if (ibudfl /= 0 .and. this%iprflow > 0) then
  !     write (this%iout, fmtiprflow) kper, kstp
  !     do n = 1, this%dis%nodes
  !       line = ''
  !       call this%dis%noder_to_string(n, tempstr)
  !       line = trim(tempstr)//':'
  !       do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
  !         m = this%dis%con%ja(ipos)
  !         call this%dis%noder_to_string(m, tempstr)
  !         line = trim(line)//' '//trim(tempstr)
  !         qnm = flowja(ipos)
  !         write (tempstr, '(1pg15.6)') qnm
  !         line = trim(line)//' '//trim(adjustl(tempstr))
  !       end do
  !       write (this%iout, '(a)') trim(line)
  !     end do
  !   end if

  ! end subroutine dfw_print_model_flows

  !> @brief deallocate memory
  !<
  subroutine jnc_da(this)
    ! modules
    use MemoryManagerModule, only: mem_deallocate
    use MemoryManagerExtModule, only: memorystore_remove
    use SimVariablesModule, only: idm_context
    ! dummy
    class(SwfJncType) :: this !< this instance

    ! Deallocate input memory
    call memorystore_remove(this%name_model, 'JNC', idm_context)

    ! Deallocate arrays
    call mem_deallocate(this%junc_ivert)
    call mem_deallocate(this%ivert_junc)
    call mem_deallocate(this%iajunction_cell)
    call mem_deallocate(this%jajunction_cell)
    call mem_deallocate(this%jstart)
    call mem_deallocate(this%jend)
    call mem_deallocate(this%junc_juncred)
    call mem_deallocate(this%juncred_junc)
    call mem_deallocate(this%irowqdn)
    call mem_deallocate(this%irowqup)
    call mem_deallocate(this%icellup)
    call mem_deallocate(this%icelldn)
  !   call mem_deallocate(this%manningsn)
  !   call mem_deallocate(this%idcxs)
  !   call mem_deallocate(this%icelltype)
  !   call mem_deallocate(this%nodedge)
  !   call mem_deallocate(this%ihcedge)
  !   call mem_deallocate(this%propsedge)
  !   call mem_deallocate(this%vcomp)
  !   call mem_deallocate(this%vmag)
  !   if (this%is2d == 1) then
  !     call mem_deallocate(this%grad_dhds_mag)
  !     call mem_deallocate(this%dhdsja)
  !   end if

    ! Scalars
    call mem_deallocate(this%moffset)
    call mem_deallocate(this%njunction)
    call mem_deallocate(this%njunc_red)
    call mem_deallocate(this%nvert)
    call mem_deallocate(this%nq)
  !   call mem_deallocate(this%is2d)
  !   call mem_deallocate(this%icentral)
  !   call mem_deallocate(this%iswrcond)
  !   call mem_deallocate(this%unitconv)
  !   call mem_deallocate(this%lengthconv)
  !   call mem_deallocate(this%timeconv)
  !   call mem_deallocate(this%isavvelocity)
  !   call mem_deallocate(this%icalcvelocity)
  !   call mem_deallocate(this%nedges)
  !   call mem_deallocate(this%lastedge)

  !   ! obs package
  !   call mem_deallocate(this%inobspkg)
  !   call this%obs%obs_da()
  !   deallocate (this%obs)
  !   nullify (this%obs)
  !   nullify (this%cxs)

  !   ! deallocate parent
  !   call this%NumericalPackageType%da()

  !   ! pointers
  !   this%hnew => null()

  end subroutine jnc_da

  ! !> @brief Calculate the 3 components of velocity at the cell center
  ! !!
  ! !! todo: duplicated from NPF; should consolidate
  ! !<
  ! subroutine calc_velocity(this, flowja)
  !   ! modules
  !   ! dummy
  !   class(SwfJncType) :: this !< this instance
  !   real(DP), intent(in), dimension(:) :: flowja !< vector of flows in CSR format
  !   ! local
  !   integer(I4B) :: n
  !   integer(I4B) :: m
  !   integer(I4B) :: ipos
  !   integer(I4B) :: isympos
  !   integer(I4B) :: ihc
  !   integer(I4B) :: ic
  !   integer(I4B) :: iz
  !   integer(I4B) :: nc
  !   integer(I4B) :: ncz
  !   real(DP) :: vx
  !   real(DP) :: vy
  !   real(DP) :: vz
  !   real(DP) :: xn
  !   real(DP) :: yn
  !   real(DP) :: zn
  !   real(DP) :: xc
  !   real(DP) :: yc
  !   real(DP) :: zc
  !   real(DP) :: cl1
  !   real(DP) :: cl2
  !   real(DP) :: dltot
  !   real(DP) :: ooclsum
  !   real(DP) :: dsumx
  !   real(DP) :: dsumy
  !   real(DP) :: dsumz
  !   real(DP) :: denom
  !   real(DP) :: area
  !   real(DP) :: axy
  !   real(DP) :: ayx
  !   real(DP), allocatable, dimension(:) :: vi
  !   real(DP), allocatable, dimension(:) :: di
  !   real(DP), allocatable, dimension(:) :: viz
  !   real(DP), allocatable, dimension(:) :: diz
  !   real(DP), allocatable, dimension(:) :: nix
  !   real(DP), allocatable, dimension(:) :: niy
  !   real(DP), allocatable, dimension(:) :: wix
  !   real(DP), allocatable, dimension(:) :: wiy
  !   real(DP), allocatable, dimension(:) :: wiz
  !   real(DP), allocatable, dimension(:) :: bix
  !   real(DP), allocatable, dimension(:) :: biy
  !   logical :: nozee = .true.

  !   ! Ensure dis has necessary information
  !   ! todo: do we need this for SWF?
  !   if (this%icalcvelocity /= 0 .and. this%dis%con%ianglex == 0) then
  !     call store_error('Error.  ANGLDEGX not provided in '// &
  !                      'discretization file.  ANGLDEGX required for '// &
  !                      'calculation of velocity.', terminate=.TRUE.)
  !   end if

  !   ! Find max number of connections and allocate weight arrays
  !   nc = 0
  !   do n = 1, this%dis%nodes

  !     ! Count internal model connections
  !     ic = this%dis%con%ia(n + 1) - this%dis%con%ia(n) - 1

  !     ! Count edge connections
  !     do m = 1, this%nedges
  !       if (this%nodedge(m) == n) then
  !         ic = ic + 1
  !       end if
  !     end do

  !     ! Set max number of connections for any cell
  !     if (ic > nc) nc = ic
  !   end do

  !   ! Allocate storage arrays needed for cell-centered calculation
  !   allocate (vi(nc))
  !   allocate (di(nc))
  !   allocate (viz(nc))
  !   allocate (diz(nc))
  !   allocate (nix(nc))
  !   allocate (niy(nc))
  !   allocate (wix(nc))
  !   allocate (wiy(nc))
  !   allocate (wiz(nc))
  !   allocate (bix(nc))
  !   allocate (biy(nc))

  !   ! Go through each cell and calculate specific discharge
  !   do n = 1, this%dis%nodes

  !     ! first calculate geometric properties for x and y directions and
  !     !    the specific discharge at a face (vi)
  !     ic = 0
  !     iz = 0
  !     vi(:) = DZERO
  !     di(:) = DZERO
  !     viz(:) = DZERO
  !     diz(:) = DZERO
  !     nix(:) = DZERO
  !     niy(:) = DZERO
  !     do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
  !       m = this%dis%con%ja(ipos)
  !       isympos = this%dis%con%jas(ipos)
  !       ihc = this%dis%con%ihc(isympos)
  !       ic = ic + 1
  !       call this%dis%connection_normal(n, m, ihc, xn, yn, zn, ipos)
  !       call this%dis%connection_vector(n, m, nozee, DONE, DONE, &
  !                                       ihc, xc, yc, zc, dltot)
  !       cl1 = this%dis%con%cl1(isympos)
  !       cl2 = this%dis%con%cl2(isympos)
  !       if (m < n) then
  !         cl1 = this%dis%con%cl2(isympos)
  !         cl2 = this%dis%con%cl1(isympos)
  !       end if
  !       ooclsum = DONE / (cl1 + cl2)
  !       nix(ic) = -xn
  !       niy(ic) = -yn
  !       di(ic) = dltot * cl1 * ooclsum
  !       area = this%get_flow_area_nm(n, m, this%hnew(n), this%hnew(m), &
  !                                    cl1, cl2, ipos)
  !       if (area > DZERO) then
  !         vi(ic) = flowja(ipos) / area
  !       else
  !         vi(ic) = DZERO
  !       end if

  !     end do

  !     ! Look through edge flows that may have been provided by an exchange
  !     ! and incorporate them into the averaging arrays
  !     do m = 1, this%nedges
  !       if (this%nodedge(m) == n) then

  !         ! propsedge: (Q, area, nx, ny, distance)
  !         ihc = this%ihcedge(m)
  !         area = this%propsedge(2, m)

  !         ic = ic + 1
  !         nix(ic) = -this%propsedge(3, m)
  !         niy(ic) = -this%propsedge(4, m)
  !         di(ic) = this%propsedge(5, m)
  !         if (area > DZERO) then
  !           vi(ic) = this%propsedge(1, m) / area
  !         else
  !           vi(ic) = DZERO
  !         end if

  !       end if
  !     end do

  !     ! Assign number of vertical and horizontal connections
  !     ncz = iz
  !     nc = ic

  !     ! calculate z weight (wiz) and z velocity
  !     if (ncz == 1) then
  !       wiz(1) = DONE
  !     else
  !       dsumz = DZERO
  !       do iz = 1, ncz
  !         dsumz = dsumz + diz(iz)
  !       end do
  !       denom = (ncz - DONE)
  !       if (denom < DZERO) denom = DZERO
  !       dsumz = dsumz + DEM10 * dsumz
  !       do iz = 1, ncz
  !         if (dsumz > DZERO) wiz(iz) = DONE - diz(iz) / dsumz
  !         if (denom > 0) then
  !           wiz(iz) = wiz(iz) / denom
  !         else
  !           wiz(iz) = DZERO
  !         end if
  !       end do
  !     end if
  !     vz = DZERO
  !     do iz = 1, ncz
  !       vz = vz + wiz(iz) * viz(iz)
  !     end do

  !     ! distance-based weighting
  !     nc = ic
  !     dsumx = DZERO
  !     dsumy = DZERO
  !     dsumz = DZERO
  !     do ic = 1, nc
  !       wix(ic) = di(ic) * abs(nix(ic))
  !       wiy(ic) = di(ic) * abs(niy(ic))
  !       dsumx = dsumx + wix(ic)
  !       dsumy = dsumy + wiy(ic)
  !     end do

  !     ! Finish computing omega weights.  Add a tiny bit
  !     ! to dsum so that the normalized omega weight later
  !     ! evaluates to (essentially) 1 in the case of a single
  !     ! relevant connection, avoiding 0/0.
  !     dsumx = dsumx + DEM10 * dsumx
  !     dsumy = dsumy + DEM10 * dsumy
  !     do ic = 1, nc
  !       wix(ic) = (dsumx - wix(ic)) * abs(nix(ic))
  !       wiy(ic) = (dsumy - wiy(ic)) * abs(niy(ic))
  !     end do

  !     ! compute B weights
  !     dsumx = DZERO
  !     dsumy = DZERO
  !     do ic = 1, nc
  !       bix(ic) = wix(ic) * sign(DONE, nix(ic))
  !       biy(ic) = wiy(ic) * sign(DONE, niy(ic))
  !       dsumx = dsumx + wix(ic) * abs(nix(ic))
  !       dsumy = dsumy + wiy(ic) * abs(niy(ic))
  !     end do
  !     if (dsumx > DZERO) dsumx = DONE / dsumx
  !     if (dsumy > DZERO) dsumy = DONE / dsumy
  !     axy = DZERO
  !     ayx = DZERO
  !     do ic = 1, nc
  !       bix(ic) = bix(ic) * dsumx
  !       biy(ic) = biy(ic) * dsumy
  !       axy = axy + bix(ic) * niy(ic)
  !       ayx = ayx + biy(ic) * nix(ic)
  !     end do

  !     ! Calculate specific discharge.  The divide by zero checking below
  !     ! is problematic for cells with only one flow, such as can happen
  !     ! with triangular cells in corners.  In this case, the resulting
  !     ! cell velocity will be calculated as zero.  The method should be
  !     ! improved so that edge flows of zero are included in these
  !     ! calculations.  But this needs to be done with consideration for LGR
  !     ! cases in which flows are submitted from an exchange.
  !     vx = DZERO
  !     vy = DZERO
  !     do ic = 1, nc
  !       vx = vx + (bix(ic) - axy * biy(ic)) * vi(ic)
  !       vy = vy + (biy(ic) - ayx * bix(ic)) * vi(ic)
  !     end do
  !     denom = DONE - axy * ayx
  !     if (denom /= DZERO) then
  !       vx = vx / denom
  !       vy = vy / denom
  !     end if

  !     this%vcomp(1, n) = vx
  !     this%vcomp(2, n) = vy
  !     this%vcomp(3, n) = vz
  !     this%vmag(n) = sqrt(vx**2 + vy**2 + vz**2)

  !   end do

  !   ! cleanup
  !   deallocate (vi)
  !   deallocate (di)
  !   deallocate (nix)
  !   deallocate (niy)
  !   deallocate (wix)
  !   deallocate (wiy)
  !   deallocate (wiz)
  !   deallocate (bix)
  !   deallocate (biy)

  ! end subroutine calc_velocity

  ! !> @brief Reserve space for nedges cells that have an edge on them.
  ! !!
  ! !! todo: duplicated from NPF; should consolidate
  ! !! This must be called before the swf%allocate_arrays routine, which is
  ! !< called from swf%ar.
  ! subroutine increase_edge_count(this, nedges)
  !   ! dummy
  !   class(SwfJncType) :: this !< this instance
  !   integer(I4B), intent(in) :: nedges

  !   this%nedges = this%nedges + nedges

  ! end subroutine increase_edge_count

  ! !> @brief Provide the swf package with edge properties
  ! !!
  ! !! todo: duplicated from NPF; should consolidate
  ! !<
  ! subroutine set_edge_properties(this, nodedge, ihcedge, q, area, nx, ny, &
  !                                distance)
  !   ! dummy
  !   class(SwfJncType) :: this !< this instance
  !   integer(I4B), intent(in) :: nodedge
  !   integer(I4B), intent(in) :: ihcedge
  !   real(DP), intent(in) :: q
  !   real(DP), intent(in) :: area
  !   real(DP), intent(in) :: nx
  !   real(DP), intent(in) :: ny
  !   real(DP), intent(in) :: distance
  !   ! local
  !   integer(I4B) :: lastedge

  !   this%lastedge = this%lastedge + 1
  !   lastedge = this%lastedge
  !   this%nodedge(lastedge) = nodedge
  !   this%ihcedge(lastedge) = ihcedge
  !   this%propsedge(1, lastedge) = q
  !   this%propsedge(2, lastedge) = area
  !   this%propsedge(3, lastedge) = nx
  !   this%propsedge(4, lastedge) = ny
  !   this%propsedge(5, lastedge) = distance

  !   ! If this is the last edge, then the next call must be starting a new
  !   ! edge properties assignment loop, so need to reset lastedge to 0
  !   if (this%lastedge == this%nedges) this%lastedge = 0

  ! end subroutine set_edge_properties

  ! !> @brief Save specific discharge in binary format to ibinun
  ! !!
  ! !! todo: should write 2d velocity; what about for 1D channel?
  ! !<
  ! subroutine sav_velocity(this, ibinun)
  !   ! dummy
  !   class(SwfJncType) :: this !< this instance
  !   integer(I4B), intent(in) :: ibinun
  !   ! local
  !   character(len=16) :: text
  !   character(len=16), dimension(3) :: auxtxt
  !   integer(I4B) :: n
  !   integer(I4B) :: naux

  !   ! Write the header
  !   text = '      DATA-VCOMP'
  !   naux = 3
  !   auxtxt(:) = ['              vx', '              vy', '              vz']
  !   call this%dis%record_srcdst_list_header(text, this%name_model, &
  !                                           this%packName, this%name_model, &
  !                                           this%packName, naux, auxtxt, ibinun, &
  !                                           this%dis%nodes, this%iout)

  !   ! Write a zero for Q, and then write qx, qy, qz as aux variables
  !   do n = 1, this%dis%nodes
  !     call this%dis%record_mf6_list_entry(ibinun, n, n, DZERO, naux, &
  !                                         this%vcomp(:, n))
  !   end do

  ! end subroutine sav_velocity

  ! !> @brief Define the observation types available in the package
  ! !!
  ! !< Method to define the observation types available in the package.
  ! subroutine dfw_df_obs(this)
  !   ! dummy variables
  !   class(SwfJncType) :: this !< this instance
  !   ! local variables
  !   integer(I4B) :: indx

  !   ! Store obs type and assign procedure pointer
  !   !    for ext-outflow observation type.
  !   call this%obs%StoreObsType('ext-outflow', .true., indx)
  !   this%obs%obsData(indx)%ProcessIdPtr => dfwobsidprocessor

  ! end subroutine dfw_df_obs

  ! subroutine dfwobsidprocessor(obsrv, dis, inunitobs, iout)
  !   ! dummy
  !   type(ObserveType), intent(inout) :: obsrv
  !   class(DisBaseType), intent(in) :: dis
  !   integer(I4B), intent(in) :: inunitobs
  !   integer(I4B), intent(in) :: iout
  !   ! local
  !   integer(I4B) :: n
  !   character(len=LINELENGTH) :: string

  !   ! Initialize variables
  !   string = obsrv%IDstring
  !   read (string, *) n

  !   if (n > 0) then
  !     obsrv%NodeNumber = n
  !   else
  !     errmsg = 'Error reading data from ID string'
  !     call store_error(errmsg)
  !     call store_error_unit(inunitobs)
  !   end if

  ! end subroutine dfwobsidprocessor

  ! !> @brief Save observations for the package
  ! !!
  ! !< Method to save simulated values for the package.
  ! subroutine dfw_bd_obs(this)
  !   ! dummy variables
  !   class(SwfJncType) :: this !< this instance
  !   ! local variables
  !   integer(I4B) :: i
  !   integer(I4B) :: j
  !   integer(I4B) :: n
  !   real(DP) :: v
  !   character(len=100) :: msg
  !   type(ObserveType), pointer :: obsrv => null()

  !   ! Write simulated values for all observations
  !   if (this%obs%npakobs > 0) then
  !     call this%obs%obs_bd_clear()
  !     do i = 1, this%obs%npakobs
  !       obsrv => this%obs%pakobs(i)%obsrv
  !       do j = 1, obsrv%indxbnds_count
  !         n = obsrv%indxbnds(j)
  !         v = DZERO
  !         select case (obsrv%ObsTypeId)
  !         case default
  !           msg = 'Unrecognized observation type: '//trim(obsrv%ObsTypeId)
  !           call store_error(msg)
  !         end select
  !         call this%obs%SaveOneSimval(obsrv, v)
  !       end do
  !     end do

  !     ! write summary of package error messages
  !     if (count_errors() > 0) then
  !       call store_error_filename(this%input_fname)
  !     end if
  !   end if

  ! end subroutine dfw_bd_obs

  ! !> @brief Read and prepare observations for a package
  ! !!
  ! !< Method to read and prepare observations for a package.
  ! subroutine dfw_rp_obs(this)
  !   ! modules
  !   use TdisModule, only: kper
  !   ! dummy
  !   class(SwfJncType), intent(inout) :: this !< this instance
  !   ! local
  !   integer(I4B) :: i
  !   integer(I4B) :: j
  !   integer(I4B) :: nn1
  !   class(ObserveType), pointer :: obsrv => null()

  !   ! process each package observation
  !   ! only done the first stress period since boundaries are fixed
  !   ! for the simulation
  !   if (kper == 1) then
  !     do i = 1, this%obs%npakobs
  !       obsrv => this%obs%pakobs(i)%obsrv

  !       ! get node number 1
  !       nn1 = obsrv%NodeNumber
  !       if (nn1 < 1 .or. nn1 > this%dis%nodes) then
  !         write (errmsg, '(a,1x,a,1x,i0,1x,a,1x,i0,a)') &
  !           trim(adjustl(obsrv%ObsTypeId)), &
  !           'reach must be greater than 0 and less than or equal to', &
  !           this%dis%nodes, '(specified value is ', nn1, ')'
  !         call store_error(errmsg)
  !       else
  !         if (obsrv%indxbnds_count == 0) then
  !           call obsrv%AddObsIndex(nn1)
  !         else
  !           errmsg = 'Programming error in dfw_rp_obs'
  !           call store_error(errmsg)
  !         end if
  !       end if

  !       ! check that node number 1 is valid; call store_error if not
  !       do j = 1, obsrv%indxbnds_count
  !         nn1 = obsrv%indxbnds(j)
  !         if (nn1 < 1 .or. nn1 > this%dis%nodes) then
  !           write (errmsg, '(a,1x,a,1x,i0,1x,a,1x,i0,a)') &
  !             trim(adjustl(obsrv%ObsTypeId)), &
  !             'reach must be greater than 0 and less than or equal to', &
  !             this%dis%nodes, '(specified value is ', nn1, ')'
  !           call store_error(errmsg)
  !         end if
  !       end do
  !     end do

  !     ! evaluate if there are any observation errors
  !     if (count_errors() > 0) then
  !       call store_error_filename(this%input_fname)
  !     end if

  !   end if

  ! end subroutine dfw_rp_obs

  !> @ brief Count number of junctions
  !!
  !! Count the number of unique vertices that are starting or ending
  !! vertices for a cell.
  subroutine count_junctions(nvert, iavert, javert, njunc_tot, njunc_red)
    ! dummy
    integer(I4B), intent(in) :: nvert
    integer(I4B), dimension(:), intent(in) :: iavert !< vertex index array
    integer(I4B), dimension(:), intent(in) :: javert !< vertex csr array
    integer(I4B), intent(inout) :: njunc_tot !< total number of junctions
    integer(I4B), intent(inout) :: njunc_red !< number of junctions without 2 connected cells
    ! local
    integer(I4B) :: ic
    integer(I4B) :: iv
    integer(I4B), dimension(:), allocatable :: icount

    ! allocate memory for counting
    allocate(icount(nvert))
    do iv = 1, nvert
      icount(iv) = 0
    end do

    ! count the number of unique junctions
    njunc_tot = 0
    do ic = 1, size(iavert) - 1
      ! starting cell vertex
      iv = javert(iavert(ic))
      if (icount(iv) == 0) njunc_tot = njunc_tot + 1
      icount(iv) = icount(iv) + 1
 
      ! ending cell vertex
      iv = javert(iavert(ic + 1) - 1)
      if (icount(iv) == 0) njunc_tot = njunc_tot + 1
      icount(iv) = icount(iv) + 1
    end do

    ! count number of junctions without 2 connected cells
    njunc_red = 0
    do iv = 1, nvert
      if (icount(iv) /= 2) njunc_red = njunc_red + 1
    end do

    ! deallocate memory
    deallocate(icount)
    
  end subroutine count_junctions

  subroutine fill_junc_ivert(nvert, iavert, javert, junc_ivert, ivert_junc)
    ! dummy
    integer(I4B), intent(in) :: nvert
    integer(I4B), dimension(:), intent(in) :: iavert
    integer(I4B), dimension(:), intent(in) :: javert
    integer(I4B), dimension(:), intent(inout) :: junc_ivert
    integer(I4B), dimension(:), intent(inout) :: ivert_junc
    ! local
    integer(I4B), dimension(:), allocatable :: icount
    integer(I4B) :: ic
    integer(I4B) :: iv
    integer(I4B) :: ijunction

    ! allocate memory for counting
    allocate(icount(nvert))
    do iv = 1, nvert
      icount(iv) = 0
    end do

    ! Save the vertex number for each junction in junc_ivert
    ijunction = 1
    do ic = 1, size(iavert) - 1
      ! starting cell vertex
      iv = javert(iavert(ic))
      if (icount(iv) == 0) then
        junc_ivert(ijunction) = iv
        ivert_junc(iv) = ijunction
        ijunction = ijunction + 1
      end if
      icount(iv) = icount(iv) + 1
 
      ! ending cell vertex
      iv = javert(iavert(ic + 1) - 1)
      if (icount(iv) == 0) then
        junc_ivert(ijunction) = iv
        ivert_junc(iv) = ijunction
        ijunction = ijunction + 1
      end if
      icount(iv) = icount(iv) + 1
    end do
    deallocate(icount)

  end subroutine fill_junc_ivert

  !> @brief 
  !<
  subroutine fill_junction_cell(njunction, iavert, javert, ivert_junc, &
                                iajunction_cell, jajunction_cell, mempath)
    ! modules
    use ConstantsModule, only: LENMEMPATH
    use SparseModule, only: sparsematrix
    ! dummy
    integer(I4B), intent(in) :: njunction
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: iavert
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: javert
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: ivert_junc
    integer(I4B), dimension(:), pointer, contiguous, intent(inout) :: iajunction_cell
    integer(I4B), dimension(:), pointer, contiguous, intent(inout) :: jajunction_cell
    character(len=LENMEMPATH), intent(in) :: mempath
    ! locals
    type(sparsematrix) :: spm
    integer(I4B) :: nodesuser
    integer(I4B) :: ijunction
    integer(I4B) :: n
    integer(I4B) :: ivert
    integer(I4B) :: ipos
    integer(I4B) :: ierr
    integer(I4B) :: maxnnz = 2

    ! Create sparse matrix with rows of junctions and columns of user nodes
    nodesuser = size(iavert) - 1
    call spm%init(njunction, nodesuser, maxnnz)
    do n = 1, nodesuser
      do ipos = iavert(n), iavert(n + 1) - 1
        ivert = javert(ipos)
        ijunction = ivert_junc(ivert)
        call spm%addconnection(ijunction, n, 1)
      end do
    end do

    ! allocate memory and fill csr vectors
    call mem_allocate(iajunction_cell, njunction + 1, 'IAJUNCTION_CELL', &
                      mempath)
    call mem_allocate(jajunction_cell, spm%nnz, 'JAJUNCTION_CELL', mempath)
    call spm%filliaja(iajunction_cell, jajunction_cell, ierr)
    call spm%destroy()

  end subroutine fill_junction_cell

  subroutine fill_jstart_jend(nvert, iavert, javert, ivert_junc, jstart, jend)
    ! dummy
    integer(I4B), intent(in) :: nvert
    integer(I4B), dimension(:), intent(in) :: iavert
    integer(I4B), dimension(:), intent(in) :: javert
    integer(I4B), dimension(:), intent(in) :: ivert_junc
    integer(I4B), dimension(:), intent(inout) :: jstart
    integer(I4B), dimension(:), intent(inout) :: jend
    ! local
    integer(I4B) :: ic
    integer(I4B) :: iv
    integer(I4B) :: ijunction

    ! Save the starting and ending junction number for each cell
    do ic = 1, size(iavert) - 1
      ! starting cell vertex
      iv = javert(iavert(ic))
      ijunction = ivert_junc(iv)
      jstart(ic) = ijunction
 
      ! ending cell vertex
      iv = javert(iavert(ic + 1) - 1)
      ijunction = ivert_junc(iv)
      jend(ic) = ijunction
    end do

  end subroutine fill_jstart_jend

  subroutine fill_junc_juncred(iajunction_cell, jajunction_cell, junc_juncred, &
                               juncred_junc)
    ! dummy
    integer(I4B), dimension(:), intent(in) :: iajunction_cell !< csr junction to cell index array
    integer(I4B), dimension(:), intent(in) :: jajunction_cell !< csr junction to cell array
    integer(I4B), dimension(:), intent(inout) :: junc_juncred !< map from junction to reduced junction
    integer(I4B), dimension(:), intent(inout) :: juncred_junc !< map from reduced junction to junction
    ! local
    integer(I4B) :: j, jred
    integer(I4B) :: n_connected_cells

    ! Loop through junctions
    jred = 0
    do j = 1, size(iajunction_cell) - 1

      ! if the number of connected cells is not 2 (the common case) then
      ! this is not an active junction, and jred will be zero
      n_connected_cells = iajunction_cell(j + 1) - iajunction_cell(j)
      if (n_connected_cells /= 2) then
        jred = jred + 1

        ! store maps between the two
        junc_juncred(j) = jred
        juncred_junc(jred) = j
      end if

    end do

  end subroutine fill_junc_juncred

  subroutine fill_irow_qupdn(jstart, jend, junc_jred, icellup, nq, irowqup, &
                             irowqdn)
    ! dummy
    integer(I4B), dimension(:), intent(in) :: jstart !< junction for up side
    integer(I4B), dimension(:), intent(in) :: jend !< junction for down side
    integer(I4B), dimension(:), intent(in) :: junc_jred !< junction for down side
    integer(I4B), dimension(:), intent(in) :: icellup !< cell number of up side cell
    integer(I4B), intent(inout) :: nq !< number of unique q equations
    integer(I4B), dimension(:), intent(inout) :: irowqup !< map from cell number to upstream q row in equations
    integer(I4B), dimension(:), intent(inout) :: irowqdn !< map from reduced junction to junction
    ! local
    integer(I4B) :: j, jred
    integer(I4B) :: icell
    integer(I4B) :: iq

    ! Loop through all cells and reserve a Q row for down side
    iq = 1
    do icell = 1, size(jstart)
      ! dn side
      irowqdn(icell) = iq
      iq = iq + 1
    end do

    ! Loop through all cells and find the row for the q up side
    do icell = 1, size(jstart)
      ! up side
      j = jstart(icell)
      jred = junc_jred(j)
      if (jred /= 0) then
        ! there is an active junction on the up side so there will be a row
        ! storing the upside q
        irowqup(icell) = iq
        iq = iq + 1
      else
        ! there is no active junction on upside, so qup(icell) = -qdn(icellup(icell))
        irowqup(icell) = irowqdn(icellup(icell))
      end if
    end do
    nq = iq - 1

  end subroutine fill_irow_qupdn

  subroutine fill_icellupdn(jend, iajunction_cell, jajunction_cell, icellup, &
                            icelldn)
    ! dummy
    integer(I4B), dimension(:), intent(in) :: jend !< down side junction for each cell
    integer(I4B), dimension(:), intent(in) :: iajunction_cell !< csr index for junction to cell
    integer(I4B), dimension(:), intent(in) :: jajunction_cell !< csr junction to cell
    integer(I4B), dimension(:), intent(inout) :: icellup !< index of up side cell; 0 if active junction
    integer(I4B), dimension(:), intent(inout) :: icelldn !< index of dn side cell; 0 if active junction
    ! local
    integer(I4B) :: ijunc, ncon, ic1, ic2
    do ijunc = 1, size(iajunction_cell) - 1
      ncon = iajunction_cell(ijunc + 1) - iajunction_cell(ijunc)
      if (ncon == 2) then
        ! evaluate two cells and assign icellup
        ic1 = jajunction_cell(iajunction_cell(ijunc))
        ic2 = jajunction_cell(iajunction_cell(ijunc) + 1)

        ! ic2 is down from ic1, so icellup(ic2) = ic1
        if(jend(ic1) == ijunc) then
          icellup(ic2) = ic1
          icelldn(ic1) = ic2
        end if

        ! ic1 is down from ic2, so icellup(ic1) = ic2
        if (jend(ic2) == ijunc) then
          icellup(ic1) = ic2
          icelldn(ic2) = ic1
        end if
        
      end if
    end do
  end subroutine fill_icellupdn

end module SwfJncModule
