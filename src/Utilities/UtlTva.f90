!> @brief This module contains the UtlTva Module
!!
!! This module contains the code for reading and storing a
!! generic input file of time varying grid data
!<
module UtlTvaModule

  use KindModule, only: DP, LGP, I4B
  use ConstantsModule, only: LENPACKAGENAME, LENMODELNAME, &
                             LENMEMPATH, DZERO, LENFTYPE, &
                             LINELENGTH, TABLEFT, TABCENTER, &
                             LENAUXNAME
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error, count_errors
  use MemoryHelperModule, only: create_mem_path
  use MemoryManagerModule, only: mem_allocate, mem_reallocate, mem_deallocate
  use BlockParserModule, only: BlockParserType
  use BaseDisModule, only: DisBaseType
  use TimeArraySeriesManagerModule, only: TimeArraySeriesManagerType, &
                                          tasmanager_cr
  use InputOutputModule, only: str_pad_left, urdaux

  implicit none
  private
  public :: UtlTvaType

  character(len=LENFTYPE) :: ftype = 'TVA'
  character(len=LENPACKAGENAME) :: text = 'TIME VARYING ARR'

  !> @brief Derived type for managing TVA input
  !!
  !! This derived type will read and process an TVA input file,
  !! make time series interpolations, and provide concentrations or
  !! temperatures to the SSM package that correspond to an individual
  !! GWF stress package.
  !<
  type :: VariableType
    real(DP), dimension(:), pointer, contiguous :: dblvec => null() !< vector of floats
  end type VariableType

  type :: UtlTvaType

    character(len=LENMODELNAME) :: component = '' !< the name of the parent component (typically model name)
    character(len=LENPACKAGENAME) :: subcomponent = '' !< name of the subcomponent (typically package name)
    character(len=LENAUXNAME), dimension(:), pointer, &
      contiguous :: auxname => null() !< vector of auxiliary variable names
    character(len=LENMEMPATH) :: memoryPath = '' !< the location in the memory manager where the variables are stored
    integer(I4B), pointer :: inunit => null() !< unit number for input
    integer(I4B), pointer :: iout => null() !< unit number for output
    integer(I4B), pointer :: nvar => null() !< number of variables
    integer(I4B), pointer :: isize => null() !< length of dblvec
    integer(I4B), pointer :: ionper => null() !< stress period for next data
    integer(I4B), pointer :: lastonper => null() !< last value of ionper (for checking)
    integer(I4B), pointer :: iprpak => null() !< flag for printing input
    class(DisBaseType), pointer :: dis => null() !< model discretization object
    type(BlockParserType) :: parser !< parser object for reading blocks of information
    type(TimeArraySeriesManagerType), pointer :: TasManager => null() !< time array series manager

    type(VariableType), dimension(:), allocatable :: variables

  contains

    procedure :: initialize
    procedure :: has
    procedure :: fill
    procedure :: tva_rp
    procedure :: tva_ad
    procedure :: tva_da
    procedure, private :: tva_rp_array
    procedure, private :: allocate_scalars
    procedure, private :: read_options
    procedure, private :: allocate_arrays
    procedure, private :: read_check_ionper

  end type UtlTvaType

contains

  !> @ brief Initialize the TVA type
  !!
  !! Initialize the TVA object by setting up the parser,
  !! and time series manager, reading options and dimensions,
  !< and allocating memory.
  subroutine initialize(this, dis, isize, inunit, iout, component, &
                        subcomponent)
    ! dummy
    class(UtlTvaType) :: this !<  UtlTvaType
    class(DisBaseType), pointer, intent(in) :: dis !<  discretization package
    integer(I4B), intent(in) :: isize !<  size of array
    integer(I4B), intent(in) :: inunit !<  unit number for input
    integer(I4B), intent(in) :: iout !<  unit number for output
    character(len=*), intent(in) :: component !<  character string containing component
    character(len=*), intent(in) :: subcomponent !<  character string containing subcomponent
    ! local

    ! construct the memory path
    this%component = component
    this%memoryPath = create_mem_path(this%component, this%subcomponent)

    ! allocate scalar variables
    call this%allocate_scalars()

    ! assign member values
    this%inunit = inunit
    this%iout = iout

    ! set pointers
    this%dis => dis

    ! Setup the parser
    call this%parser%Initialize(this%inunit, this%iout)

    ! Setup the time series manager
    call tasmanager_cr(this%TasManager, dis, component, this%iout)

    ! read options
    call this%read_options()

    ! set dimensions
    this%isize = isize
    
    ! allocate arrays
    call this%allocate_arrays()

    ! Now that time series are read, call define
    call this%tasmanager%tasmanager_df()
  end subroutine initialize

  !> @ Check to see if variable was included by user
  function has(this, auxname) result(found)
    ! dummy variables
    class(UtlTvaType) :: this !< UtlTvaType object
    character(len=*), intent(in) :: auxname
    ! return
    logical(LGP) :: found
    ! local
    integer(I4B) :: i

    found = .false.
    do i = 1, this%nvar
      if (this%auxname(i) == auxname) then
        found = .true.
        exit
      end if
    end do
  end function has

  !> @ brief Allocate package scalars
  !!
  !<  Allocate and initialize package scalars.
  subroutine allocate_scalars(this)
    ! modules
    ! dummy variables
    class(UtlTvaType) :: this !< UtlTvaType object

    ! allocate scalars in memory manager
    call mem_allocate(this%inunit, 'INUNIT', this%memoryPath)
    call mem_allocate(this%iout, 'IOUT', this%memoryPath)
    call mem_allocate(this%nvar, 'NVAR', this%memoryPath)
    call mem_allocate(this%isize, 'ISIZE', this%memoryPath)
    call mem_allocate(this%ionper, 'IONPER', this%memoryPath)
    call mem_allocate(this%lastonper, 'LASTONPER', this%memoryPath)
    call mem_allocate(this%iprpak, 'IPRPAK', this%memoryPath)
    call mem_allocate(this%auxname, LENAUXNAME, 0, 'AUXNAME', this%memoryPath)

    ! allocate special derived types
    allocate (this%TasManager)

    ! initialize
    this%inunit = 0
    this%iout = 0
    this%nvar = 0
    this%isize = 0
    this%ionper = 0
    this%lastonper = 0
    this%iprpak = 0
  end subroutine allocate_scalars

  !> @ brief Read options for package
  !!
  !<  Read options for this package.
  subroutine read_options(this)
    ! modules
    ! dummy
    class(UtlTvaType) :: this
    ! local
    character(len=:), allocatable :: line
    character(len=LINELENGTH) :: keyword, fname
    integer(I4B) :: ierr
    logical :: isfound, endOfBlock
    character(len=LENAUXNAME), dimension(:), allocatable :: vname
    integer(I4B) :: n
    integer(I4B) :: lloc
    integer(I4B) :: istart
    integer(I4B) :: istop
    ! formats
    character(len=*), parameter :: fmtiprpak = &
      &"(4x,'TVA INFORMATION WILL BE PRINTED TO LISTING FILE.')"
    character(len=*), parameter :: fmtts = &
      &"(4x, 'TIME-SERIES DATA WILL BE READ FROM FILE: ', a)"
    character(len=*), parameter :: fmttas = &
      &"(4x, 'TIME-ARRAY SERIES DATA WILL BE READ FROM FILE: ', a)"

    ! get options block
    call this%parser%GetBlock('OPTIONS', isfound, ierr, blockRequired=.false., &
                              supportOpenClose=.true.)

    ! parse options block if detected
    if (isfound) then
      write (this%iout, '(1x,a)') 'PROCESSING TVA OPTIONS'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('AUXILIARY')
          call this%parser%GetRemainingLine(line)
          lloc = 1
          call urdaux(this%nvar, this%parser%iuactive, this%iout, lloc, &
                      istart, istop, vname, line, text)
          call mem_reallocate(this%auxname, LENAUXNAME, this%nvar, &
                              'AUXNAME', this%memoryPath)
          do n = 1, this%nvar
            this%auxname(n) = vname(n)
          end do
          deallocate (vname)
        case ('PRINT_INPUT')
          this%iprpak = 1
          write (this%iout, fmtiprpak)
        case ('TAS6')
          call this%parser%GetStringCaps(keyword)
          if (trim(adjustl(keyword)) /= 'FILEIN') then
            errmsg = 'TAS6 keyword must be followed by "FILEIN" '// &
                     'then by filename.'
            call store_error(errmsg)
            call this%parser%StoreErrorUnit()
          end if
          call this%parser%GetString(fname)
          write (this%iout, fmttas) trim(fname)
          call this%TasManager%add_tasfile(fname)
        case default
          write (errmsg, '(a,a)') 'Unknown TVA option: ', trim(keyword)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end select
      end do
      write (this%iout, '(1x,a)') 'END OF TVA OPTIONS'
    end if
  end subroutine read_options

  !> @ brief Allocate package arrays
  !!
  !<  Allocate and initialize package arrays.
  subroutine allocate_arrays(this)
    ! modules
    ! dummy variables
    class(UtlTvaType) :: this !< UtlTvaType object
    ! local
    integer(I4B) :: i
    integer(I4B) :: j

    ! allocate arrays
    allocate(this%variables(this%nvar))
    do i = 1, this%nvar
      call mem_allocate(this%variables(i)%dblvec, this%isize, &
                        this%auxname(i), this%memoryPath)
    end do

    ! initialize each dblvec to zero
    do i = 1, this%nvar
      do j = 1, this%isize
        this%variables(i)%dblvec(j) = DZERO
      end do
    end do
  end subroutine allocate_arrays

  !> @ brief Fill the array with values from dblvec
  !<
  subroutine fill(this, dblvec, auxname)
    ! dummy
    class(UtlTvaType) :: this !< UtlTvaType object
    character(len=*), intent(in) :: auxname
    real(DP), dimension(:), intent(inout) :: dblvec
    ! local
    integer(I4B) :: i
    integer(I4B) :: ivar
    integer(I4B) :: nu
    integer(I4B) :: nr
    logical(LGP) :: found

    ivar = 0
    found = .false.
    do i = 1, this%nvar
      if (this%auxname(i) == auxname) then
        found = .true.
        ivar = i
      end if
    end do

    ! copy values from this object into the passed array
    do nr = 1, this%isize
      nu = this%dis%get_nodeuser(nr)
      dblvec(nr) = this%variables(ivar)%dblvec(nu)
    end do

  end subroutine fill

  !> @ brief Read and prepare
  !!
  !!  Read and prepare the period data block and fill dblvec
  !!  if the next period block corresponds to this time step.
  !!
  !<
  subroutine tva_rp(this)
    ! modules
    use TdisModule, only: kper, nper
    ! dummy
    class(UtlTvaType), intent(inout) :: this !< UtlTvaType object
    ! local
    character(len=LINELENGTH) :: line
    logical :: isfound
    integer(I4B) :: ierr
    ! formats
    character(len=*), parameter :: fmtblkerr = &
      &"('Looking for BEGIN PERIOD iper.  Found ', a, ' instead.')"
    character(len=*), parameter :: fmtlsp = &
      &"(1X,/1X,'REUSING ',A,'S FROM LAST STRESS PERIOD')"

    ! Set ionper to the stress period number for which a new block of data
    !    will be read.
    if (this%inunit == 0) return

    ! get stress period data
    if (this%ionper < kper) then

      ! get period block
      call this%parser%GetBlock('PERIOD', isfound, ierr, &
                                supportOpenClose=.true., &
                                blockRequired=.false.)
      if (isfound) then

        ! read ionper and check for increasing period numbers
        call this%read_check_ionper()
      else

        ! PERIOD block not found
        if (ierr < 0) then
          ! End of file found; data applies for remainder of simulation.
          this%ionper = nper + 1
        else
          ! Found invalid block
          call this%parser%GetCurrentLine(line)
          write (errmsg, fmtblkerr) adjustl(trim(line))
          call store_error(errmsg, terminate=.TRUE.)
        end if
      end if
    end if

    ! Read data if ionper == kper
    if (this%ionper == kper) then

      ! Remove all time-series and time-array-series links associated with
      !    this package.
      !    Do not reset as we are using a "settings" approach here in which the
      !    settings remain the same until the user changes them.
      call this%TasManager%Reset(this%subcomponent)
      call this%tva_rp_array(line)

    else
      ! using data from the last stress period
      write (this%iout, fmtlsp) trim(ftype)
    end if

    ! write summary of maw well stress period error messages
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
  end subroutine tva_rp

  !> @ brief tva_rp_array
  !!
  !<  Read the stress period data in array format
  subroutine tva_rp_array(this, line)
    use ConstantsModule, only: LENTIMESERIESNAME, LENANAME
    use SimModule, only: store_error
    use ArrayHandlersModule, only: ifind
    use Double1dReaderModule, only: read_dbl1d
    use IdmLoggerModule, only: idm_log_var
    ! dummy
    class(UtlTvaType), intent(inout) :: this !< UtlTvaType object
    character(len=LINELENGTH), intent(inout) :: line
    ! local
    integer(I4B) :: n
    integer(I4B) :: ncolbnd
    integer(I4B) :: ivar
    integer(I4B) :: jauxcol, ivarsread
    integer(I4B), dimension(:), allocatable, target :: nodelist
    character(len=LENTIMESERIESNAME) :: tasName
    ! character(len=LENANAME) :: aname
    character(len=LINELENGTH) :: keyword
    logical :: endOfBlock
    logical :: convertFlux

    ! these time array series pointers need to be non-contiguous
    !    because a slice of bound is passed
    real(DP), dimension(:), pointer :: bndArrayPtr => null()

    ! Initialize
    jauxcol = 0
    ivarsread = 0
    ncolbnd = 1
    allocate (nodelist(this%isize))
    do n = 1, size(nodelist)
      nodelist(n) = n
    end do

    ! Read variable
    do
      call this%parser%GetNextLine(endOfBlock)
      if (endOfBlock) exit
      call this%parser%GetStringCaps(keyword)

      ivar = ifind(this%auxname, keyword)

      ! Parse the keywords
      if (ivar > 0) then
        !
        ! Look for keyword TIMEARRAYSERIES and time-array series
        !    name on line, following RECHARGE
        call this%parser%GetStringCaps(keyword)
        if (keyword == 'TIMEARRAYSERIES') then
          ! Get time-array series name
          call this%parser%GetStringCaps(tasName)
          bndArrayPtr => this%variables(ivar)%dblvec(:)
          ! Make a time-array-series link and add it to the list of links
          ! contained in the TimeArraySeriesManagerType object.
          convertflux = .false.
          call this%TasManager%MakeTasLink(this%subcomponent, bndArrayPtr, &
                                           this%iprpak, tasName, &
                                           this%auxname(ivar), &
                                           convertFlux, nodelist, &
                                           this%parser%iuactive)
        else

          call read_dbl1d(this%parser, this%variables(ivar)%dblvec(:), &
                          this%auxname(ivar))
          call idm_log_var(this%variables(ivar)%dblvec(:), &
                           this%auxname(ivar), this%memoryPath, &
                           this%iout)
        end if

      else
        call store_error('Looking for variable name ' // trim(keyword) // &
                         &'Found: '//trim(line))
        call this%parser%StoreErrorUnit()
      end if
    end do
  end subroutine tva_rp_array

  !> @ brief Advance
  !!
  !!  Call the advance method on the time series so that new values
  !<  are interpolated and entered into dblvec
  subroutine tva_ad(this)
    ! modules
    ! dummy
    class(UtlTvaType), intent(inout) :: this !< UtlTvaType object
    ! local

    ! Advance the time series
    call this%TasManager%ad()

  end subroutine tva_ad

  !> @ brief Deallocate variables
  !!
  !<  Deallocate and nullify package variables.
  subroutine tva_da(this)
    ! modules
    ! dummy variables
    class(UtlTvaType) :: this !< UtlTvaType object
    integer(I4B) :: i
    !
    ! deallocate arrays in memory manager
    do i = 1, this%nvar
      call mem_deallocate(this%variables(i)%dblvec)
    end do
    deallocate(this%variables)
    !
    ! deallocate scalars in memory manager
    call mem_deallocate(this%inunit)
    call mem_deallocate(this%iout)
    call mem_deallocate(this%nvar)
    call mem_deallocate(this%isize)
    call mem_deallocate(this%ionper)
    call mem_deallocate(this%lastonper)
    call mem_deallocate(this%iprpak)
    !
    ! deallocate derived types
    call this%TasManager%da()
    deallocate (this%TasManager)
    nullify (this%TasManager)
  end subroutine tva_da

  !> @ brief Check ionper
  !!
  !!  Generic method to read and check ionperiod, which is used to determine
  !!  if new period data should be read from the input file. The check of
  !!  ionperiod also makes sure periods are increasing in subsequent period
  !<  data blocks.  Copied from NumericalPackage
  subroutine read_check_ionper(this)
    ! modules
    use TdisModule, only: kper
    ! dummy variables
    class(UtlTvaType), intent(inout) :: this !< UtlTvaType object
    !
    ! save last value and read period number
    this%lastonper = this%ionper
    this%ionper = this%parser%GetInteger()
    !
    ! make check
    if (this%ionper <= this%lastonper) then
      write (errmsg, '(a, i0, a, i0, a, i0, a)') &
        'Error in stress period ', kper, &
        '. Period numbers not increasing.  Found ', this%ionper, &
        ' but last period block was assigned ', this%lastonper, '.'
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
    end if
  end subroutine read_check_ionper

end module UtlTvaModule
