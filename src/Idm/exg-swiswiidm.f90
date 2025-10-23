! ** Do Not Modify! MODFLOW 6 system generated file. **
module ExgSwiswiInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public exg_swiswi_param_definitions
  public exg_swiswi_aggregate_definitions
  public exg_swiswi_block_definitions
  public ExgSwiswiParamFoundType
  public exg_swiswi_multi_package
  public exg_swiswi_subpackages

  type ExgSwiswiParamFoundType
    logical :: ipr_input = .false.
    logical :: ipr_flow = .false.
  end type ExgSwiswiParamFoundType

  logical :: exg_swiswi_multi_package = .true.

  character(len=16), parameter :: &
    exg_swiswi_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    exgswiswi_ipr_input = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'SWISWI', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_INPUT', & ! tag name
    'IPR_INPUT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to print input to list file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exgswiswi_ipr_flow = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'SWISWI', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'IPR_FLOW', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to print swiswi flows to list file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exg_swiswi_param_definitions(*) = &
    [ &
    exgswiswi_ipr_input, &
    exgswiswi_ipr_flow &
    ]

  type(InputParamDefinitionType), parameter :: &
    exg_swiswi_aggregate_definitions(*) = &
    [ &
    InputParamDefinitionType &
    ( &
    '', & ! component
    '', & ! subcomponent
    '', & ! block
    '', & ! tag name
    '', & ! fortran variable
    '', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    ) &
    ]

  type(InputBlockDefinitionType), parameter :: &
    exg_swiswi_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module ExgSwiswiInputModule
