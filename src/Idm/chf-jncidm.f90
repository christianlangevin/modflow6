! ** Do Not Modify! MODFLOW 6 system generated file. **
module ChfJncInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public chf_jnc_param_definitions
  public chf_jnc_aggregate_definitions
  public chf_jnc_block_definitions
  public ChfJncParamFoundType
  public chf_jnc_multi_package
  public chf_jnc_subpackages

  type ChfJncParamFoundType
    logical :: iprflow = .false.
  end type ChfJncParamFoundType

  logical :: chf_jnc_multi_package = .false.

  character(len=16), parameter :: &
    chf_jnc_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    chfjnc_iprflow = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'JNC', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'IPRFLOW', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to print JNC flows to listing file', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chf_jnc_param_definitions(*) = &
    [ &
    chfjnc_iprflow &
    ]

  type(InputParamDefinitionType), parameter :: &
    chf_jnc_aggregate_definitions(*) = &
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
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    ) &
    ]

  type(InputBlockDefinitionType), parameter :: &
    chf_jnc_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module ChfJncInputModule
