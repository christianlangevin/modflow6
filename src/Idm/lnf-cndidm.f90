! ** Do Not Modify! MODFLOW 6 system generated file. **
module LnfCndInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public lnf_cnd_param_definitions
  public lnf_cnd_aggregate_definitions
  public lnf_cnd_block_definitions
  public LnfCndParamFoundType
  public lnf_cnd_multi_package
  public lnf_cnd_subpackages

  type LnfCndParamFoundType
    logical :: print_input = .false.
    logical :: print_flows = .false.
    logical :: save_flows = .false.
    logical :: nja = .false.
    logical :: conductance = .false.
  end type LnfCndParamFoundType

  logical :: lnf_cnd_multi_package = .false.

  character(len=16), parameter :: &
    lnf_cnd_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    lnfcnd_print_input = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'CND', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_INPUT', & ! tag name
    'PRINT_INPUT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print input to listing file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfcnd_print_flows = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'CND', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'PRINT_FLOWS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print calculated flows to listing file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfcnd_save_flows = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'CND', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'SAVE_FLOWS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save flows for package to budget file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfcnd_nja = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'CND', & ! subcomponent
    'DIMENSIONS', & ! block
    'NJA', & ! tag name
    'NJA', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'total number of connections', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfcnd_conductance = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'CND', & ! subcomponent
    'GRIDDATA', & ! block
    'CONDUCTANCE', & ! tag name
    'CONDUCTANCE', & ! fortran variable
    'DOUBLE1D', & ! type
    'NJA', & ! shape
    'conductance for connections', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnf_cnd_param_definitions(*) = &
    [ &
    lnfcnd_print_input, &
    lnfcnd_print_flows, &
    lnfcnd_save_flows, &
    lnfcnd_nja, &
    lnfcnd_conductance &
    ]

  type(InputParamDefinitionType), parameter :: &
    lnf_cnd_aggregate_definitions(*) = &
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
    lnf_cnd_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'DIMENSIONS', & ! blockname
    .true., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'GRIDDATA', & ! blockname
    .true., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module LnfCndInputModule
