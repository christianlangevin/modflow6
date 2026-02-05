! ** Do Not Modify! MODFLOW 6 system generated file. **
module LnfDisuInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public lnf_disu_param_definitions
  public lnf_disu_aggregate_definitions
  public lnf_disu_block_definitions
  public LnfDisuParamFoundType
  public lnf_disu_multi_package
  public lnf_disu_subpackages

  type LnfDisuParamFoundType
    logical :: length_units = .false.
    logical :: nogrb = .false.
    logical :: export_ascii = .false.
    logical :: nodes = .false.
    logical :: nja = .false.
    logical :: nvert = .false.
    logical :: top = .false.
    logical :: bot = .false.
    logical :: area = .false.
    logical :: idomain = .false.
    logical :: iac = .false.
    logical :: ja = .false.
    logical :: ihc = .false.
    logical :: cl12 = .false.
    logical :: hwva = .false.
    logical :: angldegx = .false.
    logical :: iv = .false.
    logical :: xv = .false.
    logical :: yv = .false.
    logical :: icell2d = .false.
    logical :: xc = .false.
    logical :: yc = .false.
    logical :: ncvert = .false.
    logical :: icvert = .false.
  end type LnfDisuParamFoundType

  logical :: lnf_disu_multi_package = .false.

  character(len=16), parameter :: &
    lnf_disu_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    lnfdisu_length_units = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'DISU', & ! subcomponent
    'OPTIONS', & ! block
    'LENGTH_UNITS', & ! tag name
    'LENGTH_UNITS', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'model length units', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfdisu_nogrb = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'DISU', & ! subcomponent
    'OPTIONS', & ! block
    'NOGRB', & ! tag name
    'NOGRB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'do not write binary grid file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfdisu_export_ascii = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'DISU', & ! subcomponent
    'OPTIONS', & ! block
    'EXPORT_ARRAY_ASCII', & ! tag name
    'EXPORT_ASCII', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'export array variables to layered ascii files.', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfdisu_nodes = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'DISU', & ! subcomponent
    'DIMENSIONS', & ! block
    'NODES', & ! tag name
    'NODES', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of layers', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfdisu_nja = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'DISU', & ! subcomponent
    'DIMENSIONS', & ! block
    'NJA', & ! tag name
    'NJA', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of columns', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfdisu_nvert = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'DISU', & ! subcomponent
    'DIMENSIONS', & ! block
    'NVERT', & ! tag name
    'NVERT', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of vertices', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfdisu_top = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'DISU', & ! subcomponent
    'GRIDDATA', & ! block
    'TOP', & ! tag name
    'TOP', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'cell top elevation', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfdisu_bot = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'DISU', & ! subcomponent
    'GRIDDATA', & ! block
    'BOT', & ! tag name
    'BOT', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'cell bottom elevation', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfdisu_area = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'DISU', & ! subcomponent
    'GRIDDATA', & ! block
    'AREA', & ! tag name
    'AREA', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'cell surface area', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfdisu_idomain = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'DISU', & ! subcomponent
    'GRIDDATA', & ! block
    'IDOMAIN', & ! tag name
    'IDOMAIN', & ! fortran variable
    'INTEGER1D', & ! type
    'NODES', & ! shape
    'idomain existence array', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfdisu_iac = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'DISU', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'IAC', & ! tag name
    'IAC', & ! fortran variable
    'INTEGER1D', & ! type
    'NODES', & ! shape
    'number of cell connections', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfdisu_ja = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'DISU', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'JA', & ! tag name
    'JA', & ! fortran variable
    'INTEGER1D', & ! type
    'NJA', & ! shape
    'grid connectivity', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfdisu_ihc = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'DISU', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'IHC', & ! tag name
    'IHC', & ! fortran variable
    'INTEGER1D', & ! type
    'NJA', & ! shape
    'connection type', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfdisu_cl12 = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'DISU', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'CL12', & ! tag name
    'CL12', & ! fortran variable
    'DOUBLE1D', & ! type
    'NJA', & ! shape
    'connection lengths', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfdisu_hwva = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'DISU', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'HWVA', & ! tag name
    'HWVA', & ! fortran variable
    'DOUBLE1D', & ! type
    'NJA', & ! shape
    'connection lengths', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfdisu_angldegx = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'DISU', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'ANGLDEGX', & ! tag name
    'ANGLDEGX', & ! fortran variable
    'DOUBLE1D', & ! type
    'NJA', & ! shape
    'angle of face normal to connection', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfdisu_iv = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'DISU', & ! subcomponent
    'VERTICES', & ! block
    'IV', & ! tag name
    'IV', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'vertex number', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfdisu_xv = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'DISU', & ! subcomponent
    'VERTICES', & ! block
    'XV', & ! tag name
    'XV', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'x-coordinate for vertex', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfdisu_yv = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'DISU', & ! subcomponent
    'VERTICES', & ! block
    'YV', & ! tag name
    'YV', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'y-coordinate for vertex', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfdisu_icell2d = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'DISU', & ! subcomponent
    'CELL2D', & ! block
    'ICELL2D', & ! tag name
    'ICELL2D', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'cell2d number', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfdisu_xc = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'DISU', & ! subcomponent
    'CELL2D', & ! block
    'XC', & ! tag name
    'XC', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'x-coordinate for cell center', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfdisu_yc = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'DISU', & ! subcomponent
    'CELL2D', & ! block
    'YC', & ! tag name
    'YC', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'y-coordinate for cell center', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfdisu_ncvert = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'DISU', & ! subcomponent
    'CELL2D', & ! block
    'NCVERT', & ! tag name
    'NCVERT', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of cell vertices', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfdisu_icvert = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'DISU', & ! subcomponent
    'CELL2D', & ! block
    'ICVERT', & ! tag name
    'ICVERT', & ! fortran variable
    'INTEGER1D', & ! type
    'NCVERT', & ! shape
    'array of vertex numbers', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnf_disu_param_definitions(*) = &
    [ &
    lnfdisu_length_units, &
    lnfdisu_nogrb, &
    lnfdisu_export_ascii, &
    lnfdisu_nodes, &
    lnfdisu_nja, &
    lnfdisu_nvert, &
    lnfdisu_top, &
    lnfdisu_bot, &
    lnfdisu_area, &
    lnfdisu_idomain, &
    lnfdisu_iac, &
    lnfdisu_ja, &
    lnfdisu_ihc, &
    lnfdisu_cl12, &
    lnfdisu_hwva, &
    lnfdisu_angldegx, &
    lnfdisu_iv, &
    lnfdisu_xv, &
    lnfdisu_yv, &
    lnfdisu_icell2d, &
    lnfdisu_xc, &
    lnfdisu_yc, &
    lnfdisu_ncvert, &
    lnfdisu_icvert &
    ]

  type(InputParamDefinitionType), parameter :: &
    lnfdisu_vertices = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'DISU', & ! subcomponent
    'VERTICES', & ! block
    'VERTICES', & ! tag name
    'VERTICES', & ! fortran variable
    'RECARRAY IV XV YV', & ! type
    'NVERT', & ! shape
    'vertices data', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnfdisu_cell2d = InputParamDefinitionType &
    ( &
    'LNF', & ! component
    'DISU', & ! subcomponent
    'CELL2D', & ! block
    'CELL2D', & ! tag name
    'CELL2D', & ! fortran variable
    'RECARRAY ICELL2D XC YC NCVERT ICVERT', & ! type
    'NODES', & ! shape
    'cell2d data', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    lnf_disu_aggregate_definitions(*) = &
    [ &
    lnfdisu_vertices, &
    lnfdisu_cell2d &
    ]

  type(InputBlockDefinitionType), parameter :: &
    lnf_disu_block_definitions(*) = &
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
    ), &
    InputBlockDefinitionType( &
    'CONNECTIONDATA', & ! blockname
    .true., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'VERTICES', & ! blockname
    .false., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'CELL2D', & ! blockname
    .false., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module LnfDisuInputModule
