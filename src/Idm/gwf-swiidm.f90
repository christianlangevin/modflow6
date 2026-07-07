! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfSwiInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_swi_param_definitions
  public gwf_swi_aggregate_definitions
  public gwf_swi_block_definitions
  public GwfSwiParamFoundType
  public gwf_swi_multi_package
  public gwf_swi_subpackages

  type GwfSwiParamFoundType
    logical :: zeta_filerecord = .false.
    logical :: zeta = .false.
    logical :: fileout = .false.
    logical :: zetafile = .false.
    logical :: tva6_filerecord = .false.
    logical :: tva6 = .false.
    logical :: filein = .false.
    logical :: tva6_filename = .false.
    logical :: hsalt_user = .false.
  end type GwfSwiParamFoundType

  logical :: gwf_swi_multi_package = .false.

  character(len=16), parameter :: &
    gwf_swi_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfswi_zeta_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SWI', & ! subcomponent
    'OPTIONS', & ! block
    'ZETA_FILERECORD', & ! tag name
    'ZETA_FILERECORD', & ! fortran variable
    'RECORD ZETA FILEOUT ZETAFILE', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfswi_zeta = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SWI', & ! subcomponent
    'OPTIONS', & ! block
    'ZETA', & ! tag name
    'ZETA', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'zeta keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfswi_fileout = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SWI', & ! subcomponent
    'OPTIONS', & ! block
    'FILEOUT', & ! tag name
    'FILEOUT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfswi_zetafile = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SWI', & ! subcomponent
    'OPTIONS', & ! block
    'ZETAFILE', & ! tag name
    'ZETAFILE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfswi_tva6_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SWI', & ! subcomponent
    'OPTIONS', & ! block
    'TVA6_FILERECORD', & ! tag name
    'TVA6_FILERECORD', & ! fortran variable
    'RECORD TVA6 FILEIN TVA6_FILENAME', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfswi_tva6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SWI', & ! subcomponent
    'OPTIONS', & ! block
    'TVA6', & ! tag name
    'TVA6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'grid package component', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfswi_filein = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SWI', & ! subcomponent
    'OPTIONS', & ! block
    'FILEIN', & ! tag name
    'FILEIN', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfswi_tva6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SWI', & ! subcomponent
    'OPTIONS', & ! block
    'TVA6_FILENAME', & ! tag name
    'TVA6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfswi_hsalt_user = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SWI', & ! subcomponent
    'OPTIONS', & ! block
    'SALTWATER_HEAD', & ! tag name
    'HSALT_USER', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'saltwater head', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_swi_param_definitions(*) = &
    [ &
    gwfswi_zeta_filerecord, &
    gwfswi_zeta, &
    gwfswi_fileout, &
    gwfswi_zetafile, &
    gwfswi_tva6_filerecord, &
    gwfswi_tva6, &
    gwfswi_filein, &
    gwfswi_tva6_filename, &
    gwfswi_hsalt_user &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwf_swi_aggregate_definitions(*) = &
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
    gwf_swi_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module GwfSwiInputModule
