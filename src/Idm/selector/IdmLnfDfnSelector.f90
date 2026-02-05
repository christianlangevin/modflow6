! ** Do Not Modify! MODFLOW 6 system generated file. **
module IdmLnfDfnSelectorModule

  use ConstantsModule, only: LENVARNAME
  use SimModule, only: store_error
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  use LnfChdInputModule
  use LnfCndInputModule
  use LnfDisuInputModule
  use LnfNamInputModule
  use LnfOcInputModule

  implicit none
  private
  public :: lnf_param_definitions
  public :: lnf_aggregate_definitions
  public :: lnf_block_definitions
  public :: lnf_idm_multi_package
  public :: lnf_idm_subpackages
  public :: lnf_idm_integrated

contains

  subroutine set_param_pointer(input_dfn, input_dfn_target)
    type(InputParamDefinitionType), dimension(:), pointer :: input_dfn
    type(InputParamDefinitionType), dimension(:), target :: input_dfn_target
    input_dfn => input_dfn_target
  end subroutine set_param_pointer

  subroutine set_block_pointer(input_dfn, input_dfn_target)
    type(InputBlockDefinitionType), dimension(:), pointer :: input_dfn
    type(InputBlockDefinitionType), dimension(:), target :: input_dfn_target
    input_dfn => input_dfn_target
  end subroutine set_block_pointer

  subroutine set_subpkg_pointer(subpkg_list, subpkg_list_target)
    character(len=16), dimension(:), pointer :: subpkg_list
    character(len=16), dimension(:), target :: subpkg_list_target
    subpkg_list => subpkg_list_target
  end subroutine set_subpkg_pointer

  function lnf_param_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputParamDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('CHD')
      call set_param_pointer(input_definition, lnf_chd_param_definitions)
    case ('CND')
      call set_param_pointer(input_definition, lnf_cnd_param_definitions)
    case ('DISU')
      call set_param_pointer(input_definition, lnf_disu_param_definitions)
    case ('NAM')
      call set_param_pointer(input_definition, lnf_nam_param_definitions)
    case ('OC')
      call set_param_pointer(input_definition, lnf_oc_param_definitions)
    case default
    end select
    return
  end function lnf_param_definitions

  function lnf_aggregate_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputParamDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('CHD')
      call set_param_pointer(input_definition, lnf_chd_aggregate_definitions)
    case ('CND')
      call set_param_pointer(input_definition, lnf_cnd_aggregate_definitions)
    case ('DISU')
      call set_param_pointer(input_definition, lnf_disu_aggregate_definitions)
    case ('NAM')
      call set_param_pointer(input_definition, lnf_nam_aggregate_definitions)
    case ('OC')
      call set_param_pointer(input_definition, lnf_oc_aggregate_definitions)
    case default
    end select
    return
  end function lnf_aggregate_definitions

  function lnf_block_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputBlockDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('CHD')
      call set_block_pointer(input_definition, lnf_chd_block_definitions)
    case ('CND')
      call set_block_pointer(input_definition, lnf_cnd_block_definitions)
    case ('DISU')
      call set_block_pointer(input_definition, lnf_disu_block_definitions)
    case ('NAM')
      call set_block_pointer(input_definition, lnf_nam_block_definitions)
    case ('OC')
      call set_block_pointer(input_definition, lnf_oc_block_definitions)
    case default
    end select
    return
  end function lnf_block_definitions

  function lnf_idm_multi_package(subcomponent) result(multi_package)
    character(len=*), intent(in) :: subcomponent
    logical :: multi_package
    select case (subcomponent)
    case ('CHD')
      multi_package = lnf_chd_multi_package
    case ('CND')
      multi_package = lnf_cnd_multi_package
    case ('DISU')
      multi_package = lnf_disu_multi_package
    case ('NAM')
      multi_package = lnf_nam_multi_package
    case ('OC')
      multi_package = lnf_oc_multi_package
    case default
      call store_error('Idm selector subcomponent not found; '//&
                       &'component="LNF"'//&
                       &', subcomponent="'//trim(subcomponent)//'".', .true.)
    end select
    return
  end function lnf_idm_multi_package

  function lnf_idm_subpackages(subcomponent) result(subpackages)
    character(len=*), intent(in) :: subcomponent
    character(len=16), dimension(:), pointer :: subpackages
    select case (subcomponent)
    case ('CHD')
      call set_subpkg_pointer(subpackages, lnf_chd_subpackages)
    case ('CND')
      call set_subpkg_pointer(subpackages, lnf_cnd_subpackages)
    case ('DISU')
      call set_subpkg_pointer(subpackages, lnf_disu_subpackages)
    case ('NAM')
      call set_subpkg_pointer(subpackages, lnf_nam_subpackages)
    case ('OC')
      call set_subpkg_pointer(subpackages, lnf_oc_subpackages)
    case default
    end select
    return
  end function lnf_idm_subpackages

  function lnf_idm_integrated(subcomponent) result(integrated)
    character(len=*), intent(in) :: subcomponent
    logical :: integrated
    integrated = .false.
    select case (subcomponent)
    case ('CHD')
      integrated = .true.
    case ('CND')
      integrated = .true.
    case ('DISU')
      integrated = .true.
    case ('NAM')
      integrated = .true.
    case ('OC')
      integrated = .true.
    case default
    end select
    return
  end function lnf_idm_integrated

end module IdmLnfDfnSelectorModule
