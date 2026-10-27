!> @brief Flow reduction for stress packages
!!
!! A stress package that reduces a withdrawal as the fluid runs out needs to
!! know how much fluid is actually available in a cell. The default answer is
!! the whole cell: the reduction is tapered over an interval measured up from
!! the cell bottom. An alternative flow formulation can answer differently,
!! and this type is how it does so without the stress package naming it.
!!
!! Two questions can be answered, and a formulation declares which of them it
!! can answer. They are not interchangeable:
!!
!!   - defines_interval: the fluid occupies a known elevation range within the
!!     cell. A sharp interface can say this, and a consumer is free to use the
!!     range as geometry, not only as a taper datum.
!!   - defines_limit: there is no such range, but the formulation can give the
!!     reduction factor and its derivative directly. A continuous unsaturated
!!     profile has no elevation at which the fluid runs out, so this is the
!!     only question it can answer.
!!
!! A formulation that declares neither leaves the stress package to behave as
!! it always has.
!<
module FlowReductionInterfaceModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DONE

  implicit none
  private
  public :: FlowReductionType

  !> @brief Provider of flow-reduction information for a stress package
  !<
  type :: FlowReductionType
    logical(LGP) :: defines_interval = .false. !< can report the fluid interval
    logical(LGP) :: defines_limit = .false. !< can report the factor directly
  contains
    procedure :: interval => fr_interval
    procedure :: limit => fr_limit
  end type FlowReductionType

contains

  !> @brief Narrow an interval of cell elevations to the fluid it holds
  !!
  !! The bounds arrive holding the interval the consumer would otherwise use,
  !! normally the cell bottom and top. An implementation narrows them to the
  !! part of the cell its fluid occupies. The default leaves them alone.
  !<
  subroutine fr_interval(this, n, bt, tp)
    class(FlowReductionType), intent(in) :: this !< this provider
    integer(I4B), intent(in) :: n !< reduced node number
    real(DP), intent(inout) :: bt !< bottom of the interval
    real(DP), intent(inout) :: tp !< top of the interval
  end subroutine fr_interval

  !> @brief Reduction factor for a withdrawal, and its derivative
  !!
  !! frac is the dimensionless fraction of the available fluid over which the
  !! withdrawal is tapered, as asked for by the stress package input. An
  !! implementation decides what the fluid available in this cell is and
  !! returns the factor to apply to the requested rate together with its
  !! derivative with respect to head. The default reduces nothing.
  !<
  subroutine fr_limit(this, n, h, frac, f, dfdh)
    class(FlowReductionType), intent(in) :: this !< this provider
    integer(I4B), intent(in) :: n !< reduced node number
    real(DP), intent(in) :: h !< head in the cell
    real(DP), intent(in) :: frac !< fraction to taper over
    real(DP), intent(out) :: f !< factor to apply to the rate
    real(DP), intent(out) :: dfdh !< derivative of f with respect to head

    f = DONE
    dfdh = DZERO

  end subroutine fr_limit

end module FlowReductionInterfaceModule
