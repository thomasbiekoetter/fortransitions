module cosmotransitions__potentials
!! Abstract interfaces for the potentials used throughout the package.
!!
!! The Python package passes plain callables (V, dV, d2V) around. In
!! Fortran, an object with type-bound procedures is used instead so that
!! stateful potentials (like the spline path potential) work naturally.

  use cosmotransitions__config, only : wp

  implicit none

  private

  public :: potential_1d
  public :: potential_nd
  public :: pot1d_func

  type, abstract :: potential_1d
    !! A one-dimensional potential V(phi). The first and second derivatives
    !! default to fourth-order finite differences with step `fd_eps`
    !! (mirroring tunneling1D.SingleFieldInstanton.dV/d2V), but can be
    !! overridden with analytic expressions.
    real(wp) :: fd_eps = 1.0e-3_wp
      !! Step size for the default finite-difference derivatives. It is
      !! overwritten by single_field_instanton%init with
      !! phi_eps*|phi_absMin - phi_metaMin|, like in the Python class.
  contains
    procedure(pot1d_f), deferred :: v
    procedure :: dv => potential_1d_dv
    procedure :: d2v => potential_1d_d2v
  end type potential_1d

  abstract interface
    function pot1d_f(self, phi) result(y)
      import :: wp, potential_1d
      implicit none
      class(potential_1d), intent(inout) :: self
      real(wp), intent(in) :: phi
      real(wp) :: y
    end function pot1d_f
  end interface

  type, abstract :: potential_nd
    !! A multi-field potential V(x) with its gradient, x having ndim
    !! components. This plays the role of the (V, dV) callable pair passed
    !! to pathDeformation.fullTunneling.
  contains
    procedure(potnd_v), deferred :: v
    procedure(potnd_dv), deferred :: grad
  end type potential_nd

  abstract interface
    function potnd_v(self, x) result(y)
      import :: wp, potential_nd
      implicit none
      class(potential_nd), intent(inout) :: self
      real(wp), intent(in) :: x(:)
      real(wp) :: y
    end function potnd_v

    function potnd_dv(self, x) result(dv)
      import :: wp, potential_nd
      implicit none
      class(potential_nd), intent(inout) :: self
      real(wp), intent(in) :: x(:)
      real(wp) :: dv(size(x))
    end function potnd_dv
  end interface

  abstract interface
    function plain_func(phi) result(y)
      import :: wp
      implicit none
      real(wp), intent(in) :: phi
      real(wp) :: y
    end function plain_func
  end interface

  type, extends(potential_1d) :: pot1d_func
    !! Convenience wrapper turning plain functions into a potential_1d.
    !! `dvf` and `d2vf` are optional; if not associated, the inherited
    !! finite-difference derivatives are used.
    procedure(plain_func), pointer, nopass :: vf => null()
    procedure(plain_func), pointer, nopass :: dvf => null()
    procedure(plain_func), pointer, nopass :: d2vf => null()
  contains
    procedure :: v => pot1d_func_v
    procedure :: dv => pot1d_func_dv
    procedure :: d2v => pot1d_func_d2v
  end type pot1d_func

contains

  function potential_1d_dv(self, phi) result(y)
    !! dV/dphi via fourth-order finite differences with step `fd_eps`.
    !! Port of tunneling1D.SingleFieldInstanton.dV.

    class(potential_1d), intent(inout) :: self
    real(wp), intent(in) :: phi
    real(wp) :: y

    real(wp) :: eps

    eps = self%fd_eps
    y = (self%v(phi - 2.0_wp*eps) - 8.0_wp*self%v(phi - eps)  &
      + 8.0_wp*self%v(phi + eps) - self%v(phi + 2.0_wp*eps))/(12.0_wp*eps)

  end function potential_1d_dv

  function potential_1d_d2v(self, phi) result(y)
    !! d2V/dphi2 via fourth-order finite differences with step `fd_eps`.
    !! Port of tunneling1D.SingleFieldInstanton.d2V.

    class(potential_1d), intent(inout) :: self
    real(wp), intent(in) :: phi
    real(wp) :: y

    real(wp) :: eps

    eps = self%fd_eps
    y = (-self%v(phi - 2.0_wp*eps) + 16.0_wp*self%v(phi - eps)  &
      - 30.0_wp*self%v(phi) + 16.0_wp*self%v(phi + eps)  &
      - self%v(phi + 2.0_wp*eps))/(12.0_wp*eps*eps)

  end function potential_1d_d2v

  function pot1d_func_v(self, phi) result(y)

    class(pot1d_func), intent(inout) :: self
    real(wp), intent(in) :: phi
    real(wp) :: y

    y = self%vf(phi)

  end function pot1d_func_v

  function pot1d_func_dv(self, phi) result(y)

    class(pot1d_func), intent(inout) :: self
    real(wp), intent(in) :: phi
    real(wp) :: y

    if (associated(self%dvf)) then
      y = self%dvf(phi)
    else
      y = potential_1d_dv(self, phi)
    end if

  end function pot1d_func_dv

  function pot1d_func_d2v(self, phi) result(y)

    class(pot1d_func), intent(inout) :: self
    real(wp), intent(in) :: phi
    real(wp) :: y

    if (associated(self%d2vf)) then
      y = self%d2vf(phi)
    else
      y = potential_1d_d2v(self, phi)
    end if

  end function pot1d_func_d2v

end module cosmotransitions__potentials
