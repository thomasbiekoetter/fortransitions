module fortransitions__c_potentials
!! Potentials whose values are supplied by C function pointers, e.g. by
!! Python callables wrapped with ctypes. They are used by the bind(C)
!! interface (modules fortransitions__c_*) for the Python package
!! pyfortransitions.

  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_c_binding, only : c_double
  use, intrinsic :: iso_c_binding, only : c_funptr
  use, intrinsic :: iso_c_binding, only : c_associated
  use, intrinsic :: iso_c_binding, only : c_f_procpointer

  use fortransitions__config, only : wp
  use fortransitions__potentials, only : potential_1d
  use fortransitions__potentials, only : potential_nd
  use fortransitions__potentials, only : potential_1d_dv
  use fortransitions__potentials, only : potential_1d_d2v

  implicit none

  private

  public :: v_nd_cb
  public :: grad_nd_cb
  public :: f_1d_cb
  public :: c_potential_nd
  public :: c_potential_1d
  public :: new_c_potential_nd
  public :: new_c_potential_1d

  abstract interface
    function v_nd_cb(n, x) result(y) bind(C)
      !! C signature: double v(int n, const double *x).
      import :: c_int, c_double
      implicit none
      integer(c_int), value :: n
      real(c_double), intent(in) :: x(n)
      real(c_double) :: y
    end function v_nd_cb

    subroutine grad_nd_cb(n, x, g) bind(C)
      !! C signature: void grad(int n, const double *x, double *g).
      import :: c_int, c_double
      implicit none
      integer(c_int), value :: n
      real(c_double), intent(in) :: x(n)
      real(c_double), intent(out) :: g(n)
    end subroutine grad_nd_cb

    function f_1d_cb(phi) result(y) bind(C)
      !! C signature: double f(double phi). Used for V, dV and d2V of a
      !! one-dimensional potential.
      import :: c_double
      implicit none
      real(c_double), value :: phi
      real(c_double) :: y
    end function f_1d_cb
  end interface

  type, extends(potential_nd) :: c_potential_nd
    !! Multi-field potential V(x) and its gradient given by C callbacks.
    procedure(v_nd_cb), pointer, nopass :: vf => null()
    procedure(grad_nd_cb), pointer, nopass :: gradf => null()
  contains
    procedure :: v => c_potential_nd_v
    procedure :: grad => c_potential_nd_grad
  end type c_potential_nd

  type, extends(potential_1d) :: c_potential_1d
    !! One-dimensional potential given by C callbacks. `dvf` and `d2vf`
    !! are optional; if not associated, the inherited finite-difference
    !! derivatives of potential_1d are used.
    procedure(f_1d_cb), pointer, nopass :: vf => null()
    procedure(f_1d_cb), pointer, nopass :: dvf => null()
    procedure(f_1d_cb), pointer, nopass :: d2vf => null()
  contains
    procedure :: v => c_potential_1d_v
    procedure :: dv => c_potential_1d_dv
    procedure :: d2v => c_potential_1d_d2v
  end type c_potential_1d

contains

  function new_c_potential_nd(v_ptr, grad_ptr) result(pot)
    !! Build a c_potential_nd from C function pointers.

    type(c_funptr), intent(in) :: v_ptr
    type(c_funptr), intent(in) :: grad_ptr
    type(c_potential_nd) :: pot

    call c_f_procpointer(v_ptr, pot%vf)
    call c_f_procpointer(grad_ptr, pot%gradf)

  end function new_c_potential_nd

  function new_c_potential_1d(v_ptr, dv_ptr, d2v_ptr) result(pot)
    !! Build a c_potential_1d from C function pointers. A null `dv_ptr` or
    !! `d2v_ptr` selects the finite-difference derivative.

    type(c_funptr), intent(in) :: v_ptr
    type(c_funptr), intent(in) :: dv_ptr
    type(c_funptr), intent(in) :: d2v_ptr
    type(c_potential_1d) :: pot

    call c_f_procpointer(v_ptr, pot%vf)
    if (c_associated(dv_ptr)) call c_f_procpointer(dv_ptr, pot%dvf)
    if (c_associated(d2v_ptr)) call c_f_procpointer(d2v_ptr, pot%d2vf)

  end function new_c_potential_1d

  function c_potential_nd_v(self, x) result(y)

    class(c_potential_nd), intent(inout) :: self
    real(wp), intent(in) :: x(:)
    real(wp) :: y

    real(c_double) :: xc(size(x))

    xc = real(x, c_double)
    y = real(self%vf(int(size(x), c_int), xc), wp)

  end function c_potential_nd_v

  function c_potential_nd_grad(self, x) result(dv)

    class(c_potential_nd), intent(inout) :: self
    real(wp), intent(in) :: x(:)
    real(wp) :: dv(size(x))

    real(c_double) :: xc(size(x))
    real(c_double) :: gc(size(x))

    xc = real(x, c_double)
    call self%gradf(int(size(x), c_int), xc, gc)
    dv = real(gc, wp)

  end function c_potential_nd_grad

  function c_potential_1d_v(self, phi) result(y)

    class(c_potential_1d), intent(inout) :: self
    real(wp), intent(in) :: phi
    real(wp) :: y

    y = real(self%vf(real(phi, c_double)), wp)

  end function c_potential_1d_v

  function c_potential_1d_dv(self, phi) result(y)

    class(c_potential_1d), intent(inout) :: self
    real(wp), intent(in) :: phi
    real(wp) :: y

    if (associated(self%dvf)) then
      y = real(self%dvf(real(phi, c_double)), wp)
    else
      y = potential_1d_dv(self, phi)
    end if

  end function c_potential_1d_dv

  function c_potential_1d_d2v(self, phi) result(y)

    class(c_potential_1d), intent(inout) :: self
    real(wp), intent(in) :: phi
    real(wp) :: y

    if (associated(self%d2vf)) then
      y = real(self%d2vf(real(phi, c_double)), wp)
    else
      y = potential_1d_d2v(self, phi)
    end if

  end function c_potential_1d_d2v

end module fortransitions__c_potentials
