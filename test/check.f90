module check__potentials
!! Test potentials: the two single-field potentials from the tunneling1D
!! docstring and the two-field potential from examples/fullTunneling.py of
!! the Python package.

  use cosmotransitions, only : wp
  use cosmotransitions, only : potential_nd

  implicit none

  private

  public :: v1_thin
  public :: dv1_thin
  public :: v2_thick
  public :: dv2_thick
  public :: example_potential

  type, extends(potential_nd) :: example_potential
    !! The two-dimensional potential of examples/fullTunneling.py.
    real(wp) :: c = 5.0_wp
    real(wp) :: fx = 10.0_wp
    real(wp) :: fy = 10.0_wp
  contains
    procedure :: v => example_v
    procedure :: grad => example_grad
  end type example_potential

contains

  function v1_thin(phi) result(y)
    real(wp), intent(in) :: phi
    real(wp) :: y
    y = 0.25_wp*phi**4 - 0.49_wp*phi**3 + 0.235_wp*phi**2
  end function v1_thin

  function dv1_thin(phi) result(y)
    real(wp), intent(in) :: phi
    real(wp) :: y
    y = phi*(phi - 0.47_wp)*(phi - 1.0_wp)
  end function dv1_thin

  function v2_thick(phi) result(y)
    real(wp), intent(in) :: phi
    real(wp) :: y
    y = 0.25_wp*phi**4 - 0.4_wp*phi**3 + 0.1_wp*phi**2
  end function v2_thick

  function dv2_thick(phi) result(y)
    real(wp), intent(in) :: phi
    real(wp) :: y
    y = phi*(phi - 0.2_wp)*(phi - 1.0_wp)
  end function dv2_thick

  function example_v(self, x) result(y)
    class(example_potential), intent(inout) :: self
    real(wp), intent(in) :: x(:)
    real(wp) :: y
    real(wp) :: r1
    real(wp) :: r2
    real(wp) :: r3
    r1 = x(1)*x(1) + self%c*x(2)*x(2)
    r2 = self%c*(x(1) - 1.0_wp)**2 + (x(2) - 1.0_wp)**2
    r3 = self%fx*(0.25_wp*x(1)**4 - x(1)**3/3.0_wp)
    r3 = r3 + self%fy*(0.25_wp*x(2)**4 - x(2)**3/3.0_wp)
    y = r1*r2 + r3
  end function example_v

  function example_grad(self, x) result(dv)
    class(example_potential), intent(inout) :: self
    real(wp), intent(in) :: x(:)
    real(wp) :: dv(size(x))
    real(wp) :: r1
    real(wp) :: r2
    r1 = x(1)*x(1) + self%c*x(2)*x(2)
    r2 = self%c*(x(1) - 1.0_wp)**2 + (x(2) - 1.0_wp)**2
    dv(1) = r1*2.0_wp*self%c*(x(1) - 1.0_wp) + 2.0_wp*x(1)*r2  &
      + self%fx*x(1)*x(1)*(x(1) - 1.0_wp)
    dv(2) = r1*2.0_wp*(x(2) - 1.0_wp) + 2.0_wp*self%c*x(2)*r2  &
      + self%fy*x(2)*x(2)*(x(2) - 1.0_wp)
  end function example_grad

end module check__potentials

program check
!! Validation of the Fortran port against reference values produced with
!! the Python package (numpy 2.2.6 / scipy 1.16.1):
!!
!!   1d thin-walled:  action = 1.0927555625e+03
!!   1d thick-walled: action = 6.6489883968e+00
!!   2d thin-walled:  action = 1.7669376756e+03
!!   2d thick-walled: action = 4.5036661308e+00
!!
!! The tolerances are loose-ish (0.5%) because the Fortran port uses
!! different spline knots (bspline-fortran vs. FITPACK) and a slightly
!! different Simpson rule, so exact agreement is not expected.

  use cosmotransitions, only : wp
  use cosmotransitions, only : status_ok
  use cosmotransitions, only : pot1d_func
  use cosmotransitions, only : single_field_instanton
  use cosmotransitions, only : profile1d
  use cosmotransitions, only : full_tunneling
  use cosmotransitions, only : full_tunneling_result

  use check__potentials, only : v1_thin
  use check__potentials, only : dv1_thin
  use check__potentials, only : v2_thick
  use check__potentials, only : dv2_thick
  use check__potentials, only : example_potential

  implicit none

  integer :: nfail

  nfail = 0

  call test_1d()
  call test_2d()

  if (nfail > 0) then
    print "(a,i0,a)", "FAILED: ", nfail, " test(s) failed."
    error stop 1
  else
    print "(a)", "All tests passed."
  end if

contains

  subroutine test_1d()

    type(pot1d_func), target :: pot
    type(single_field_instanton) :: inst
    type(profile1d) :: prof
    real(wp) :: action
    integer :: status

    ! Thin-walled
    pot%vf => v1_thin
    pot%dvf => dv1_thin
    call inst%init(1.0_wp, 0.0_wp, pot, status)
    call assert_status("1d thin: init", status)
    call assert_close("1d thin: phi_bar", inst%phi_bar,  &
      8.37171431e-1_wp, 1.0e-4_wp)
    call assert_close("1d thin: rscale", inst%rscale,  &
      1.66770894_wp, 1.0e-4_wp)
    call inst%find_profile(prof, status)
    call assert_status("1d thin: find_profile", status)
    action = inst%find_action(prof)
    call assert_close("1d thin: action", action, 1.0927555625e3_wp, 5.0e-3_wp)

    ! Thick-walled
    pot%vf => v2_thick
    pot%dvf => dv2_thick
    call inst%init(1.0_wp, 0.0_wp, pot, status)
    call assert_status("1d thick: init", status)
    call inst%find_profile(prof, status)
    call assert_status("1d thick: find_profile", status)
    action = inst%find_action(prof)
    call assert_close("1d thick: action", action, 6.6489883968_wp, 5.0e-3_wp)
    call assert_close("1d thick: phi(0)", prof%phi(1),  &
      7.41882997e-1_wp, 5.0e-3_wp)

  end subroutine test_1d

  subroutine test_2d()

    type(example_potential), target :: pot
    type(full_tunneling_result) :: res
    real(wp) :: path_pts(2, 2)
    integer :: status
    integer :: n

    path_pts(1, :) = [1.0_wp, 1.0_wp]
    path_pts(2, :) = [0.0_wp, 0.0_wp]

    ! Thin-walled instanton
    pot = example_potential(c=5.0_wp, fx=0.0_wp, fy=2.0_wp)
    call full_tunneling(path_pts, pot, res, status, verbose=.true.)
    call assert_status("2d thin: full_tunneling", status)
    call assert_close("2d thin: action", res%action,  &
      1.7669376756e3_wp, 5.0e-3_wp)
    n = size(res%phi, 1)
    ! The path should start at the absolute minimum near (1, 1) and end at
    ! the metastable minimum near (0, 0).
    call assert_close("2d thin: path start x", res%phi(1, 1),  &
      1.0_wp, 1.0e-2_wp)
    call assert_close("2d thin: path start y", res%phi(1, 2),  &
      1.0_wp, 1.0e-2_wp)
    call assert_true("2d thin: path end x", abs(res%phi(n, 1)) < 1.0e-2_wp)
    call assert_true("2d thin: path end y", abs(res%phi(n, 2)) < 1.0e-2_wp)

    ! Thick-walled instanton
    pot = example_potential(c=5.0_wp, fx=0.0_wp, fy=80.0_wp)
    call full_tunneling(path_pts, pot, res, status, verbose=.true.)
    call assert_status("2d thick: full_tunneling", status)
    call assert_close("2d thick: action", res%action,  &
      4.5036661308_wp, 5.0e-3_wp)

  end subroutine test_2d

  subroutine assert_status(name, status)

    character(len=*), intent(in) :: name
    integer, intent(in) :: status

    if (status /= status_ok) then
      print "(a,a,a,i0)", "FAIL ", name, ": status = ", status
      nfail = nfail + 1
    end if

  end subroutine assert_status

  subroutine assert_close(name, val, ref, rtol)

    character(len=*), intent(in) :: name
    real(wp), intent(in) :: val
    real(wp), intent(in) :: ref
    real(wp), intent(in) :: rtol

    if (abs(val - ref) > rtol*abs(ref)) then
      print "(a,a,a,es16.8,a,es16.8,a,es9.2)", "FAIL ", name, ": got ",  &
        val, ", expected ", ref, ", rtol ", rtol
      nfail = nfail + 1
    else
      print "(a,a,a,es16.8,a,es16.8,a)", "ok   ", name, ": got ", val,  &
        " (ref ", ref, ")"
    end if

  end subroutine assert_close

  subroutine assert_true(name, cond)

    character(len=*), intent(in) :: name
    logical, intent(in) :: cond

    if (.not. cond) then
      print "(a,a)", "FAIL ", name
      nfail = nfail + 1
    else
      print "(a,a)", "ok   ", name
    end if

  end subroutine assert_true

end program check
