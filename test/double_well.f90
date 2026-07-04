module double_well__potential
!! Two-field test potential for the finite-temperature bounce action:
!!
!!   V(x, y) = (x^2 + y^2) * ( a*(x-1)^2 + b*(y-1)^2 - c )
!!
!! The false (metastable) minimum sits at (0, 0) with V = 0; the true
!! minimum lies close to (1, 1) and has to be located numerically. The
!! coefficient c controls the depth of the true vacuum and is set at run
!! time, so the same potential can be scanned over a range of c values.

  use cosmotransitions, only : wp
  use cosmotransitions, only : potential_nd

  implicit none

  private

  public :: bounce_potential
  public :: v_plain

  real(wp), parameter, public :: a = 1.8e0_wp
  real(wp), parameter, public :: b = 0.2e0_wp
  real(wp), public :: c = 0.1e0_wp

  type, extends(potential_nd) :: bounce_potential
  contains
    procedure :: v => pot_v
    procedure :: grad => pot_grad
  end type bounce_potential

contains

  function pot_v(self, x) result(y)

    class(bounce_potential), intent(inout) :: self
    real(wp), intent(in) :: x(:)
    real(wp) :: y

    y = v_plain(x)

  end function pot_v

  function pot_grad(self, x) result(dv)

    class(bounce_potential), intent(inout) :: self
    real(wp), intent(in) :: x(:)
    real(wp) :: dv(size(x))

    real(wp) :: f
    real(wp) :: r

    f = a*(x(1) - 1.0e0_wp)**2 + b*(x(2) - 1.0e0_wp)**2 - c
    r = x(1)**2 + x(2)**2
    dv(1) = 2.0e0_wp*x(1)*f + r*2.0e0_wp*a*(x(1) - 1.0e0_wp)
    dv(2) = 2.0e0_wp*x(2)*f + r*2.0e0_wp*b*(x(2) - 1.0e0_wp)

  end function pot_grad

  function v_plain(phi) result(y)
    !! Plain-function version of the potential, as needed by the gradmin
    !! minimizer.

    real(wp), intent(in) :: phi(:)
    real(wp) :: y

    y = (phi(1)**2 + phi(2)**2)*(  &
      a*(phi(1) - 1.0e0_wp)**2 +  &
      b*(phi(2) - 1.0e0_wp)**2 - c)

  end function v_plain

end module double_well__potential

program double_well
!! Computes the finite-temperature bounce action S_3 for the two-field
!! potential in double_well__potential over a range of the coefficient c
!! between 0.01 and 0.3. The instanton is O(3) symmetric, i.e. alpha = 2,
!! which is the default of full_tunneling.
!!
!! For each c, the true minimum is first located with the gradmin Newton
!! minimizer. Reference values from the Python package (numpy 2.2.6/
!! scipy 1.16.1). Small c means a nearly degenerate (thin-walled) vacuum
!! with a large action; larger c gives a deeper true vacuum and easier
!! tunneling.

  use gradmin__newton, only : minimize

  use cosmotransitions, only : wp
  use cosmotransitions, only : status_ok
  use cosmotransitions, only : full_tunneling
  use cosmotransitions, only : full_tunneling_result

  use double_well__potential, only : bounce_potential
  use double_well__potential, only : v_plain
  use double_well__potential, only : c

  implicit none

  integer, parameter :: nc = 7

  real(wp), parameter :: c_vals(nc) = [  &
    0.01e0_wp, 0.05e0_wp, 0.10e0_wp, 0.15e0_wp,  &
    0.20e0_wp, 0.25e0_wp, 0.30e0_wp]
  real(wp), parameter :: ref_xmin(nc) = [  &
    1.002674189573_wp, 1.011663435734_wp, 1.020199200731_wp,  &
    1.026807410988_wp, 1.032131164050_wp, 1.036548270078_wp,  &
    1.040296717098_wp]
  real(wp), parameter :: ref_ymin(nc) = [  &
    1.024593855202_wp, 1.115773451273_wp, 1.216831369462_wp,  &
    1.307134557354_wp, 1.389232396274_wp, 1.464850194264_wp,  &
    1.535207198375_wp]
  ! Actions from reference/ref_double_well.py, with the findProfile
  ! tolerances tightened to xtol = phitol = 1e-9 in both codes (the
  ! defaults of 1e-4 leave a few-1e-3 relative error in the action).
  ! The residual differences of up to ~1e-3 come from slightly
  ! different path-deformation trajectories: the deformation stops at a
  ! path with fRatio of a few percent, and which such path is reached
  ! depends on the implementation, the compiler and the optimization
  ! level (the action varies at order fRatio**2 among them). Hence the
  ! relative tolerance of 2e-3 in the assertion.
  real(wp), parameter :: ref_action(nc) = [  &
    7.8805837752e3_wp, 3.8913237834e2_wp, 1.1963573348e2_wp,  &
    6.2230671791e1_wp, 3.9636887635e1_wp, 2.7922475045e1_wp,  &
    2.0887209824e1_wp]

  type(bounce_potential), target :: pot
  type(full_tunneling_result) :: res
  real(wp) :: x0(2)
  real(wp) :: xtrue(2)
  real(wp) :: vtrue
  real(wp) :: path_pts(2, 2)
  integer :: status
  integer :: nfail
  integer :: i
  character(len=16) :: label

  nfail = 0

  do i = 1, nc
    c = c_vals(i)
    print "(a,f4.2,a)", "--- c = ", c, " ---"

    ! Find the precise location of the true minimum near (1, 1).
    x0 = [1.0e0_wp, 1.0e0_wp]
    status = 0
    call minimize(v_plain, x0, xtrue, fmin=vtrue, status=status,  &
      xtol=1.0e-12_wp)
    call assert_true("minimization status", status == 0)
    print "(a,2f16.12)", "true minimum at ", xtrue
    print "(a,es16.8)", "V(true minimum) = ", vtrue
    call assert_close("true minimum x", xtrue(1), ref_xmin(i), 1.0e-6_wp)
    call assert_close("true minimum y", xtrue(2), ref_ymin(i), 1.0e-6_wp)
    call assert_true("V(true) < V(false)", vtrue < 0.0e0_wp)

    ! Finite-temperature bounce: the path starts at the true minimum and
    ! ends at the false minimum (0, 0).
    path_pts(1, :) = xtrue
    path_pts(2, :) = [0.0e0_wp, 0.0e0_wp]
    call full_tunneling(path_pts, pot, res, status, verbose=.false.,  &
      xtol=1.0e-9_wp, phitol=1.0e-9_wp)
    call assert_true("full_tunneling status", status == status_ok)
    print "(a,es16.8)", "bounce action S_3 = ", res%action
    print "(a,es16.8,a,es16.8)", "S_3 from potential/kinetic term only: ",  &
      res%action_pot, " /", res%action_kin
    print "(a,es10.2)", "fRatio = ", res%fratio
    call assert_close("bounce action S_3", res%action,  &
      ref_action(i), 2.0e-3_wp)
  end do

  if (nfail > 0) then
    print "(a,i0,a)", "FAILED: ", nfail, " test(s) failed."
    error stop 1
  else
    print "(a)", "All tests passed."
  end if

contains

  subroutine assert_close(name, val, ref, rtol)

    character(len=*), intent(in) :: name
    real(wp), intent(in) :: val
    real(wp), intent(in) :: ref
    real(wp), intent(in) :: rtol

    write (label, "(a,f4.2,a)") "[c=", c, "] "
    if (abs(val - ref) > rtol*abs(ref)) then
      print "(a,a,a,a,es16.8,a,es16.8,a,es9.2)", "FAIL ", trim(label),  &
        name, ": got ", val, ", expected ", ref, ", rtol ", rtol
      nfail = nfail + 1
    else
      print "(a,a,a,a,es16.8,a,es16.8,a)", "ok   ", trim(label), name,  &
        ": got ", val, " (ref ", ref, ")"
    end if

  end subroutine assert_close

  subroutine assert_true(name, cond)

    character(len=*), intent(in) :: name
    logical, intent(in) :: cond

    write (label, "(a,f4.2,a)") "[c=", c, "] "
    if (.not. cond) then
      print "(a,a,a)", "FAIL ", trim(label), name
      nfail = nfail + 1
    else
      print "(a,a,a)", "ok   ", trim(label), name
    end if

  end subroutine assert_true

end program double_well
