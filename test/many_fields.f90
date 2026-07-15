module many_fields__potential
!! Test potential with a variable number of field dimensions ndim:
!!
!!   V(x) = 0.4*phi^2 - 1.6*phi^3 + phi^4
!!          + 1/2 * sum_{i=2}^{ndim} m2_i * (x_i - a_i*phi^2)^2
!!
!! with phi = x(1), a_i = 0.04*(i-1) and m2_i = 2 + 0.2*(i-1). The
!! minima are known exactly by construction: the metastable minimum is
!! the origin (V = 0) and the true minimum is x = (1, a_2, ..., a_ndim)
!! (V = -0.2), with the barrier top along phi at phi = 0.2. The
!! transverse fields want to sit on the curve x_i = a_i*phi^2, so the
!! tunneling path bends in every transverse direction and the path
!! deformation is genuinely exercised.

  use cosmotransitions, only : wp
  use cosmotransitions, only : potential_nd

  implicit none

  private

  public :: many_fields_potential

  type, extends(potential_nd) :: many_fields_potential
    real(wp), allocatable :: a(:)
      !! Positions a_i of the transverse valley, i = 2..ndim.
    real(wp), allocatable :: m2(:)
      !! Transverse curvatures m2_i, i = 2..ndim.
  contains
    procedure :: init => pot_init
    procedure :: v => pot_v
    procedure :: grad => pot_grad
  end type many_fields_potential

contains

  subroutine pot_init(self, ndim)

    class(many_fields_potential), intent(inout) :: self
    integer, intent(in) :: ndim

    integer :: i

    self%a = 0.04e0_wp*real([(i, i = 1, ndim - 1)], wp)
    self%m2 = 2.0e0_wp + 0.2e0_wp*real([(i, i = 1, ndim - 1)], wp)

  end subroutine pot_init

  function pot_v(self, x) result(y)

    class(many_fields_potential), intent(inout) :: self
    real(wp), intent(in) :: x(:)
    real(wp) :: y

    real(wp) :: phi
    real(wp), allocatable :: d(:)

    phi = x(1)
    d = x(2:) - self%a*phi*phi
    y = 0.4e0_wp*phi**2 - 1.6e0_wp*phi**3 + phi**4  &
      + 0.5e0_wp*sum(self%m2*d*d)

  end function pot_v

  function pot_grad(self, x) result(dv)

    class(many_fields_potential), intent(inout) :: self
    real(wp), intent(in) :: x(:)
    real(wp) :: dv(size(x))

    real(wp) :: phi
    real(wp), allocatable :: d(:)

    phi = x(1)
    d = x(2:) - self%a*phi*phi
    dv(1) = 0.8e0_wp*phi - 4.8e0_wp*phi**2 + 4.0e0_wp*phi**3  &
      - 2.0e0_wp*phi*sum(self%m2*self%a*d)
    dv(2:) = self%m2*d

  end function pot_grad

end module many_fields__potential

program many_fields
!! Computes the O(3)-symmetric bounce action for the potential in
!! many_fields__potential with ndim = 2..20 field dimensions, starting
!! from the straight line between the two (exactly known) minima.
!!
!! For ndim = 2..10 the action is asserted against reference values from
!! the Python package (reference/ref_many_fields.py, numpy 2.2.6/
!! scipy 1.16.1). The Python package cannot go beyond ndim = 10: its
!! SplinePath fits the path with scipy.interpolate.splprep, whose
!! underlying FITPACK parcur routine is hard-limited to 10 dimensions
!! ("0 < idim < 11 must hold"). This port fits each field dimension as
!! its own 1d B-spline and has no such limit, which this test
!! demonstrates by continuing the scan to ndim = 20. For ndim > 10 the
!! solution is validated through internal consistency instead: the
!! deformation must converge, the Derrick-rescaled actions computed from
!! the potential and kinetic terms alone must agree with the full
!! action, and the action must keep growing with ndim (each additional
!! transverse field adds potential energy along the curved path, and the
!! observed growth per dimension of a few percent is far above the
!! deformation noise of order 1e-3).

  use cosmotransitions, only : wp
  use cosmotransitions, only : status_ok
  use cosmotransitions, only : full_tunneling
  use cosmotransitions, only : full_tunneling_result

  use many_fields__potential, only : many_fields_potential

  implicit none

  integer, parameter :: ndim_max = 20
  integer, parameter :: npath = 41

  ! Actions from reference/ref_many_fields.py, with the findProfile
  ! tolerances tightened to xtol = phitol = 1e-9 in both codes (the
  ! defaults of 1e-4 leave a few-1e-3 relative error in the action).
  ! The residual differences of up to ~1e-3 come from slightly
  ! different path-deformation trajectories (see test/double_well.f90);
  ! hence the relative tolerance of 2e-3 in the assertion.
  real(wp), parameter :: ref_action(2:10) = [  &
    3.3399851500e0_wp, 3.3265446619e0_wp, 3.3478390266e0_wp,  &
    3.3864753646e0_wp, 3.4488214698e0_wp, 3.5359964603e0_wp,  &
    3.6595890944e0_wp, 3.8232885005e0_wp, 4.0330848003e0_wp]

  type(many_fields_potential), target :: pot
  type(full_tunneling_result) :: res
  real(wp), allocatable :: x_true(:)
  real(wp), allocatable :: x_false(:)
  real(wp), allocatable :: path_pts(:, :)
  real(wp) :: t
  real(wp) :: action_prev
  integer :: ndim
  integer :: i
  integer :: status
  integer :: nfail
  character(len=16) :: label

  nfail = 0
  action_prev = 0.0e0_wp

  do ndim = 2, ndim_max
    print "(a,i0,a)", "--- ndim = ", ndim, " ---"
    call pot%init(ndim)

    ! The minima are known exactly by construction; verify.
    x_true = [1.0e0_wp, pot%a]
    x_false = 0.0e0_wp*x_true
    call assert_true("V(false minimum) = 0",  &
      abs(pot%v(x_false)) < 1.0e-13_wp)
    call assert_true("V(true minimum) = -0.2",  &
      abs(pot%v(x_true) + 0.2e0_wp) < 1.0e-13_wp)
    call assert_true("grad V(false minimum) = 0",  &
      maxval(abs(pot%grad(x_false))) < 1.0e-13_wp)
    call assert_true("grad V(true minimum) = 0",  &
      maxval(abs(pot%grad(x_true))) < 1.0e-13_wp)

    ! Straight-line initial guess from the true to the false minimum.
    if (allocated(path_pts)) deallocate(path_pts)
    allocate(path_pts(npath, ndim))
    do i = 1, npath
      t = real(i - 1, wp)/real(npath - 1, wp)
      path_pts(i, :) = x_true*(1.0e0_wp - t)
    end do

    call full_tunneling(path_pts, pot, res, status, verbose=.false.,  &
      xtol=1.0e-9_wp, phitol=1.0e-9_wp)
    call assert_true("full_tunneling status", status == status_ok)
    if (status /= status_ok) exit
    print "(a,es16.8,a,i0,a)", "bounce action S_3 = ", res%action,  &
      "   (", res%num_iters, " iterations)"
    print "(a,es16.8,a,es16.8)", "S_3 from potential/kinetic term only: ",  &
      res%action_pot, " /", res%action_kin
    print "(a,es10.2)", "fRatio = ", res%fratio
    call assert_true("deformation converged", res%converged)
    if (ndim <= 10) then
      call assert_close("bounce action S_3", res%action,  &
        ref_action(ndim), 2.0e-3_wp)
    end if
    ! Derrick consistency: for an exact solution of the bounce equation
    ! both rescaled actions equal the full action; at these tolerances
    ! they agree to ~1e-6.
    call assert_close("Derrick action (potential term)",  &
      res%action_pot, res%action, 1.0e-3_wp)
    call assert_close("Derrick action (kinetic term)",  &
      res%action_kin, res%action, 1.0e-3_wp)
    call assert_true("action grows with ndim",  &
      ndim <= 3 .or. res%action > action_prev)
    action_prev = res%action
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

    write (label, "(a,i0,a)") "[ndim=", ndim, "] "
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

    write (label, "(a,i0,a)") "[ndim=", ndim, "] "
    if (.not. cond) then
      print "(a,a,a)", "FAIL ", trim(label), name
      nfail = nfail + 1
    else
      print "(a,a,a)", "ok   ", trim(label), name
    end if

  end subroutine assert_true

end program many_fields
