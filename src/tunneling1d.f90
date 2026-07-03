module cosmotransitions__tunneling1d
!! Port of cosmoTransitions/tunneling1D.py (class SingleFieldInstanton):
!! everything needed to calculate instantons in one field dimension using
!! the overshoot/undershoot method.

  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

  use cosmotransitions__config, only : wp
  use cosmotransitions__config, only : pi
  use cosmotransitions__config, only : status_ok
  use cosmotransitions__config, only : err_integration
  use cosmotransitions__config, only : err_no_barrier
  use cosmotransitions__config, only : err_stable
  use cosmotransitions__config, only : err_numerical
  use cosmotransitions__helpers, only : rkqs
  use cosmotransitions__helpers, only : cubic_interp
  use cosmotransitions__helpers, only : linspace
  use cosmotransitions__helpers, only : simpson
  use cosmotransitions__helpers, only : interp_linear
  use cosmotransitions__helpers, only : monotonic_indices
  use cosmotransitions__optimize, only : brentq
  use cosmotransitions__optimize, only : minimize_bounded
  use cosmotransitions__special, only : besseli_nu
  use cosmotransitions__special, only : besselj_nu
  use cosmotransitions__potentials, only : potential_1d

  implicit none

  private

  public :: profile1d
  public :: single_field_instanton
  public :: ctype_converged
  public :: ctype_undershoot
  public :: ctype_overshoot

  integer, parameter :: ctype_converged = 0
  integer, parameter :: ctype_undershoot = 1
  integer, parameter :: ctype_overshoot = 2

  type :: profile1d
    !! The bubble profile, i.e. the namedtuple "Profile1D" of the Python
    !! package.
    real(wp), allocatable :: r(:)
    real(wp), allocatable :: phi(:)
    real(wp), allocatable :: dphi(:)
    real(wp) :: rerr = -1.0_wp
      !! First value of r at which dr < drmin, or negative if that never
      !! happened (None in Python).
  end type profile1d

  type :: single_field_instanton
    !! Calculates the properties of an instanton with a single scalar field
    !! (without gravity) using the overshoot/undershoot method.
    !! Port of tunneling1D.SingleFieldInstanton.
    real(wp) :: phi_absmin
      !! Field value of the stable vacuum to which the field tunnels.
    real(wp) :: phi_metamin
      !! Field value of the metastable vacuum.
    real(wp) :: phi_bar
      !! Field value at the edge of the barrier.
    real(wp) :: rscale
      !! Approximate radial scale of the instanton.
    real(wp) :: alpha = 2.0_wp
      !! Friction coefficient in the ODE (spacetime dimensions minus 1).
    real(wp) :: phi_eps
      !! Absolute step used in dv_from_absmin (phi_eps_rel rescaled by
      !! |phi_absmin - phi_metamin|).
    class(potential_1d), pointer :: pot => null()
      !! The potential V(phi) with its derivatives.
  contains
    procedure :: init => sfi_init
    procedure :: dv_from_absmin => sfi_dv_from_absmin
    procedure :: find_barrier_location => sfi_find_barrier_location
    procedure :: find_rscale => sfi_find_rscale
    procedure :: exact_solution => sfi_exact_solution
    procedure :: initial_conditions => sfi_initial_conditions
    procedure :: equation_of_motion => sfi_equation_of_motion
    procedure :: integrate_profile => sfi_integrate_profile
    procedure :: integrate_and_save_profile => sfi_integrate_and_save_profile
    procedure :: find_profile => sfi_find_profile
    procedure :: find_action => sfi_find_action
    procedure :: evenly_spaced_phi => sfi_evenly_spaced_phi
  end type single_field_instanton

contains

  subroutine sfi_init(self, phi_absmin, phi_metamin, pot, status,  &
      alpha, phi_eps_rel, phi_bar, rscale)
    !! Port of SingleFieldInstanton.__init__.
    !!
    !! Fails with `err_stable` when the presumably stable minimum has a
    !! higher energy than the metastable one, and with `err_no_barrier`
    !! when there is no barrier between the minima.

    class(single_field_instanton), intent(inout) :: self
    real(wp), intent(in) :: phi_absmin
    real(wp), intent(in) :: phi_metamin
    class(potential_1d), intent(inout), target :: pot
      !! The potential. The actual argument must have the `target` (or
      !! `pointer`) attribute and outlive this object.
    integer, intent(out) :: status
    real(wp), intent(in), optional :: alpha
      !! Friction coefficient (default 2).
    real(wp), intent(in), optional :: phi_eps_rel
      !! Small unitless value used for derivatives near the minimum
      !! (default 1e-3); rescaled by |phi_absmin - phi_metamin|.
    real(wp), intent(in), optional :: phi_bar
      !! Field value at the barrier edge; found automatically if absent.
    real(wp), intent(in), optional :: rscale
      !! Radial scale of the instanton; found automatically if absent.

    real(wp) :: phi_eps_rel_

    status = status_ok
    self%phi_absmin = phi_absmin
    self%phi_metamin = phi_metamin
    self%pot => pot
    self%alpha = 2.0_wp
    if (present(alpha)) self%alpha = alpha
    phi_eps_rel_ = 1.0e-3_wp
    if (present(phi_eps_rel)) phi_eps_rel_ = phi_eps_rel

    if (pot%v(phi_metamin) <= pot%v(phi_absmin)) then
      ! "V(phi_metaMin) <= V(phi_absMin); tunneling cannot occur."
      status = err_stable
      return
    end if

    self%phi_eps = phi_eps_rel_*abs(phi_absmin - phi_metamin)
    ! Make the default finite-difference derivatives of the potential use
    ! the same step as the Python class does.
    pot%fd_eps = self%phi_eps

    if (present(phi_bar)) then
      self%phi_bar = phi_bar
    else
      self%phi_bar = self%find_barrier_location()
    end if
    if (present(rscale)) then
      self%rscale = rscale
    else
      self%rscale = self%find_rscale(status)
      if (status /= status_ok) return
    end if

  end subroutine sfi_init

  function sfi_dv_from_absmin(self, delta_phi) result(dv)
    !! dV/dphi at phi = phi_absmin + delta_phi, blending the direct
    !! derivative with d2V*delta_phi close to the minimum for numerical
    !! stability. Port of SingleFieldInstanton.dV_from_absMin.

    class(single_field_instanton), intent(inout) :: self
    real(wp), intent(in) :: delta_phi
    real(wp) :: dv

    real(wp) :: phi
    real(wp) :: dv_
    real(wp) :: blend_factor

    phi = self%phi_absmin + delta_phi
    dv = self%pot%dv(phi)
    if (self%phi_eps > 0.0_wp) then
      dv_ = self%pot%d2v(phi)*delta_phi
      blend_factor = exp(-(delta_phi/self%phi_eps)**2)
      dv = dv_*blend_factor + dv*(1.0_wp - blend_factor)
    end if

  end function sfi_dv_from_absmin

  function sfi_find_barrier_location(self) result(phi0)
    !! The field value phi0 between the minima such that
    !! V(phi0) = V(phi_metamin), found by binary search.
    !! Port of SingleFieldInstanton.findBarrierLocation.

    class(single_field_instanton), intent(inout) :: self
    real(wp) :: phi0

    real(wp) :: phi_tol
    real(wp) :: v_phimeta
    real(wp) :: phi1
    real(wp) :: phi2
    real(wp) :: v0

    phi_tol = abs(self%phi_metamin - self%phi_absmin)*1.0e-12_wp
    v_phimeta = self%pot%v(self%phi_metamin)
    phi1 = self%phi_metamin
    phi2 = self%phi_absmin
    phi0 = 0.5_wp*(phi1 + phi2)

    do while (abs(phi1 - phi2) > phi_tol)
      v0 = self%pot%v(phi0)
      if (v0 > v_phimeta) then
        phi1 = phi0
      else
        phi2 = phi0
      end if
      phi0 = 0.5_wp*(phi1 + phi2)
    end do

  end function sfi_find_barrier_location

  function sfi_find_rscale(self, status) result(rscale)
    !! Characteristic length scale for tunneling over the barrier: the
    !! period of oscillations about the top of a cubic fitted to the
    !! barrier. Port of SingleFieldInstanton.findRScale.

    class(single_field_instanton), intent(inout) :: self
    integer, intent(out) :: status
    real(wp) :: rscale

    real(wp) :: phi_tol
    real(wp) :: x1
    real(wp) :: x2
    real(wp) :: phi_bar_top
    real(wp) :: vtop
    real(wp) :: xtop
    integer :: st

    status = status_ok
    rscale = 0.0_wp
    phi_tol = abs(self%phi_bar - self%phi_metamin)*1.0e-6_wp
    x1 = min(self%phi_bar, self%phi_metamin)
    x2 = max(self%phi_bar, self%phi_metamin)
    phi_bar_top = minimize_bounded(neg_v, x1, x2, phi_tol, st)
    if (st /= status_ok) then
      status = err_numerical
      return
    end if
    if (phi_bar_top + phi_tol > x2 .or. phi_bar_top - phi_tol < x1) then
      ! "Minimization is placing the top of the potential barrier outside
      ! of the interval defined by phi_bar and phi_metaMin. Assume that the
      ! barrier does not exist."
      status = err_no_barrier
      return
    end if

    vtop = self%pot%v(phi_bar_top) - self%pot%v(self%phi_metamin)
    xtop = phi_bar_top - self%phi_metamin
    if (vtop <= 0.0_wp) then
      ! "Barrier height is not positive, does not exist."
      status = err_no_barrier
      return
    end if
    rscale = abs(xtop)/sqrt(abs(6.0_wp*vtop))

  contains

    function neg_v(x) result(y)
      real(wp), intent(in) :: x
      real(wp) :: y
      y = -self%pot%v(x)
    end function neg_v

  end function sfi_find_rscale

  subroutine sfi_exact_solution(self, r, phi0, dv, d2v, phi, dphi)
    !! phi(r) and dphi(r) at radius `r` given phi(r=0) = phi0, assuming a
    !! quadratic potential with derivatives dv and d2v at phi0.
    !! Port of SingleFieldInstanton.exactSolution. The solution involves
    !! (modified) Bessel functions of order nu = (alpha-1)/2.

    class(single_field_instanton), intent(inout) :: self
    real(wp), intent(in) :: r
    real(wp), intent(in) :: phi0
    real(wp), intent(in) :: dv
    real(wp), intent(in) :: d2v
    real(wp), intent(out) :: phi
    real(wp), intent(out) :: dphi

    real(wp) :: beta
    real(wp) :: beta_r
    real(wp) :: nu
    real(wp) :: s
    real(wp) :: term
    real(wp) :: iv
    integer :: k

    beta = sqrt(abs(d2v))
    beta_r = beta*r
    nu = 0.5_wp*(self%alpha - 1.0_wp)

    if (beta_r < 1.0e-2_wp) then
      ! Use a small-r approximation for the Bessel function.
      s = merge(1.0_wp, -1.0_wp, d2v > 0.0_wp)
      phi = 0.0_wp
      dphi = 0.0_wp
      do k = 1, 3
        term = (0.5_wp*beta_r)**(2*k - 2)*s**k  &
          /(gamma(real(k, wp) + 1.0_wp)*gamma(real(k, wp) + 1.0_wp + nu))
        phi = phi + term
        dphi = dphi + term*real(2*k, wp)
      end do
      phi = phi*0.25_wp*gamma(nu + 1.0_wp)*r*r*dv*s
      dphi = dphi*0.25_wp*gamma(nu + 1.0_wp)*r*dv*s
      phi = phi + phi0
    else if (d2v > 0.0_wp) then
      iv = besseli_nu(nu, beta_r)
      phi = (gamma(nu + 1.0_wp)*(0.5_wp*beta_r)**(-nu)*iv - 1.0_wp)*dv/d2v
      dphi = -nu*((0.5_wp*beta_r)**(-nu)/r)*iv
      dphi = dphi + (0.5_wp*beta_r)**(-nu)*0.5_wp*beta  &
        *(besseli_nu(nu - 1.0_wp, beta_r) + besseli_nu(nu + 1.0_wp, beta_r))
      dphi = dphi*gamma(nu + 1.0_wp)*dv/d2v
      phi = phi + phi0
    else
      iv = besselj_nu(nu, beta_r)
      phi = (gamma(nu + 1.0_wp)*(0.5_wp*beta_r)**(-nu)*iv - 1.0_wp)*dv/d2v
      dphi = -nu*((0.5_wp*beta_r)**(-nu)/r)*iv
      dphi = dphi + (0.5_wp*beta_r)**(-nu)*0.5_wp*beta  &
        *(besselj_nu(nu - 1.0_wp, beta_r) - besselj_nu(nu + 1.0_wp, beta_r))
      dphi = dphi*gamma(nu + 1.0_wp)*dv/d2v
      phi = phi + phi0
    end if

  end subroutine sfi_exact_solution

  subroutine sfi_initial_conditions(self, delta_phi0, rmin,  &
      delta_phi_cutoff, r0, phi_r0, dphi_r0, status)
    !! Finds the initial conditions for integration: the value r0 such that
    !! phi(r0) = phi_absmin + delta_phi_cutoff (or the conditions at rmin
    !! if no such value exists). Port of
    !! SingleFieldInstanton.initialConditions.

    class(single_field_instanton), intent(inout) :: self
    real(wp), intent(in) :: delta_phi0
      !! phi(r=0) - phi_absmin.
    real(wp), intent(in) :: rmin
      !! The smallest acceptable radius at which to start integration.
    real(wp), intent(in) :: delta_phi_cutoff
      !! The desired phi(r0) - phi_absmin.
    real(wp), intent(out) :: r0
    real(wp), intent(out) :: phi_r0
    real(wp), intent(out) :: dphi_r0
    integer, intent(out) :: status

    real(wp) :: phi0
    real(wp) :: dv
    real(wp) :: d2v
    real(wp) :: r
    real(wp) :: rlast
    real(wp) :: phi
    real(wp) :: dphi
    integer :: st

    status = status_ok
    phi0 = self%phi_absmin + delta_phi0
    dv = self%dv_from_absmin(delta_phi0)
    d2v = self%pot%d2v(phi0)
    call self%exact_solution(rmin, phi0, dv, d2v, phi_r0, dphi_r0)
    r0 = rmin
    if (abs(phi_r0 - self%phi_absmin) > abs(delta_phi_cutoff)) then
      ! The initial conditions at rmin work. Stop here.
      return
    end if
    if (sgn(dphi_r0) /= sgn(delta_phi0)) then
      ! The field is evolving in the wrong direction. Increasing r0 won't
      ! increase |delta_phi_r0|.
      return
    end if

    ! Find the smallest r0 such that delta_phi_r0 > delta_phi_cutoff.
    r = rmin
    do
      rlast = r
      r = r*10.0_wp
      if (r > 1.0e250_wp) then
        ! In Python r runs off to infinity here; treat it as a failure to
        ! find initial conditions.
        status = err_numerical
        return
      end if
      call self%exact_solution(r, phi0, dv, d2v, phi, dphi)
      if (abs(phi - self%phi_absmin) > abs(delta_phi_cutoff)) exit
    end do

    ! Now find where phi - phi_absmin = delta_phi_cutoff exactly.
    r0 = brentq(delta_phi_diff, rlast, r, st)
    if (st /= status_ok) then
      status = err_numerical
      return
    end if
    call self%exact_solution(r0, phi0, dv, d2v, phi_r0, dphi_r0)

  contains

    function delta_phi_diff(r_) result(y)
      real(wp), intent(in) :: r_
      real(wp) :: y
      real(wp) :: p
      real(wp) :: dp
      call self%exact_solution(r_, phi0, dv, d2v, p, dp)
      y = abs(p - self%phi_absmin) - abs(delta_phi_cutoff)
    end function delta_phi_diff

  end subroutine sfi_initial_conditions

  subroutine sfi_equation_of_motion(self, y, r, dydr)
    !! The bubble wall equation of motion:
    !! d2phi/dr2 + (alpha/r) dphi/dr = dV/dphi.

    class(single_field_instanton), intent(inout) :: self
    real(wp), intent(in) :: y(:)
    real(wp), intent(in) :: r
    real(wp), intent(out) :: dydr(:)

    dydr(1) = y(2)
    dydr(2) = self%pot%dv(y(1)) - self%alpha*y(2)/r

  end subroutine sfi_equation_of_motion

  subroutine sfi_integrate_profile(self, r0_in, y0_in, dr0, epsfrac,  &
      epsabs, drmin, rmax, r, y, convergence_type, status)
    !! Integrate the bubble wall equation until the field either overshoots
    !! or undershoots the false vacuum, or converges on it.
    !! Port of SingleFieldInstanton.integrateProfile.

    class(single_field_instanton), intent(inout) :: self
    real(wp), intent(in) :: r0_in
    real(wp), intent(in) :: y0_in(2)
      !! Starting values [phi(r0), dphi(r0)].
    real(wp), intent(in) :: dr0
    real(wp), intent(in) :: epsfrac(2)
    real(wp), intent(in) :: epsabs(2)
    real(wp), intent(in) :: drmin
    real(wp), intent(in) :: rmax
      !! Maximum allowed value of r - r0 before raising an error.
    real(wp), intent(out) :: r
    real(wp), intent(out) :: y(2)
    integer, intent(out) :: convergence_type
      !! One of ctype_converged, ctype_undershoot, ctype_overshoot.
    integer, intent(out) :: status

    real(wp) :: r0
    real(wp) :: r1
    real(wp) :: dr
    real(wp) :: drnext
    real(wp) :: dr_did
    real(wp) :: ysign
    real(wp) :: rmax_tot
    real(wp) :: x
    real(wp), dimension(2) :: y0
    real(wp), dimension(2) :: y1
    real(wp), dimension(2) :: dy
    real(wp), dimension(2) :: dydr0
    real(wp), dimension(2) :: dydr1
    integer :: st

    status = status_ok
    convergence_type = -1
    r0 = r0_in
    y0 = y0_in
    dr = dr0
    call self%equation_of_motion(y0, r0, dydr0)
    ysign = sign(1.0_wp, y0(1) - self%phi_metamin)
      ! Positive means we're heading down, negative means heading up.
    rmax_tot = rmax + r0

    do
      call rkqs(y0, dydr0, r0, deriv, dr, epsfrac, epsabs, dy, dr_did,  &
        drnext, st)
      if (st /= status_ok) then
        status = st
        return
      end if
      r1 = r0 + dr_did
      y1 = y0 + dy
      call self%equation_of_motion(y1, r1, dydr1)

      ! Check for completion
      if (r1 > rmax_tot) then
        ! "r > rmax"
        status = err_integration
        return
      else if (dr_did < drmin) then
        ! "dr < drmin"
        status = err_integration
        return
      else if (abs(y1(1) - self%phi_metamin) < 3.0_wp*epsabs(1) .and.  &
          abs(y1(2)) < 3.0_wp*epsabs(2)) then
        r = r1
        y = y1
        convergence_type = ctype_converged
        exit
      else if (y1(2)*ysign > 0.0_wp .or.  &
          (y1(1) - self%phi_metamin)*ysign < 0.0_wp) then
        ! Extrapolate with a cubic interpolation over the last step.
        if (y1(2)*ysign > 0.0_wp) then
          ! Extrapolate to where dphi(r) = 0.
          x = brentq(interp_dphi, 0.0_wp, 1.0_wp, st)
          convergence_type = ctype_undershoot
        else
          ! Extrapolate to where phi(r) = phi_metamin.
          x = brentq(interp_phi_diff, 0.0_wp, 1.0_wp, st)
          convergence_type = ctype_overshoot
        end if
        if (st /= status_ok) then
          status = err_numerical
          return
        end if
        r = r0 + dr_did*x
        y = cubic_interp(x, y0, dr_did*dydr0, y1, dr_did*dydr1)
        exit
      end if
      ! Advance the integration variables.
      r0 = r1
      y0 = y1
      dydr0 = dydr1
      dr = drnext
    end do

    ! Check convergence for a second time. The extrapolation in
    ! overshoot/undershoot might have gotten us within the acceptable error.
    if (abs(y(1) - self%phi_metamin) < 3.0_wp*epsabs(1) .and.  &
        abs(y(2)) < 3.0_wp*epsabs(2)) then
      convergence_type = ctype_converged
    end if

  contains

    subroutine deriv(y_, r_, dydr_)
      real(wp), intent(in) :: y_(:)
      real(wp), intent(in) :: r_
      real(wp), intent(out) :: dydr_(:)
      call self%equation_of_motion(y_, r_, dydr_)
    end subroutine deriv

    function interp_dphi(x_) result(f_)
      real(wp), intent(in) :: x_
      real(wp) :: f_
      real(wp) :: yi(2)
      yi = cubic_interp(x_, y0, dr_did*dydr0, y1, dr_did*dydr1)
      f_ = yi(2)
    end function interp_dphi

    function interp_phi_diff(x_) result(f_)
      real(wp), intent(in) :: x_
      real(wp) :: f_
      real(wp) :: yi(2)
      yi = cubic_interp(x_, y0, dr_did*dydr0, y1, dr_did*dydr1)
      f_ = yi(1) - self%phi_metamin
    end function interp_phi_diff

  end subroutine sfi_integrate_profile

  subroutine sfi_integrate_and_save_profile(self, r_array, y0_in, dr0,  &
      epsfrac, epsabs, drmin, profile, status)
    !! Integrate the bubble profile, saving the output at the radii given
    !! in `r_array`. Port of SingleFieldInstanton.integrateAndSaveProfile.

    class(single_field_instanton), intent(inout) :: self
    real(wp), intent(in) :: r_array(:)
    real(wp), intent(in) :: y0_in(2)
    real(wp), intent(in) :: dr0
    real(wp), intent(in) :: epsfrac(2)
    real(wp), intent(in) :: epsabs(2)
    real(wp), intent(in) :: drmin
    type(profile1d), intent(out) :: profile
    integer, intent(out) :: status

    integer :: n
    integer :: i
    real(wp) :: r0
    real(wp) :: r1
    real(wp) :: dr
    real(wp) :: dr_did
    real(wp) :: drnext
    real(wp) :: x
    real(wp), dimension(2) :: y0
    real(wp), dimension(2) :: y1
    real(wp), dimension(2) :: dy
    real(wp), dimension(2) :: dydr0
    real(wp), dimension(2) :: dydr1
    real(wp), dimension(2) :: yx
    real(wp), allocatable :: yout(:, :)
    integer :: st

    status = status_ok
    n = size(r_array)
    r0 = r_array(1)
    allocate(yout(n, 2))
    yout(1, :) = y0_in
    y0 = y0_in
    dr = dr0
    call self%equation_of_motion(y0, r0, dydr0)
    ! Note: the Python version contains "if Rerr is not None: Rerr = r1"
    ! inside the small-step branch, which never fires because Rerr starts
    ! as None. We reproduce that behaviour by leaving profile%rerr at its
    ! sentinel value.

    i = 2
    do while (i <= n)
      call rkqs(y0, dydr0, r0, deriv, dr, epsfrac, epsabs, dy, dr_did,  &
        drnext, st)
      if (st /= status_ok) then
        status = st
        return
      end if
      if (dr_did >= drmin) then
        r1 = r0 + dr_did
        y1 = y0 + dy
      else
        y1 = y0 + dy*drmin/dr_did
        dr_did = drmin
        drnext = drmin
        r1 = r0 + dr_did
      end if
      call self%equation_of_motion(y1, r1, dydr1)
      ! Fill the arrays, if necessary.
      do while (i <= n)
        if (r0 < r_array(i) .and. r_array(i) <= r1) then
          x = (r_array(i) - r0)/dr_did
          yx = cubic_interp(x, y0, dr_did*dydr0, y1, dr_did*dydr1)
          yout(i, :) = yx
          i = i + 1
        else
          exit
        end if
      end do
      ! Advance the integration variables.
      r0 = r1
      y0 = y1
      dydr0 = dydr1
      dr = drnext
    end do

    profile%r = r_array
    profile%phi = yout(:, 1)
    profile%dphi = yout(:, 2)
    profile%rerr = -1.0_wp

  contains

    subroutine deriv(y_, r_, dydr_)
      real(wp), intent(in) :: y_(:)
      real(wp), intent(in) :: r_
      real(wp), intent(out) :: dydr_(:)
      call self%equation_of_motion(y_, r_, dydr_)
    end subroutine deriv

  end subroutine sfi_integrate_and_save_profile

  subroutine sfi_find_profile(self, profile, status, xguess, xtol,  &
      phitol, thin_cutoff, npoints, rmin_rel, rmax_rel, max_interior_pts)
    !! Calculate the bubble profile by iteratively over/undershooting.
    !! Port of SingleFieldInstanton.findProfile. Rather than varying
    !! phi(r=0) directly, the parameter x defined by
    !! phi(r=0) = phi_absmin + exp(-x)*(phi_metamin - phi_absmin)
    !! is varied.

    class(single_field_instanton), intent(inout) :: self
    type(profile1d), intent(out) :: profile
    integer, intent(out) :: status
    real(wp), intent(in), optional :: xguess
      !! Initial guess for x. By default set from phi_bar.
    real(wp), intent(in), optional :: xtol
      !! Target accuracy in x (default 1e-4).
    real(wp), intent(in), optional :: phitol
      !! Fractional error tolerance in integration (default 1e-4).
    real(wp), intent(in), optional :: thin_cutoff
      !! delta_phi_cutoff/(phi_metamin - phi_absmin) used in
      !! initial_conditions (default 0.01).
    integer, intent(in), optional :: npoints
      !! Number of points to return in the profile (default 500).
    real(wp), intent(in), optional :: rmin_rel
      !! Smallest starting radius relative to rscale (default 1e-4).
    real(wp), intent(in), optional :: rmax_rel
      !! Maximum integration distance relative to rscale (default 1e4).
    integer, intent(in), optional :: max_interior_pts
      !! Maximum number of points to place between r=0 and the start of
      !! integration. Defaults to npoints/2; 0 disables interior points.

    real(wp) :: xtol_
    real(wp) :: phitol_
    real(wp) :: thin_cutoff_
    integer :: npoints_
    real(wp) :: rmin
    real(wp) :: rmax
    integer :: max_interior_
    real(wp) :: xmin
    real(wp) :: xmax
    logical :: has_xmax
    real(wp) :: x
    real(wp), parameter :: xincrease = 5.0_wp
    real(wp) :: dr0
    real(wp) :: drmin
    real(wp) :: delta_phi
    real(wp), dimension(2) :: epsabs
    real(wp), dimension(2) :: epsfrac
    real(wp) :: delta_phi_cutoff
    real(wp) :: delta_phi0
    real(wp) :: r0_
    real(wp) :: r0
    real(wp) :: phi0
    real(wp) :: dphi0
    real(wp) :: rf
    logical :: have_rf
    real(wp), dimension(2) :: y0
    real(wp), dimension(2) :: yf
    integer :: ctype
    integer :: st
    integer :: n_int
    integer :: i
    real(wp) :: dx0
    real(wp) :: a
    real(wp) :: nk
    real(wp), allocatable :: r_arr(:)
    real(wp), allocatable :: r_int(:)
    real(wp), allocatable :: phi_int(:)
    real(wp), allocatable :: dphi_int(:)
    real(wp) :: dv
    real(wp) :: d2v
    type(profile1d) :: prof_out

    status = status_ok
    xtol_ = 1.0e-4_wp
    phitol_ = 1.0e-4_wp
    thin_cutoff_ = 0.01_wp
    npoints_ = 500
    if (present(xtol)) xtol_ = xtol
    if (present(phitol)) phitol_ = phitol
    if (present(thin_cutoff)) thin_cutoff_ = thin_cutoff
    if (present(npoints)) npoints_ = npoints

    ! Set x parameters
    xmin = xtol_*10.0_wp
    xmax = 0.0_wp
    has_xmax = .false.
    if (present(xguess)) then
      x = xguess
    else
      x = -log(abs((self%phi_bar - self%phi_absmin)  &
        /(self%phi_metamin - self%phi_absmin)))
    end if
    ! Set r parameters
    rmin = 1.0e-4_wp
    if (present(rmin_rel)) rmin = rmin_rel
    rmin = rmin*self%rscale
    dr0 = rmin
    drmin = 0.01_wp*rmin
    rmax = 1.0e4_wp
    if (present(rmax_rel)) rmax = rmax_rel
    rmax = rmax*self%rscale
    ! Set phi parameters
    delta_phi = self%phi_metamin - self%phi_absmin
    epsabs(1) = abs(delta_phi*phitol_)
    epsabs(2) = abs(delta_phi*phitol_/self%rscale)
    epsfrac = phitol_
    delta_phi_cutoff = thin_cutoff_*delta_phi

    have_rf = .false.
    rf = 0.0_wp
    r0 = rmin
    y0 = 0.0_wp
    delta_phi0 = 0.0_wp
    do
      delta_phi0 = exp(-x)*delta_phi
      call self%initial_conditions(delta_phi0, rmin, delta_phi_cutoff,  &
        r0_, phi0, dphi0, st)
      if (st /= status_ok .or. .not. ieee_is_finite(r0_)  &
          .or. .not. ieee_is_finite(x)) then
        ! Use the last finite values instead (assuming there are such
        ! values).
        if (.not. have_rf) then
          ! "Failed to retrieve initial conditions on the first try."
          status = err_numerical
          return
        end if
        exit
      end if
      r0 = r0_
      y0(1) = phi0
      y0(2) = dphi0
      call self%integrate_profile(r0, y0, dr0, epsfrac, epsabs, drmin,  &
        rmax, rf, yf, ctype, st)
      if (st /= status_ok) then
        status = st
        return
      end if
      have_rf = .true.
      ! Check for overshoot, undershoot
      if (ctype == ctype_converged) then
        exit
      else if (ctype == ctype_undershoot) then
        ! x is too low
        xmin = x
        if (has_xmax) then
          x = 0.5_wp*(xmin + xmax)
        else
          x = x*xincrease
        end if
      else if (ctype == ctype_overshoot) then
        ! x is too high
        xmax = x
        has_xmax = .true.
        x = 0.5_wp*(xmin + xmax)
      end if
      ! Check if we've reached xtol
      if (has_xmax) then
        if (xmax - xmin < xtol_) exit
      end if
    end do

    ! Integrate a second time, this time getting the points along the way.
    r_arr = linspace(r0, rf, npoints_)
    call self%integrate_and_save_profile(r_arr, y0, dr0, epsfrac, epsabs,  &
      drmin, prof_out, status)
    if (status /= status_ok) return

    ! Make points interior to the bubble.
    max_interior_ = npoints_/2
    if (present(max_interior_pts)) max_interior_ = max_interior_pts
    if (max_interior_ > 0) then
      dx0 = r_arr(2) - r_arr(1)
      if (r_arr(1)/dx0 <= real(max_interior_, wp)) then
        n_int = int(ceiling(r_arr(1)/dx0))
        block
          real(wp), allocatable :: r_lin(:)
          r_lin = linspace(0.0_wp, r_arr(1), n_int + 1)
          r_int = r_lin(1:n_int)
        end block
      else
        n_int = max_interior_
        ! r_arr(1) = dx0*(n + a*n*(n+1)/2)
        a = (r_arr(1)/dx0 - real(n_int, wp))  &
          *2.0_wp/(real(n_int, wp)*real(n_int + 1, wp))
        allocate(r_int(n_int))
        do i = 1, n_int
          nk = real(n_int + 1 - i, wp)
          r_int(i) = r_arr(1) - dx0*(nk + 0.5_wp*a*nk*(nk + 1.0_wp))
        end do
        r_int(1) = 0.0_wp  ! enforce this exactly
      end if
      if (size(r_int) > 0) then
        allocate(phi_int(size(r_int)))
        allocate(dphi_int(size(r_int)))
        phi_int(1) = self%phi_absmin + delta_phi0
        dphi_int(1) = 0.0_wp
        dv = self%dv_from_absmin(delta_phi0)
        d2v = self%pot%d2v(phi_int(1))
        do i = 2, size(r_int)
          call self%exact_solution(r_int(i), phi_int(1), dv, d2v,  &
            phi_int(i), dphi_int(i))
        end do
        profile%r = [r_int, prof_out%r]
        profile%phi = [phi_int, prof_out%phi]
        profile%dphi = [dphi_int, prof_out%dphi]
        profile%rerr = prof_out%rerr
      else
        profile = prof_out
      end if
    else
      profile = prof_out
    end if

  end subroutine sfi_find_profile

  function sfi_find_action(self, profile) result(action)
    !! The Euclidean action of the instanton:
    !! S = int [ (dphi/dr)^2/2 + V(phi) - V(phi_metamin) ] r^alpha dr dOmega.
    !! Port of SingleFieldInstanton.findAction.

    class(single_field_instanton), intent(inout) :: self
    type(profile1d), intent(in) :: profile
    real(wp) :: action

    real(wp) :: d
    real(wp) :: area_coef
    real(wp) :: volume
    real(wp) :: v_meta
    real(wp), allocatable :: integrand(:)
    integer :: n
    integer :: i

    n = size(profile%r)
    d = self%alpha + 1.0_wp  ! Number of dimensions in the integration
    area_coef = 2.0_wp*pi**(0.5_wp*d)/gamma(0.5_wp*d)
    v_meta = self%pot%v(self%phi_metamin)
    allocate(integrand(n))
    do i = 1, n
      integrand(i) = (0.5_wp*profile%dphi(i)**2  &
        + self%pot%v(profile%phi(i)) - v_meta)  &
        *profile%r(i)**self%alpha*area_coef
    end do
    action = simpson(integrand, profile%r)
    ! Find the bulk term in the bubble interior.
    volume = profile%r(1)**d*pi**(0.5_wp*d)/gamma(0.5_wp*d + 1.0_wp)
    action = action + volume*(self%pot%v(profile%phi(1)) - v_meta)

  end function sfi_find_action

  subroutine sfi_evenly_spaced_phi(self, phi_in, dphi_in, npoints,  &
      fix_abs, phi_out, dphi_out)
    !! Returns phi and dphi(phi) on a linearly spaced grid in phi (instead
    !! of r). Port of SingleFieldInstanton.evenlySpacedPhi with k = 1
    !! (linear interpolation, the default used by fullTunneling).

    class(single_field_instanton), intent(inout) :: self
    real(wp), intent(in) :: phi_in(:)
    real(wp), intent(in) :: dphi_in(:)
    integer, intent(in) :: npoints
    logical, intent(in) :: fix_abs
      !! If true, make phi go all the way to phi_absmin.
    real(wp), allocatable, intent(out) :: phi_out(:)
    real(wp), allocatable, intent(out) :: dphi_out(:)

    real(wp), allocatable :: phi(:)
    real(wp), allocatable :: dphi(:)
    real(wp), allocatable :: phi_m(:)
    real(wp), allocatable :: dphi_m(:)
    integer, allocatable :: idx(:)
    integer :: i

    if (fix_abs) then
      phi = [self%phi_absmin, phi_in, self%phi_metamin]
      dphi = [0.0_wp, dphi_in, 0.0_wp]
    else
      phi = [phi_in, self%phi_metamin]
      dphi = [dphi_in, 0.0_wp]
    end if
    ! Make sure that phi is increasing everywhere.
    idx = monotonic_indices(phi)
    phi_m = phi(idx)
    dphi_m = dphi(idx)

    if (fix_abs) then
      phi_out = linspace(self%phi_absmin, self%phi_metamin, npoints)
    else
      phi_out = linspace(phi_m(1), self%phi_metamin, npoints)
    end if
    allocate(dphi_out(npoints))
    do i = 1, npoints
      dphi_out(i) = interp_linear(phi_m, dphi_m, phi_out(i))
    end do

  end subroutine sfi_evenly_spaced_phi

  pure function sgn(x) result(s)
    !! Sign function with sgn(0) = 0, like numpy.sign.

    real(wp), intent(in) :: x
    real(wp) :: s

    if (x > 0.0_wp) then
      s = 1.0_wp
    else if (x < 0.0_wp) then
      s = -1.0_wp
    else
      s = 0.0_wp
    end if

  end function sgn

end module cosmotransitions__tunneling1d
