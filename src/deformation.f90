module cosmotransitions__deformation
!! Port of pathDeformation.Deformation_Spline: deform a path in the
!! presence of a potential such that the normal forces along the path
!! vanish. The path is fit to a set of spline basis functions and the
!! deformation acts on the spline coefficients.

  use cosmotransitions__config, only : wp
  use cosmotransitions__config, only : status_ok
  use cosmotransitions__config, only : err_deformation
  use cosmotransitions__helpers, only : nbspld2
  use cosmotransitions__helpers, only : linspace
  use cosmotransitions__helpers, only : lstsq
  use cosmotransitions__potentials, only : potential_nd

  implicit none

  private

  public :: deformation_spline

  type :: deformation_spline
    !! Port of the class Deformation_Spline. The `save_all_steps` option
    !! of the Python class is not implemented.
    real(wp), allocatable :: phi(:, :)
      !! Current path points, shape (n, ndim). Rewritten at each step.
    real(wp), allocatable :: v2(:)
      !! Squared 'speed' along the path at the initial points (does not
      !! change as the path deforms).
    integer :: num_steps = 0
      !! Total number of steps taken.
    logical :: fix_start = .false.
    logical :: fix_end = .false.
    class(potential_nd), pointer :: pot => null()
    ! Private-ish members (leading underscore attributes in Python):
    real(wp) :: length = 0.0_wp
      !! Total length of the path, set during initialization (`_L`).
    real(wp), allocatable :: t(:)
      !! Path parameter in (0, 1] marking the location of each point.
    real(wp), allocatable :: xb(:, :)
    real(wp), allocatable :: dxb(:, :)
    real(wp), allocatable :: d2xb(:, :)
      !! Spline basis functions and their derivatives evaluated at `t`.
    real(wp), allocatable :: beta(:, :)
      !! Spline coefficients for each dimension, shape (nb, ndim).
    real(wp), allocatable :: phi_prev(:, :)
    real(wp), allocatable :: f_prev(:, :)
    logical :: has_prev = .false.
  contains
    procedure :: init => dfs_init
    procedure :: forces => dfs_forces
    procedure :: step => dfs_step
    procedure :: deform_path => dfs_deform_path
  end type deformation_spline

contains

  subroutine dfs_init(self, phi, dphidr, pot, status, nb, kb, v2min,  &
      fix_start, fix_end)
    !! Port of Deformation_Spline.__init__.

    class(deformation_spline), intent(inout) :: self
    real(wp), intent(in) :: phi(:, :)
      !! Initial path, shape (n_points, ndim).
    real(wp), intent(in) :: dphidr(:)
      !! The 'speed' along the path at the initial points, shape
      !! (n_points,).
    class(potential_nd), intent(inout), target :: pot
      !! The potential; only its gradient is used. The actual argument
      !! must have the `target` (or `pointer`) attribute and outlive this
      !! object.
    integer, intent(out) :: status
    integer, intent(in), optional :: nb
      !! Number of basis splines (default 10).
    integer, intent(in), optional :: kb
      !! Order of the basis splines (default 3).
    real(wp), intent(in), optional :: v2min
      !! Smallest allowed square of dphidr, relative to the characteristic
      !! force exerted by the potential (default 0).
    logical, intent(in), optional :: fix_start
    logical, intent(in), optional :: fix_end
      !! If true, the force on the first/last point is set to zero, so the
      !! point will not move (default false).

    integer :: nb_
    integer :: kb_
    real(wp) :: v2min_
    integer :: n
    integer :: i
    real(wp), allocatable :: dl(:)
    real(wp), allocatable :: knots(:)
    real(wp), allocatable :: phi_lin(:, :)
    real(wp) :: dvmax
    integer :: st

    status = status_ok
    nb_ = 10
    kb_ = 3
    v2min_ = 0.0_wp
    if (present(nb)) nb_ = nb
    if (present(kb)) kb_ = kb
    if (present(v2min)) v2min_ = v2min
    self%fix_start = .false.
    self%fix_end = .false.
    if (present(fix_start)) self%fix_start = fix_start
    if (present(fix_end)) self%fix_end = fix_end
    self%pot => pot
    self%num_steps = 0
    self%has_prev = .false.
    if (allocated(self%phi_prev)) deallocate(self%phi_prev)
    if (allocated(self%f_prev)) deallocate(self%f_prev)

    ! First step: convert phi to a set of path lengths.
    n = size(phi, 1)
    self%phi = phi
    allocate(dl(n - 1))
    do i = 1, n - 1
      dl(i) = norm2(phi(i+1, :) - phi(i, :))
    end do
    if (allocated(self%t)) deallocate(self%t)
    allocate(self%t(n))
    self%t(1) = 0.0_wp
    do i = 2, n
      self%t(i) = self%t(i-1) + dl(i-1)
    end do
    self%length = self%t(n)
    self%t = self%t/self%length
    self%t(1) = 1.0e-100_wp  ! Without this, the first data point isn't in
                             ! any bin (this matters for dxb).

    ! Create the starting spline: make the knots and then the spline
    ! matrices at each point t.
    allocate(knots(2*(kb_ - 1) + nb_ + 3 - kb_))
    knots(1:kb_-1) = 0.0_wp
    knots(kb_:kb_+nb_+2-kb_) = linspace(0.0_wp, 1.0_wp, nb_ + 3 - kb_)
    knots(kb_+nb_+3-kb_:) = 1.0_wp
    call nbspld2(knots, self%t, kb_, self%xb, self%dxb, self%d2xb)

    ! Subtract off the linear component and fit the spline coefficients.
    allocate(phi_lin(n, size(phi, 2)))
    do i = 1, n
      phi_lin(i, :) = phi(1, :) + (phi(n, :) - phi(1, :))*self%t(i)
    end do
    call lstsq(self%xb, phi - phi_lin, self%beta, st)
    if (st /= status_ok) then
      status = st
      return
    end if

    ! Ensure that v2 isn't too small.
    self%v2 = dphidr**2
    if (v2min_ > 0.0_wp) then
      block
        ! Contiguous copy of the point: passing a strided section directly
        ! to the deferred type-bound grad() can miscompile with gfortran.
        real(wp), allocatable :: xpt(:)
        allocate(xpt(size(phi, 2)))
        dvmax = 0.0_wp
        do i = 1, n
          xpt = phi(i, :)
          dvmax = max(dvmax, norm2(pot%grad(xpt)))
        end do
      end block
      v2min_ = v2min_*dvmax*self%length/real(nb_, wp)
      do i = 1, n
        if (self%v2(i) < v2min_) self%v2(i) = v2min_
      end do
    end if

  end subroutine dfs_init

  subroutine dfs_forces(self, f_norm, dv)
    !! Calculate the normal force and the potential gradient on the path.
    !! Port of Deformation_Spline.forces.

    class(deformation_spline), intent(inout) :: self
    real(wp), allocatable, intent(out) :: f_norm(:, :)
    real(wp), allocatable, intent(out) :: dv(:, :)

    integer :: n
    integer :: ndim
    integer :: i
    real(wp), allocatable :: dphi(:, :)
    real(wp), allocatable :: d2phi(:, :)
    real(wp), allocatable :: dphids(:, :)
    real(wp), allocatable :: d2phids2(:, :)
    real(wp), allocatable :: xpt(:)
    real(wp) :: dphi_sq
    real(wp) :: proj

    n = size(self%phi, 1)
    ndim = size(self%phi, 2)
    allocate(xpt(ndim))

    ! First find dphi and d2phi with respect to the path parameter t. Note
    ! that dphi needs a linear component added in, while d2phi does not.
    ! Nota bene: the Python code adds (phi[-1] - phi[1]), i.e. it uses the
    ! *second* point instead of the first one. We reproduce that behaviour
    ! faithfully; the converged path satisfies F_norm = 0 either way.
    dphi = matmul(self%dxb, self%beta)
    d2phi = matmul(self%d2xb, self%beta)
    do i = 1, n
      dphi(i, :) = dphi(i, :) + (self%phi(n, :) - self%phi(2, :))
    end do

    ! Compute dphi/ds, where s is the path length instead of the path
    ! parameter t. This is just the direction along the path. Then find
    ! the acceleration along the path, d2phi/ds2.
    allocate(dphids(n, ndim))
    allocate(d2phids2(n, ndim))
    allocate(f_norm(n, ndim))
    allocate(dv(n, ndim))
    do i = 1, n
      dphi_sq = sum(dphi(i, :)*dphi(i, :))
      dphids(i, :) = dphi(i, :)/sqrt(dphi_sq)
      d2phids2(i, :) = (d2phi(i, :)  &
        - dphi(i, :)*sum(dphi(i, :)*d2phi(i, :))/dphi_sq)/dphi_sq
      ! Now get the normal force acting on the path.
      xpt = self%phi(i, :)
      dv(i, :) = self%pot%grad(xpt)
      proj = sum(dv(i, :)*dphids(i, :))
      f_norm(i, :) = d2phids2(i, :)*self%v2(i)  &
        - (dv(i, :) - proj*dphids(i, :))
    end do
    if (self%fix_start) f_norm(1, :) = 0.0_wp
    if (self%fix_end) f_norm(n, :) = 0.0_wp

  end subroutine dfs_forces

  subroutine dfs_step(self, last_step, stepsize, step_reversed, fratio,  &
      status, maxstep, minstep, reverse_check, step_increase,  &
      step_decrease, check_after_fit, verbose)
    !! Deform the path one step: push each point in the direction of the
    !! normal force. Port of Deformation_Spline.step.

    class(deformation_spline), intent(inout) :: self
    real(wp), intent(in) :: last_step
      !! Size of the last step.
    real(wp), intent(out) :: stepsize
      !! The stepsize used for this step.
    logical, intent(out) :: step_reversed
      !! True if this step was reversed.
    real(wp), intent(out) :: fratio
      !! Ratio of the maximum normal force to the maximum potential
      !! gradient; goes to zero when the path is a perfect fit.
    integer, intent(out) :: status
    real(wp), intent(in), optional :: maxstep
    real(wp), intent(in), optional :: minstep
    real(wp), intent(in), optional :: reverse_check
      !! Fraction of points for which the force can reverse direction
      !! (relative to the last step) before the stepsize is decreased.
    real(wp), intent(in), optional :: step_increase
    real(wp), intent(in), optional :: step_decrease
    logical, intent(in), optional :: check_after_fit
      !! If true, the convergence test is performed after the points are
      !! fit to the spline.
    logical, intent(in), optional :: verbose

    real(wp) :: maxstep_
    real(wp) :: minstep_
    real(wp) :: reverse_check_
    real(wp) :: step_increase_
    real(wp) :: step_decrease_
    logical :: check_after_fit_
    logical :: verbose_
    integer :: n
    integer :: i
    integer :: n_reversed
    real(wp) :: f_max
    real(wp) :: dv_max
    real(wp) :: fratio1
    real(wp) :: fratio2
    real(wp), allocatable :: f(:, :)
    real(wp), allocatable :: dv(:, :)
    real(wp), allocatable :: phi(:, :)
    real(wp), allocatable :: phi_lin(:, :)
    real(wp), allocatable :: ffit(:, :)
    integer :: st

    maxstep_ = 0.1_wp
    minstep_ = 1.0e-4_wp
    reverse_check_ = 0.15_wp
    step_increase_ = 1.5_wp
    step_decrease_ = 5.0_wp
    check_after_fit_ = .true.
    verbose_ = .false.
    if (present(maxstep)) maxstep_ = maxstep
    if (present(minstep)) minstep_ = minstep
    if (present(reverse_check)) reverse_check_ = reverse_check
    if (present(step_increase)) step_increase_ = step_increase
    if (present(step_decrease)) step_decrease_ = step_decrease
    if (present(check_after_fit)) check_after_fit_ = check_after_fit
    if (present(verbose)) verbose_ = verbose

    status = status_ok
    n = size(self%phi, 1)

    ! Find out the direction of the deformation.
    call self%forces(f, dv)
    f_max = 0.0_wp
    dv_max = 0.0_wp
    do i = 1, n
      f_max = max(f_max, norm2(f(i, :)))
      dv_max = max(dv_max, norm2(dv(i, :)))
    end do
    fratio1 = f_max/dv_max
    ! Rescale the normal force so that it's relative to the path length.
    f = f*self%length/dv_max

    ! Now see how big the stepsize should be.
    stepsize = last_step
    phi = self%phi
    step_reversed = .false.
    if (reverse_check_ < 1.0_wp .and. self%has_prev) then
      n_reversed = 0
      do i = 1, n
        if (sum(f(i, :)*self%f_prev(i, :)) < 0.0_wp) then
          n_reversed = n_reversed + 1
        end if
      end do
      if (real(n_reversed, wp) > real(n, wp)*reverse_check_) then
        ! We want to reverse the last step.
        if (stepsize > minstep_) then
          step_reversed = .true.
          phi = self%phi_prev
          f = self%f_prev
          if (verbose_) print "(a)", "step reversed"
          stepsize = last_step/step_decrease_
        end if
      else
        ! No (large number of) indices reversed, just do a regular step.
        ! Increase the stepsize a bit over the last one.
        stepsize = last_step*step_increase_
      end if
    end if
    if (stepsize > maxstep_) stepsize = maxstep_
    if (stepsize < minstep_) stepsize = minstep_

    ! Save the state before the step.
    self%phi_prev = phi
    self%f_prev = f
    self%has_prev = .true.

    ! Now make the step.
    phi = phi + f*stepsize

    ! Fit to the spline.
    allocate(phi_lin(n, size(phi, 2)))
    do i = 1, n
      phi_lin(i, :) = phi(1, :) + (phi(n, :) - phi(1, :))*self%t(i)
    end do
    phi = phi - phi_lin
    call lstsq(self%xb, phi, self%beta, st)
    if (st /= status_ok) then
      status = st
      return
    end if
    phi = matmul(self%xb, self%beta) + phi_lin
    self%phi = phi

    ffit = (phi - self%phi_prev)/stepsize
    fratio2 = 0.0_wp
    do i = 1, n
      fratio2 = max(fratio2, norm2(ffit(i, :)))
    end do
    fratio2 = fratio2/self%length

    if (verbose_) then
      print "(a,i0,a,es9.2,a,es9.2,a,es9.2)", "step: ", self%num_steps,  &
        "; stepsize: ", stepsize, "; fRatio1: ", fratio1, "; fRatio2: ",  &
        fratio2
    end if

    if (check_after_fit_) then
      fratio = fratio2
    else
      fratio = fratio1
    end if

  end subroutine dfs_step

  subroutine dfs_deform_path(self, converged, status, startstep,  &
      fratio_conv, converge_0, fratio_increase, maxiter, verbose)
    !! Deform the path in many individual steps, stopping when the
    !! convergence criterion is reached, when the maximum number of
    !! iterations is reached, or when the path appears to be running away
    !! from convergence. Port of Deformation_Spline.deformPath.
    !!
    !! In the runaway case the status is set to `err_deformation` (the
    !! DeformationError of the Python code) and the path is reset to the
    !! point of best convergence.

    class(deformation_spline), intent(inout) :: self
    logical, intent(out) :: converged
      !! True if the deformation converged (as determined by fratio_conv).
    integer, intent(out) :: status
    real(wp), intent(in), optional :: startstep
      !! Starting stepsize (default 2e-3).
    real(wp), intent(in), optional :: fratio_conv
      !! Convergence criterion (default 0.02).
    real(wp), intent(in), optional :: converge_0
      !! On the first step, converge if fratio < converge_0*fratio_conv
      !! (default 5).
    real(wp), intent(in), optional :: fratio_increase
      !! Maximum fractional amount that fratio can increase before
      !! aborting (default 5).
    integer, intent(in), optional :: maxiter
      !! Maximum number of steps (default 500).
    logical, intent(in), optional :: verbose
      !! If true, print the ending condition (default true).

    real(wp) :: startstep_
    real(wp) :: fratio_conv_
    real(wp) :: converge_0_
    real(wp) :: fratio_increase_
    integer :: maxiter_
    logical :: verbose_
    real(wp) :: minfratio
    real(wp) :: stepsize
    real(wp) :: stepsize_new
    real(wp) :: fratio
    logical :: step_reversed
    real(wp), allocatable :: minfratio_beta(:, :)
    real(wp), allocatable :: minfratio_phi(:, :)
    integer :: st

    startstep_ = 2.0e-3_wp
    fratio_conv_ = 0.02_wp
    converge_0_ = 5.0_wp
    fratio_increase_ = 5.0_wp
    maxiter_ = 500
    verbose_ = .true.
    if (present(startstep)) startstep_ = startstep
    if (present(fratio_conv)) fratio_conv_ = fratio_conv
    if (present(converge_0)) converge_0_ = converge_0
    if (present(fratio_increase)) fratio_increase_ = fratio_increase
    if (present(maxiter)) maxiter_ = maxiter
    if (present(verbose)) verbose_ = verbose

    status = status_ok
    converged = .false.
    minfratio = huge(1.0_wp)
    stepsize = startstep_

    do
      self%num_steps = self%num_steps + 1
      call self%step(stepsize, stepsize_new, step_reversed, fratio, st)
      stepsize = stepsize_new
      if (st /= status_ok) then
        status = st
        return
      end if
      minfratio = min(minfratio, fratio)
      if (fratio < fratio_conv_ .or.  &
          (self%num_steps == 1 .and. fratio < converge_0_*fratio_conv_)) then
        if (verbose_) then
          print "(a,i0,a,es12.5)", "Path deformation converged. ",  &
            self%num_steps, " steps. fRatio = ", fratio
        end if
        converged = .true.
        exit
      end if
      if (minfratio == fratio) then
        minfratio_beta = self%beta
        minfratio_phi = self%phi
      end if
      if (fratio > fratio_increase_*minfratio .and. .not. step_reversed) then
        self%beta = minfratio_beta
        self%phi = minfratio_phi
        if (verbose_) then
          print "(a)", "Deformation doesn't appear to be converging. "//  &
            "Stopping at the point of best convergence."
        end if
        status = err_deformation
        return
      end if
      if (self%num_steps >= maxiter_) then
        if (verbose_) then
          print "(a)", "Maximum number of deformation iterations reached."
        end if
        exit
      end if
    end do

  end subroutine dfs_deform_path

end module cosmotransitions__deformation
