module cosmotransitions__splinepath
!! Port of pathDeformation.SplinePath: fit a spline to a path in
!! multi-dimensional field space and evaluate the potential (and its
!! derivatives) as a function of the distance along the path.
!!
!! The parametric spline (scipy.interpolate.splprep) is replaced by one
!! bspline-fortran interpolant per field dimension, and the potential
!! spline (scipy.interpolate.splrep) by another bspline-fortran
!! interpolant. Only the `V_spline_samples > 0` mode of the Python class
!! is implemented, which is what fullTunneling uses by default.

  use bspline_module, only : bspline_1d
  use, intrinsic :: iso_fortran_env, only : ip => int32

  use cosmotransitions__config, only : wp
  use cosmotransitions__config, only : status_ok
  use cosmotransitions__config, only : err_numerical
  use cosmotransitions__helpers, only : deriv14_const_dx
  use cosmotransitions__helpers, only : cumtrapz
  use cosmotransitions__helpers, only : linspace
  use cosmotransitions__helpers, only : interp_linear
  use cosmotransitions__optimize, only : minimize_unbounded
  use cosmotransitions__potentials, only : potential_1d
  use cosmotransitions__potentials, only : potential_nd

  implicit none

  private

  public :: spline_path
  public :: path_deriv

  ! 10-point Gauss-Legendre nodes and weights on [-1, 1], used to
  ! re-evaluate the distance to each knot by integrating |dp/dx| along the
  ! spline (replacing the scipy.integrate.odeint call in the Python class).
  real(wp), parameter :: gl_x(10) = [  &
    -0.9739065285171717_wp, -0.8650633666889845_wp,  &
    -0.6794095682990244_wp, -0.4333953941292472_wp,  &
    -0.1488743389816312_wp, 0.1488743389816312_wp,  &
    0.4333953941292472_wp, 0.6794095682990244_wp,  &
    0.8650633666889845_wp, 0.9739065285171717_wp]
  real(wp), parameter :: gl_w(10) = [  &
    0.0666713443086881_wp, 0.1494513491505806_wp,  &
    0.2190863625159820_wp, 0.2692667193099963_wp,  &
    0.2955242247147529_wp, 0.2955242247147529_wp,  &
    0.2692667193099963_wp, 0.2190863625159820_wp,  &
    0.1494513491505806_wp, 0.0666713443086881_wp]

  type, extends(potential_1d) :: spline_path
    !! A path in field space, parametrized by the distance x along it.
    !! As a `potential_1d`, it provides V(x), dV/dx and d2V/dx2 on the
    !! path, which is exactly what single_field_instanton needs.
    real(wp) :: length = 0.0_wp
      !! Total length of the path (attribute `L` in Python).
    integer :: ndim = 0
    type(bspline_1d), allocatable :: path_spl(:)
      !! One interpolant per field dimension, knotted by path distance.
    type(bspline_1d) :: v_spl
      !! Interpolant of the potential along the path.
  contains
    procedure :: init => spline_path_init
    procedure :: v => spline_path_v
    procedure :: dv => spline_path_dv
    procedure :: d2v => spline_path_d2v
    procedure :: pts => spline_path_pts
    procedure :: pts_many => spline_path_pts_many
  end type spline_path

contains

  subroutine spline_path_init(self, pts_in, pot, status, v_spline_samples,  &
      extend_to_minima, reeval_distances)
    !! Port of SplinePath.__init__.

    class(spline_path), intent(inout) :: self
    real(wp), intent(in) :: pts_in(:, :)
      !! The points that describe the path, shape (num_points, ndim).
    class(potential_nd), intent(inout) :: pot
      !! The full multi-field potential.
    integer, intent(out) :: status
    integer, intent(in), optional :: v_spline_samples
      !! Number of samples along the path used to build the potential
      !! spline (default 100).
    logical, intent(in), optional :: extend_to_minima
      !! If true, extend the path at both ends until it hits local minima
      !! (default false).
    logical, intent(in), optional :: reeval_distances
      !! If true, get more accurate knot distances by integrating along
      !! the spline (default true).

    integer :: nsamples
    logical :: extend_
    logical :: reeval_
    real(wp), allocatable :: pts(:, :)
    real(wp), allocatable :: dpts(:, :)
    real(wp), allocatable :: pdist(:)
    real(wp), allocatable :: speeds(:)
    real(wp), allocatable :: newdist(:)
    real(wp), allocatable :: xs(:)
    real(wp), allocatable :: x_ext(:)
    real(wp), allocatable :: xv(:)
    real(wp), allocatable :: yv(:)
    real(wp), allocatable :: p0(:)
    real(wp), allocatable :: dp0(:)
    real(wp), allocatable :: xg(:)
    real(wp), allocatable :: pt_ext(:, :)
    real(wp), allocatable :: tmp(:, :)
    real(wp) :: xmin
    real(wp) :: dx1
    real(wp) :: xm
    real(wp) :: hh
    real(wp) :: seg
    integer :: n
    integer :: nx
    integer :: nxg
    integer :: n_ext
    integer :: d
    integer :: i
    integer :: j
    integer :: k
    integer(ip) :: iflag
    integer :: st

    status = status_ok
    nsamples = 100
    if (present(v_spline_samples)) nsamples = v_spline_samples
    extend_ = .false.
    if (present(extend_to_minima)) extend_ = extend_to_minima
    reeval_ = .true.
    if (present(reeval_distances)) reeval_ = reeval_distances

    pts = pts_in
    self%ndim = size(pts, 2)

    ! 1. Find derivs
    dpts = path_deriv(pts)

    ! 2. Extend the path
    if (extend_) then
      ! Extend at the front of the path.
      n = size(pts, 1)
      p0 = pts(1, :)
      dp0 = dpts(1, :)
      xmin = minimize_unbounded(v_lin, 0.0_wp, 1.0e-6_wp, st)
      if (st /= status_ok) xmin = 0.0_wp
      if (xmin > 0.0_wp) xmin = 0.0_wp
      nxg = int(ceiling(abs(xmin) - 0.5_wp)) + 1
      xg = linspace(xmin, 0.0_wp, nxg)
      allocate(pt_ext(nxg, self%ndim))
      do i = 1, nxg
        pt_ext(i, :) = p0 + xg(i)*dp0
      end do
      allocate(tmp(nxg + n - 1, self%ndim))
      tmp(1:nxg, :) = pt_ext
      tmp(nxg+1:, :) = pts(2:n, :)
      call move_alloc(tmp, pts)
      deallocate(pt_ext)
      deallocate(xg)

      ! Extend at the end of the path. Like the Python code, this uses the
      ! derivative of the *original* path at its last point (the last point
      ! itself is unchanged by the front extension).
      n = size(pts, 1)
      p0 = pts(n, :)
      dp0 = dpts(size(dpts, 1), :)
      xmin = minimize_unbounded(v_lin, 0.0_wp, 1.0e-6_wp, st)
      if (st /= status_ok) xmin = 0.0_wp
      if (xmin < 0.0_wp) xmin = 0.0_wp
      nxg = int(ceiling(abs(xmin) - 0.5_wp)) + 1
      xg = linspace(xmin, 0.0_wp, nxg)
      allocate(pt_ext(nxg, self%ndim))
      do i = 1, nxg
        ! Reversed: the extension starts at the old end point (x = 0) and
        ! goes outward to xmin.
        pt_ext(i, :) = p0 + xg(nxg + 1 - i)*dp0
      end do
      allocate(tmp(n - 1 + nxg, self%ndim))
      tmp(1:n-1, :) = pts(1:n-1, :)
      tmp(n:, :) = pt_ext
      call move_alloc(tmp, pts)
      deallocate(pt_ext)
      deallocate(xg)

      ! Recalculate the derivative.
      dpts = path_deriv(pts)
    end if

    ! The bspline interpolation below needs at least order = k+1 = 4
    ! points for a cubic. Densify very short paths by linear interpolation
    ! along the polyline; for the typical straight two-point starting
    ! guess this represents the same path exactly.
    n = size(pts, 1)
    if (n < 5) then
      call densify_path(pts, 5)
      dpts = path_deriv(pts)
      n = size(pts, 1)
    end if

    ! 3. Find knot positions and fit the spline.
    allocate(speeds(n))
    do i = 1, n
      speeds(i) = norm2(dpts(i, :))
    end do
    pdist = cumtrapz(speeds)
    self%length = pdist(n)
    k = min(n - 1, 3)  ! degree of the spline

    if (allocated(self%path_spl)) deallocate(self%path_spl)
    allocate(self%path_spl(self%ndim))
    do d = 1, self%ndim
      call self%path_spl(d)%initialize(pdist, pts(:, d), int(k + 1, ip),  &
        iflag, extrap=.true.)
      if (iflag /= 0_ip) then
        status = err_numerical
        return
      end if
    end do

    ! 4. Re-evaluate the distance to each point.
    if (reeval_) then
      allocate(newdist(n))
      newdist(1) = pdist(1)
      do i = 2, n
        ! Arc length over [pdist(i-1), pdist(i)] by Gauss-Legendre
        ! quadrature of |dp/dx|.
        xm = 0.5_wp*(pdist(i) + pdist(i-1))
        hh = 0.5_wp*(pdist(i) - pdist(i-1))
        seg = 0.0_wp
        do j = 1, size(gl_x)
          seg = seg + gl_w(j)*path_speed(self, xm + hh*gl_x(j))
        end do
        newdist(i) = newdist(i-1) + hh*seg
      end do
      pdist = newdist
      self%length = pdist(n)
      do d = 1, self%ndim
        call self%path_spl(d)%initialize(pdist, pts(:, d), int(k + 1, ip),  &
          iflag, extrap=.true.)
        if (iflag /= 0_ip) then
          status = err_numerical
          return
        end if
      end do
    end if

    ! Now make the potential spline.
    if (nsamples <= 0) then
      error stop "spline_path_init: v_spline_samples must be positive "//  &
        "(the V_spline_samples=None mode of the Python class is not ported)"
    end if
    xs = linspace(0.0_wp, self%length, nsamples)
    ! Extend 20% beyond the path so that we more accurately model the
    ! path end points.
    dx1 = xs(2)
    n_ext = 0
    do while (dx1*real(n_ext + 1, wp) < self%length*0.2_wp)
      n_ext = n_ext + 1
    end do
    allocate(x_ext(n_ext))
    do i = 1, n_ext
      x_ext(i) = dx1*real(i, wp)
    end do
    nx = nsamples + 2*n_ext
    allocate(xv(nx))
    xv(1:n_ext) = -x_ext(n_ext:1:-1)
    xv(n_ext+1:n_ext+nsamples) = xs
    xv(n_ext+nsamples+1:nx) = self%length + x_ext
    allocate(yv(nx))
    do i = 1, nx
      yv(i) = pot%v(self%pts(xv(i)))
    end do
    call self%v_spl%initialize(xv, yv, 4_ip, iflag, extrap=.true.)
    if (iflag /= 0_ip) then
      status = err_numerical
      return
    end if

  contains

    function v_lin(x_) result(y_)
      !! The potential along the linear extension of the path,
      !! V(p0 + x*dp0), used to find the nearest minimum beyond an end
      !! point.
      real(wp), intent(in) :: x_
      real(wp) :: y_
      y_ = pot%v(p0 + x_*dp0)
    end function v_lin

  end subroutine spline_path_init

  subroutine densify_path(pts, ntarget)
    !! Resample the polyline `pts` at `ntarget` parameter values uniformly
    !! spaced in cumulative chord length, using linear interpolation.

    real(wp), allocatable, intent(inout) :: pts(:, :)
    integer, intent(in) :: ntarget

    real(wp), allocatable :: s(:)
    real(wp), allocatable :: snew(:)
    real(wp), allocatable :: out(:, :)
    integer :: n
    integer :: ndim
    integer :: i
    integer :: d

    n = size(pts, 1)
    ndim = size(pts, 2)
    if (n >= ntarget) return
    allocate(s(n))
    s(1) = 0.0_wp
    do i = 2, n
      s(i) = s(i-1) + norm2(pts(i, :) - pts(i-1, :))
    end do
    snew = linspace(0.0_wp, s(n), ntarget)
    allocate(out(ntarget, ndim))
    do i = 1, ntarget
      do d = 1, ndim
        out(i, d) = interp_linear(s, pts(:, d), snew(i))
      end do
    end do
    call move_alloc(out, pts)

  end subroutine densify_path

  function path_deriv(phi) result(dphi)
    !! Derivative of the path points with respect to the point index:
    !! 4th order if there are at least 5 points, otherwise 1st/2nd order.
    !! Port of pathDeformation._pathDeriv.

    real(wp), intent(in) :: phi(:, :)
    real(wp), allocatable :: dphi(:, :)

    integer :: n

    n = size(phi, 1)
    allocate(dphi(n, size(phi, 2)))
    if (n >= 5) then
      dphi = deriv14_const_dx(phi)
    else if (n > 2) then
      dphi(2:n-1, :) = 0.5_wp*(phi(3:n, :) - phi(1:n-2, :))
      dphi(1, :) = -1.5_wp*phi(1, :) + 2.0_wp*phi(2, :) - 0.5_wp*phi(3, :)
      dphi(n, :) = 1.5_wp*phi(n, :) - 2.0_wp*phi(n-1, :) + 0.5_wp*phi(n-2, :)
    else
      dphi(1, :) = phi(2, :) - phi(1, :)
      dphi(2, :) = phi(2, :) - phi(1, :)
    end if

  end function path_deriv

  function path_speed(self, x) result(speed)
    !! |dp/dx| along the path spline.

    class(spline_path), intent(inout) :: self
    real(wp), intent(in) :: x
    real(wp) :: speed

    real(wp) :: dp
    integer :: d
    integer(ip) :: iflag

    speed = 0.0_wp
    do d = 1, self%ndim
      call self%path_spl(d)%evaluate(x, 1_ip, dp, iflag)
      speed = speed + dp*dp
    end do
    speed = sqrt(speed)

  end function path_speed

  function spline_path_v(self, phi) result(y)
    !! The potential as a function of the distance `phi` along the path.

    class(spline_path), intent(inout) :: self
    real(wp), intent(in) :: phi
    real(wp) :: y

    integer(ip) :: iflag

    call self%v_spl%evaluate(phi, 0_ip, y, iflag)

  end function spline_path_v

  function spline_path_dv(self, phi) result(y)
    !! dV/dx as a function of the distance along the path.

    class(spline_path), intent(inout) :: self
    real(wp), intent(in) :: phi
    real(wp) :: y

    integer(ip) :: iflag

    call self%v_spl%evaluate(phi, 1_ip, y, iflag)

  end function spline_path_dv

  function spline_path_d2v(self, phi) result(y)
    !! d2V/dx2 as a function of the distance along the path.

    class(spline_path), intent(inout) :: self
    real(wp), intent(in) :: phi
    real(wp) :: y

    integer(ip) :: iflag

    call self%v_spl%evaluate(phi, 2_ip, y, iflag)

  end function spline_path_d2v

  function spline_path_pts(self, x) result(p)
    !! The field-space point at distance `x` along the path.

    class(spline_path), intent(inout) :: self
    real(wp), intent(in) :: x
    real(wp) :: p(self%ndim)

    integer :: d
    integer(ip) :: iflag

    do d = 1, self%ndim
      call self%path_spl(d)%evaluate(x, 0_ip, p(d), iflag)
    end do

  end function spline_path_pts

  function spline_path_pts_many(self, x) result(p)
    !! The field-space points at the distances `x` along the path,
    !! shape (size(x), ndim).

    class(spline_path), intent(inout) :: self
    real(wp), intent(in) :: x(:)
    real(wp) :: p(size(x), self%ndim)

    integer :: i

    do i = 1, size(x)
      p(i, :) = self%pts(x(i))
    end do

  end function spline_path_pts_many

end module cosmotransitions__splinepath
