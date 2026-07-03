module cosmotransitions__helpers
!! Port of the parts of cosmoTransitions/helper_functions.py that are needed
!! by pathDeformation.fullTunneling, plus a few small numerical utilities
!! (trapezoid/Simpson integration, linear interpolation, least squares).

  use cosmotransitions__config, only : wp
  use cosmotransitions__config, only : status_ok
  use cosmotransitions__config, only : err_integration
  use cosmotransitions__config, only : err_numerical
  use gradmin__linearsolve, only : gauss_jordan

  implicit none

  private

  public :: dydt_iface
  public :: rkqs
  public :: rkck
  public :: cubic_interp
  public :: deriv14_const_dx
  public :: nbspld2
  public :: monotonic_indices
  public :: cumtrapz
  public :: simpson
  public :: linspace
  public :: interp_linear
  public :: lstsq

  abstract interface
    subroutine dydt_iface(y, t, dydt)
      !! Derivative function for the ODE integrators.
      import :: wp
      implicit none
      real(wp), intent(in) :: y(:)
      real(wp), intent(in) :: t
      real(wp), intent(out) :: dydt(:)
    end subroutine dydt_iface
  end interface

contains

  subroutine rkqs(y, dydt, t, f, dt_try, epsfrac, epsabs, dy_out, dt_out,  &
      dtnext, status)
    !! Take a single 5th order Runge-Kutta step with error monitoring.
    !! Port of helper_functions.rkqs (adapted from Numerical Recipes).
    !! The step size dynamically changes such that the error in `y` is
    !! smaller than the larger of `epsfrac` and `epsabs`.

    real(wp), intent(in) :: y(:)
    real(wp), intent(in) :: dydt(:)
      !! Derivative at the start of the step; must satisfy dydt = f(y, t).
    real(wp), intent(in) :: t
    procedure(dydt_iface) :: f
    real(wp), intent(in) :: dt_try
      !! Initial guess for the step size.
    real(wp), intent(in) :: epsfrac(:)
    real(wp), intent(in) :: epsabs(:)
      !! Maximum fractional and absolute errors, same size as `y`.
    real(wp), intent(out) :: dy_out(:)
      !! Change in `y` during this step.
    real(wp), intent(out) :: dt_out
      !! Change in `t` during this step.
    real(wp), intent(out) :: dtnext
      !! Best guess for the next step size.
    integer, intent(out) :: status
      !! `err_integration` if the step size rounds down to zero.

    real(wp) :: dt
    real(wp) :: errmax
    real(wp) :: dttemp
    real(wp) :: e1
    real(wp) :: e2
    real(wp), dimension(size(y)) :: dy
    real(wp), dimension(size(y)) :: yerr
    integer :: i

    status = status_ok
    dt = dt_try
    do
      call rkck(y, dydt, t, f, dt, dy, yerr)
      errmax = 0.0_wp
      do i = 1, size(y)
        e1 = abs(yerr(i)/epsabs(i))
        e2 = abs(yerr(i))/((abs(y(i)) + 1.0e-300_wp)*epsfrac(i))
        errmax = max(errmax, min(e1, e2))
      end do
      if (errmax < 1.0_wp) exit  ! Step succeeded
      dttemp = 0.9_wp*dt*errmax**(-0.25_wp)
      if (dt > 0.0_wp) then
        dt = max(dttemp, 0.1_wp*dt)
      else
        dt = min(dttemp, 0.1_wp*dt)
      end if
      if (t + dt == t) then
        ! "Stepsize rounds down to zero."
        status = err_integration
        return
      end if
    end do
    if (errmax > 1.89e-4_wp) then
      dtnext = 0.9_wp*dt*errmax**(-0.2_wp)
    else
      dtnext = 5.0_wp*dt
    end if
    dy_out = dy
    dt_out = dt

  end subroutine rkqs

  subroutine rkck(y, dydt, t, f, dt, dyout, yerr)
    !! Take one 5th-order Cash-Karp Runge-Kutta step.
    !! Port of helper_functions._rkck.

    real(wp), intent(in) :: y(:)
    real(wp), intent(in) :: dydt(:)
    real(wp), intent(in) :: t
    procedure(dydt_iface) :: f
    real(wp), intent(in) :: dt
    real(wp), intent(out) :: dyout(:)
      !! The change in `y` during this step.
    real(wp), intent(out) :: yerr(:)
      !! An error estimate for `y`.

    real(wp), parameter :: a2 = 0.2_wp, a3 = 0.3_wp, a4 = 0.6_wp,  &
      a5 = 1.0_wp, a6 = 0.875_wp
    real(wp), parameter :: b21 = 0.2_wp
    real(wp), parameter :: b31 = 3.0_wp/40.0_wp, b32 = 9.0_wp/40.0_wp
    real(wp), parameter :: b41 = 0.3_wp, b42 = -0.9_wp, b43 = 1.2_wp
    real(wp), parameter :: b51 = -11.0_wp/54.0_wp, b52 = 2.5_wp,  &
      b53 = -70.0_wp/27.0_wp, b54 = 35.0_wp/27.0_wp
    real(wp), parameter :: b61 = 1631.0_wp/55296.0_wp,  &
      b62 = 175.0_wp/512.0_wp, b63 = 575.0_wp/13824.0_wp,  &
      b64 = 44275.0_wp/110592.0_wp, b65 = 253.0_wp/4096.0_wp
    real(wp), parameter :: c1 = 37.0_wp/378.0_wp, c3 = 250.0_wp/621.0_wp,  &
      c4 = 125.0_wp/594.0_wp, c6 = 512.0_wp/1771.0_wp
    real(wp), parameter :: dc5 = -277.0_wp/14336.0_wp
    real(wp), parameter :: dc1 = c1 - 2825.0_wp/27648.0_wp,  &
      dc3 = c3 - 18575.0_wp/48384.0_wp, dc4 = c4 - 13525.0_wp/55296.0_wp,  &
      dc6 = c6 - 0.25_wp

    real(wp), dimension(size(y)) :: ytemp
    real(wp), dimension(size(y)) :: ak2
    real(wp), dimension(size(y)) :: ak3
    real(wp), dimension(size(y)) :: ak4
    real(wp), dimension(size(y)) :: ak5
    real(wp), dimension(size(y)) :: ak6

    ytemp = y + b21*dt*dydt
    call f(ytemp, t + a2*dt, ak2)
    ytemp = y + dt*(b31*dydt + b32*ak2)
    call f(ytemp, t + a3*dt, ak3)
    ytemp = y + dt*(b41*dydt + b42*ak2 + b43*ak3)
    call f(ytemp, t + a4*dt, ak4)
    ytemp = y + dt*(b51*dydt + b52*ak2 + b53*ak3 + b54*ak4)
    call f(ytemp, t + a5*dt, ak5)
    ytemp = y + dt*(b61*dydt + b62*ak2 + b63*ak3 + b64*ak4 + b65*ak5)
    call f(ytemp, t + a6*dt, ak6)
    dyout = dt*(c1*dydt + c3*ak3 + c4*ak4 + c6*ak6)
    yerr = dt*(dc1*dydt + dc3*ak3 + dc4*ak4 + dc5*ak5 + dc6*ak6)

  end subroutine rkck

  pure function cubic_interp(t, y0, dy0, y1, dy1) result(y)
    !! Cubic interpolation between two points, given the values and the
    !! derivatives (with respect to `t` in [0, 1]) at the endpoints.
    !! Port of helper_functions.cubicInterpFunction (Bezier form).

    real(wp), intent(in) :: t
    real(wp), intent(in) :: y0(:)
    real(wp), intent(in) :: dy0(:)
    real(wp), intent(in) :: y1(:)
    real(wp), intent(in) :: dy1(:)
    real(wp) :: y(size(y0))

    real(wp) :: mt

    mt = 1.0_wp - t
    y = y0*mt**3 + 3.0_wp*(y0 + dy0/3.0_wp)*mt*mt*t  &
      + 3.0_wp*(y1 - dy1/3.0_wp)*mt*t*t + y1*t**3

  end function cubic_interp

  pure function deriv14_const_dx(y) result(dy)
    !! dy/dx to fourth order in dx = 1, where the derivative is taken along
    !! the first dimension of `y` (one row per point).
    !! Port of helper_functions.deriv14_const_dx (which acts on the last
    !! axis of the transposed array, i.e. the same thing).
    !! Requires size(y, 1) >= 5.

    real(wp), intent(in) :: y(:, :)
    real(wp) :: dy(size(y, 1), size(y, 2))

    integer :: n
    integer :: i

    n = size(y, 1)
    do i = 3, n - 2
      dy(i, :) = y(i-2, :) - 8.0_wp*y(i-1, :) + 8.0_wp*y(i+1, :) - y(i+2, :)
    end do
    dy(1, :) = -25.0_wp*y(1, :) + 48.0_wp*y(2, :) - 36.0_wp*y(3, :)  &
      + 16.0_wp*y(4, :) - 3.0_wp*y(5, :)
    dy(2, :) = -3.0_wp*y(1, :) - 10.0_wp*y(2, :) + 18.0_wp*y(3, :)  &
      - 6.0_wp*y(4, :) + y(5, :)
    dy(n-1, :) = 3.0_wp*y(n, :) + 10.0_wp*y(n-1, :) - 18.0_wp*y(n-2, :)  &
      + 6.0_wp*y(n-3, :) - y(n-4, :)
    dy(n, :) = 25.0_wp*y(n, :) - 48.0_wp*y(n-1, :) + 36.0_wp*y(n-2, :)  &
      - 16.0_wp*y(n-3, :) + 3.0_wp*y(n-4, :)
    dy = dy/12.0_wp

  end function deriv14_const_dx

  subroutine nbspld2(tk, x, k, nmat, dnmat, d2nmat)
    !! B-spline basis functions for the knots `tk` evaluated at the points
    !! `x`, together with their first and second derivatives.
    !! Port of helper_functions.Nbspld2.
    !! The output matrices have shape (size(x), size(tk)-k-1).

    real(wp), intent(in) :: tk(:)
      !! Knots which define the basis functions.
    real(wp), intent(in) :: x(:)
      !! Values at which to calculate the functions.
    integer, intent(in) :: k
      !! Order of the spline. Must satisfy k <= size(tk) - 2.
    real(wp), allocatable, intent(out) :: nmat(:, :)
    real(wp), allocatable, intent(out) :: dnmat(:, :)
    real(wp), allocatable, intent(out) :: d2nmat(:, :)

    integer :: nt
    integer :: nx
    integer :: kk
    integer :: j
    integer :: i
    integer :: ncols
    real(wp), allocatable :: nb(:, :)
    real(wp), allocatable :: dnb(:, :)
    real(wp), allocatable :: d2nb(:, :)
    real(wp), allocatable :: idt(:)
    real(wp), dimension(size(x)) :: nj
    real(wp), dimension(size(x)) :: dnj
    real(wp), dimension(size(x)) :: d2nj

    nt = size(tk)
    nx = size(x)
    if (k > nt - 2) then
      error stop "nbspld2: require that k <= size(tk) - 2"
    end if

    allocate(nb(nx, nt-1))
    allocate(dnb(nx, nt-1))
    allocate(d2nb(nx, nt-1))
    dnb = 0.0_wp
    d2nb = 0.0_wp
    do j = 1, nt - 1
      do i = 1, nx
        nb(i, j) = merge(1.0_wp, 0.0_wp,  &
          x(i) > tk(j) .and. x(i) <= tk(j+1))
      end do
    end do

    do kk = 1, k
      allocate(idt(nt-kk))
      do j = 1, nt - kk
        if (tk(j+kk) /= tk(j)) then
          idt(j) = 1.0_wp/(tk(j+kk) - tk(j))
        else
          idt(j) = 0.0_wp
        end if
      end do
      ! Update columns in ascending order; new column j only depends on the
      ! old columns j and j+1, which are still untouched when we get there.
      do j = 1, nt - kk - 1
        d2nj = d2nb(:, j)*(x - tk(j))*idt(j)  &
          - d2nb(:, j+1)*(x - tk(j+kk+1))*idt(j+1)  &
          + 2.0_wp*dnb(:, j)*idt(j) - 2.0_wp*dnb(:, j+1)*idt(j+1)
        dnj = dnb(:, j)*(x - tk(j))*idt(j)  &
          - dnb(:, j+1)*(x - tk(j+kk+1))*idt(j+1)  &
          + nb(:, j)*idt(j) - nb(:, j+1)*idt(j+1)
        nj = nb(:, j)*(x - tk(j))*idt(j)  &
          - nb(:, j+1)*(x - tk(j+kk+1))*idt(j+1)
        d2nb(:, j) = d2nj
        dnb(:, j) = dnj
        nb(:, j) = nj
      end do
      deallocate(idt)
    end do

    ncols = nt - k - 1
    nmat = nb(:, 1:ncols)
    dnmat = dnb(:, 1:ncols)
    d2nmat = d2nb(:, 1:ncols)

  end subroutine nbspld2

  function monotonic_indices(x) result(idx)
    !! Indices of `x` such that x(idx) is purely increasing.
    !! Port of helper_functions.monotonicIndices.

    real(wp), intent(in) :: x(:)
    integer, allocatable :: idx(:)

    real(wp), allocatable :: xw(:)
    integer, allocatable :: iw(:)
    logical :: is_reversed
    integer :: n
    integer :: i
    integer :: m

    n = size(x)
    if (x(1) > x(n)) then
      xw = x(n:1:-1)
      is_reversed = .true.
    else
      xw = x
      is_reversed = .false.
    end if

    allocate(iw(n))
    iw(1) = 1
    m = 1
    do i = 2, n - 1
      if (xw(i) > xw(iw(m)) .and. xw(i) < xw(n)) then
        m = m + 1
        iw(m) = i
      end if
    end do
    m = m + 1
    iw(m) = n

    if (is_reversed) then
      idx = n + 1 - iw(1:m)
    else
      idx = iw(1:m)
    end if

  end function monotonic_indices

  pure function cumtrapz(y) result(c)
    !! Cumulative trapezoidal integral of `y` with unit spacing, starting
    !! at zero (equivalent to scipy cumulative_trapezoid(y, initial=0)).

    real(wp), intent(in) :: y(:)
    real(wp) :: c(size(y))

    integer :: i

    c(1) = 0.0_wp
    do i = 2, size(y)
      c(i) = c(i-1) + 0.5_wp*(y(i-1) + y(i))
    end do

  end function cumtrapz

  pure function simpson(y, x) result(s)
    !! Composite Simpson integration of the samples `y` at the (possibly
    !! non-uniformly spaced) points `x`. For an even number of samples the
    !! last interval is handled with a parabolic (Cartwright) correction,
    !! matching scipy.integrate.simpson.

    real(wp), intent(in) :: y(:)
    real(wp), intent(in) :: x(:)
    real(wp) :: s

    integer :: n
    integer :: m
    integer :: i
    real(wp) :: h1
    real(wp) :: h2
    real(wp) :: hsum

    n = size(y)
    s = 0.0_wp
    if (n < 2) return
    if (n == 2) then
      s = 0.5_wp*(y(1) + y(2))*(x(2) - x(1))
      return
    end if

    m = n
    if (mod(n, 2) == 0) m = n - 1
    do i = 1, m - 2, 2
      h1 = x(i+1) - x(i)
      h2 = x(i+2) - x(i+1)
      hsum = h1 + h2
      s = s + hsum/6.0_wp*((2.0_wp - h2/h1)*y(i)  &
        + hsum**2/(h1*h2)*y(i+1) + (2.0_wp - h1/h2)*y(i+2))
    end do
    if (mod(n, 2) == 0) then
      ! Parabola through the last three points, integrated over the last
      ! interval only.
      h1 = x(n-1) - x(n-2)
      h2 = x(n) - x(n-1)
      s = s + y(n)*(2.0_wp*h2**2 + 3.0_wp*h1*h2)/(6.0_wp*(h1 + h2))
      s = s + y(n-1)*(h2**2 + 3.0_wp*h1*h2)/(6.0_wp*h1)
      s = s - y(n-2)*h2**3/(6.0_wp*h1*(h1 + h2))
    end if

  end function simpson

  pure function linspace(a, b, n) result(x)
    !! `n` evenly spaced values from `a` to `b` (inclusive), like
    !! numpy.linspace.

    real(wp), intent(in) :: a
    real(wp), intent(in) :: b
    integer, intent(in) :: n
    real(wp) :: x(n)

    integer :: i

    if (n == 1) then
      x(1) = a
      return
    end if
    do i = 1, n
      x(i) = a + (b - a)*real(i - 1, wp)/real(n - 1, wp)
    end do
    x(n) = b

  end function linspace

  pure function interp_linear(xp, yp, x) result(y)
    !! Linear interpolation of the data (xp, yp) at the point `x`.
    !! `xp` must be strictly increasing. Points outside the range are
    !! extrapolated from the end segments (like a k=1 scipy spline).

    real(wp), intent(in) :: xp(:)
    real(wp), intent(in) :: yp(:)
    real(wp), intent(in) :: x
    real(wp) :: y

    integer :: lo
    integer :: hi
    integer :: mid
    real(wp) :: w

    lo = 1
    hi = size(xp)
    do while (hi - lo > 1)
      mid = (lo + hi)/2
      if (xp(mid) > x) then
        hi = mid
      else
        lo = mid
      end if
    end do
    w = (x - xp(lo))/(xp(hi) - xp(lo))
    y = yp(lo)*(1.0_wp - w) + yp(hi)*w

  end function interp_linear

  subroutine lstsq(a, b, x, status)
    !! Least-squares solution of the overdetermined system a x = b via the
    !! normal equations (a^T a) x = a^T b. Fortran replacement for the
    !! numpy.linalg.lstsq calls in pathDeformation.py. The spline basis
    !! matrices used there are well conditioned, so the normal equations
    !! are adequate.

    real(wp), intent(in) :: a(:, :)
      !! Shape (m, n) with m >= n.
    real(wp), intent(in) :: b(:, :)
      !! Shape (m, nrhs).
    real(wp), allocatable, intent(out) :: x(:, :)
      !! Shape (n, nrhs).
    integer, intent(out) :: status

    real(wp), allocatable :: ata(:, :)
    real(wp), allocatable :: atb(:)
    real(wp), allocatable :: sol(:)
    integer :: j
    integer :: st

    status = status_ok
    ata = matmul(transpose(a), a)
    allocate(x(size(a, 2), size(b, 2)))
    do j = 1, size(b, 2)
      atb = matmul(b(:, j), a)
      st = 0
      call gauss_jordan(ata, atb, sol, st)
      if (st /= 0) then
        status = err_numerical
        return
      end if
      x(:, j) = sol
      deallocate(sol)
    end do

  end subroutine lstsq

end module cosmotransitions__helpers
