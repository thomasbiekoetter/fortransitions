module cosmotransitions__optimize
!! One-dimensional root finding and minimization, replacing the
!! scipy.optimize routines used by the Python package:
!!   brentq            -> brentq
!!   fminbound         -> minimize_bounded
!!   fmin (1d usage)   -> minimize_unbounded (bracket + bounded Brent)

  use cosmotransitions__config, only : wp
  use cosmotransitions__config, only : status_ok
  use cosmotransitions__config, only : err_numerical

  implicit none

  private

  public :: scalar_func
  public :: brentq
  public :: minimize_bounded
  public :: minimize_unbounded

  abstract interface
    function scalar_func(x) result(y)
      import :: wp
      implicit none
      real(wp), intent(in) :: x
      real(wp) :: y
    end function scalar_func
  end interface

contains

  function brentq(f, xa, xb, status, xtol, rtol, maxiter) result(root)
    !! Find a root of f in the bracketing interval [xa, xb] using Brent's
    !! method. Port of scipy.optimize.brentq (defaults match scipy).

    procedure(scalar_func) :: f
    real(wp), intent(in) :: xa
    real(wp), intent(in) :: xb
    integer, intent(out) :: status
    real(wp), intent(in), optional :: xtol
    real(wp), intent(in), optional :: rtol
    integer, intent(in), optional :: maxiter
    real(wp) :: root

    real(wp) :: xtol_
    real(wp) :: rtol_
    integer :: maxiter_
    real(wp) :: xpre
    real(wp) :: xcur
    real(wp) :: xblk
    real(wp) :: fpre
    real(wp) :: fcur
    real(wp) :: fblk
    real(wp) :: spre
    real(wp) :: scur
    real(wp) :: sbis
    real(wp) :: stry
    real(wp) :: dpre
    real(wp) :: dblk
    real(wp) :: delta
    integer :: i

    xtol_ = 2.0e-12_wp
    rtol_ = 4.0_wp*epsilon(1.0_wp)
    maxiter_ = 100
    if (present(xtol)) xtol_ = xtol
    if (present(rtol)) rtol_ = rtol
    if (present(maxiter)) maxiter_ = maxiter

    status = status_ok
    xpre = xa
    xcur = xb
    xblk = 0.0_wp
    fblk = 0.0_wp
    spre = 0.0_wp
    scur = 0.0_wp
    root = xcur

    fpre = f(xpre)
    fcur = f(xcur)
    if (fpre == 0.0_wp) then
      root = xpre
      return
    end if
    if (fcur == 0.0_wp) then
      root = xcur
      return
    end if
    if (sign(1.0_wp, fpre) == sign(1.0_wp, fcur)) then
      status = err_numerical  ! f(xa) and f(xb) must have different signs
      return
    end if

    do i = 1, maxiter_
      if (fpre /= 0.0_wp .and. fcur /= 0.0_wp .and.  &
          sign(1.0_wp, fpre) /= sign(1.0_wp, fcur)) then
        xblk = xpre
        fblk = fpre
        spre = xcur - xpre
        scur = xcur - xpre
      end if
      if (abs(fblk) < abs(fcur)) then
        xpre = xcur
        xcur = xblk
        xblk = xpre
        fpre = fcur
        fcur = fblk
        fblk = fpre
      end if

      delta = 0.5_wp*(xtol_ + rtol_*abs(xcur))
      sbis = 0.5_wp*(xblk - xcur)
      if (fcur == 0.0_wp .or. abs(sbis) < delta) then
        root = xcur
        return
      end if

      if (abs(spre) > delta .and. abs(fcur) < abs(fpre)) then
        if (xpre == xblk) then
          ! Interpolate (secant)
          stry = -fcur*(xcur - xpre)/(fcur - fpre)
        else
          ! Extrapolate (inverse quadratic)
          dpre = (fpre - fcur)/(xpre - xcur)
          dblk = (fblk - fcur)/(xblk - xcur)
          stry = -fcur*(fblk*dblk - fpre*dpre)/(dblk*dpre*(fblk - fpre))
        end if
        if (2.0_wp*abs(stry) < min(abs(spre), 3.0_wp*abs(sbis) - delta)) then
          ! Good short step
          spre = scur
          scur = stry
        else
          ! Bisect
          spre = sbis
          scur = sbis
        end if
      else
        ! Bisect
        spre = sbis
        scur = sbis
      end if

      xpre = xcur
      fpre = fcur
      if (abs(scur) > delta) then
        xcur = xcur + scur
      else
        xcur = xcur + sign(delta, sbis)
      end if
      fcur = f(xcur)
    end do

    root = xcur

  end function brentq

  function minimize_bounded(f, x1, x2, xatol, status) result(xf)
    !! Bounded minimization of f on the interval [x1, x2].
    !! Port of scipy.optimize.fminbound (golden section with parabolic
    !! interpolation; maxfun = 500).

    procedure(scalar_func) :: f
    real(wp), intent(in) :: x1
    real(wp), intent(in) :: x2
    real(wp), intent(in) :: xatol
    integer, intent(out) :: status
    real(wp) :: xf

    real(wp), parameter :: golden_mean = 0.5_wp*(3.0_wp - sqrt(5.0_wp))
    integer, parameter :: maxfun = 500

    real(wp) :: sqrt_eps
    real(wp) :: a
    real(wp) :: b
    real(wp) :: fulc
    real(wp) :: nfc
    real(wp) :: rat
    real(wp) :: e
    real(wp) :: x
    real(wp) :: fx
    real(wp) :: fu
    real(wp) :: ffulc
    real(wp) :: fnfc
    real(wp) :: xm
    real(wp) :: tol1
    real(wp) :: tol2
    real(wp) :: r
    real(wp) :: q
    real(wp) :: p
    real(wp) :: si
    integer :: num
    logical :: golden

    status = status_ok
    sqrt_eps = sqrt(2.2e-16_wp)
    a = x1
    b = x2
    fulc = a + golden_mean*(b - a)
    nfc = fulc
    xf = fulc
    rat = 0.0_wp
    e = 0.0_wp
    x = xf
    fx = f(x)
    num = 1
    ffulc = fx
    fnfc = fx
    xm = 0.5_wp*(a + b)
    tol1 = sqrt_eps*abs(xf) + xatol/3.0_wp
    tol2 = 2.0_wp*tol1

    do while (abs(xf - xm) > tol2 - 0.5_wp*(b - a))
      golden = .true.
      ! Check for parabolic fit
      if (abs(e) > tol1) then
        golden = .false.
        r = (xf - nfc)*(fx - ffulc)
        q = (xf - fulc)*(fx - fnfc)
        p = (xf - fulc)*q - (xf - nfc)*r
        q = 2.0_wp*(q - r)
        if (q > 0.0_wp) p = -p
        q = abs(q)
        r = e
        e = rat
        ! Check for acceptability of parabola
        if (abs(p) < abs(0.5_wp*q*r) .and. p > q*(a - xf) .and.  &
            p < q*(b - xf)) then
          rat = p/q
          x = xf + rat
          if ((x - a) < tol2 .or. (b - x) < tol2) then
            si = sgn(xm - xf)
            if (xm - xf == 0.0_wp) si = 1.0_wp
            rat = tol1*si
          end if
        else
          golden = .true.
        end if
      end if
      if (golden) then
        ! Do a golden-section step
        if (xf >= xm) then
          e = a - xf
        else
          e = b - xf
        end if
        rat = golden_mean*e
      end if

      si = sgn(rat)
      if (rat == 0.0_wp) si = 1.0_wp
      x = xf + si*max(abs(rat), tol1)
      fu = f(x)
      num = num + 1

      if (fu <= fx) then
        if (x >= xf) then
          a = xf
        else
          b = xf
        end if
        fulc = nfc
        ffulc = fnfc
        nfc = xf
        fnfc = fx
        xf = x
        fx = fu
      else
        if (x < xf) then
          a = x
        else
          b = x
        end if
        if (fu <= fnfc .or. nfc == xf) then
          fulc = nfc
          ffulc = fnfc
          nfc = x
          fnfc = fu
        else if (fu <= ffulc .or. fulc == xf .or. fulc == nfc) then
          fulc = x
          ffulc = fu
        end if
      end if

      xm = 0.5_wp*(a + b)
      tol1 = sqrt_eps*abs(xf) + xatol/3.0_wp
      tol2 = 2.0_wp*tol1

      if (num >= maxfun) then
        status = err_numerical
        return
      end if
    end do

  end function minimize_bounded

  function minimize_unbounded(f, x0, xtol, status) result(xmin)
    !! Local minimization of f starting from x0 (no bounds). Replacement for
    !! the one-dimensional scipy.optimize.fmin calls in the Python package:
    !! the minimum is first bracketed (like scipy.optimize.bracket) and then
    !! refined with the bounded Brent routine.

    procedure(scalar_func) :: f
    real(wp), intent(in) :: x0
    real(wp), intent(in) :: xtol
    integer, intent(out) :: status
    real(wp) :: xmin

    real(wp) :: xa
    real(wp) :: xb
    real(wp) :: xc

    xmin = x0
    call bracket_min(f, x0, xa, xb, xc, status)
    if (status /= status_ok) return
    xmin = minimize_bounded(f, min(xa, xc), max(xa, xc), xtol, status)

  end function minimize_unbounded

  subroutine bracket_min(f, x0, xa, xb, xc, status)
    !! Bracket a local minimum of f downhill from x0, such that
    !! f(xb) <= min(f(xa), f(xc)) with xb between xa and xc.
    !! Port of scipy.optimize.bracket.

    procedure(scalar_func) :: f
    real(wp), intent(in) :: x0
    real(wp), intent(out) :: xa
    real(wp), intent(out) :: xb
    real(wp), intent(out) :: xc
    integer, intent(out) :: status

    real(wp), parameter :: gold = 1.618034_wp
    real(wp), parameter :: verysmall = 1.0e-21_wp
    real(wp), parameter :: grow_limit = 110.0_wp
    integer, parameter :: maxiter = 1000

    real(wp) :: fa
    real(wp) :: fb
    real(wp) :: fc
    real(wp) :: fw
    real(wp) :: tmp1
    real(wp) :: tmp2
    real(wp) :: val
    real(wp) :: denom
    real(wp) :: w
    real(wp) :: wlim
    real(wp) :: swap
    integer :: iter

    status = status_ok
    xa = x0
    ! Initial step comparable to the one scipy.optimize.fmin uses for its
    ! starting simplex.
    if (x0 /= 0.0_wp) then
      xb = x0*1.05_wp
    else
      xb = 0.00025_wp
    end if
    fa = f(xa)
    fb = f(xb)
    if (fa < fb) then
      swap = xa
      xa = xb
      xb = swap
      swap = fa
      fa = fb
      fb = swap
    end if
    xc = xb + gold*(xb - xa)
    fc = f(xc)
    iter = 0

    do while (fc < fb)
      tmp1 = (xb - xa)*(fb - fc)
      tmp2 = (xb - xc)*(fb - fa)
      val = tmp2 - tmp1
      if (abs(val) < verysmall) then
        denom = 2.0_wp*verysmall
      else
        denom = 2.0_wp*val
      end if
      w = xb - ((xb - xc)*tmp2 - (xb - xa)*tmp1)/denom
      wlim = xb + grow_limit*(xc - xb)
      if (iter > maxiter) then
        status = err_numerical
        return
      end if
      iter = iter + 1
      if ((w - xc)*(xb - w) > 0.0_wp) then
        fw = f(w)
        if (fw < fc) then
          xa = xb
          xb = w
          fa = fb
          fb = fw
          return
        else if (fw > fb) then
          xc = w
          fc = fw
          return
        end if
        w = xc + gold*(xc - xb)
        fw = f(w)
      else if ((w - wlim)*(wlim - xc) >= 0.0_wp) then
        w = wlim
        fw = f(w)
      else if ((w - wlim)*(xc - w) > 0.0_wp) then
        fw = f(w)
        if (fw < fc) then
          xb = xc
          xc = w
          w = xc + gold*(xc - xb)
          fb = fc
          fc = fw
          fw = f(w)
        end if
      else
        w = xc + gold*(xc - xb)
        fw = f(w)
      end if
      xa = xb
      xb = xc
      xc = w
      fa = fb
      fb = fc
      fc = fw
    end do

  end subroutine bracket_min

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

end module cosmotransitions__optimize
