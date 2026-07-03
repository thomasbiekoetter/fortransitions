module cosmotransitions__special
!! Bessel functions of (possibly) fractional order, replacing
!! scipy.special.iv and scipy.special.jv in tunneling1D.exactSolution.
!! For the instanton problem the order is nu = (alpha - 1)/2 with alpha the
!! friction coefficient, so nu and nu +- 1 are typically half-integers
!! (alpha = 2) or integers (alpha = 1, 3).

  use cosmotransitions__config, only : wp
  use cosmotransitions__config, only : pi

  implicit none

  private

  public :: besseli_nu
  public :: besselj_nu

  real(wp), parameter :: bessel_sat = 1.0e150_wp
    !! Saturation value for the exponentially growing I_nu. scipy lets iv()
    !! overflow to inf, which the Python code deliberately ignores because
    !! only the comparison "|phi - phi_absMin| > cutoff" matters. Saturating
    !! instead of overflowing keeps the debug builds (-ffpe-trap=overflow)
    !! usable while preserving that behaviour.

contains

  function besseli_nu(nu, x) result(iv)
    !! Modified Bessel function of the first kind I_nu(x) for x >= 0 and
    !! real order nu >= -1/2. Uses the ascending series for small x and the
    !! large-argument asymptotic expansion otherwise.

    real(wp), intent(in) :: nu
    real(wp), intent(in) :: x
    real(wp) :: iv

    real(wp) :: term
    real(wp) :: s
    real(wp) :: xh
    real(wp) :: ak
    integer :: k

    if (x <= 0.0_wp) then
      iv = merge(1.0_wp, 0.0_wp, nu == 0.0_wp)
      return
    end if

    if (x <= 15.0_wp) then
      ! Ascending series: I_nu(x) = sum_k (x/2)^(nu+2k) / (k! Gamma(nu+k+1))
      xh = 0.5_wp*x
      term = exp(nu*log(xh) - log_gamma(nu + 1.0_wp))
      s = term
      do k = 1, 200
        term = term*xh*xh/(real(k, wp)*(real(k, wp) + nu))
        s = s + term
        if (term < s*1.0e-17_wp) exit
      end do
      iv = s
    else if (x > 690.0_wp) then
      iv = bessel_sat
    else
      ! Asymptotic expansion:
      ! I_nu(x) ~ e^x/sqrt(2 pi x) * sum_k (-1)^k a_k(nu)/x^k with
      ! a_k = (4 nu^2 - 1)(4 nu^2 - 9)...(4 nu^2 - (2k-1)^2)/(k! 8^k).
      ak = 1.0_wp
      s = 1.0_wp
      do k = 1, 10
        ak = -ak*(4.0_wp*nu*nu - real(2*k - 1, wp)**2)/(real(k, wp)*8.0_wp*x)
        s = s + ak
      end do
      iv = exp(x)/sqrt(2.0_wp*pi*x)*s
      if (iv > bessel_sat) iv = bessel_sat
    end if

  end function besseli_nu

  function besselj_nu(nu, x) result(jv)
    !! Bessel function of the first kind J_nu(x) for x > 0. Handles integer
    !! orders (via the intrinsic bessel_jn), half-integer orders (closed
    !! forms plus the standard recurrence) and, for moderate x, arbitrary
    !! real orders via the ascending series.

    real(wp), intent(in) :: nu
    real(wp), intent(in) :: x
    real(wp) :: jv

    real(wp) :: jm
    real(wp) :: j0_
    real(wp) :: jp
    real(wp) :: prefac
    real(wp) :: term
    real(wp) :: s
    real(wp) :: xh
    real(wp) :: mu
    integer :: n
    integer :: k

    if (is_near_integer(nu)) then
      n = nint(nu)
      if (n >= 0) then
        jv = bessel_jn(n, x)
      else
        jv = real((-1)**(-n), wp)*bessel_jn(-n, x)
      end if
      return
    end if

    if (is_near_integer(nu - 0.5_wp)) then
      ! Half-integer order: start from the closed forms for J_{-1/2} and
      ! J_{1/2} and recurse upward, J_{mu+1} = (2 mu/x) J_mu - J_{mu-1}.
      if (nu < -0.6_wp) then
        error stop "besselj_nu: half-integer orders below -1/2 not supported"
      end if
      prefac = sqrt(2.0_wp/(pi*x))
      jm = prefac*cos(x)   ! J_{-1/2}
      j0_ = prefac*sin(x)  ! J_{+1/2}
      if (is_near_integer(nu + 0.5_wp) .and. nint(nu + 0.5_wp) == 0) then
        jv = jm
        return
      end if
      mu = 0.5_wp
      jv = j0_
      do while (mu < nu - 0.25_wp)
        jp = (2.0_wp*mu/x)*jv - jm
        jm = jv
        jv = jp
        mu = mu + 1.0_wp
      end do
      return
    end if

    ! Generic real order: ascending series
    ! J_nu(x) = sum_k (-1)^k (x/2)^(nu+2k) / (k! Gamma(nu+k+1)).
    ! This suffers from cancellation for large x, so it is restricted to
    ! moderate arguments (which is all the instanton solver ever needs).
    if (x > 30.0_wp) then
      error stop "besselj_nu: series evaluation not supported for x > 30 "//  &
        "at non-(half-)integer order"
    end if
    xh = 0.5_wp*x
    term = exp(nu*log(xh) - log_gamma(nu + 1.0_wp))
    s = term
    do k = 1, 200
      term = -term*xh*xh/(real(k, wp)*(real(k, wp) + nu))
      s = s + term
      if (abs(term) < abs(s)*1.0e-17_wp + 1.0e-300_wp) exit
    end do
    jv = s

  end function besselj_nu

  pure function is_near_integer(x) result(res)

    real(wp), intent(in) :: x
    logical :: res

    res = abs(x - real(nint(x), wp)) < 1.0e-12_wp

  end function is_near_integer

end module cosmotransitions__special
