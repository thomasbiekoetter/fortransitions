module c_interface__callbacks
!! bind(C) versions of test potentials, passed to the C interface as
!! function pointers (as ctypes does for the Python package), together
!! with the same potentials for the native Fortran API.

  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_c_binding, only : c_double

  use fortransitions, only : wp
  use fortransitions, only : potential_nd

  implicit none

  private

  public :: v1_c
  public :: dv1_c
  public :: v1
  public :: dv1
  public :: v2d_c
  public :: grad2d_c
  public :: grad2d_many_c
  public :: nbatch
  public :: pot2d
  public :: ncalls

  integer, save :: ncalls = 0
    !! Number of potential evaluations through the C callbacks.
  integer, save :: nbatch = 0
    !! Number of calls of the batched gradient callback.

  type, extends(potential_nd) :: pot2d
  contains
    procedure :: v => pot2d_v
    procedure :: grad => pot2d_grad
  end type pot2d

contains

  ! One-dimensional thin-walled potential of test/check.f90.

  function v1(phi) result(y)
    real(wp), intent(in) :: phi
    real(wp) :: y
    y = 0.25_wp*phi**4 - 0.49_wp*phi**3 + 0.235_wp*phi**2
  end function v1

  function dv1(phi) result(y)
    real(wp), intent(in) :: phi
    real(wp) :: y
    y = phi*(phi - 0.47_wp)*(phi - 1.0_wp)
  end function dv1

  function v1_c(phi) result(y) bind(C)
    real(c_double), value :: phi
    real(c_double) :: y
    ncalls = ncalls + 1
    y = v1(phi)
  end function v1_c

  function dv1_c(phi) result(y) bind(C)
    real(c_double), value :: phi
    real(c_double) :: y
    y = dv1(phi)
  end function dv1_c

  ! Two-field potential of examples/fullTunneling.py with c = 5, fx = 0,
  ! fy = 80 (thick-walled case of test/check.f90).

  function v2d(x) result(y)
    real(wp), intent(in) :: x(:)
    real(wp) :: y
    real(wp) :: r1, r2, r3
    r1 = x(1)*x(1) + 5.0_wp*x(2)*x(2)
    r2 = 5.0_wp*(x(1) - 1.0_wp)**2 + (x(2) - 1.0_wp)**2
    r3 = 80.0_wp*(0.25_wp*x(2)**4 - x(2)**3/3.0_wp)
    y = r1*r2 + r3
  end function v2d

  function grad2d(x) result(dv)
    real(wp), intent(in) :: x(:)
    real(wp) :: dv(size(x))
    real(wp) :: r1, r2
    r1 = x(1)*x(1) + 5.0_wp*x(2)*x(2)
    r2 = 5.0_wp*(x(1) - 1.0_wp)**2 + (x(2) - 1.0_wp)**2
    dv(1) = r1*2.0_wp*5.0_wp*(x(1) - 1.0_wp) + 2.0_wp*x(1)*r2
    dv(2) = r1*2.0_wp*(x(2) - 1.0_wp) + 2.0_wp*5.0_wp*x(2)*r2  &
      + 80.0_wp*x(2)*x(2)*(x(2) - 1.0_wp)
  end function grad2d

  function v2d_c(n, x) result(y) bind(C)
    integer(c_int), value :: n
    real(c_double), intent(in) :: x(n)
    real(c_double) :: y
    ncalls = ncalls + 1
    y = v2d(x)
  end function v2d_c

  subroutine grad2d_c(n, x, g) bind(C)
    integer(c_int), value :: n
    real(c_double), intent(in) :: x(n)
    real(c_double), intent(out) :: g(n)
    g = grad2d(x)
  end subroutine grad2d_c

  subroutine grad2d_many_c(npts, n, x, g) bind(C)
    integer(c_int), value :: npts
    integer(c_int), value :: n
    real(c_double), intent(in) :: x(npts, n)
    real(c_double), intent(out) :: g(npts, n)
    integer :: i
    nbatch = nbatch + 1
    do i = 1, npts
      g(i, :) = grad2d(x(i, :))
    end do
  end subroutine grad2d_many_c

  function pot2d_v(self, x) result(y)
    class(pot2d), intent(inout) :: self
    real(wp), intent(in) :: x(:)
    real(wp) :: y
    y = v2d(x)
  end function pot2d_v

  function pot2d_grad(self, x) result(dv)
    class(pot2d), intent(inout) :: self
    real(wp), intent(in) :: x(:)
    real(wp) :: dv(size(x))
    dv = grad2d(x)
  end function pot2d_grad

end module c_interface__callbacks

program c_interface
!! Tests the bind(C) interface used by the Python package
!! fortransitions (modules fortransitions__c_*): the same problems are
!! solved through the C interface, with the potential given by C function
!! pointers, and through the native Fortran API. Both must give identical
!! results, since they run the same code on the same potential.

  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_c_binding, only : c_double
  use, intrinsic :: iso_c_binding, only : c_ptr
  use, intrinsic :: iso_c_binding, only : c_funloc
  use, intrinsic :: iso_c_binding, only : c_null_funptr
  use, intrinsic :: iso_c_binding, only : c_funptr
  use, intrinsic :: iso_c_binding, only : c_associated

  use fortransitions, only : wp
  use fortransitions, only : status_ok
  use fortransitions, only : err_stable
  use fortransitions, only : pot1d_func
  use fortransitions, only : single_field_instanton
  use fortransitions, only : profile1d
  use fortransitions, only : full_tunneling
  use fortransitions, only : full_tunneling_result

  use fortransitions__c_tunneling1d, only : fortransitions_sfi_new
  use fortransitions__c_tunneling1d, only : fortransitions_sfi_params
  use fortransitions__c_tunneling1d, only : fortransitions_sfi_find_profile
  use fortransitions__c_tunneling1d, only : fortransitions_sfi_find_action
  use fortransitions__c_tunneling1d, only : fortransitions_sfi_free
  use fortransitions__c_tunneling1d, only : fortransitions_profile_size
  use fortransitions__c_tunneling1d, only : fortransitions_profile_arrays
  use fortransitions__c_tunneling1d, only : fortransitions_profile_free
  use fortransitions__c_pathdeformation, only : fortransitions_full_tunneling
  use fortransitions__c_pathdeformation, only : fortransitions_ft_result_sizes
  use fortransitions__c_pathdeformation, only :  &
    fortransitions_ft_result_scalars
  use fortransitions__c_pathdeformation, only :  &
    fortransitions_ft_result_arrays
  use fortransitions__c_pathdeformation, only : fortransitions_ft_result_free

  use c_interface__callbacks, only : v1_c
  use c_interface__callbacks, only : dv1_c
  use c_interface__callbacks, only : v1
  use c_interface__callbacks, only : dv1
  use c_interface__callbacks, only : v2d_c
  use c_interface__callbacks, only : grad2d_c
  use c_interface__callbacks, only : grad2d_many_c
  use c_interface__callbacks, only : nbatch
  use c_interface__callbacks, only : pot2d
  use c_interface__callbacks, only : ncalls

  implicit none

  integer :: nfail

  nfail = 0

  call test_sfi(.true.)
  call test_sfi(.false.)
  call test_sfi_stable()
  call test_full_tunneling()

  if (nfail > 0) then
    print "(a,i0,a)", "FAILED: ", nfail, " test(s) failed."
    error stop 1
  else
    print "(a)", "All tests passed."
  end if

contains

  subroutine test_sfi(analytic_dv)
    !! single_field_instanton through the C interface vs. natively, with
    !! the analytic derivative or with finite differences.

    logical, intent(in) :: analytic_dv

    type(pot1d_func), target :: pot
    type(single_field_instanton) :: inst
    type(profile1d) :: prof
    real(wp) :: action
    real(wp) :: action_pot
    real(wp) :: action_kin
    type(c_ptr) :: h
    type(c_ptr) :: hp
    integer(c_int) :: status
    integer(c_int) :: n
    real(c_double) :: phi_bar
    real(c_double) :: rscale
    real(c_double) :: rerr
    real(c_double) :: action_c
    real(c_double) :: action_pot_c
    real(c_double) :: action_kin_c
    real(c_double), allocatable :: r(:)
    real(c_double), allocatable :: phi(:)
    real(c_double), allocatable :: dphi(:)
    character(len=8) :: tag
    integer :: st

    tag = merge("dV      ", "FD dV   ", analytic_dv)

    ! Native
    pot%vf => v1
    if (analytic_dv) pot%dvf => dv1
    call inst%init(1.0_wp, 0.0_wp, pot, st)
    call inst%find_profile(prof, st)
    action = inst%find_action(prof, action_pot=action_pot,  &
      action_kin=action_kin)

    ! C interface
    ncalls = 0
    if (analytic_dv) then
      call fortransitions_sfi_new(1.0_c_double, 0.0_c_double,  &
        c_funloc(v1_c), c_funloc(dv1_c), c_null_funptr, 2.0_c_double,  &
        1.0e-3_c_double, 0_c_int, 0.0_c_double, 0_c_int, 0.0_c_double,  &
        h, status)
    else
      call fortransitions_sfi_new(1.0_c_double, 0.0_c_double,  &
        c_funloc(v1_c), c_null_funptr, c_null_funptr, 2.0_c_double,  &
        1.0e-3_c_double, 0_c_int, 0.0_c_double, 0_c_int, 0.0_c_double,  &
        h, status)
    end if
    call check_true(trim(tag)//" sfi_new status", status == status_ok)
    call check_true(trim(tag)//" sfi_new handle", c_associated(h))
    call fortransitions_sfi_params(h, phi_bar, rscale)
    call check_equal(trim(tag)//" phi_bar", real(phi_bar, wp), inst%phi_bar)
    call check_equal(trim(tag)//" rscale", real(rscale, wp), inst%rscale)

    call fortransitions_sfi_find_profile(h, 0_c_int, 0.0_c_double,  &
      1.0e-4_c_double, 1.0e-4_c_double, 1.0e-2_c_double, 500_c_int,  &
      1.0e-4_c_double, 1.0e4_c_double, 0_c_int, 0_c_int, hp, status)
    call check_true(trim(tag)//" find_profile status", status == status_ok)
    call fortransitions_profile_size(hp, n)
    call check_true(trim(tag)//" profile size", n == size(prof%r))
    allocate(r(n), phi(n), dphi(n))
    call fortransitions_profile_arrays(hp, n, r, phi, dphi, rerr)
    call check_true(trim(tag)//" profile arrays",  &
      all(r == prof%r) .and. all(phi == prof%phi)  &
      .and. all(dphi == prof%dphi) .and. rerr == prof%rerr)
    call fortransitions_sfi_find_action(h, n, r, phi, dphi, rerr,  &
      action_c, action_pot_c, action_kin_c)
    call check_equal(trim(tag)//" action", real(action_c, wp), action)
    call check_equal(trim(tag)//" action_pot", real(action_pot_c, wp),  &
      action_pot)
    call check_equal(trim(tag)//" action_kin", real(action_kin_c, wp),  &
      action_kin)
    call check_true(trim(tag)//" callbacks used", ncalls > 0)
    print "(a,a,es20.12,a,i0,a)", trim(tag), " action = ", action_c,  &
      " (", ncalls, " callback evaluations)"
    call fortransitions_profile_free(hp)
    call fortransitions_sfi_free(h)

  end subroutine test_sfi

  subroutine test_sfi_stable()
    !! Error path: swapped minima must give err_stable and a null handle.

    type(c_ptr) :: h
    integer(c_int) :: status

    call fortransitions_sfi_new(0.0_c_double, 1.0_c_double,  &
      c_funloc(v1_c), c_null_funptr, c_null_funptr, 2.0_c_double,  &
      1.0e-3_c_double, 0_c_int, 0.0_c_double, 0_c_int, 0.0_c_double,  &
      h, status)
    call check_true("stable: status err_stable", status == err_stable)
    call check_true("stable: null handle", .not. c_associated(h))
    call fortransitions_sfi_free(h)  ! must be a no-op

  end subroutine test_sfi_stable

  subroutine test_full_tunneling()
    !! full_tunneling through the C interface vs. natively, with the
    !! gradient on the path evaluated point by point (null grad_many_cb)
    !! and in batches (grad_many_cb).

    type(pot2d), target :: pot
    type(full_tunneling_result) :: res
    real(wp) :: path_pts(2, 2)
    integer :: st
    type(c_ptr) :: h
    integer(c_int) :: status
    integer(c_int) :: nprof
    integer(c_int) :: ndim
    integer(c_int) :: num_iters
    integer(c_int) :: converged
    real(c_double) :: action
    real(c_double) :: action_pot
    real(c_double) :: action_kin
    real(c_double) :: fratio
    real(c_double) :: rerr
    real(c_double), allocatable :: r(:)
    real(c_double), allocatable :: phi(:)
    real(c_double), allocatable :: dphi(:)
    real(c_double), allocatable :: phi_path(:, :)
    integer(c_int) :: v_mode
    integer :: batched
    type(c_funptr) :: grad_many_ptr
    character(len=24) :: tag

    path_pts(1, :) = [1.0_wp, 1.0_wp]
    path_pts(2, :) = [0.0_wp, 0.0_wp]

    do v_mode = 0, 1
    do batched = 0, 1
      if (v_mode == 0) then
        call full_tunneling(path_pts, pot, res, st)
      else
        call full_tunneling(path_pts, pot, res, st, v_mode="direct")
      end if
      call check_true("ft native status", st == status_ok)
      if (batched == 1) then
        grad_many_ptr = c_funloc(grad2d_many_c)
      else
        grad_many_ptr = c_null_funptr
      end if
      write (tag, "(a,i0,a,i0)") "v_mode ", v_mode, ", batched ", batched
      print "(a,a)", "--- ", trim(tag)

      ncalls = 0
      nbatch = 0
      call fortransitions_full_tunneling(2_c_int, 2_c_int,  &
        real(path_pts, c_double), c_funloc(v2d_c), c_funloc(grad2d_c),  &
        grad_many_ptr, 20_c_int, 0_c_int, 0_c_int, 100_c_int, v_mode, 2.0_c_double,  &
        500_c_int, 10_c_int, 3_c_int, 0_c_int, 0_c_int, 1.0e-4_c_double,  &
        1.0e-4_c_double, h, status)
      call check_true("ft status", status == status_ok)
      call fortransitions_ft_result_sizes(h, nprof, ndim)
      call check_true("ft sizes",  &
        nprof == size(res%phi, 1) .and. ndim == size(res%phi, 2))
      call fortransitions_ft_result_scalars(h, action, action_pot,  &
        action_kin, fratio, rerr, num_iters, converged)
      call check_equal("ft action", real(action, wp), res%action)
      call check_equal("ft action_pot", real(action_pot, wp), res%action_pot)
      call check_equal("ft action_kin", real(action_kin, wp), res%action_kin)
      call check_equal("ft fratio", real(fratio, wp), res%fratio)
      call check_true("ft num_iters/converged/rerr",  &
        num_iters == res%num_iters  &
        .and. ((converged /= 0) .eqv. res%converged)  &
        .and. rerr == res%profile%rerr)
      allocate(r(nprof), phi(nprof), dphi(nprof), phi_path(nprof, ndim))
      call fortransitions_ft_result_arrays(h, nprof, ndim, r, phi, dphi,  &
        phi_path)
      call check_true("ft arrays",  &
        all(r == res%profile%r) .and. all(phi == res%profile%phi)  &
        .and. all(dphi == res%profile%dphi) .and. all(phi_path == res%phi))
      call check_true("ft callbacks used", ncalls > 0)
      call check_true("ft batched callback used iff given",  &
        (nbatch > 0) .eqv. (batched == 1))
      print "(a,a,es20.12,a,i0,a,i0,a)", trim(tag), ": action = ",  &
        action, " (", ncalls, " V evaluations, ", nbatch,  &
        " batched gradient calls)"
      deallocate(r, phi, dphi, phi_path)
      call fortransitions_ft_result_free(h)
    end do
    end do

  end subroutine test_full_tunneling

  subroutine check_equal(name, val, ref)
    !! Exact equality: both paths run identical code.

    character(len=*), intent(in) :: name
    real(wp), intent(in) :: val
    real(wp), intent(in) :: ref

    if (val /= ref) then
      print "(a,a,a,es22.14,a,es22.14)", "FAIL ", name, ": got ", val,  &
        ", native ", ref
      nfail = nfail + 1
    else
      print "(a,a,a,es22.14)", "ok   ", name, ": ", val
    end if

  end subroutine check_equal

  subroutine check_true(name, cond)

    character(len=*), intent(in) :: name
    logical, intent(in) :: cond

    if (.not. cond) then
      print "(a,a)", "FAIL ", name
      nfail = nfail + 1
    else
      print "(a,a)", "ok   ", name
    end if

  end subroutine check_true

end program c_interface
