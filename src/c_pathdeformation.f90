module fortransitions__c_pathdeformation
!! bind(C) interface to full_tunneling for the Python package
!! pyfortransitions.
!!
!! The result is returned as an opaque handle (a C pointer to a
!! heap-allocated full_tunneling_result), because the size of the profile
!! is only known after the calculation. Query the sizes with
!! fortransitions_ft_result_sizes, copy the data into caller-owned buffers
!! with fortransitions_ft_result_scalars/_arrays, and release the handle
!! with fortransitions_ft_result_free.

  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_c_binding, only : c_double
  use, intrinsic :: iso_c_binding, only : c_ptr
  use, intrinsic :: iso_c_binding, only : c_funptr
  use, intrinsic :: iso_c_binding, only : c_null_ptr
  use, intrinsic :: iso_c_binding, only : c_loc
  use, intrinsic :: iso_c_binding, only : c_f_pointer
  use, intrinsic :: iso_c_binding, only : c_associated
  use, intrinsic :: iso_fortran_env, only : output_unit

  use fortransitions__config, only : wp
  use fortransitions__config, only : status_ok
  use fortransitions__pathdeformation, only : full_tunneling
  use fortransitions__pathdeformation, only : full_tunneling_result
  use fortransitions__c_potentials, only : c_potential_nd
  use fortransitions__c_potentials, only : new_c_potential_nd

  implicit none

  private

  public :: fortransitions_full_tunneling
  public :: fortransitions_ft_result_sizes
  public :: fortransitions_ft_result_scalars
  public :: fortransitions_ft_result_arrays
  public :: fortransitions_ft_result_free

contains

  subroutine fortransitions_full_tunneling(npts, ndim, path_pts, v_cb,  &
      grad_cb, maxiter, verbose, deform_verbose, v_spline_samples,  &
      v_mode, alpha, npoints, nb, kb, fix_start, fix_end, xtol, phitol,  &
      handle, status) bind(C, name="fortransitions_full_tunneling")
    !! Calls full_tunneling with a potential given by C callbacks. All
    !! options are passed explicitly (no optional arguments); logical
    !! flags are 0/1, and v_mode is 0 ("spline") or 1 ("direct"). On
    !! success, `handle` points to the result; on failure, it is null and
    !! `status` holds the error code.

    integer(c_int), value, intent(in) :: npts
    integer(c_int), value, intent(in) :: ndim
    real(c_double), intent(in) :: path_pts(npts, ndim)
      !! Initial path, column-major shape (npts, ndim).
    type(c_funptr), value, intent(in) :: v_cb
    type(c_funptr), value, intent(in) :: grad_cb
    integer(c_int), value, intent(in) :: maxiter
    integer(c_int), value, intent(in) :: verbose
    integer(c_int), value, intent(in) :: deform_verbose
    integer(c_int), value, intent(in) :: v_spline_samples
    integer(c_int), value, intent(in) :: v_mode
    real(c_double), value, intent(in) :: alpha
    integer(c_int), value, intent(in) :: npoints
    integer(c_int), value, intent(in) :: nb
    integer(c_int), value, intent(in) :: kb
    integer(c_int), value, intent(in) :: fix_start
    integer(c_int), value, intent(in) :: fix_end
    real(c_double), value, intent(in) :: xtol
    real(c_double), value, intent(in) :: phitol
    type(c_ptr), intent(out) :: handle
    integer(c_int), intent(out) :: status

    type(c_potential_nd), target :: pot
    type(full_tunneling_result), pointer :: res
    character(len=6) :: mode
    integer :: st

    handle = c_null_ptr
    pot = new_c_potential_nd(v_cb, grad_cb)
    if (v_mode == 1) then
      mode = "direct"
    else
      mode = "spline"
    end if

    allocate(res)
    call full_tunneling(real(path_pts, wp), pot, res, st,  &
      maxiter=int(maxiter), verbose=(verbose /= 0),  &
      deform_verbose=(deform_verbose /= 0),  &
      v_spline_samples=int(v_spline_samples), alpha=real(alpha, wp),  &
      npoints=int(npoints), nb=int(nb), kb=int(kb),  &
      fix_start=(fix_start /= 0), fix_end=(fix_end /= 0),  &
      xtol=real(xtol, wp), phitol=real(phitol, wp), v_mode=trim(mode))
    flush(output_unit)

    status = int(st, c_int)
    if (st == status_ok) then
      handle = c_loc(res)
    else
      deallocate(res)
    end if

  end subroutine fortransitions_full_tunneling

  subroutine fortransitions_ft_result_sizes(handle, nprof, ndim)  &
      bind(C, name="fortransitions_ft_result_sizes")
    !! Number of profile points (= number of path points) and field
    !! dimensions of a result.

    type(c_ptr), value, intent(in) :: handle
    integer(c_int), intent(out) :: nprof
    integer(c_int), intent(out) :: ndim

    type(full_tunneling_result), pointer :: res

    call c_f_pointer(handle, res)
    nprof = int(size(res%phi, 1), c_int)
    ndim = int(size(res%phi, 2), c_int)

  end subroutine fortransitions_ft_result_sizes

  subroutine fortransitions_ft_result_scalars(handle, action, action_pot,  &
      action_kin, fratio, rerr, num_iters, converged)  &
      bind(C, name="fortransitions_ft_result_scalars")

    type(c_ptr), value, intent(in) :: handle
    real(c_double), intent(out) :: action
    real(c_double), intent(out) :: action_pot
    real(c_double), intent(out) :: action_kin
    real(c_double), intent(out) :: fratio
    real(c_double), intent(out) :: rerr
      !! Negative if no step size error occurred (None in Python).
    integer(c_int), intent(out) :: num_iters
    integer(c_int), intent(out) :: converged
      !! 1 if the deformation converged, else 0.

    type(full_tunneling_result), pointer :: res

    call c_f_pointer(handle, res)
    action = real(res%action, c_double)
    action_pot = real(res%action_pot, c_double)
    action_kin = real(res%action_kin, c_double)
    fratio = real(res%fratio, c_double)
    rerr = real(res%profile%rerr, c_double)
    num_iters = int(res%num_iters, c_int)
    converged = merge(1_c_int, 0_c_int, res%converged)

  end subroutine fortransitions_ft_result_scalars

  subroutine fortransitions_ft_result_arrays(handle, nprof, ndim, r, phi,  &
      dphi, phi_path) bind(C, name="fortransitions_ft_result_arrays")
    !! Copy the profile (r, phi, dphi) and the points of the final path
    !! (column-major shape (nprof, ndim)) into caller-owned buffers whose
    !! sizes must match fortransitions_ft_result_sizes.

    type(c_ptr), value, intent(in) :: handle
    integer(c_int), value, intent(in) :: nprof
    integer(c_int), value, intent(in) :: ndim
    real(c_double), intent(out) :: r(nprof)
    real(c_double), intent(out) :: phi(nprof)
    real(c_double), intent(out) :: dphi(nprof)
    real(c_double), intent(out) :: phi_path(nprof, ndim)

    type(full_tunneling_result), pointer :: res

    call c_f_pointer(handle, res)
    r = real(res%profile%r, c_double)
    phi = real(res%profile%phi, c_double)
    dphi = real(res%profile%dphi, c_double)
    phi_path = real(res%phi, c_double)

  end subroutine fortransitions_ft_result_arrays

  subroutine fortransitions_ft_result_free(handle)  &
      bind(C, name="fortransitions_ft_result_free")

    type(c_ptr), value, intent(in) :: handle

    type(full_tunneling_result), pointer :: res

    if (.not. c_associated(handle)) return
    call c_f_pointer(handle, res)
    deallocate(res)

  end subroutine fortransitions_ft_result_free

end module fortransitions__c_pathdeformation
