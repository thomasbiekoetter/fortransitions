module fortransitions__c_tunneling1d
!! bind(C) interface to single_field_instanton for the Python package
!! pyfortransitions.
!!
!! An instanton is an opaque handle (C pointer) to a heap-allocated
!! container holding the potential and the single_field_instanton object.
!! find_profile returns a second handle to the profile, because its size
!! is only known after the calculation. Both handles must be released
!! with their *_free routines.

  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_c_binding, only : c_double
  use, intrinsic :: iso_c_binding, only : c_ptr
  use, intrinsic :: iso_c_binding, only : c_funptr
  use, intrinsic :: iso_c_binding, only : c_null_ptr
  use, intrinsic :: iso_c_binding, only : c_loc
  use, intrinsic :: iso_c_binding, only : c_f_pointer
  use, intrinsic :: iso_c_binding, only : c_associated

  use fortransitions__config, only : wp
  use fortransitions__config, only : status_ok
  use fortransitions__tunneling1d, only : single_field_instanton
  use fortransitions__tunneling1d, only : profile1d
  use fortransitions__c_potentials, only : c_potential_1d
  use fortransitions__c_potentials, only : new_c_potential_1d

  implicit none

  private

  public :: fortransitions_sfi_new
  public :: fortransitions_sfi_params
  public :: fortransitions_sfi_find_profile
  public :: fortransitions_sfi_find_action
  public :: fortransitions_sfi_free
  public :: fortransitions_profile_size
  public :: fortransitions_profile_arrays
  public :: fortransitions_profile_free

  type :: sfi_container
    !! Keeps the potential alive next to the instanton, which only holds
    !! a pointer to it.
    type(c_potential_1d) :: pot
    type(single_field_instanton) :: inst
  end type sfi_container

contains

  subroutine fortransitions_sfi_new(phi_absmin, phi_metamin, v_cb, dv_cb,  &
      d2v_cb, alpha, phi_eps_rel, has_phi_bar, phi_bar, has_rscale,  &
      rscale, handle, status) bind(C, name="fortransitions_sfi_new")
    !! Create and initialize a single_field_instanton. `dv_cb` and
    !! `d2v_cb` may be null pointers, which selects finite-difference
    !! derivatives. `phi_bar` and `rscale` are only used if the
    !! corresponding has_* flag is nonzero; otherwise they are found
    !! automatically. On failure, `handle` is null.

    real(c_double), value, intent(in) :: phi_absmin
    real(c_double), value, intent(in) :: phi_metamin
    type(c_funptr), value, intent(in) :: v_cb
    type(c_funptr), value, intent(in) :: dv_cb
    type(c_funptr), value, intent(in) :: d2v_cb
    real(c_double), value, intent(in) :: alpha
    real(c_double), value, intent(in) :: phi_eps_rel
    integer(c_int), value, intent(in) :: has_phi_bar
    real(c_double), value, intent(in) :: phi_bar
    integer(c_int), value, intent(in) :: has_rscale
    real(c_double), value, intent(in) :: rscale
    type(c_ptr), intent(out) :: handle
    integer(c_int), intent(out) :: status

    type(sfi_container), pointer :: c
    integer :: st

    handle = c_null_ptr
    allocate(c)
    c%pot = new_c_potential_1d(v_cb, dv_cb, d2v_cb)

    if (has_phi_bar /= 0 .and. has_rscale /= 0) then
      call c%inst%init(real(phi_absmin, wp), real(phi_metamin, wp), c%pot,  &
        st, alpha=real(alpha, wp), phi_eps_rel=real(phi_eps_rel, wp),  &
        phi_bar=real(phi_bar, wp), rscale=real(rscale, wp))
    else if (has_phi_bar /= 0) then
      call c%inst%init(real(phi_absmin, wp), real(phi_metamin, wp), c%pot,  &
        st, alpha=real(alpha, wp), phi_eps_rel=real(phi_eps_rel, wp),  &
        phi_bar=real(phi_bar, wp))
    else if (has_rscale /= 0) then
      call c%inst%init(real(phi_absmin, wp), real(phi_metamin, wp), c%pot,  &
        st, alpha=real(alpha, wp), phi_eps_rel=real(phi_eps_rel, wp),  &
        rscale=real(rscale, wp))
    else
      call c%inst%init(real(phi_absmin, wp), real(phi_metamin, wp), c%pot,  &
        st, alpha=real(alpha, wp), phi_eps_rel=real(phi_eps_rel, wp))
    end if

    status = int(st, c_int)
    if (st == status_ok) then
      handle = c_loc(c)
    else
      deallocate(c)
    end if

  end subroutine fortransitions_sfi_new

  subroutine fortransitions_sfi_params(handle, phi_bar, rscale)  &
      bind(C, name="fortransitions_sfi_params")
    !! The barrier edge and radial scale (given or found by init).

    type(c_ptr), value, intent(in) :: handle
    real(c_double), intent(out) :: phi_bar
    real(c_double), intent(out) :: rscale

    type(sfi_container), pointer :: c

    call c_f_pointer(handle, c)
    phi_bar = real(c%inst%phi_bar, c_double)
    rscale = real(c%inst%rscale, c_double)

  end subroutine fortransitions_sfi_params

  subroutine fortransitions_sfi_find_profile(handle, has_xguess, xguess,  &
      xtol, phitol, thin_cutoff, npoints, rmin_rel, rmax_rel,  &
      has_max_interior, max_interior_pts, prof_handle, status)  &
      bind(C, name="fortransitions_sfi_find_profile")
    !! Calculate the bubble profile. `xguess` and `max_interior_pts` are
    !! only used if the corresponding has_* flag is nonzero. On success,
    !! `prof_handle` points to the profile; on failure, it is null.

    type(c_ptr), value, intent(in) :: handle
    integer(c_int), value, intent(in) :: has_xguess
    real(c_double), value, intent(in) :: xguess
    real(c_double), value, intent(in) :: xtol
    real(c_double), value, intent(in) :: phitol
    real(c_double), value, intent(in) :: thin_cutoff
    integer(c_int), value, intent(in) :: npoints
    real(c_double), value, intent(in) :: rmin_rel
    real(c_double), value, intent(in) :: rmax_rel
    integer(c_int), value, intent(in) :: has_max_interior
    integer(c_int), value, intent(in) :: max_interior_pts
    type(c_ptr), intent(out) :: prof_handle
    integer(c_int), intent(out) :: status

    type(sfi_container), pointer :: c
    type(profile1d), pointer :: prof
    integer :: st

    prof_handle = c_null_ptr
    call c_f_pointer(handle, c)
    allocate(prof)

    if (has_xguess /= 0 .and. has_max_interior /= 0) then
      call c%inst%find_profile(prof, st, xguess=real(xguess, wp),  &
        xtol=real(xtol, wp), phitol=real(phitol, wp),  &
        thin_cutoff=real(thin_cutoff, wp), npoints=int(npoints),  &
        rmin_rel=real(rmin_rel, wp), rmax_rel=real(rmax_rel, wp),  &
        max_interior_pts=int(max_interior_pts))
    else if (has_xguess /= 0) then
      call c%inst%find_profile(prof, st, xguess=real(xguess, wp),  &
        xtol=real(xtol, wp), phitol=real(phitol, wp),  &
        thin_cutoff=real(thin_cutoff, wp), npoints=int(npoints),  &
        rmin_rel=real(rmin_rel, wp), rmax_rel=real(rmax_rel, wp))
    else if (has_max_interior /= 0) then
      call c%inst%find_profile(prof, st,  &
        xtol=real(xtol, wp), phitol=real(phitol, wp),  &
        thin_cutoff=real(thin_cutoff, wp), npoints=int(npoints),  &
        rmin_rel=real(rmin_rel, wp), rmax_rel=real(rmax_rel, wp),  &
        max_interior_pts=int(max_interior_pts))
    else
      call c%inst%find_profile(prof, st,  &
        xtol=real(xtol, wp), phitol=real(phitol, wp),  &
        thin_cutoff=real(thin_cutoff, wp), npoints=int(npoints),  &
        rmin_rel=real(rmin_rel, wp), rmax_rel=real(rmax_rel, wp))
    end if

    status = int(st, c_int)
    if (st == status_ok) then
      prof_handle = c_loc(prof)
    else
      deallocate(prof)
    end if

  end subroutine fortransitions_sfi_find_profile

  subroutine fortransitions_sfi_find_action(handle, n, r, phi, dphi, rerr,  &
      action, action_pot, action_kin)  &
      bind(C, name="fortransitions_sfi_find_action")
    !! The Euclidean action of a profile given as arrays (not necessarily
    !! one returned by find_profile), with the Derrick-rescaled estimates
    !! from the potential and kinetic terms alone.

    type(c_ptr), value, intent(in) :: handle
    integer(c_int), value, intent(in) :: n
    real(c_double), intent(in) :: r(n)
    real(c_double), intent(in) :: phi(n)
    real(c_double), intent(in) :: dphi(n)
    real(c_double), value, intent(in) :: rerr
    real(c_double), intent(out) :: action
    real(c_double), intent(out) :: action_pot
    real(c_double), intent(out) :: action_kin

    type(sfi_container), pointer :: c
    type(profile1d) :: prof
    real(wp) :: s_pot
    real(wp) :: s_kin

    call c_f_pointer(handle, c)
    prof%r = real(r, wp)
    prof%phi = real(phi, wp)
    prof%dphi = real(dphi, wp)
    prof%rerr = real(rerr, wp)
    action = real(c%inst%find_action(prof, action_pot=s_pot,  &
      action_kin=s_kin), c_double)
    action_pot = real(s_pot, c_double)
    action_kin = real(s_kin, c_double)

  end subroutine fortransitions_sfi_find_action

  subroutine fortransitions_sfi_free(handle)  &
      bind(C, name="fortransitions_sfi_free")

    type(c_ptr), value, intent(in) :: handle

    type(sfi_container), pointer :: c

    if (.not. c_associated(handle)) return
    call c_f_pointer(handle, c)
    deallocate(c)

  end subroutine fortransitions_sfi_free

  subroutine fortransitions_profile_size(prof_handle, n)  &
      bind(C, name="fortransitions_profile_size")

    type(c_ptr), value, intent(in) :: prof_handle
    integer(c_int), intent(out) :: n

    type(profile1d), pointer :: prof

    call c_f_pointer(prof_handle, prof)
    n = int(size(prof%r), c_int)

  end subroutine fortransitions_profile_size

  subroutine fortransitions_profile_arrays(prof_handle, n, r, phi, dphi,  &
      rerr) bind(C, name="fortransitions_profile_arrays")
    !! Copy a profile into caller-owned buffers of size
    !! fortransitions_profile_size. `rerr` is negative if no step size
    !! error occurred (None in Python).

    type(c_ptr), value, intent(in) :: prof_handle
    integer(c_int), value, intent(in) :: n
    real(c_double), intent(out) :: r(n)
    real(c_double), intent(out) :: phi(n)
    real(c_double), intent(out) :: dphi(n)
    real(c_double), intent(out) :: rerr

    type(profile1d), pointer :: prof

    call c_f_pointer(prof_handle, prof)
    r = real(prof%r, c_double)
    phi = real(prof%phi, c_double)
    dphi = real(prof%dphi, c_double)
    rerr = real(prof%rerr, c_double)

  end subroutine fortransitions_profile_arrays

  subroutine fortransitions_profile_free(prof_handle)  &
      bind(C, name="fortransitions_profile_free")

    type(c_ptr), value, intent(in) :: prof_handle

    type(profile1d), pointer :: prof

    if (.not. c_associated(prof_handle)) return
    call c_f_pointer(prof_handle, prof)
    deallocate(prof)

  end subroutine fortransitions_profile_free

end module fortransitions__c_tunneling1d
