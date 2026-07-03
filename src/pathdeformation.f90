module cosmotransitions__pathdeformation
!! Port of pathDeformation.fullTunneling: find the instanton solution in
!! multiple field dimensions by looping over
!!   1. fitting a spline to the path,
!!   2. solving the one-dimensional tunneling problem along the path,
!!   3. deforming the path to satisfy the transverse equations of motion,
!!   4. checking for convergence.

  use cosmotransitions__config, only : wp
  use cosmotransitions__config, only : status_ok
  use cosmotransitions__config, only : err_deformation
  use cosmotransitions__potentials, only : potential_nd
  use cosmotransitions__tunneling1d, only : single_field_instanton
  use cosmotransitions__tunneling1d, only : profile1d
  use cosmotransitions__splinepath, only : spline_path
  use cosmotransitions__deformation, only : deformation_spline

  implicit none

  private

  public :: full_tunneling_result
  public :: full_tunneling

  type :: full_tunneling_result
    !! The namedtuple "fullTunneling_rval" of the Python package.
    type(profile1d) :: profile
      !! The one-dimensional bubble profile along the path. `profile%phi`
      !! is the distance along the path, in one-to-one correspondence with
      !! the rows of `phi`.
    real(wp), allocatable :: phi(:, :)
      !! Points along the final deformed path, shape (npoints, ndim).
    real(wp) :: action = 0.0_wp
      !! Euclidean action of the instanton.
    real(wp) :: fratio = 0.0_wp
      !! Largest transverse force on the final path relative to the
      !! largest potential gradient (zero for a perfect solution).
    integer :: num_iters = 0
      !! Number of tunneling/deformation iterations performed.
    logical :: converged = .false.
      !! True if the deformation converged before maxiter was reached.
  end type full_tunneling_result

contains

  subroutine full_tunneling(path_pts, pot, res, status, maxiter, verbose,  &
      v_spline_samples, alpha, npoints, nb, kb, fix_start, fix_end,  &
      deform_verbose)
    !! Calculate the instanton solution in multiple field dimensions.
    !! Port of pathDeformation.fullTunneling (with the default
    !! tunneling_class = SingleFieldInstanton and
    !! deformation_class = Deformation_Spline).

    real(wp), intent(in) :: path_pts(:, :)
      !! Initial guess for the tunneling path, shape (num_points, ndim).
      !! The first point should be at (or near) the stable minimum (the
      !! minimum to which the field tunnels), and the last point at the
      !! metastable minimum.
    class(potential_nd), intent(inout), target :: pot
      !! The potential V and its gradient.
    type(full_tunneling_result), intent(out) :: res
    integer, intent(out) :: status
    integer, intent(in), optional :: maxiter
      !! Maximum number of allowed deformation/tunneling iterations
      !! (default 20).
    logical, intent(in), optional :: verbose
      !! If true, print a message at the start of each iteration
      !! (default false).
    integer, intent(in), optional :: v_spline_samples
      !! Number of potential samples for the path spline (default 100).
    real(wp), intent(in), optional :: alpha
      !! Friction coefficient for the 1d tunneling (default 2).
    integer, intent(in), optional :: npoints
      !! Number of points in the 1d profile (default 500).
    integer, intent(in), optional :: nb
      !! Number of deformation basis splines (default 10).
    integer, intent(in), optional :: kb
      !! Order of the deformation basis splines (default 3).
    logical, intent(in), optional :: fix_start
    logical, intent(in), optional :: fix_end
      !! Hold the first/last path point fixed during deformation
      !! (default false).
    logical, intent(in), optional :: deform_verbose
      !! Verbosity of the deformation loop (default true, like Python).

    integer :: maxiter_
    logical :: verbose_
    type(spline_path), target :: path
    type(single_field_instanton) :: tobj
    type(deformation_spline) :: deform_obj
    type(profile1d) :: prof
    real(wp), allocatable :: pts(:, :)
    real(wp), allocatable :: phi(:)
    real(wp), allocatable :: dphi(:)
    real(wp), allocatable :: f(:, :)
    real(wp), allocatable :: dv(:, :)
    real(wp) :: f_max
    real(wp) :: dv_max
    integer :: num_iter
    integer :: i
    logical :: converged
    integer :: st

    maxiter_ = 20
    verbose_ = .false.
    if (present(maxiter)) maxiter_ = maxiter
    if (present(verbose)) verbose_ = verbose
    if (maxiter_ <= 0) then
      error stop "full_tunneling: maxiter must be positive"
    end if

    status = status_ok
    res%converged = .false.
    pts = path_pts

    do num_iter = 1, maxiter_
      res%num_iters = num_iter
      if (verbose_) print "(a,i0)", "Starting tunneling step ", num_iter
      ! 1. Fit the spline to the path.
      call path%init(pts, pot, st, v_spline_samples=v_spline_samples,  &
        extend_to_minima=.true., reeval_distances=.true.)
      if (st /= status_ok) then
        status = st
        return
      end if
      ! 2. Do 1d tunneling along the path.
      call tobj%init(0.0_wp, path%length, path, st, alpha=alpha)
      if (st /= status_ok) then
        status = st
        return
      end if
      call tobj%find_profile(prof, st, npoints=npoints)
      if (st /= status_ok) then
        status = st
        return
      end if
      call tobj%evenly_spaced_phi(prof%phi, prof%dphi, size(prof%phi),  &
        .false., phi, dphi)
      dphi(1) = 0.0_wp
      dphi(size(dphi)) = 0.0_wp  ! enforce this
      ! 3. Deform the path.
      pts = path%pts_many(phi)  ! multi-dimensional points
      call deform_obj%init(pts, dphi, pot, st, nb=nb, kb=kb,  &
        fix_start=fix_start, fix_end=fix_end)
      if (st /= status_ok) then
        status = st
        return
      end if
      call deform_obj%deform_path(converged, st, verbose=deform_verbose)
      if (st == err_deformation) then
        ! Mirrors the Python code, which catches DeformationError, prints
        ! it and continues with converged = False.
        converged = .false.
      else if (st /= status_ok) then
        status = st
        return
      end if
      pts = deform_obj%phi
      ! 4. Check convergence. If the deformation converged after one step,
      ! then assume that the path is a good solution.
      if (converged .and. deform_obj%num_steps < 2) then
        res%converged = .true.
        exit
      end if
    end do
    if (.not. res%converged .and. verbose_) then
      print "(a)", "Reached maxiter in full_tunneling. No convergence."
    end if

    ! Calculate the ratio of the maximum perpendicular force to the
    ! maximum gradient. Make sure to go back a step and use the forces on
    ! the path, not the most recently deformed path.
    call deform_obj%init(pts, dphi, pot, st, nb=nb, kb=kb,  &
      fix_start=fix_start, fix_end=fix_end)
    if (st /= status_ok) then
      status = st
      return
    end if
    call deform_obj%forces(f, dv)
    f_max = 0.0_wp
    dv_max = 0.0_wp
    do i = 1, size(f, 1)
      f_max = max(f_max, norm2(f(i, :)))
      dv_max = max(dv_max, norm2(dv(i, :)))
    end do
    res%fratio = f_max/dv_max

    ! Assemble the output.
    res%profile = prof
    res%phi = path%pts_many(prof%phi)
    res%action = tobj%find_action(prof)

  end subroutine full_tunneling

end module cosmotransitions__pathdeformation
