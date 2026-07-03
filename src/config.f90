module cosmotransitions__config
!! Working precision and status codes shared by all cosmotransitions modules.

  use, intrinsic :: iso_fortran_env, only : real64

  implicit none

  private

  integer, parameter, public :: wp = real64
    !! Working precision. Must match the `wp` of bspline-fortran (real64 by
    !! default).

  real(wp), parameter, public :: pi = 3.14159265358979323846_wp

  ! Status codes. `status_ok` (0) means success; everything else mirrors one
  ! of the exceptions raised by the Python package.
  integer, parameter, public :: status_ok = 0
  integer, parameter, public :: err_integration = 1
    !! helper_functions.IntegrationError
  integer, parameter, public :: err_no_barrier = 2
    !! tunneling1D.PotentialError ("no barrier")
  integer, parameter, public :: err_stable = 3
    !! tunneling1D.PotentialError ("stable, not metastable")
  integer, parameter, public :: err_deformation = 4
    !! pathDeformation.DeformationError
  integer, parameter, public :: err_numerical = 5
    !! Generic numerical failure (root finding, linear solve, spline setup...)

end module cosmotransitions__config
