module cosmotransitions
!! Top-level module of the cosmotransitions package: re-exports the public
!! API of the individual modules.
!!
!! This is a Fortran port of the Python package CosmoTransitions
!! (Comput. Phys. Commun. 183 (2012) 2006, arXiv:1109.4189). So far the
!! functionality of pathDeformation.fullTunneling has been ported,
!! together with everything it needs (tunneling1D.SingleFieldInstanton,
!! pathDeformation.SplinePath, pathDeformation.Deformation_Spline and the
!! required helper functions).

  use cosmotransitions__config, only : wp
  use cosmotransitions__config, only : status_ok
  use cosmotransitions__config, only : err_integration
  use cosmotransitions__config, only : err_no_barrier
  use cosmotransitions__config, only : err_stable
  use cosmotransitions__config, only : err_deformation
  use cosmotransitions__config, only : err_numerical
  use cosmotransitions__potentials, only : potential_1d
  use cosmotransitions__potentials, only : potential_nd
  use cosmotransitions__potentials, only : pot1d_func
  use cosmotransitions__tunneling1d, only : single_field_instanton
  use cosmotransitions__tunneling1d, only : profile1d
  use cosmotransitions__splinepath, only : spline_path
  use cosmotransitions__deformation, only : deformation_spline
  use cosmotransitions__pathdeformation, only : full_tunneling
  use cosmotransitions__pathdeformation, only : full_tunneling_result

  implicit none

  private

  ! Kind and status codes
  public :: wp
  public :: status_ok
  public :: err_integration
  public :: err_no_barrier
  public :: err_stable
  public :: err_deformation
  public :: err_numerical
  ! Potential interfaces
  public :: potential_1d
  public :: potential_nd
  public :: pot1d_func
  ! One-dimensional tunneling
  public :: single_field_instanton
  public :: profile1d
  ! Path deformation
  public :: spline_path
  public :: deformation_spline
  public :: full_tunneling
  public :: full_tunneling_result

end module cosmotransitions
