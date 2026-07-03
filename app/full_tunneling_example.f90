module full_tunneling_example__potential
!! The two-dimensional example potential from examples/fullTunneling.py of
!! the Python package.

  use cosmotransitions, only : wp
  use cosmotransitions, only : potential_nd

  implicit none

  private

  public :: example_potential

  type, extends(potential_nd) :: example_potential
    real(wp) :: c = 5.0_wp
    real(wp) :: fx = 10.0_wp
    real(wp) :: fy = 10.0_wp
  contains
    procedure :: v => example_v
    procedure :: grad => example_grad
  end type example_potential

contains

  function example_v(self, x) result(y)
    class(example_potential), intent(inout) :: self
    real(wp), intent(in) :: x(:)
    real(wp) :: y
    real(wp) :: r1
    real(wp) :: r2
    real(wp) :: r3
    r1 = x(1)*x(1) + self%c*x(2)*x(2)
    r2 = self%c*(x(1) - 1.0_wp)**2 + (x(2) - 1.0_wp)**2
    r3 = self%fx*(0.25_wp*x(1)**4 - x(1)**3/3.0_wp)
    r3 = r3 + self%fy*(0.25_wp*x(2)**4 - x(2)**3/3.0_wp)
    y = r1*r2 + r3
  end function example_v

  function example_grad(self, x) result(dv)
    class(example_potential), intent(inout) :: self
    real(wp), intent(in) :: x(:)
    real(wp) :: dv(size(x))
    real(wp) :: r1
    real(wp) :: r2
    r1 = x(1)*x(1) + self%c*x(2)*x(2)
    r2 = self%c*(x(1) - 1.0_wp)**2 + (x(2) - 1.0_wp)**2
    dv(1) = r1*2.0_wp*self%c*(x(1) - 1.0_wp) + 2.0_wp*x(1)*r2  &
      + self%fx*x(1)*x(1)*(x(1) - 1.0_wp)
    dv(2) = r1*2.0_wp*(x(2) - 1.0_wp) + 2.0_wp*self%c*x(2)*r2  &
      + self%fy*x(2)*x(2)*(x(2) - 1.0_wp)
  end function example_grad

end module full_tunneling_example__potential

program full_tunneling_example
!! Reproduces examples/fullTunneling.py of the Python package: calculates
!! a thin-walled and a thick-walled two-field instanton and writes the
!! resulting paths and one-dimensional profiles to CSV files.

  use cosmotransitions, only : wp
  use cosmotransitions, only : status_ok
  use cosmotransitions, only : full_tunneling
  use cosmotransitions, only : full_tunneling_result

  use full_tunneling_example__potential, only : example_potential

  implicit none

  type(example_potential), target :: pot
  type(full_tunneling_result) :: res
  real(wp) :: path_pts(2, 2)
  integer :: status

  path_pts(1, :) = [1.0_wp, 1.0_wp]
  path_pts(2, :) = [0.0_wp, 0.0_wp]

  print "(a)", "=== Thin-walled instanton ==="
  pot = example_potential(c=5.0_wp, fx=0.0_wp, fy=2.0_wp)
  call full_tunneling(path_pts, pot, res, status, verbose=.true.)
  call report(res, status, "thin")

  print "(a)", ""
  print "(a)", "=== Thick-walled instanton ==="
  pot = example_potential(c=5.0_wp, fx=0.0_wp, fy=80.0_wp)
  call full_tunneling(path_pts, pot, res, status, verbose=.true.)
  call report(res, status, "thick")

contains

  subroutine report(res, status, label)

    type(full_tunneling_result), intent(in) :: res
    integer, intent(in) :: status
    character(len=*), intent(in) :: label

    if (status /= status_ok) then
      print "(a,i0)", "full_tunneling failed with status ", status
      return
    end if
    print "(a,es16.8)", "action = ", res%action
    print "(a,es10.2)", "fRatio = ", res%fratio
    print "(a,i0)", "iterations = ", res%num_iters
    call write_csv(res, label)

  end subroutine report

  subroutine write_csv(res, label)
    !! Write the tunneling path and the 1d profile to CSV files.

    use csv_module, only : csv_file

    type(full_tunneling_result), intent(in) :: res
    character(len=*), intent(in) :: label

    type(csv_file) :: f
    logical :: ok
    integer :: i
    integer :: n

    n = size(res%profile%r)
    call f%initialize()
    call f%open("instanton_"//label//".csv", n_cols=5, status_ok=ok)
    if (.not. ok) then
      print "(a)", "Could not open CSV file for writing."
      return
    end if
    call f%add(["r   ", "phi ", "dphi", "x   ", "y   "])
    call f%next_row()
    do i = 1, n
      call f%add(res%profile%r(i))
      call f%add(res%profile%phi(i))
      call f%add(res%profile%dphi(i))
      call f%add(res%phi(i, 1))
      call f%add(res%phi(i, 2))
      call f%next_row()
    end do
    call f%close(ok)
    print "(a)", "Wrote instanton_"//label//".csv"

  end subroutine write_csv

end program full_tunneling_example
