! temperature_solve2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-06-27 (YYYY-MM-DD)
!
! This subroutine solves the temperature equation of the SIMPLER algorithm
!
subroutine temperature_solve2d

  ! Include variable header
  include "var2d.dec"

  ! Define internal variables
  integer :: i, j

  ! Update source terms
  call temperature_source2d

  ! Print coefficients
  !print *, "Aw_T:", Aw_T
  !print *, "Ae_T:", Ae_T
  !print *, "An_T:", As_T
  !print *, "As_T:", An_T
  !print *, "Ap_T:", Ap_T
  !print *, "b_T:", b_T

  ! Update source terms
  do i = 1, m-1
    do j = 1, n-1

      Ap_T(i,j) = Ap_T(i,j)/alpha_t
      b_T(i,j) = Su_T(i,j)*dx*dy+(1.0-alpha_t)*Ap_T(i,j)*T(i,j)

    end do
  end do

  ! Solve velocity Equations
  if (solver .eq. "b") then
    call solver2d_bicgstab2(As_T, Aw_T, Ap_T, Ae_T, An_T, b_T, T, m-1, n-1, solver_tol, maxit)
  else
    call solver2d_tdma(Aw_T, Ae_T, As_T, An_T, Ap_T, b_T, T, m-1, n-1, solver_tol, maxit)
  end if

  return

end subroutine temperature_solve2d
