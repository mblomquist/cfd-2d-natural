! temperature_solve2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-02-26 (YYYY-MM-DD)
!
! This subroutine solves the temperature equation of the SIMPLER algorithm
!
! These definitions are defined for a 2D natural convection problem that 
! simulates a vertical plate at i = 0.

subroutine temperature_solve2d

  ! Include variable header
  include "var2d.dec"

  ! Update source terms
  call temperature_source2d

  ! Solve velocity Equations
  call solver2d_bicgstab(As_T(1:m-1,1:n-1), Aw_T(1:m-1,1:n-1), Ap_T(1:m-1,1:n-1), Ae_T(1:m-1,1:n-1), An_T(1:m-1,1:n-1), b_T(1:m-1,1:n-1), T(1:m-1,1:n-1), m-1, n-1, tol, maxit)

  ! Update boundaries
  T(m, 2:n) = T(m-1, 2:n-1)
  T(2:m-1, n) = T(2:m-1, n)

  return

end subroutine temperature_solve2d
