! temperature_solve2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-03-09 (YYYY-MM-DD)
!
! This subroutine solves the temperature equation of the SIMPLER algorithm
!
subroutine temperature_solve2d

  ! Include variable header
  include "var2d.dec"

  ! Update source terms
  call temperature_source2d

  ! Solve velocity Equations
  call solver2d_bicgstab(As_T, Aw_T, Ap_T, Ae_T, An_T, b_T, T, m-1, n-1, tol, maxit)

  return

end subroutine temperature_solve2d
