! velocity_solve2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-03-08 (YYYY-MM-DD)
!
! This subroutine updates the source terms for the solution of the momentum
! equations in the SIMPLER algorithm.
!
subroutine velocity_solve2d

  ! Include variable header
  include "var2d.dec"

  ! Update source terms
  call velocity_source2d("u")

  ! Solve u-velocity equation
  call solver2d_bicgstab(As_u, Aw_u, Ap_u, Ae_u, An_u, b_u, u_star, m, n-1, tol, maxit)

  ! Update source terms
  call velocity_source2d("v")

  ! Solve v-velocity equation
  call solver2d_bicgstab(As_v, Aw_v, Ap_v, Ae_v, An_v, b_v, v_star, m, n, tol, maxit)

  return

end subroutine velocity_solve2d
