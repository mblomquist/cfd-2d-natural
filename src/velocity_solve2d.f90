! velocity_solve2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-02-26 (YYYY-MM-DD)
!
! This subroutine updates the source terms for the solution of the momentum
! equations in the SIMPLER algorithm.
!
! These definitions are defined for a 2D natural convection problem that 
! simulates a vertical plate at i = 0.

subroutine velocity_solve2d

  ! Include variable header
  include "var2d.dec"

  ! Update source terms
  call velocity_source2d("u")

  ! Solve u-velocity equation
  call solver2d_bicgstab(As_u(1:m-1,1:n-1), Aw_u(1:m-1,1:n-1), Ap_u(1:m-1,1:n-1), Ae_u(1:m-1,1:n-1), An_u(1:m-1,1:n-1), b_u(1:m-1,1:n-1), u_star(1:m-1,1:n-1), m-1, n-1, tol, maxit)

  ! Outlet outlet boundaries
  u_star(m,1:n-1) = u_star(m-1,1:n-1)
  u_star(2:m,n) = u_star(2:m,n-1)

  ! Update source terms
  call velocity_source2d("v")

  ! Solve v-velocity equation
  call solver2d_bicgstab(As_v(1:m-1,1:n-1), Aw_v(1:m-1,1:n-1), Ap_v(1:m-1,1:n-1), Ae_v(1:m-1,1:n-1), An_v(1:m-1,1:n-1), b_v(1:m-1,1:n-1), v_star(1:m-1,1:n-1), m-1, n-2, tol, maxit)

  ! Outlet outlet boundaries
  v_star(m,1:n-1) = v_star(m-1,1:n-1)
  v_star(2:m,n) = v_star(2:m,n-1)

  return

end subroutine velocity_solve2d
