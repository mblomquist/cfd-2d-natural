! velocity_solve2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-04-28 (YYYY-MM-DD)
!
! This subroutine updates the source terms for the solution of the momentum
! equations in the SIMPLER algorithm.
!
subroutine velocity_solve2d

  ! Include variable header
  include "var2d.dec"

  ! Update source terms
  call velocity_source2d("u")

  ! Print coefficients
  !print *, "Aw_u:", Aw_u
  !print *, "Ae_u:", Ae_u
  !print *, "An_u:", As_u
  !print *, "As_u:", An_u
  !print *, "Ap_u:", Ap_u
  !print *, "b_u:", b_u

  ! Solve u-velocity equation
  call solver2d_bicgstab2(As_u, Aw_u, Ap_u, Ae_u, An_u, b_u, u_star, m, n-1, solver_tol, maxit)
  !call solver2d_tdma(Aw_u, Ae_u, As_u, An_u, Ap_u, b_u, u_star, m, n-1, solver_tol, maxit)

  ! Update source terms
  call velocity_source2d("v")

  ! Print coefficients
  !print *, "Aw_v:", Aw_v
  !print *, "Ae_v:", Ae_v
  !print *, "An_v:", As_v
  !print *, "As_v:", An_v
  !print *, "Ap_v:", Ap_v
  !print *, "b_v:", b_v

  ! Solve v-velocity equation
  call solver2d_bicgstab2(As_v, Aw_v, Ap_v, Ae_v, An_v, b_v, v_star, m-1, n, solver_tol, maxit)
  !call solver2d_tdma(Aw_v, Ae_v, As_v, An_v, Ap_v, b_v, v_star, m-1, n, solver_tol, maxit)

  return

end subroutine velocity_solve2d
