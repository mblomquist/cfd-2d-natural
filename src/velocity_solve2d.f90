! velocity_solve2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-06-27 (YYYY-MM-DD)
!
! This subroutine updates the source terms for the solution of the momentum
! equations in the SIMPLER algorithm.
!
subroutine velocity_solve2d

  implicit none

  ! Include variable header
  include "var2d.dec"

  ! Define internal variables
  integer :: i, j

  ! Update source terms :: u
  do i = 2,m-1
    do j = 2,n-1

      Ap_u(i,j) = Ap_u(i,j)/alpha_v
      b_u(i,j) = b_v(i,j)+dy*(P_star(i-1,j)-P_star(i,j))+(1.0-alpha_v)*Ap_u(i,j)*u_hat(i,j)

    end do
  end do

  ! Solve u-velocity equation
  call solver2d_bicgstab2(As_u, Aw_u, Ap_u, Ae_u, An_u, b_u, u_star, m, n-1, solver_tol, maxit)
  !call solver2d_tdma(Aw_u, Ae_u, As_u, An_u, Ap_u, b_u, u_star, m, n-1, solver_tol, maxit)

  ! Update source terms :: v
  do i = 2, m-2
    do j = 2, n-1

      Ap_v(i,j) = Ap_v(i,j)/alpha_v
      b_v(i,j) = b_v(i,j)+dx*(P_star(i,j-1)-P_star(i,j))+(1.0-alpha_v)*Ap_u(i,j)*v_hat(i,j)

    end do
  end do

  ! Solve v-velocity equation
  if (solver .eq. 1) then
    call solver2d_bicgstab2(As_v, Aw_v, Ap_v, Ae_v, An_v, b_v, v_star, m-1, n, solver_tol, maxit)
  else
    call solver2d_tdma(Aw_v, Ae_v, As_v, An_v, Ap_v, b_v, v_star, m-1, n, solver_tol, maxit)
  end if

  return

end subroutine velocity_solve2d
