! pressure_correct2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-02-26 (YYYY-MM-DD)
!
! This subroutine solves the pressure equation for the SIMPLER algorithm
! using u_star and v_star.
!
! These definitions are defined for a 2D natural convection problem that 
! simulates a vertical plate at i = 0.

subroutine pressure_correct2d

  implicit none

  ! Include variable header
  include "var2d.dec"

  ! Define internal variables
  integer :: i, j

  ! Update west boundary
  do j = 1,n-1
    
	! Set i
	i = 1

    ! Update coefficients
	Ae_p(i,j) = rho*A_x/Ap_u(i+1,j)*alpha
	Aw_p(i,j) = 0
	An_p(i,j) = rho*A_y/Ap_v(i,j+1)*alpha
	As_p(i,j) = rho*A_y/Ap_v(i,j)*alpha

    Ap_p(i,j) = As_p(i,j)+Aw_p(i,j)+Ae_p(i,j)+An_p(i,j)-Sp_p(i,j)

	! Update b values
	b_p(i,j) = rho*(u_star(i,j)-u_star(i+1,j))*A_x+rho*(v_star(i,j)-v_star(i,j+1))*A_y+Su_p(i,j)

  end do

  ! Update east boundary :: constant pressure
  do j = 1,n-1
    
	! Set i
	i = m-1
	
	! Update coefficients
	Ae_p(i,j) = 0
	Aw_p(i,j) = rho*A_x/Ap_u(i,j)*alpha
	An_p(i,j) = rho*A_y/Ap_v(i,j+1)*alpha
	As_p(i,j) = rho*A_y/Ap_v(i,j)*alpha

    Ap_p(i,j) = As_p(i,j)+Aw_p(i,j)+Ae_p(i,j)+An_p(i,j)-Sp_p(i,j)

	! Update b values
	b_p(i,j) = rho*(u_star(i,j)-u_star(i+1,j))*A_x+rho*(v_star(i,j)-v_star(i,j+1))*A_y+Su_p(i,j)

  end do

  ! Update internal coefficients
  do i = 2,m-1
    do j = 2,n-1

      ! Update coefficients
	  Ae_p(i,j) = rho*A_x/Ap_u(i+1,j)*alpha
	  Aw_p(i,j) = rho*A_x/Ap_u(i,j)*alpha
	  An_p(i,j) = rho*A_y/Ap_v(i,j+1)*alpha
	  As_p(i,j) = rho*A_y/Ap_v(i,j)*alpha

      Ap_p(i,j) = As_p(i,j)+Aw_p(i,j)+Ae_p(i,j)+An_p(i,j)-Sp_p(i,j)

	  ! Update b values
	  b_p(i,j) = rho*(u_star(i,j)-u_star(i+1,j))*A_x+rho*(v_star(i,j)-v_star(i,j+1))*A_y+Su_p(i,j)

    end do
  end do

  P(1:m-1,:) = 0

  ! Solve pressure equation
  call solver2d_bicgstab(As_p(1:m-1,:), Aw_p(1:m-1,:), Ap_p(1:m-1,:), Ae_p(1:m-1,:), An_p(1:m-1,:), b_p(1:m-1,:), P(1:m-1,:), m-1, n-1, tol, maxit)

  return

end subroutine pressure_correct2d
