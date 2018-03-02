! pseudo_solve2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-02-26 (YYYY-MM-DD)
!
! This subroutine computes the pseudo velocities (u_hat and v_hat) for a 2D
! CFD problem.
!
! These definitions are defined for a 2D natural convection problem that 
! simulates a vertical plate at i = 0.

subroutine pseudo_solve2d

  implicit none

  ! Include standard variable header
  include "var2d.dec"

  ! Define internal variables
  integer :: i, j

  ! ========================== u_hat ========================== !

  ! Update stiffness coefficients
  call pseudo_source2d("u")

  ! Calculate interior nodes
  do i = 2,m-1
    do j = 2,n-1

      u_hat(i,j) = (Aw_u(i,j)*u_star(i-1,j)+Ae_u(i,j)*u_star(i+1,j)+As_u(i,j)*u_star(i,j-1)+An_u(i,j)*u_star(i,j+1)+b_u(i,j))/Ap_u(i,j)

	  end do
  end do

  ! Copy boundary nodes
  u_hat(m,:) = u_hat(m,:)
  u_hat(2:m-1,n) = u_hat(2:m-1,n-1)
  u_hat(2:m-1,1) = u_hat(2:m-1,2)

  ! Hard set west boundary to 0
  u_hat(1,:) = 0

  ! ========================== v_hat ========================== !

  ! Update stiffness coefficients
  call pseudo_source2d("v")

  ! Calculate interior nodes
  do j = 2,n-1
    do i = 2,m-1

      v_hat(i,j) = (Aw_v(i,j)*v_star(i-1,j)+Ae_v(i,j)*v_star(i+1,j)+As_v(i,j)*v_star(i,j-1)+An_v(i,j)*v_star(i,j+1)+b_v(i,j))/Ap_v(i,j)

	  end do
  end do


  ! Calculate west nodes
  do j = 1, n

    ! Set i
	i = 1
	v_hat(i,j) = (Ae_v(i,j)*v_star(i+1,j)+As_v(i,j)*v_star(i,j-1)+An_v(i,j)*v_star(i,j+1)+b_v(i,j))/Ap_v(i,j)

  end do

  ! Copy boundary nodes
  v_hat(2:m-1,1) = v_hat(2:m-1,2)
  v_hat(2:m-1,n) = v_hat(2:m-1,n-1)
  v_hat(m,:) = v_hat(m-1,:)

  return

end subroutine pseudo_solve2d
