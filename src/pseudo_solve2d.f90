! pseudo_solve2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-03-07 (YYYY-MM-DD)
!
! This subroutine computes the pseudo velocities (u_hat and v_hat) for a 2D
! CFD problem.
!
! These definitions are defined for a 2D cfd problem with an inlet boundary
! condition (west i=0), an outlet boundary condition (east i=m), a wall
! (north j=1), and a wall (south j=n).

subroutine pseudo_solve2d

  implicit none

  ! Include standard variable header
  include "var2d.dec"

  ! Define internal variables
  integer :: i, j

  ! ========================== v_hat ========================== !

  ! Update stiffness coefficients
  call pseudo_source2d("v")

  ! Calculate west nodes
  ! Set i
  i = 1
  do j = 2, n-1

	  v_hat(i,j) = (Ae_v(i,j)*v_star(i+1,j)+As_v(i,j)*v_star(i,j-1)+An_v(i,j)*v_star(i,j+1)+b_v(i,j))/Ap_v(i,j)

  end do

  ! Calulate south nodes :: fixed value
  v_hat(1:m-1, 1) = v_star(1:m-1, 1)

  ! Calculate north nodes :: fixed value
  v_hat(1:m-1, n) = v_star(1:m-1, n)

  ! Calculate interior nodes
  do j = 2,n-1
    do i = 2,m-1

      v_hat(i,j) = (Aw_v(i,j)*v_star(i-1,j)+Ae_v(i,j)*v_star(i+1,j)+As_v(i,j)*v_star(i,j-1)+An_v(i,j)*v_star(i,j+1)+b_v(i,j))/Ap_v(i,j)

	  end do
  end do

  v_hat(m,:) = v_hat(m-1,:)

  ! ========================== u_hat ========================== !

  ! Update stiffness coefficients
  call pseudo_source2d("u")

  ! Calculate west nodes :: fixed value
  u_hat(1, :) = u_star(1, :)

  ! Calculate north and south nodes
  do i = 2,m-1

    ! Calculate south nodes :: no input from south
    j = 1
    u_hat(i,j) = (Aw_u(i,j)*u_star(i-1,j)+Ae_u(i,j)*u_star(i+1,j)+An_u(i,j)*u_star(i,j+1)+b_u(i,j))/Ap_u(i,j)

	  ! Calculate north nodes :: no input from north
	  j = n-1
	  u_hat(i,j) = (Aw_u(i,j)*u_star(i-1,j)+Ae_u(i,j)*u_star(i+1,j)+As_u(i,j)*u_star(i,j-1)+b_u(i,j))/Ap_u(i,j)

  end do

  ! Calculate interior nodes
  do i = 2,m-1
    do j = 2,n-2

      u_hat(i,j) = (Aw_u(i,j)*u_star(i-1,j)+Ae_u(i,j)*u_star(i+1,j)+As_u(i,j)*u_star(i,j-1)+An_u(i,j)*u_star(i,j+1)+b_u(i,j))/Ap_u(i,j)

	  end do
  end do

  ! Update east nodes
  u_hat(m,:) = u_hat(m-1,:)

  return

end subroutine pseudo_solve2d
