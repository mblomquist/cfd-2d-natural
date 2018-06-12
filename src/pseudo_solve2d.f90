! pseudo_solve2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-03=5-15 (YYYY-MM-DD)
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

  ! ========================== u_hat ========================== !

  ! Update stiffness coefficients
  call pseudo_source2d("u")

  ! Print coefficients
  !print *, "Aw_u:", Aw_u
  !print *, "Ae_u:", Ae_u
  !print *, "An_u:", As_u
  !print *, "As_u:", An_u
  !print *, "Ap_u:", Ap_u
  !print *, "b_u:", b_u

  u_hat = 0

  ! Calculate nodes
  do i = 1,m
    do j = 1,n-1

      if (i .eq. 1) then

        if (j .eq. 1) then
          u_hat(i,j) = (Ae_u(i,j)*u_star(i+1,j)+An_u(i,j)*u_star(i,j+1)+b_u(i,j))/Ap_u(i,j)
        elseif (j .eq. n-1) then
          u_hat(i,j) = (Ae_u(i,j)*u_star(i+1,j)+As_u(i,j)*u_star(i,j-1)+b_u(i,j))/Ap_u(i,j)
        else
          u_hat(i,j) = (Ae_u(i,j)*u_star(i+1,j)+As_u(i,j)*u_star(i,j-1)+An_u(i,j)*u_star(i,j+1)+b_u(i,j))/Ap_u(i,j)
        end if

      elseif (i .eq. m) then

        if (j .eq. 1) then
          u_hat(i,j) = (Aw_u(i,j)*u_star(i-1,j)+An_u(i,j)*u_star(i,j+1)+b_u(i,j))/Ap_u(i,j)
        elseif (j .eq. n-1) then
          u_hat(i,j) = (Aw_u(i,j)*u_star(i-1,j)+As_u(i,j)*u_star(i,j-1)+b_u(i,j))/Ap_u(i,j)
        else
          u_hat(i,j) = (Aw_u(i,j)*u_star(i-1,j)+As_u(i,j)*u_star(i,j-1)+An_u(i,j)*u_star(i,j+1)+b_u(i,j))/Ap_u(i,j)
        end if

      else

        if (j .eq. 1) then
          u_hat(i,j) = (Aw_u(i,j)*u_star(i-1,j)+Ae_u(i,j)*u_star(i+1,j)+An_u(i,j)*u_star(i,j+1)+b_u(i,j))/Ap_u(i,j)
        elseif (j .eq. n-1) then
          u_hat(i,j) = (Aw_u(i,j)*u_star(i-1,j)+Ae_u(i,j)*u_star(i+1,j)+As_u(i,j)*u_star(i,j-1)+b_u(i,j))/Ap_u(i,j)
        else
          u_hat(i,j) = (Aw_u(i,j)*u_star(i-1,j)+Ae_u(i,j)*u_star(i+1,j)+As_u(i,j)*u_star(i,j-1)+An_u(i,j)*u_star(i,j+1)+b_u(i,j))/Ap_u(i,j)
        end if

      end if

	  end do
  end do

  ! ========================== v_hat ========================== !

  ! Update stiffness coefficients
  call pseudo_source2d("v")

  ! Print coefficients
  !print *, "Aw_v:", Aw_v
  !print *, "Ae_v:", Ae_v
  !print *, "An_v:", As_v
  !print *, "As_v:", An_v
  !print *, "Ap_v:", Ap_v
  !print *, "b_v:", b_v

  v_hat = 0

  ! Calculate nodes
  do j = 1,n
    do i = 1,m-1

      if (j .eq. 1) then

        if (i .eq. 1) then
          v_hat(i,j) = (Ae_v(i,j)*v_star(i+1,j)+An_v(i,j)*v_star(i,j+1)+b_v(i,j))/Ap_v(i,j)
        elseif (i .eq. m-1) then
          v_hat(i,j) = (Aw_v(i,j)*v_star(i-1,j)+An_v(i,j)*v_star(i,j+1)+b_v(i,j))/Ap_v(i,j)
        else
          v_hat(i,j) = (Aw_v(i,j)*v_star(i-1,j)+Ae_v(i,j)*v_star(i+1,j)+An_v(i,j)*v_star(i,j+1)+b_v(i,j))/Ap_v(i,j)
        end if

      elseif (j .eq. n) then

        if (i .eq. 1) then
          v_hat(i,j) = (Ae_v(i,j)*v_star(i+1,j)+As_v(i,j)*v_star(i,j-1)+b_v(i,j))/Ap_v(i,j)
        elseif (i .eq. m-1) then
          v_hat(i,j) = (Aw_v(i,j)*v_star(i-1,j)+As_v(i,j)*v_star(i,j-1)+b_v(i,j))/Ap_v(i,j)
        else
          v_hat(i,j) = (Aw_v(i,j)*v_star(i-1,j)+Ae_v(i,j)*v_star(i+1,j)+As_v(i,j)*v_star(i,j-1)+b_v(i,j))/Ap_v(i,j)
        end if

      else

        if (i .eq. 1) then
          v_hat(i,j) = (Ae_v(i,j)*v_star(i+1,j)+As_v(i,j)*v_star(i,j-1)+An_v(i,j)*v_star(i,j+1)+b_v(i,j))/Ap_v(i,j)
        elseif (i .eq. m-1) then
          v_hat(i,j) = (Aw_v(i,j)*v_star(i-1,j)+As_v(i,j)*v_star(i,j-1)+An_v(i,j)*v_star(i,j+1)+b_v(i,j))/Ap_v(i,j)
        else
          v_hat(i,j) = (Aw_v(i,j)*v_star(i-1,j)+Ae_v(i,j)*v_star(i+1,j)+As_v(i,j)*v_star(i,j-1)+An_v(i,j)*v_star(i,j+1)+b_v(i,j))/Ap_v(i,j)
        end if

      end if

	  end do
  end do
  
  return

end subroutine pseudo_solve2d
