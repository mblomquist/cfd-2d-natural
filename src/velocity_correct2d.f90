! velocity_correct2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-05-15 (YYYY-MM-DD)
!
! This subroutine corrects the velocity values for a 2D
! CFD problem
!
subroutine velocity_correct2d

  ! Include standard variable header
  include "var2d.dec"

  !u_hat = 0

  ! Print coefficients
  !print *, "Ap_v:", Ap_v


  ! Correct velocity values
  do i = 2,m-1
    do j = 1,n-1
      u(i,j) = u_star(i,j)+dy/Ap_u(i,j)*(P_prime(i-1,j)-P_prime(i,j))
    end do
  end do

  ! Update u-velocity :: Update with relaxation
  !do i = 2,m-1
  !  do j = 1,n-1

  !    if (j .eq. 1) then
  !      u(i,j) = u_hat(i,j)+relax*((Aw_u(i,j)*u_star(i-1,j)+Ae_u(i,j)*u_star(i+1,j)+An_u(i,j)*u_star(i,j+1)+b_u(i,j))/Ap_u(i,j)-u_hat(i,j))
  !    elseif (j .eq. n-1) then
  !      u(i,j) = u_hat(i,j)+relax*((Aw_u(i,j)*u_star(i-1,j)+Ae_u(i,j)*u_star(i+1,j)+As_u(i,j)*u_star(i,j-1)+b_u(i,j))/Ap_u(i,j)-u_hat(i,j))
  !    else
  !      u(i,j) = u_hat(i,j)+relax*((Aw_u(i,j)*u_star(i-1,j)+Ae_u(i,j)*u_star(i+1,j)+As_u(i,j)*u_star(i,j-1)+An_u(i,j)*u_star(i,j+1)+b_u(i,j))/Ap_u(i,j)-u_hat(i,j))
  !    end if

	!  end do
  !end do

  v_hat = 0

  ! Correct velocity values
  do i = 1,m-1
    do j = 2,n-1
      v(i,j) = v_star(i,j)+dx/Ap_v(i,j)*(P_prime(i,j-1)-P_prime(i,j))
    end do
  end do

  ! Update v-velocity :: Update with relaxation
  !do j = 2,n-1
  !  do i = 1,m-1

  !    if (i .eq. 1) then
  !      v(i,j) = v_hat(i,j)+relax*((Ae_v(i,j)*v_star(i+1,j)+As_v(i,j)*v_star(i,j-1)+An_v(i,j)*v_star(i,j+1)+b_v(i,j))/Ap_v(i,j)-v_hat(i,j))
  !    elseif (i .eq. m-1) then
  !      v(i,j) = v_hat(i,j)+relax*((Aw_v(i,j)*v_star(i-1,j)+As_v(i,j)*v_star(i,j-1)+An_v(i,j)*v_star(i,j+1)+b_v(i,j))/Ap_v(i,j)-v_hat(i,j))
  !    else
  !      v(i,j) = v_hat(i,j)+relax*((Aw_v(i,j)*v_star(i-1,j)+Ae_v(i,j)*v_star(i+1,j)+As_v(i,j)*v_star(i,j-1)+An_v(i,j)*v_star(i,j+1)+b_v(i,j))/Ap_v(i,j)-v_hat(i,j))
  !    end if

	!  end do
  !end do

  return

end subroutine velocity_correct2d
