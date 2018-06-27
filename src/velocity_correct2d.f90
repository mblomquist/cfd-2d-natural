! velocity_correct2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-06-27 (YYYY-MM-DD)
!
! This subroutine corrects the velocity values for a 2D
! CFD problem
!
subroutine velocity_correct2d

  ! Include standard variable header
  include "var2d.dec"

  ! Correct velocity values
  do i = 2,m-1
    do j = 1,n-1
      u(i,j) = u_star(i,j)+dy/Ap_u(i,j)*(P_prime(i-1,j)-P_prime(i,j))*alpha_v
    end do
  end do

  ! Correct wall values
  u(1, :) = 0
  u(m, :) = 0

  ! Correct velocity values
  do i = 1,m-1
    do j = 2,n-1
      v(i,j) = v_star(i,j)+dx/Ap_v(i,j)*(P_prime(i,j-1)-P_prime(i,j))*alpha_v
    end do
  end do

  ! Correct wall values
  v(:, 1) = 0
  v(:, n) = 0

  return

end subroutine velocity_correct2d
