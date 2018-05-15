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

  ! Correct velocity values
  do i = 2,m-1
    do j = 1,n-1
      u(i,j) = u_star(i,j)+dy*dz/Ap_u(i,j)*(P_prime(i-1,j)-P_prime(i,j))
    end do
  end do

  ! Correct velocity values
  do i = 1,m-1
    do j = 2,n-1
      v(i,j) = v_star(i,j)+dx*dz/Ap_v(i,j)*(P_prime(i,j-1)-P_prime(i,j))
    end do
  end do

  return

end subroutine velocity_correct2d
