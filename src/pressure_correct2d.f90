! pressure_correct2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-03-08 (YYYY-MM-DD)
!
! This subroutine solves the pressure equation for the SIMPLER algorithm
! using u_star and v_star.

subroutine pressure_correct2d

  implicit none

  ! Include variable header
  include "var2d.dec"

  ! Define internal variables
  integer :: i, j

  ! Update reference point :: set to 0
  Ae_p(m-1,1) = 0
  Aw_p(m-1,1) = 0
  An_p(m-1,1) = 0
  As_p(m-1,1) = 0

  Ap_p(m-1,1) = 1
  b_p(m-1,1) = 0

  ! Update east boundary :: outlet boundary // no link to Ae_p
  ! Set i
  i = m-1

  do j = 2,n-1

    ! Update coefficients
    Ae_p(i,j) = 0
    Aw_p(i,j) = 0 !rho*A_x/Ap_u(i,j)*alpha
    An_p(i,j) = 0 !rho*A_y/Ap_v(i,j+1)*alpha
    As_p(i,j) = 0 !rho*A_y/Ap_v(i,j)*alpha

    if (j .eq. 1) then
      As_p(i,j) = 0
    end if

    if (j .eq. n-1) then
      An_p(i,j) = 0
    end if

    Ap_p(i,j) = 1 !As_p(i,j)+Aw_p(i,j)+Ae_p(i,j)+An_p(i,j)-Sp_p(i,j)

    ! Update b values
    b_p(i,j) = 0 !rho*(v_star(i,j)-v_star(i,j+1))*A_y+Su_p(i,j)

  end do

  ! Update internal coefficients
  do i = 1,m-2
    do j = 1,n-1

      ! Update coefficients
	    Ae_p(i,j) = rho*A_x/Ap_u(i+1,j)*alpha
	    Aw_p(i,j) = rho*A_x/Ap_u(i,j)*alpha
	    An_p(i,j) = rho*A_y/Ap_v(i,j+1)*alpha
	    As_p(i,j) = rho*A_y/Ap_v(i,j)*alpha

      if (i .eq. 1) then
        Aw_p(i,j) = 0
      end if

      if (j .eq. 1) then
        As_p(i,j) = 0
      end if

      if (j .eq. n-1) then
        An_p(i,j) = 0
      end if

      Ap_p(i,j) = As_p(i,j)+Aw_p(i,j)+Ae_p(i,j)+An_p(i,j)-Sp_p(i,j)

	  ! Update b values
	  b_p(i,j) = rho*(u_star(i,j)-u_star(i+1,j))*A_x+rho*(v_star(i,j)-v_star(i,j+1))*A_y+Su_p(i,j)

    end do
  end do

  ! Reset P to 0
  !P_prime = 0

  ! Solve pressure equation
  call solver2d_bicgstab(As_p, Aw_p, Ap_p, Ae_p, An_p, b_p, P_prime, m-1, n-1, tol, maxit)

  return

end subroutine pressure_correct2d
