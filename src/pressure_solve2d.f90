! pressure_solve2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-05-15 (YYYY-MM-DD)
!
! This subroutine solves the pressure equation for the SIMPLER algorithm
! using u_hat and v_hat.

subroutine pressure_solve2d

  implicit none

  ! Include variable header
  include "var2d.dec"

  ! Define internal variables
  integer :: i, j

  ! Update coefficients
  do i = 1,m-1
    do j = 1,n-1

      ! Update coefficients
	    Ae_p(i,j) = rho*u0*dy*dy/Ap_u(i+1,j)
	    Aw_p(i,j) = rho*u0*dy*dy/Ap_u(i,j)
	    An_p(i,j) = rho*u0*dx*dx/Ap_v(i,j+1)
	    As_p(i,j) = rho*u0*dx*dx/Ap_v(i,j)

	    ! Check west node
      if (i .eq. 1) then
        Aw_p(i,j) = 0
      end if

	    ! Check east node
	    if (i .eq. m-1) then
	      Ae_p(i,j) = 0
	    end if

	    ! Check south node
      if (j .eq. 1) then
        As_p(i,j) = 0
      end if

	    ! Check north node
      if (j .eq. n-1) then
        An_p(i,j) = 0
      end if

	    ! Ap coefficient
      Ap_p(i,j) = As_p(i,j)+Aw_p(i,j)+Ae_p(i,j)+An_p(i,j)-Sp_p(i,j)

	    ! Update b values
	    b_p(i,j) = rho*dx*(u_hat(i,j)-u_hat(i+1,j))+rho*dy*(v_hat(i,j)-v_hat(i,j+1))

    end do
  end do

  ! Set reference pressure node (east-north corner)
  Aw_p(m-1,n-1) = 0
  Ae_p(m-1,n-1) = 0
  As_p(m-1,n-1) = 0
  An_p(m-1,n-1) = 0

  Ap_p(m-1,n-1) = 1
  b_p(m-1,n-1) = 0

  ! Solve pressure equation
  !call solver2d_bicgstab2(As_p, Aw_p, Ap_p, Ae_p, An_p, b_p, P, m-1, n-1, solver_tol, maxit)
  call solver2d_tdma(Aw_p, Ae_p, As_p, An_p, Ap_p, b_p, P, m-1, n-1, solver_tol, maxit)

  return

end subroutine pressure_solve2d
