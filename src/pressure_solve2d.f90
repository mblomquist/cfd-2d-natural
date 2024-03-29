! pressure_solve2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-06-27 (YYYY-MM-DD)
!
! This subroutine solves the pressure equation for the SIMPLER algorithm
! using u_hat and v_hat.

subroutine pressure_solve2d

  implicit none

  ! Include variable header
  include "var2d.dec"

  ! Define internal variables
  integer :: i, j, fault

  Aw_p = 0.
  Ae_p = 0.
  As_p = 0.
  An_p = 0.
  Ap_p = 0.

  ! Update coefficients
  do i = 1,m-1
    do j = 1,n-1

      ! Update coefficients
	    Aw_p(i,j) = dy/Ap_u(i,j)*alpha_v
      Ae_p(i,j) = dy/Ap_u(i+1,j)*alpha_v
	    As_p(i,j) = dx/Ap_v(i,j)*alpha_v
      An_p(i,j) = dx/Ap_v(i,j+1)*alpha_v

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
      Ap_p(i,j) = As_p(i,j)+Aw_p(i,j)+Ae_p(i,j)+An_p(i,j)

	    ! Update b values
	    b_p(i,j) = ((u_hat(i,j)-u_hat(i+1,j))*dy+(v_hat(i,j)-v_hat(i,j+1))*dx)

    end do
  end do

  ! Set reference pressure node (east-north corner)
  Aw_p(m-1,1) = 0
  Ae_p(m-1,1) = 0
  As_p(m-1,1) = 0
  An_p(m-1,1) = 0

  Ap_p(m-1,1) = 1
  b_p(m-1,1) = 0

  ! Print coefficients
  !print *, "Aw_p:", Aw_p
  !print *, "Ae_p:", Ae_p
  !print *, "An_p:", As_p
  !print *, "As_p:", An_p
  !print *, "Ap_p:", Ap_p
  !print *, "b_p:", b_p

  P = 0

  ! Solve pressure equation
  if (solver .eq. 0) then
    call solver2d_bicgstab(As_p, Aw_p, Ap_p, Ae_p, An_p, b_p, P, m-1, n-1, solver_tol, maxit)
  elseif (solver .eq. 1) then
    call solver2d_bicgstab2(As_p, Aw_p, Ap_p, Ae_p, An_p, b_p, P, m-1, n-1, solver_tol, maxit)
  elseif (solver .eq. 2) then
    fault = 0
    do i = 3,maxit
      if (fault .eq. 0) then
        call solver2d_gmres(As_p, Aw_p, Ap_p, Ae_p, An_p, b_p, P, m-1, n-1, solver_tol, i, fault)
        if (fault .eq. 0) then
          !print *, "Restarting GMRES."
        end if
      end if
    end do
  elseif (solver .eq. 3) then
    call solver2d_bicg(As_p, Aw_p, Ap_p, Ae_p, An_p, b_p, P, m-1, n-1, solver_tol, maxit)
  else
    call solver2d_tdma(As_p, Aw_p, Ap_p, Ae_p, An_p, b_p, P, m-1, n-1, solver_tol, maxit)
  end if

  return

end subroutine pressure_solve2d
