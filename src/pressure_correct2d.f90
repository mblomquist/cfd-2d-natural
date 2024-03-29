! pressure_correct2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-06-27 (YYYY-MM-DD)
!
! This subroutine solves the pressure equation for the SIMPLER algorithm
! using u_star and v_star.

subroutine pressure_correct2d

  implicit none

  ! Include variable header
  include "var2d.dec"

  ! Define internal variables
  integer :: i, j, fault

  ! Update coefficients
  do i = 1,m-1
    do j = 1,n-1

	    ! Update b values
	    b_p(i,j) = ((u_star(i,j)-u_star(i+1,j))*dy+(v_star(i,j)-v_star(i,j+1))*dx)

    end do
  end do

  ! Set reference pressure node (east-south corner)
  Aw_p(m-1,1) = 0
  Ae_p(m-1,1) = 0
  As_p(m-1,1) = 0
  An_p(m-1,1) = 0

  Ap_p(m-1,1) = 1
  b_p(m-1,1) = 0

  ! Print coefficients
  !print *, "b_p:", b_p

  P_prime = 0

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

end subroutine pressure_correct2d
