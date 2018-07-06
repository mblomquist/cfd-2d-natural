! simpler2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-06-27 (YYYY-MM-DD)
!
! This subroutine runs the SIMPLER algorithm for a 2D CFD problem.
!
subroutine simpler2d

  implicit none

  ! Pull in standard variable header
  include "var2d.dec"

  ! Define Internal Variables
  integer :: i, j

  print *, 'Start SIMPLER Algorithm.'

  ! Set initial guesses
  P = 0
  u = 0
  v = 0
  T = 0

  P_star = P
  u_star = u
  v_star = v

  P_prime = P
  u_hat = u
  v_hat = v

  ! Solve Temperature for Natural Convection First
  !print *, "Step 0: Solve Temperature Equation"
  call temperature_solve2d
  !print *, ".............."
  !print *, "T:", T
  !print *, ".............."

  do i = 1,itrmax

    ! Update stiffness coefficients :: u
    call velocity_source2d("u")

    ! Print coefficients
    !print *, "Aw_u:", Aw_u
    !print *, "Ae_u:", Ae_u
    !print *, "An_u:", As_u
    !print *, "As_u:", An_u
    !print *, "Ap_u:", Ap_u
    !print *, "b_u:", b_u

    ! Update stiffness coefficients :: v
    call velocity_source2d("v")

    ! Print coefficients
    !print *, "Aw_v:", Aw_v
    !print *, "Ae_v:", Ae_v
    !print *, "An_v:", As_v
    !print *, "As_v:", An_v
    !print *, "Ap_v:", Ap_v
    !print *, "b_v:", b_v

    ! Step 2: Calculate Pseudo-Velocities
    !print *, "Step 1: Solve Pseudo-Velocities"
    call pseudo_solve2d
    !print *, ".............."
    !print *, "u_hat:", u_hat
    !print *, "v_hat:", v_hat
    !print *, ".............."

    ! Step 3: Solve Pressure Equation
    !print *, "Step 2: Solve Pressure Equation"
    call pressure_solve2d
    !print *, ".............."
    !print *, "P:", P
    !print *, ".............."

	  ! Set p_star := P
	  P_star = P

    ! Step 8: Check Convergence
    print *, "Step 7: Check Convergence"
    call convergence2d(i)

    if (i .eq. 1) then
      R_e = 1.0
      R_t = 1.0
    end if

    print *, "Iteration:", i
    print *, "Relative Momentum Error: ", R_e
    print *, "Relative Energy Error:", R_t
    !print *, "Maximum Error: ", e_max_new
    !print *, "Previous Maximum Error:", e_max_old

    if ((R_e .le. simpler_tol) .and. (R_t .le. simpler_tol)) then

      print *, "Simpler completed in: ", i
      exit

    end if

    ! Step 4: Solve Momentum Equations
    !print *, "Step 3: Solve Momentum Equations"
    call velocity_solve2d
    !print *, ".............."
    !print *, "v_star:", v_star
    !print *, "u_star:", u_star
    !print *, ".............."

    ! Step 5: Solve Pressure Equation
    !print *, "Step 4: Solve Pressure Correction"
    call pressure_correct2d
    !print *, ".............."
    !print *, "P_prime:", P_prime
    !print *, ".............."

    ! Step 6: Correct Velocities
    !print *, "Step 5: Correct Velocities"
    call velocity_correct2d
    !print *, ".............."
    !print *, "u:", u
    !print *, "v:", v
    !print *, ".............."

    ! Step 7: Solve Temperature Equation
    !print *, "Step 6: Solve Temperature Equation"
    call temperature_solve2d
    !print *, ".............."
    !print *, "T:", T
    !print *, ".............."

    P_star = P
	  u_star = u
	  v_star = v

  end do

  return

end subroutine simpler2d
