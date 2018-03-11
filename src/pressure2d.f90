! pressure2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-03-07 (YYYY-MM-DD)
!
! This subroutine calculates the boundary source terms, sets interior source
! terms to 0, sets boundary values, and initializes the pressure grid.
!
! These definitions are defined for a 2D cfd problem with an inlet boundary
! condition (west i=0), an outlet boundary condition (east i=m), a wall
! (north j=1), and a wall (south j=n).

subroutine pressure2d

  ! Pull in standard variable header
  include "var2d.dec"

  ! Initialize Pressure Field
  P = 0
  P_star = 0
  P_prime = 0

  ! Set pressure source terms
  Su_p = 0
  Sp_p = 0

  ! Set pressure reference
  !Sp_p(m-1,1) = -1e30
  !Su_p(m-1,1) = 0

  return

end subroutine pressure2d
