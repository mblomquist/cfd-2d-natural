! pressure2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-04-02 (YYYY-MM-DD)
!
! This subroutine calculates the boundary source terms, sets interior source
! terms to 0, sets boundary values, and initializes the pressure grid.
!

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

  return

end subroutine pressure2d
