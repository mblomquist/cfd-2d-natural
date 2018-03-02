! pressure2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-02-26 (YYYY-MM-DD)
!
! This subroutine calculates the boundary source terms, sets interior source
! terms to 0, sets boundary values, and initializes the pressure grid.
!
! These definitions are defined for a 2D natural convection problem that 
! simulates a vertical plate at i = 0.

subroutine pressure2d

  ! Pull in standard variable header
  include "var2d.dec"

  ! West boundary source terms :: none
  Su_p(0, :) = 0
  Sp_p(0, :) = 0

  ! South boundary source terms :: none
  Su_p(1:m, 1) = 0
  Sp_p(1:m, 1) = 0

  ! North boundary source terms :: none
  Su_p(1:m, n-1) = 0
  Sp_p(1:m, n-1) = 0

  ! East boundary source terms :: reference pressure
  Su_p(m, :) = 1e30*p_east
  Sp_p(m, :) = -1e30

  ! Interior node source terms
  Su_p(1:m-1, 2:n-2) = 0
  Sp_p(1:m-1, 2:n-2) = 0

  ! Create initial pressure field (set to 0 i>0,:) in future := P_in
  P = p_east

  return

end subroutine pressure2d
