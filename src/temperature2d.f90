! temperature2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-02-26 (YYYY-MM-DD)
!
! This subroutine calculates the boundary source terms, sets interior source
! terms to 0, sets boundary values, and initializes the temperature grid.
!
! These definitions are defined for a 2D natural convection problem that 
! simulates a vertical plate at i = 0.

subroutine temperature2d

  ! Pull in standard variable header
  include "var2d.dec"

  ! Vertical Plate :: Fixed Temperature 
  T(0, :) = T_west
  Su_t(1, :) = 1e30*T_west
  Sp_t(1, :) = -1e30

  ! West boundary source terms ::
  Su_t(1, :) = -2*k_const*A_y/dy
  Sp_t(1, :) = 2*k_const*A_y/dy*T_west

  ! South boundary source terms :: fixed temperature (T_inf)
  Su_t(2:m, 1) = 1e30*T_inf
  Sp_t(2:m, 1) = -1e30

  ! North boundary source terms :: adiabatic
  Su_t(2:m, n) = 0
  Sp_t(2:m, n) = 0

  ! East boundary source terms :: adiabatic
  Su_t(m, 2:n) = 0
  Sp_t(m, 2:n) = 0

  ! Interior nodes source terms :: none
  Su_t(2:m-1, 2:n-1) = 0
  Sp_t(2:m-1, 2:n-1) = 0

  ! Create initial temperature field (set i>0, :)
  T(1:m, :) = T_inf

  return

end subroutine temperature2d
