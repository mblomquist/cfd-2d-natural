! velocity2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-02-26 (YYYY-MM-DD)
!
! This subroutine calculates the boundary source terms, sets interior source
! terms, sets boundary values, and initializes the velocity grids.
!
! These definitions are defined for a 2D natural convection problem that 
! simulates a vertical plate at i = 0.

subroutine velocity2d

  ! Pull in standard variable header
  include "var2d.dec"

  ! Define internal variables
  integer :: i, j

  ! ====================== U-Velocity ====================== !
  ! West boundary source terms :: fixed velocity (no-slip)
  Su_u(1, :) = 1e30*u_west
  Sp_u(1, :) = -1e30
  u(1, :) = u_west

  ! South boundary source terms :: none
  Su_u(2:m, 1) = 0
  Sp_u(2:m, 1) = 0

  ! North boundary source terms :: none
  Su_u(2:m, n) = 0
  Sp_u(2:m, n) = 0

  ! East boundary source terms :: none
  Su_u(m, 2:n) = 0
  Sp_u(m, 2:n) = 0

  ! Interior node source terms :: wall shear stress
  Su_u(2:m-1, 2:n-1) = 0
  Sp_u(2:m-1, 2:n-1) = 0

  ! Create initial u-velocity field
  u = 0

  ! ====================== V-Velocity ====================== !
  ! West boundary source terms :: wall shear stress
  Su_v(1, :) = 0
  Sp_v(1, :) = -mu*(dx*dy)/(dx/2)

  ! South boundary source terms :: none
  Su_v(2:m, 1) = 0
  Sp_v(2:m, 1) = 0

  ! North boundary source terms :: none
  Su_v(2:m, n) = 0
  Sp_v(2:m, n) = 0

  ! East boundary source terms :: none
  Su_v(m, 2:n-1) = 0
  Sp_v(m, 2:n-1) = 0

  ! Interior node source terms :: wall shear stress
  Su_v(1:m-1, 2:n-2) = 0
  Sp_v(1:m-1, 2:n-2) = 0

  ! Create initial v-velocity field (set to i>0, :)
  v = 0

  return

end subroutine velocity2d
