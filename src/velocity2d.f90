! velocity2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-04-03 (YYYY-MM-DD)
!
! This subroutine calculates the boundary source terms, sets interior source
! terms, sets boundary values, and initializes the velocity grids.

subroutine velocity2d

  ! Pull in standard variable header
  include "var2d.dec"

  ! Define internal variables
  integer :: i, j

  ! ====================== U-Velocity ====================== !
  ! Initialize u-velocity field
  u = 0
  u_hat = u
  u_star = u

  ! Set u-velocity source terms
  Su_u = 0
  Sp_u = 0

  ! South boundary source terms :: none
  Su_u(:, 1) = 0
  Sp_u(:, 1) = -2*mu*dy/dx

  ! North boundary source terms :: none
  Su_u(:, n-1) = 0
  Sp_u(:, n-1) = -2*mu*dy/dx


  ! ====================== V-Velocity ====================== !
  ! Initialize v-velocity field
  v = 0
  v_hat = v
  v_star = v

  ! Set v-velocity source terms
  Su_v = 0
  Sp_v = 0

  ! West boundary source terms :: no slip
  Su_v(1, :) = 0
  Sp_v(1, :) = -2*mu*dx/dy

  ! East boundary source terms :: no slip
  Su_v(m-1, :) = 0
  Sp_v(m-1, :) = -2*mu*dx/dy

  return

end subroutine velocity2d
