! velocity2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-03-07 (YYYY-MM-DD)
!
! This subroutine calculates the boundary source terms, sets interior source
! terms, sets boundary values, and initializes the velocity grids.
!
! These definitions are defined for a 2D cfd problem with an inlet boundary
! condition (west i=0), an outlet boundary condition (east i=m), a wall
! (north j=1), and a wall (south j=n).

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
  Su_u(2:m-1, 1) = 0
  Sp_u(2:m-1, 1) = -mu*(dx*dy)/(dy/2)

  ! North boundary source terms :: none
  Su_u(2:m-1, n-1) = 0
  Sp_u(2:m-1, n-1) = -mu*(dx*dy)/(dy/2)


  ! ====================== V-Velocity ====================== !
  ! Initialize v-velocity field
  v = 0
  v_hat = v
  v_star = v

  ! Set v-velocity source terms
  Su_v = 0
  Sp_v = 0

  ! South boundary source terms :: no slip
  Su_v(:, 1) = 0
  Sp_v(:, 1) = -1e30

  ! North boundary source terms :: no slip
  Su_v(:, n) = 0
  Sp_v(:, n) = -1e30

  return

end subroutine velocity2d
