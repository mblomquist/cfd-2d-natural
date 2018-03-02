! geometry2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-02-19 (YYYY-MM-DD)
!
! This subroutine calculates the geometry properties and associated
! non-dimensional quantities for use in the SIMPLER method.

subroutine geometry2d

  ! Pull in standard variable header
  include "var2d.dec"

  ! Calculate non-dimensional quantities for SIMPLER
  dx = length/m
  dy = width/n
  dz = depth/l

  A_x = dy*dz
  A_y = dz*dx
  A_z = dx*dy

  return

end subroutine geometry2d
