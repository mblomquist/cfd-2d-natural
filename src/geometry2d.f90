! geometry2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-04-03 (YYYY-MM-DD)
!
! This subroutine calculates the geometry properties and associated
! non-dimensional quantities for use in the SIMPLER method.

subroutine geometry2d

  ! Pull in standard variable header
  include "var2d.dec"

  ! Calculate non-dimensional quantities for SIMPLER
  dx = 1.0/m
  dy = 1.0/n
  dz = 1

  return

end subroutine geometry2d
