! initialize2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-02-26 (YYYY-MM-DD)
!
! This subroutine runs the static calculations for geometry properties,
! pressure properties, velocity, properties, temperature properties, and
! sets the boundary conditions for pressure, velocity, and temperature.
!
! These definitions are defined for a 2D natural convection problem that 
! simulates a vertical plate at i = 0.

subroutine initialize2d

  ! Pull in standard variable header
  include "var2d.dec"

  ! Read Input File ..........
  ! ..........................
  ! ..........................

  ! Define geometry variables
  length = 1   ! 1 meter long
  width = 1    ! 1 meter wide
  depth = 1    ! 1 meter deep
  g = 9.81

  ! Define media variables
  ReL = 0
  Pr = 3.56
  Pe = ReL*Pr
  rho = 1.093
  mu = 13.49
  k_const = 24.35
  Cp = 1.006
  beta = 0.0031

  ! Define Pressure boundary condition
  P_east = 1 ! 202650     ! Pa, 2 atm

  ! Define U-Velocity boundary condition (no slip)
  u_west = 0

  ! Define Temperature Boundary conditions
  T_west = 100       ! Inlet fluid temperature
  T_south = 50     ! Heat flux from south

  ! Define solution parameters
  itrmax = 1
  maxit = 1000
  tol = 1e-4
  alpha = 1.0

  ! Calculate geometry properties.
  call geometry2d

  ! Calculate pressure source terms and initialize grid
  call pressure2d

  ! Calculate velocity source terms and initialize grid
  call velocity2d

  ! Calculate temperature source terms and initialize grid
  call temperature2d

  return

end subroutine initialize2d
