! initialize2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-03-08 (YYYY-MM-DD)
!
! This subroutine runs the static calculations for geometry properties,
! pressure properties, velocity, properties, temperature properties, and
! sets the boundary conditions for pressure, velocity, and temperature.
!
! These definitions are defined for a 2D cfd problem with an inlet boundary
! condition (left i=0), an outlet boundary condition (right i=m), a wall
! (north j=1), and a wall (south j=n).

subroutine initialize2d

  ! Pull in standard variable header
  include "var2d.dec"

  ! Read Input File ..........
  ! ..........................
  ! ..........................

  ! Define geometry variables
  length = 1   ! 8 inches long
  width = 1      ! 1 inch wide
  depth = 1      ! 1 inch deep

  ! Define media variables
  ReL = 0
  Pr = 3.56
  Pe = ReL*Pr
  rho = 1000
  mu = 547.4
  k_const = 640.6
  Cp = 4.187
  beta = 0.000214

  ! Define Pressure boundary condition
  P_east = 0

  ! Define U-Velocity boundary condition
  u_west = ReL*mu/rho/length

  ! Define V-Velocity boundary conditions
  v_west = 0
  v_north = 0
  v_south = 0

  ! Define Temperature Boundary conditions
  T_north = 50       ! Inlet fluid temperature
  T_south = 100     ! Heat flux from south

  ! Define solution parameters
  itrmax = 4
  maxit = 1000
  tol = 1e-6
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
