! initialize2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-03-19 (YYYY-MM-DD)
!

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
  Re = 7.397e4
  Pr = 7.127e-1
  Ra = 2.417e5

  rho = 1.276e0
  mu = 1.725e-5
  k_const = 2.435e-2
  Cp = 1.006e3
  beta = 3.663e-3

  ! Define V-Velocity boundary conditions
  v_west = 0
  v_north = 0
  v_south = 0

  ! Define Temperature Boundary conditions
  T_h = 1       ! Inlet fluid temperature
  T_c = 0     ! Heat flux from south

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
