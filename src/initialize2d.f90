! initialize2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-04-03 (YYYY-MM-DD)
!

subroutine initialize2d

  ! Pull in standard variable header
  include "var2d.dec"

  ! Read Input File ..........
  ! ..........................
  ! ..........................

  ! Define geometry variables
  length = 1	! 1 meter long
  width = 1     ! 1 meter wide
  depth = 1     ! 1 meter deep

  ! Define media variables
  Re = 7.397e4
  Pr = 7.127e-1
  Ra = 2.417e5

  rho = 1.276e0
  mu = 1.725e-5
  k_const = 2.435e-2
  Cp = 1.006e3

  alpha = k_const/Cp/rho
  beta = 3.663e-3

  ! Define high and low temperature
  T_h = 373     ! High temperature wall
  T_c = 100     ! Low temperature wall
  delta_T = T_h - T_c

  ! Define dimensionless temperature at boundaries
  T_w = T_c
  T_e = T_c
  T_s = T_h
  T_n = T_c

  ! Define solution parameters
  itrmax = 10
  maxit = 1000
  solver_tol = 1e-9
  simpler_tol = 1e-6

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
