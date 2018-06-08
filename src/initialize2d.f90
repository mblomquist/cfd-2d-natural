! initialize2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-05-14 (YYYY-MM-DD)
!

subroutine initialize2d

  ! Pull in standard variable header
  include "var2d.dec"

  ! Read Input File ..........
  ! ..........................
  ! ..........................

  ! Define Dimensionless Inputs
  Ra = 1.4e2
  Pr = 7.0e0

  ! Define Dimensionless Temperatures
  T_h = 1
  T_c = 0
  delta_T = 1

  ! Define media variables
  g = 9.81e0
  rho = 9.97e2
  mu = 1.07e-3
  k_const = 6.40e-1
  Cp = 4.19e3
  beta = 6.90e-5

  ! Calculate variables
  nu = mu / rho
  alpha = k_const / rho / Cp

  ! Calculate Dimensionless Values
  Gr = Ra / Pr
  Re = (Gr/10.0)**(0.5)

  ! Calculate Scale Values
  length = ((Ra*nu**2)/(Pr*g*beta*delta_T))**(0.333)
  width = length
  depth = 1

  u0 = Re*nu/length

  ! Define dimensionless temperature at boundaries
  T_w = T_h
  T_e = T_c
  T_s = T_c
  T_n = T_c

  ! Define solution parameters
  itrmax = 7
  maxit = 1e6
  solver_tol = 1e-9
  simpler_tol = 1e-2
  relax = 1.0

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
