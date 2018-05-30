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

  ! Define geometry variables
  length = 1  ! 1 meter long
  width = 1     ! 1 meter wide
  depth = 1     ! 1 meter deep

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

  ! Define Non-Dimensional parameters
  Ra = 1.4e2

  ! Calculate Non-Dimensional parameters
  Pr = nu / alpha
  Gr = Ra / Pr
  Re = (Gr/10)**(0.5)

  ! Calculate delta temperature
  delta_T = Ra * alpha * nu / g / beta

  ! Calculate u0
  u0 = (g*beta*delta_T)**(0.5)

  ! Define dimensionless temperature at boundaries
  T_w = 1
  T_e = 0
  T_s = 0
  T_n = 0

  ! Define solution parameters
  itrmax = 100
  maxit = 1000
  solver_tol = 1e-6
  simpler_tol = 1e-6
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
