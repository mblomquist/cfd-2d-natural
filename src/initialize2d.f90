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

  ! Define length Scale
  length = 0.01     ! (m)
  width = 0.01      ! (m)
  depth = 1.0         ! (m)

  ! Define velocity scale
  u0 = 0.01         ! (m/s)

  ! Define temperature sclae
  T_h = 274
  T_c = 273

  delta_T = T_h - T_c

  ! Define media parameters
  g = 9.81e0
  rho = 1.28e0
  mu = 1.73e-5
  k_const = 2.44e-2
  Cp = 1.01e3
  beta = 3.66e-3

  ! Calculate parameters
  alpha = k_const / Cp / rho
  nu = mu / rho

  ! Calculate dimensionless numbers
  Ra = g*beta*delta_T*length**3.0/alpha/nu
  Pr = nu/alpha
  Gr = Ra/Pr
  Re = u0*length/nu


  ! Define dimensionless temperature at boundaries
  T_w = 1
  T_e = 0
  T_s = 0
  T_n = 0

  ! Define solution parameters
  itrmax = 4
  maxit = 1e8
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
