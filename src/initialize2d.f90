! initialize2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-06-27 (YYYY-MM-DD)
!

subroutine initialize2d

  ! Pull in standard variable header
  include "var2d.dec"

  ! Read Input File ..........
  open(unit = 2, file = "input2d.txt")
  read(2,*)
  read(2,*)
  read(2,*) length, width, depth
  read(2,*)
  read(2,*) g, rho, mu, k_const, Cp, beta
  read(2,*)
  read(2,*) u0, T_h, T_c
  read(2,*)
  read(2,*) itrmax, maxit, solver_tol, simpler_tol, alpha_v, alpha_t, solver
  close(2)

  ! Calculate parameters
  alpha = k_const / Cp / rho
  nu = mu / rho
  delta_T = T_h - T_c

  ! Calculate dimensionless numbers
  Ra = g*beta*delta_T*length**3.0/alpha/nu
  Pr = nu/alpha
  Gr = Ra/Pr
  Ga = g*length**3.0/nu**2.0
  Re = u0*length/nu


  ! Define dimensionless temperature at boundaries
  T_w = (T_c-T_c)/delta_T
  T_e = (T_c-T_c)/delta_T
  T_s = (T_h-T_c)/delta_T
  T_n = (T_c-T_c)/delta_T

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
