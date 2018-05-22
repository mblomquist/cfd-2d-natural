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
  length = 1	  ! 1 meter long
  width = 1     ! 1 meter wide
  depth = 1     ! 1 meter deep

  ! Define media variables
  Pr = 7.000e0  ! Ref :: Kimura, Bejan 1983
  Ra = 1.400e5  ! Ref :: Kimura, Bejan 1983
  Re = Ra*Pr    ! Grashoff Number for Natural Convection

  rho = 9.970e2
  mu = 1.070e-3
  k_const = 6.400e-1
  Cp = 4.186e3

  alpha = k_const/Cp/rho
  beta = 6.9e-5

  ! Define high and low temperature
  T_h = 343     ! High temperature wall
  T_c = 293     ! Low temperature wall
  delta_T = T_h - T_c

  ! Define dimensionless temperature at boundaries
  T_w = 1
  T_e = 0
  T_s = 0
  T_n = 0

  ! Define solution parameters
  itrmax = 3
  maxit = 10000
  solver_tol = 1e-6
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
