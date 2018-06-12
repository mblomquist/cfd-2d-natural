! 2D CFD Solver
!
! Written by Matt Blomquist
! Last Update: 2018-05-14 (YYYY-MM-DD)
!
! This program solves a two-dimensional, steady state CFD problem
! using the SIMPLER method.

program main2d

  implicit none

  ! Define time variables
  real(8) :: start_time, end_time

  ! Include standard variable header
  !   Standard variable header establishes parameters and global values
  !   subroutines will call these modify these values.
  include "var2d.dec"

  ! Initialize Problem
  !   The initialize subroutine modifies global vaules to incorporate
  !   boundary conditions, calculate global values, etc.
  call initialize2d
  print *, 'Problem Initialization Complete.'
  print *, ''
  print *, 'Grid size: ', m, n
  print *, 'Rayleigh Number: ', Ra
  print *, 'Prandtl Number: ', Pr
  print *, 'Grashoff Number: ', Gr
  print *, 'Reynolds Number', Re
  print *, 'Gr/Re2:', Gr/Re/Re
  print *, 'delta_T:', delta_T
  print *, 'length:', length
  print *, 'u0:', u0
  print *, ''

  call cpu_time(start_time)

  ! Start SIMPLER Algorithm
  !   The SIMPLER Algorithm solves the problem using the iterative method.
  print *, 'Starting SIMPLER Algorithm...'
  call simpler2d
  print *, 'SIMPLER Algorithm Complete.'

  call cpu_time(end_time)

  print *, "SIMPLER Algorithm Duration:", end_time-start_time

  ! Output results
  !   The output subroutine writes the final values of P, T, Uvel, Vvel,
  !   and convergance outputs.
  print *, 'Writing results to file.'
  call output_results2d
  print *, 'Program complete.'

! End program
end program main2d
