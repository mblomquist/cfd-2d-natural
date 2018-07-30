! output_results2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-04-15 (YYYY-MM-DD)
!
! This subroutine generates data files for the pressure, temperature, and velocity fields.

subroutine output_results2d

  ! Include standard variable header
  include "var2d.dec"

  ! Define internal variables
  integer :: i, j

  ! Write field data to file
  open(unit=2, file="output/output2d_results.dat")

  write(2, *), 'TITLE = "2D Natual Convection Field Data"'
  write(2, *), 'VARIABLES = "X", "Y", "Pressure", "u-velocity", "v-velocity", "Temperature"'
  write(2, *), 'ZONE I=100, J=100, DATAPACKING=POINT'

  do j = 1, n-1
    do i = 1, m-1
      write (2,*), length*i*dx, length*j*dy, P(i,j), (u(i+1,j)+u(i,j))/2., (v(i,j)+v(i,j+1))/2., T(i,j)
    end do
  end do

  close(2)

  ! Write residual vector data to file
  open(unit=5, file="output/terminal_data.dat")
  write (5, '("Grid size: ", 5I, 1X, 5I, /)'), m, n
  write (5, '("Rayleigh Number: ", E15.4, /)', advance="no"), Ra
  write (5, '("Prandtl Number: ", E15.4, /)', advance="no"), Pr
  write (5, '("delta_T: ", E15.4, /)', advance="no"), delta_T
  write (5, '("SIMPLER Algorithm Duration:", E15.4, /)', advance="no"), end_time-start_time
  close(5)


  return

end subroutine output_results2d
