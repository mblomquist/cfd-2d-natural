! output_results2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-04-03 (YYYY-MM-DD)
!
! This subroutine generates data files for the pressure, temperature, and velocity fields.

subroutine output_results2d

  ! Include standard variable header
  include "var2d.dec"

  ! Define internal variables
  integer :: i, j

  ! Write field data to file
  open(unit=2, file="output/field_data.dat")

  write (2, '("# vtk DataFile Version 2.0", /)', advance="no")
  write (2, '("field data: pressure, temperature", /)', advance="no")
  write (2, '("ASCII", /)', advance="no")
  write (2, '("DATASET STRUCTURED_POINT", /)', advance="no")
  write (2, '("DIMENSIONS ", 5I, 1X, 5I, /)', advance="no"), M-1, N-1
  write (2, '("ORIGIN ", 5D, 1X, 5D, /")', advance="no"), DX/2, DY/2
  write (2, '("SPACING ", 5D, 1X, 5D, /)', advance="no"), DX, DY

  do j = 1, n-1
    do i = 1, m-1
      write (2, '(E15.4, 1x, E15.4, /)', advance="no"), P(i,j), T(i,j)
    end do
  end do

  close(2)

  ! Write u-velocity data
  open(unit=3, file="output/u_velocity_data.dat")

  write (3, '("# vtk DataFile Version 2.0", /)', advance="no")
  write (3, '("u_velocity_data", /)', advance="no")
  write (3, '("ASCII", /)', advance="no")
  write (3, '("DATASET STRUCTURED_POINT", /)', advance="no")
  write (3, '("DIMENSIONS ", 5I, 1X, 5I, /)', advance="no"), M, N-1
  write (3, '("ORIGIN ", 5D, 1X, 5D, /")', advance="no"), 0.0, DY/2
  write (3, '("SPACING ", 5D, 1X, 5D, /)', advance="no"), DX, DY

  do j = 1, n-1
    do i = 1, m
      write (3, '(E15.4, /)', advance="no"), u(i,j)
    end do
  end do

  close(3)

  ! Write v-velocity data
  open(unit=4, file="output/v_velocity_data.dat")

  write (4, '("# vtk DataFile Version 2.0", /)', advance="no")
  write (4, '("v_velocity_data", /)', advance="no")
  write (4, '("ASCII", /)', advance="no")
  write (4, '("DATASET STRUCTURED_POINT", /)', advance="no")
  write (4, '("DIMENSIONS ", 5I, 1X, 5I, /)', advance="no"), M-1, N
  write (4, '("ORIGIN ", 5D, 1X, 5D, /")', advance="no"), DX/2, 0.0
  write (4, '("SPACING ", 5D, 1X, 5D, /)', advance="no"), DX, DY

  do j = 1, n
    do i = 1, m-1
      write (4, '(E15.4, /)', advance="no"), v(i,j)
    end do
  end do

  close(4)

  return

end subroutine output_results2d
