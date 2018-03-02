! output_results2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-01-19 (YYYY-MM-DD)
!
! This subroutine writes a data file that includes all of the outputs of the
! results.

subroutine output_results2d

  ! Include standard variable header
  include "var2d.dec"

  ! Define internal variables
  integer :: i, j

  open(unit=2, file="output/output_results.txt")

  do j = 1, n-1
    do i = 1, m-1
      write (2, '("Node: ",I4,2x,I4,2x,"Pressure: ",E15.4,2x,"U-Velocity:",E15.4,2x,"V-Velocity: ",E15.4,2x,"Temperature: ",E15.4,/)', advance="no"), i, j, P(i,j), u(i,j), v(i,j), T(i,j)
    end do
  end do

  i = m
  do j = 1, n-1
    write (2, '("Node: ",I4,2x,I4,2x,"Pressure: ",E15.4,2x,"U-Velocity:",E15.4,2x,"V-Velocity: ",E15.4,2x,"Temperature: ",15x,/)', advance="no"), i, j, P(i,j), u(i,j), v(i,j)
  end do

  j = n
  do i = 1, m
    write (2, '("Node: ",I4,2x,I4,2x,"Pressure: ",15x,2x,"U-Velocity:",E15.4,2x,"V-Velocity: ",E15.4,2x,"Temperature: ",15x,/)', advance="no"), i, j, u(i,j), v(i,j)
  end do

  close(2)

  return

end subroutine output_results2d
