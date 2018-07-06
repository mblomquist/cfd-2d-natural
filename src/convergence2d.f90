! convergence2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-06-27 (YYYY-MM-DD)
!
! This subroutine computes the error of the SIMPLER solution based on the
! residuals for the u- and v-momentum equations.
!

subroutine convergence2d(itr)

  ! Include variable header
  include "var2d.dec"

  ! Define internal variables
  integer :: i, j, itr
  real(8) :: e_temp

  R_e = 0.0

  e_5 = e_4
  e_4 = e_3
  e_3 = e_2
  e_2 = e_1

  do i = 2, m-2
    do j = 2, n-2

      if (R_e .le. b_p(i,j)) then
        e_1 = e_1+abs(b_p(i,j))
      end if

    end do
  end do

  e_1 = e_1/((m-3.0)*(n-3.0))


  R_e = abs((e_1-(e_2+e_3+e_4+e_5)/4.0)/((e_2+e_3+e_4+e_5)/4.0))

  return

end subroutine convergence2d
