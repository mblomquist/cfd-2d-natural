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

  R_e = 0.0

  do i = 2, m-2
    do j = 2, n-2

      if (R_e .le. b_p(i,j)) then
        e_max_new = abs(b_p(i,j))
      end if

    end do
  end do

  if (itr .eq. 1) then
    e_max_old = 1
  end if

  R_e = abs((e_max_new-e_max_old)/e_max_old)

  e_max_old = e_max_new

  return

end subroutine convergence2d
