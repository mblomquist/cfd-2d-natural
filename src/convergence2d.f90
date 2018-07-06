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
  R_t = 0.0

  e_1 = 0.0
  e_3 = 0.0

  do i = 2, m-2
    do j = 2, n-2

      if (e_1 .le. abs(b_p(i,j))) then
        e_1 = abs(b_p(i,j))
      end if

      e_temp = abs(Ap_T(i,j)*T(i,j)-(An_T(i,j+1)*T(i,j)+As_T(i,j)*T(i,j-1)+Aw_T(i,j)*T(i-1,j)+Ae_T(i,j)*T(i+1,j)))

      if (e_3 .le. e_temp) then
        e_3 = e_temp
      end if

    end do
  end do

  R_e = abs(e_1 - e_2)/abs(e_1)
  R_t = abs(e_3 - e_4)/abs(e_4)

  e_2 = e_1
  e_4 = e_3

  return

end subroutine convergence2d
