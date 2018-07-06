! temperature2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-04-03 (YYYY-MM-DD)
!

subroutine temperature2d

  ! Pull in standard variable header
  include "var2d.dec"

  ! Internal variables
  integer :: i

  ! Initialize Temperature Field
  T = 0

  ! Set temperature source terms
  Su_T = 0
  Sp_T = 0

  ! North boundary source terms :: Wall
  Su_T(:, n-1) = 2/Re/Pr*T_n
  Sp_T(:, n-1) = -2/Re/Pr

  ! West boundary source terms :: Wall
  Su_T(1, :) = 2/Re/Pr*T_w
  Sp_T(1, :) = -2/Re/Pr

  ! East boundary source terms :: Wall
  Su_T(m-1, :) = 2/Re/Pr*T_e
  Sp_T(m-1, :) = -2/Re/Pr

  ! South boundary source terms :: Wall
  Su_T(:, 1) = 2/Re/Pr*T_s
  Sp_T(:, 1) = -2/Re/Pr

  ! Set Initial Temperature
  do i = 1,m-1
    T(i,:) = 1.0/i
  end do

  return

end subroutine temperature2d
