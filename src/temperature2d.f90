! temperature2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-03-09 (YYYY-MM-DD)
!
subroutine temperature2d

  ! Pull in standard variable header
  include "var2d.dec"

  ! Initialize Temperature Field
  T = 0

  ! Set temperature source terms
  Su_T = 0
  Sp_T = 0

  ! North boundary source terms :: Inlet (fluid temperature)
  Su_T(:, n-1) = 2*(k_const/Cp)*dx/dy*T_north
  Sp_T(:, n-1) = -2*(k_const/Cp)*dx/dy

  ! West boundary source terms ::
  !Su_T(1, :) = 2*(k_const/Cp)*dx/dy*T_north
  !Sp_T(1, :) = -2*(k_const/Cp)*dx/dy

  !Su_T(m-1, :) = 2*(k_const/Cp)*dx/dy*T_north
  !Sp_T(m-1, :) = -2*(k_const/Cp)*dx/dy

  ! South boundary source terms :: heat flux
  Su_T(:, 1) = 2*(k_const/Cp)*dx/dy*T_south
  Sp_T(:, 1) = -2*(k_const/Cp)*dx/dy

  return

end subroutine temperature2d
