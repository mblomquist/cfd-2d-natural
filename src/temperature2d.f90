! temperature2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-04-03 (YYYY-MM-DD)
!

subroutine temperature2d

  ! Pull in standard variable header
  include "var2d.dec"

  ! Initialize Temperature Field
  T = 0

  ! Set temperature source terms
  Su_T = 0
  Sp_T = 0

  ! North boundary source terms :: Wall
  Su_T(:, n-1) = 2*mu*Cp*delta_T*dx/Pr/dy*T_n
  Sp_T(:, n-1) = -2*mu*Cp*delta_T*dx/Pr/dy

  ! West boundary source terms :: Wall
  !Su_T(1, :) = 2*mu*Cp*delta_T*dy/Pr/dx*T_w
  !Sp_T(1, :) = -2*mu*Cp*delta_T*dy/Pr/dx

  ! East boundary source terms :: Wall
  !Su_T(m-1, :) = 0 !2*mu*Cp*delta_T*dy/Pr/dx*T_e
  !Sp_T(m-1, :) = 0 !-2*mu*Cp*delta_T*dy/Pr/dx

  ! South boundary source terms :: Wall
  Su_T(:, 1) = 2*mu*Cp*delta_T*dx/Pr/dy*T_s
  Sp_T(:, 1) = -2*mu*Cp*delta_T*dx/Pr/dy

  return

end subroutine temperature2d
