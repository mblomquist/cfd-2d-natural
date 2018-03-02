! temperature_source2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-02-26 (YYYY-MM-DD)
!
! This subourtine calculates the coefficients used in the solution for
! temperature in the SIMPLER method.
!
! These definitions are defined for a 2D natural convection problem that 
! simulates a vertical plate at i = 0.

subroutine temperature_source2d

  ! Include standard variable header
  include "var2d.dec"

  ! Define internal variables
  integer :: i, j
  real(8) :: Fw, Fe, Fs, Fn, Dw, De, Dn, Ds

  ! Update west boundary
  do j = 2, n
    
	! Set i
	i = 1

	! Update convective terms
    Fw = 0
    Fe = rho*A_x*(u(i+1,j)+u(i,j))/2
    Fs = rho*A_y*(v(i,j)+v(i-1,j))/2
    Fn = rho*A_y*(v(i,j+1)+v(i-1,j+1))/2

    ! Update diffusion terms
    Dw = (k_const/Cp)*dy/dx
    De = (k_const/Cp)*dy/dx
    Ds = (k_const/Cp)*dx/dy
    Dn = (k_const/Cp)*dx/dy

    ! Update energy equation coefficients :: hybrid
	Aw_T(i,j) = 0
	Ae_T(i,j) = max(-Fe,(De-Fe/2),0.0)
	As_T(i,j) = max(Fs,(Ds+Fs/2),0.0)
	An_T(i,j) = max(-Fn,(Dn-Fn/2),0.0)

    ! Update Ap coefficient
	Ap_T(i,j) = Ae_T(i,j)+Aw_T(i,j)+An_T(i,j)+As_T(i,j)+(Fe-Fw)+(Fn-Fs)-Sp_t(i,j)

	! Update b values
	b_T(i,j) = Su_t(i,j)

  end do

  ! Update south boundary :: fixed temperature
  do i = 2, m
    
	! Set j
	j = 1

	! Update energy equation coefficients :: hybrid
	Aw_T(i,j) = 0
	Ae_T(i,j) = 0
	As_T(i,j) = 0
	An_T(i,j) = 0

    ! Update Ap coefficient
	Ap_T(i,j) = -Sp_t(i,j)

	! Update b values
	b_T(i,j) = Su_t(i,j)

  end do

  ! Solve for interior coefficients
  do i = 2,m-1
    do j = 2,n-1

      ! Update convective terms
      Fw = rho*A_x*(u(i,j)+u(i-1,j))/2
      Fe = rho*A_x*(u(i+1,j)+u(i,j))/2
      Fs = rho*A_y*(v(i,j)+v(i-1,j))/2
      Fn = rho*A_y*(v(i,j+1)+v(i-1,j+1))/2

      ! Update diffusion terms
      Dw = (k_const/Cp)*dy/dx
      De = (k_const/Cp)*dy/dx
      Ds = (k_const/Cp)*dx/dy
      Dn = (k_const/Cp)*dx/dy

      ! Update energy equation coefficients :: hybrid
	  Aw_T(i,j) = max(Fw,(Dw+Fw/2),0.0)
	  Ae_T(i,j) = max(-Fe,(De-Fe/2),0.0)
	  As_T(i,j) = max(Fs,(Ds+Fs/2),0.0)
	  An_T(i,j) = max(-Fn,(Dn-Fn/2),0.0)

      ! Update Ap coefficient
	  Ap_T(i,j) = Ae_T(i,j)+Aw_T(i,j)+An_T(i,j)+As_T(i,j)+(Fe-Fw)+(Fn-Fs)-Sp_t(i,j)

	  ! Update b values
	  b_T(i,j) = Su_t(i,j)

    end do
  end do

  return

end subroutine temperature_source2d
