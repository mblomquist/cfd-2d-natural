! temperature_source2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-05-15 (YYYY-MM-DD)
!
! This subourtine calculates the coefficients used in the solution for
! temperature in the SIMPLER method.
!


subroutine temperature_source2d

  ! Include standard variable header
  include "var2d.dec"

  ! Define internal variables
  integer :: i, j
  real(8) :: Fw, Fe, Fs, Fn, Dw, De, Dn, Ds

  ! Solve for source coefficients
  do i = 1,m-1
    do j = 1,n-1

      ! Update convective terms
      Fw = rho*dy*u(i,j)
      Fe = rho*dy*u(i+1,j)
      Fs = rho*dx*v(i,j)
      Fn = rho*dx*v(i,j+1)

      ! Update diffusion terms
      Dw = k_const/Cp*dy/dx/Re/Pr
      De = k_const/Cp*dy/dx/Re/Pr
      Ds = k_const/Cp*dx/dy/Re/Pr
      Dn = k_const/Cp*dx/dy/Re/Pr

	    ! Compute Coefficients - Power Law Differening Scheme
	    Aw_T(i,j) = Dw*max(0.0,(1-0.1*abs(Fw/Dw))**5)+max(Fw,0.0)
	    Ae_T(i,j) = De*max(0.0,(1-0.1*abs(Fe/De))**5)+max(-Fe,0.0)
	    As_T(i,j) = Ds*max(0.0,(1-0.1*abs(Fs/Ds))**5)+max(Fs,0.0)
	    An_T(i,j) = Dn*max(0.0,(1-0.1*abs(Fn/Dn))**5)+max(-Fn,0.0)

  	  ! Check sourth node
  	  if (j .eq. 1) then
  	    As_T(i,j) = 0
	    end if

	    ! Check north node
	    if (j .eq. n-1) then
  	    An_T(i,j) = 0
      end if

  	  ! Update Ap coefficient
  	  Ap_T(i,j) = Ae_T(i,j)+Aw_T(i,j)+An_T(i,j)+As_T(i,j)-Sp_T(i,j)

  	  ! Update b values
  	  b_T(i,j) = Su_T(i,j)

    end do
  end do

  ! Update west coefficients :: Symmetry
  Aw_T(1,2:n-1) = 0
  Ae_T(1,2:n-1) = 1
  As_T(1,2:n-1) = 0
  An_T(1,2:n-1) = 0

  Ap_T(1,2:n-1) = 1
  b_T(1,2:n-1) = 0

  ! Update east coefficients :: Symmetry
  Aw_T(m-1,2:n-1) = 1
  Ae_T(m-1,2:n-1) = 0
  As_T(m-1,2:n-1) = 0
  An_T(m-1,2:n-1) = 0

  Ap_T(m-1,2:n-1) = 1
  b_T(m-1,2:n-1) = 0

  return

end subroutine temperature_source2d
