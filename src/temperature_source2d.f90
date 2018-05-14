! temperature_source2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-04-02 (YYYY-MM-DD)
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
      Fw = dy*u(i,j)
      Fe = dy*u(i+1,j)
      Fs = dx*v(i,j)
      Fn = dx*v(i,j+1)

      ! Update diffusion terms
		  Dw = Pr*dy/dx
      De = Pr*dy/dx
      Ds = Pr*dx/dy
      Dn = Pr*dx/dy

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

	  ! Check west node
	  if (i .eq. 1) then
        Aw_T(i,j) = 0
      end if

	  ! Check east node
	  if (i .eq. m-1) then
        Ae_T(i,j) = 0
  	  end if

  	  ! Update Ap coefficient
  	  Ap_T(i,j) = Ae_T(i,j)+Aw_T(i,j)+An_T(i,j)+As_T(i,j)-Sp_T(i,j)

  	  ! Update b values
  	  b_T(i,j) = Su_T(i,j)

    end do
  end do

  return

end subroutine temperature_source2d
