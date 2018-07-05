! temperature_source2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-06-27 (YYYY-MM-DD)
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
  do i = 2,m-2
    do j = 1,n-1

      ! Update convective terms
      Fw = dy*u(i,j)
      Fe = dy*u(i+1,j)
      Fs = dx*v(i,j)
      Fn = dx*v(i,j+1)

      ! Update diffusion terms
      Dw = dy/dx/(Pr*Ra)**(0.5)*Pr**(3./2.)
      De = dy/dx/(Pr*Ra)**(0.5)*Pr**(3./2.)
      Ds = dx/dy/(Pr*Ra)**(0.5)*Pr**(3./2.)
      Dn = dx/dy/(Pr*Ra)**(0.5)*Pr**(3./2.)

	    ! Compute Coefficients - Power Law Differening Scheme
	    Aw_T(i,j) = Dw*max(0.0,(1-0.1*abs(Fw/Dw))**5)+max(Fw,0.0)
	    Ae_T(i,j) = De*max(0.0,(1-0.1*abs(Fe/De))**5)+max(-Fe,0.0)
	    As_T(i,j) = Ds*max(0.0,(1-0.1*abs(Fs/Ds))**5)+max(Fs,0.0)
	    An_T(i,j) = Dn*max(0.0,(1-0.1*abs(Fn/Dn))**5)+max(-Fn,0.0)

  	  ! Update Ap coefficient
  	  Ap_T(i,j) = Ae_T(i,j)+Aw_T(i,j)+An_T(i,j)+As_T(i,j)-Sp_T(i,j)*dx*dy

      if (Ap_T(i,j) .eq. 0) then
        Aw_T(i,j) = 0
        Ae_T(i,j) = 0
        As_T(i,j) = 0
        An_T(i,j) = 0

        Ap_T(i,j) = 1.0

        print *, "False Diffusion (temperature)@:", i,j
      end if

  	  ! Update b values
  	  b_T(i,j) = Su_T(i,j)*dx*dy

    end do
  end do

  ! Compute Coefficients - North Wall
  Aw_T(:,n-1) = 0
  Ae_T(:,n-1) = 0
  As_T(:,n-1) = 1
  An_T(:,n-1) = 0

  Ap_T(:,n-1) = 1
  b_T(:,n-1) = 0

  ! Compute Coefficients - South Wall
  Aw_T(:,1) = 0
  Ae_T(:,1) = 0
  As_T(:,1) = 0
  An_T(:,1) = 1

  Ap_T(:,1) = 1
  b_T(:,1) = 0

  ! Compute Coefficients - East Wall
  Aw_T(m-1,:) = 0
  Ae_T(m-1,:) = 0
  As_T(m-1,:) = 0
  An_T(m-1,:) = 0

  Ap_T(m-1,:) = 1
  b_T(m-1,:) = 0

  ! Compute Coefficients - West Wall
  Aw_T(1,:) = 0
  Ae_T(1,:) = 0
  As_T(1,:) = 0
  An_T(1,:) = 0

  Ap_T(1,:) = 1
  b_T(1,:) = 1

  return

end subroutine temperature_source2d
