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
    do j = 2,n-2

      ! Update convective terms
      Fw = rho*dy*u(i,j)
      Fe = rho*dy*u(i+1,j)
      Fs = rho*dx*v(i,j)
      Fn = rho*dx*v(i,j+1)

      ! Update diffusion terms
      Dw = dy/dx/Re/Pr
      De = dy/dx/Re/Pr
      Ds = dx/dy/Re/Pr
      Dn = dx/dy/Re/Pr

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

  ! Compute Coefficients - West Wall
  Aw_T(1,:) = 0
  Ae_T(1,:) = 0
  As_T(1,:) = 0
  An_T(1,:) = 0

  Ap_T(1,:) = Sp_T(1,:)
  b_T(1,:) = Su_T(1,:)

  ! Compute Coefficients - East Wall
  Aw_T(m-1,1) = 0
  Ae_T(m-1,1) = 0
  As_T(m-1,1) = 0
  An_T(m-1,1) = 0

  Ap_T(m-1,:) = Sp_T(m-1, 1)
  b_T(m-1,:) = Su_T(m-1, 1)

  ! Compute Coefficients - North Wall
  Aw_T(:,n-1) = 0
  Ae_T(:,n-1) = 0
  As_T(:,n-1) = 0
  An_T(:,n-1) = 0

  Ap_T(:,n-1) = Sp_T(:,n-1)
  b_T(:,n-1) = Su_T(:,n-1)

  ! Compute Coefficients - South Wall
  Aw_T(:,1) = 0
  Ae_T(:,1) = 0
  As_T(:,1) = 0
  An_T(:,1) = 0

  Ap_T(:,1) = Sp_T(:,1)
  b_T(:,1) = Su_T(:,1)


  return

end subroutine temperature_source2d
