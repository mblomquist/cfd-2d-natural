! temperature_source2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-03-09 (YYYY-MM-DD)
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

  ! Update west boundary
  ! Set i
  i = 1

  do j = 1,n-1

    ! Update energy equation coefficients :: hybrid
    Aw_T(i,j) = 0
    Ae_T(i,j) = 1
    As_T(i,j) = 0
    An_T(i,j) = 0

    ! Update Ap coefficient
    Ap_T(i,j) = 1

    ! Update b values
    b_T(i,j) = 0

  end do

  ! Update east boundary
  ! Set i
  i = m-1

  do j = 1,n-1

    ! Update energy equation coefficients :: hybrid
    Aw_T(i,j) = 1
    Ae_T(i,j) = 0
    As_T(i,j) = 0
    An_T(i,j) = 0

    ! Update Ap coefficient
    Ap_T(i,j) = 1

    ! Update b values
    b_T(i,j) = 0

  end do

  ! Solve for interior coefficients
  do i = 2,m-2
    do j = 1,n-1

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

      ! Compute Coefficients - Power Law Differening Scheme
  		Aw_T(i,j) = Dw*max(0.0,(1-0.1*abs(Fw/Dw))**5)+max(Fw,0.0)
  		Ae_T(i,j) = De*max(0.0,(1-0.1*abs(Fe/De))**5)+max(-Fe,0.0)
  		As_T(i,j) = Ds*max(0.0,(1-0.1*abs(Fs/Ds))**5)+max(Fs,0.0)
  		An_T(i,j) = Dn*max(0.0,(1-0.1*abs(Fn/Dn))**5)+max(-Fn,0.0)

  		! Check South / North Nodes
  		if (j .eq. 1) then
  		  As_T(i,j) = 0
  		elseif (j .eq. n-1) then
  		  An_T(i,j) = 0
      elseif (i .eq. 1) then
        Aw_T(i,j) = 0
      elseif (i .eq. m-1) then
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
