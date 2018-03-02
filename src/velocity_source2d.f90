! velocity_solve2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-02-26 (YYYY-MM-DD)
!
! This subroutine updates the source terms for the solution of the momentum
! equations in the SIMPLER algorithm.
!
! These definitions are defined for a 2D natural convection problem that 
! simulates a vertical plate at i = 0.

subroutine velocity_source2d(direction)

  ! Include standard variable header
  include "var2d.dec"

  ! Define input variables
  character :: direction

  ! Define internal variables
  integer :: i, j
  real(8) :: Fw, Fe, Fs, Fn, Dw, De, Dn, Ds

  ! u-velocity update loop
  if (direction .eq. "u") then

    ! Calculate west coefficients :: no slip condition
	do j = 1,n

	  ! Set i
	  i = 1

	  ! Set coefficients
	  Aw_u(i,j) = 0
	  Ae_u(i,j) = 0
	  As_u(i,j) = 0
	  An_u(i,j) = 0

	  ! Update Ap coefficient
	  Ap_u(i,j) = -Sp_u(i,j)

	  ! Update b values
	  b_u(i,j) = Su_u(i,j)

	end do

    ! Calculate interior coefficients
    do i = 2,m-1
      do j = 2,n-1

        ! Update convection terms
		Fw = rho*A_x*(u_star(i,j)+u_star(i-1,j))/2
		Fe = rho*A_x*(u_star(i+1,j)+u_star(i,j))/2
		Fs = rho*A_y*(v_star(i,j)+v_star(i-1,j))/2
		Fn = rho*A_y*(v_star(i,j+1)+v_star(i-1,j+1))/2

        ! Update diffusion terms
		Dw = mu*dy/dx
        De = mu*dy/dx
        Ds = mu*dy/dx
        Dn = mu*dy/dx

		! Compute Coefficients - Hybrid Differening Scheme
		Aw_u(i,j) = max(Fw,(Dw+Fw/2),0.0)
		Ae_u(i,j) = max(-Fe,(De-Fe/2),0.0)
		As_u(i,j) = max(Fs,(Ds+Fs/2),0.0)
		An_u(i,j) = max(-Fn,(Dn-Fn/2),0.0)

		! Update Ap coefficient
		Ap_u(i,j) = Ae_u(i,j)+Aw_u(i,j)+An_u(i,j)+As_u(i,j)+(Fe-Fw)+(Fn-Fs)-Sp_u(i,j)

		! Update b values
		b_u(i,j) = Su_u(i,j)+A_x*(P_star(i-1,j)-P_star(i,j))

      end do
    end do

  end if

  ! v-velocity update loop
  if (direction .eq. "v") then

    ! Calculate west coefficients
	do j = 1,n

	  ! Set i
	  i = 1

	  ! Update coefficient terms
      Fw = 0
      Fe = rho*dy*(u_star(i+1,j)+u_star(i+1,j-1))/2
      Fs = rho*dx*(v_star(i,j)+v_star(i,j-1))/2
      Fn = rho*dy*(v_star(i,j+1)+v_star(i,j))/2

	  ! Update diffusion terms
      Dw = mu*dy/dx
      De = mu*dy/dx
      Ds = mu*dy/dx
      Dn = mu*dy/dx

	  ! Compute Coefficients - Hybrid Differening Scheme
	  Aw_v(i,j) = 0
	  Ae_v(i,j) = max(-Fe,(De-Fe/2),0.0)
	  As_v(i,j) = max(Fs,(Ds+Fs/2),0.0)
	  An_v(i,j) = max(-Fn,(Dn-Fn/2),0.0)

      ! Update Ap coefficient
      Ap_v(i,j) = Ae_v(i,j)+Aw_v(i,j)+An_v(i,j)+As_v(i,j)+(Fe-Fw)+(Fn-Fs)-Sp_v(i,j)-2

      ! Update b values
      b_v(i,j) = Su_v(i,j)+g*beta*(T(i,j)-T_inf)-rho*g+A_y*(P_star(i,j)-P_star(i,j-1))
	
	end do

    ! Calculate interior coefficients
    do i = 2,m-1
      do j = 2,n-1

        Fw = rho*dy*(u_star(i,j)+u_star(i,j-1))/2
        Fe = rho*dy*(u_star(i+1,j)+u_star(i+1,j-1))/2
        Fs = rho*dx*(v_star(i,j)+v_star(i,j-1))/2
        Fn = rho*dy*(v_star(i,j+1)+v_star(i,j))/2

        Dw = mu*dy/dx
        De = mu*dy/dx
        Ds = mu*dy/dx
        Dn = mu*dy/dx

		! Compute Coefficients - Hybrid Differening Scheme
		Aw_v(i,j) = max(Fw,(Dw+Fw/2),0.0)
		Ae_v(i,j) = max(-Fe,(De-Fe/2),0.0)
		As_v(i,j) = max(Fs,(Ds+Fs/2),0.0)
		An_v(i,j) = max(-Fn,(Dn-Fn/2),0.0)

        ! Update Ap coefficient
        Ap_v(i,j) = Ae_v(i,j)+Aw_v(i,j)+An_v(i,j)+As_v(i,j)+(Fe-Fw)+(Fn-Fs)-Sp_v(i,j)

        ! Update b values
        b_v(i,j) = Su_v(i,j)+g*beta*(T(i,j)-T_inf)-rho*g+A_y*(P_star(i,j)-P_star(i,j-1))

      end do
    end do

  end if

  return

end subroutine velocity_source2d
