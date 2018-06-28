! velocity_solve2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-06-27 (YYYY-MM-DD)
!
! This subroutine updates the source terms for the solution of the momentum
! equations in the SIMPLER algorithm.
!

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

    ! Calculate interior coefficients
    do i = 2,m-1
      do j = 1,n-1

        ! Update convection terms
		    Fw = rho*dy*(u_star(i,j)+u_star(i-1,j))/2
		    Fe = rho*dy*(u_star(i+1,j)+u_star(i,j))/2
		    Fs = rho*dx*(v_star(i,j)+v_star(i-1,j))/2
		    Fn = rho*dx*(v_star(i,j+1)+v_star(i-1,j+1))/2

        ! Update diffusion terms
        Dw = dy/dx/Re
        De = dy/dx/Re
        Ds = dx/dy/Re
        Dn = dx/dy/Re

		    ! Compute Coefficients - Power Law Differening Scheme
		    Aw_u(i,j) = Dw*max(0.0,(1-0.1*abs(Fw/Dw))**5)+max(Fw,0.0)
		    Ae_u(i,j) = De*max(0.0,(1-0.1*abs(Fe/De))**5)+max(-Fe,0.0)
		    As_u(i,j) = Ds*max(0.0,(1-0.1*abs(Fs/Ds))**5)+max(Fs,0.0)
		    An_u(i,j) = Dn*max(0.0,(1-0.1*abs(Fn/Dn))**5)+max(-Fn,0.0)

		    ! Check South / North Nodes
		    if (j .eq. 1) then
		      As_u(i,j) = 0
		    elseif (j .eq. n-1) then
		      An_u(i,j) = 0
		    end if

		    ! Update Ap coefficient
		    Ap_u(i,j) = Ae_u(i,j)+Aw_u(i,j)+An_u(i,j)+As_u(i,j)-Sp_u(i,j)*dx*dy

        if (Ap_u(i,j) .eq. 0) then
          Aw_u(i,j) = 0
          Ae_u(i,j) = 0
          As_u(i,j) = 0
          An_u(i,j) = 0

          Ap_u(i,j) = 1

          print *, "False Diffusion (u-velocity)@:", i,j
        end if

		    ! Update b values
		    b_u(i,j) = Su_u(i,j)*dx*dy

      end do
    end do

    ! West Boundary Coefficients :: Wall
    Aw_u(1,:) = 0
    Ae_u(1,:) = 0
    As_u(1,:) = 0
    An_u(1,:) = 0
    Ap_u(1,:) = 1
    b_u(1,:) = 0

    ! East Bounday Coefficients :: Wall
    Aw_u(m,:) = 0
    Ae_u(m,:) = 0
    As_u(m,:) = 0
    An_u(m,:) = 0
    Ap_u(m,:) = 1
    b_u(m,:) = 0

  end if

  ! v-velocity update loop
  if (direction .eq. "v") then

    ! Calculate interior source terms
	  do j = 2, n-1
	    do i = 1,m-1

	      ! Update convection terms
		    Fw = rho*dy*(u_star(i,j-1)+u_star(i,j))/2
		    Fe = rho*dy*(u_star(i+1,j-1)+u_star(i+1,j))/2
		    Fs = rho*dx*(v_star(i,j-1)+v_star(i,j))/2
		    Fn = rho*dx*(v_star(i,j)+v_star(i,j+1))/2

        ! Update diffusion terms
        Dw = dy/dx/Re
        De = dy/dx/Re
        Ds = dx/dy/Re
        Dn = dx/dy/Re

		    ! Compute Coefficients - Power Law Differening Scheme
		    Aw_v(i,j) = Dw*max(0.0,(1-0.1*abs(Fw/Dw))**5)+max(Fw,0.0)
		    Ae_v(i,j) = De*max(0.0,(1-0.1*abs(Fe/De))**5)+max(-Fe,0.0)
		    As_v(i,j) = Ds*max(0.0,(1-0.1*abs(Fs/Ds))**5)+max(Fs,0.0)
		    An_v(i,j) = Dn*max(0.0,(1-0.1*abs(Fn/Dn))**5)+max(-Fn,0.0)

		    ! Check South / North Nodes
		    if (i .eq. 1) then
		      Aw_v(i,j) = 0
		    elseif (i .eq. m-1) then
		      Ae_v(i,j) = 0
		    end if

		    ! Update Ap coefficient
		    Ap_v(i,j) = Ae_v(i,j)+Aw_v(i,j)+An_v(i,j)+As_v(i,j)-Sp_v(i,j)*dx*dy

        if (Ap_v(i,j) .eq. 0) then
          Aw_v(i,j) = 0
          Ae_v(i,j) = 0
          As_v(i,j) = 0
          An_v(i,j) = 0
          Ap_v(i,j) = 1

          print *, "False Diffusion (v-velocity)@:", i,j
        end if

		    ! Update b values
		    b_v(i,j) = Su_v(i,j)*dx*dy-2.0*dx*dy+Gr/Re/Re*((T(i,j)+T(i,j+1))/2.0)*dx*dy

	    end do
	  end do

    ! South Boundary :: Wall
    Aw_v(:,1) = 0
    Ae_v(:,1) = 0
    As_v(:,1) = 0
    An_v(:,1) = 0
    Ap_v(:,1) = 1
    b_v(:,1) = 0

    ! North Boundary :: Wall
    Aw_v(:,n) = 0
    Ae_v(:,n) = 0
    As_v(:,n) = 0
    An_v(:,n) = 0
    Ap_v(:,n) = 1
    b_v(:,n) = 0

  end if

  return

end subroutine velocity_source2d
