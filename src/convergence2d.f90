! convergence2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-02-25 (YYYY-MM-DD)
!
! This subroutine computes the error of the SIMPLER solution based on the
! residuals for the u- and v-momentum equations.
!
subroutine convergence2d(itr)

  ! Include variable header
  include "var2d.dec"

  ! Define internal variables
  integer :: i, j
  real(8), dimension(m,n) :: u_momentum_residual, v_momentum_residual

   R_u = 0
   R_v = 0
   u_momentum_residual = 0

   do i = 3,m-2
     do j = 3,n-2
       u_momentum_residual = abs(rho*(u(i,j)-u(i-1,j))*A_x+rho*A_y*(v(i,j)-v(i,j-1)))
     end do
   end do

  ! Compute the momentum residuals for u
  !call velocity_source2d("u")
  !do i = 1,m-1
  !  do j = 1,n-1
      !u_momentum_residual(i,j) = abs(As_u(i,j)*u(i,j-1)+Aw_u(i,j)*u(i-1,j)+Ae_u(i,j)*u(i+1,j)+An_u(i,j)*u(i,j+1)+b_u(i,j)-Ap_u(i,j)*u(i,j))
      !u_momentum_residual(i,j) = abs(u(i,j)-u(i-1,j)+v(i,j)-v(i,j-1))
  !  end do
  !end do

  ! Compute the momentum residuals for v
  !call velocity_source2d("v")
  !do i = 3,m-2
  !  do j = 3,n-2
  !    v_momentum_residual(i,j) = abs(As_v(i,j)*v(i,j-1)+Aw_v(i,j)*v(i-1,j)+Ae_v(i,j)*v(i+1,j)+An_v(i,j)*v(i,j+1)+b_v(i,j)-Ap_v(i,j)*v(i,j))
  !  end do
  !end do

  R_u = sum(u_momentum_residual)/(m*n)
  !R_v = sum(v_momentum_residual)
  R_v = R_u

  !if (itr .eq. 1) then
  !  R_u0 = R_u
	!  R_v0 = R_v
  !else
  !  R_u = R_u/R_u0
  !  R_v = R_v/R_v0
  !end if

  return

end subroutine convergence2d
