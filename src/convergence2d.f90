! convergence2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-04-03 (YYYY-MM-DD)
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
   v_momentum_residual = 0

  ! Compute the momentum residuals for u
  !call velocity_source2d("u")
  do i = 2,m-1
    do j = 2,n-2
      u_momentum_residual(i,j) = abs(As_u(i,j)*u(i,j-1)+Aw_u(i,j)*u(i-1,j)+Ae_u(i,j)*u(i+1,j)+An_u(i,j)*u(i,j+1)+b_u(i,j)-Ap_u(i,j)*u(i,j))
    end do
  end do

  ! Compute the momentum residuals for v
  !call velocity_source2d("v")
  do i = 2,m-1
    do j = 2,n-1
      v_momentum_residual(i,j) = abs(As_v(i,j)*v(i,j-1)+Aw_v(i,j)*v(i-1,j)+Ae_v(i,j)*v(i+1,j)+An_v(i,j)*v(i,j+1)+b_v(i,j)-Ap_v(i,j)*v(i,j))
    end do
  end do

  R_u = sum(u_momentum_residual)/(m*n)
  R_v = sum(v_momentum_residual)

  if (itr .eq. 1) then
    R_u = 1
    R_v = 1
    
    if (R_u .le. 1) then
      R_u0 = 1
    else
      R_u0 = R_u
    end if

    if (R_v .le. 1) then
      R_v0 = 1
    else
      R_v0 = R_v
    end if
  else
    R_u = R_u/R_u0
    R_v = R_v/R_v0
  end if

  return

end subroutine convergence2d
