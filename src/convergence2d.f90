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
  do i = 1,m
    do j = 1,n-1

      if (i .eq. 1) then
        if (j .eq. 1) then
          u_momentum_residual(i,j) = abs(Ae_u(i,j)*u(i+1,j)+An_u(i,j)*u(i,j+1)+b_u(i,j)-Ap_u(i,j)*u(i,j))
        elseif (j .eq. n-1) then
          u_momentum_residual(i,j) = abs(As_u(i,j)*u(i,j-1)+Ae_u(i,j)*u(i+1,j)+b_u(i,j)-Ap_u(i,j)*u(i,j))
        else
          u_momentum_residual(i,j) = abs(As_u(i,j)*u(i,j-1)+Ae_u(i,j)*u(i+1,j)+An_u(i,j)*u(i,j+1)+b_u(i,j)-Ap_u(i,j)*u(i,j))
        end if
      elseif (i .eq. m) then
        if (j .eq. 1) then
          u_momentum_residual(i,j) = abs(Aw_u(i,j)*u(i-1,j)+An_u(i,j)*u(i,j+1)+b_u(i,j)-Ap_u(i,j)*u(i,j))
        elseif (j .eq. n-1) then
          u_momentum_residual(i,j) = abs(As_u(i,j)*u(i,j-1)+Aw_u(i,j)*u(i-1,j)+b_u(i,j)-Ap_u(i,j)*u(i,j))
        else
          u_momentum_residual(i,j) = abs(As_u(i,j)*u(i,j-1)+Aw_u(i,j)*u(i-1,j)+An_u(i,j)*u(i,j+1)+b_u(i,j)-Ap_u(i,j)*u(i,j))
        end if
      else
        if (j .eq. 1) then
          u_momentum_residual(i,j) = abs(Aw_u(i,j)*u(i-1,j)+Ae_u(i,j)*u(i+1,j)+An_u(i,j)*u(i,j+1)+b_u(i,j)-Ap_u(i,j)*u(i,j))
        elseif (j .eq. n-1) then
          u_momentum_residual(i,j) = abs(As_u(i,j)*u(i,j-1)+Aw_u(i,j)*u(i-1,j)+Ae_u(i,j)*u(i+1,j)+b_u(i,j)-Ap_u(i,j)*u(i,j))
        else
          u_momentum_residual(i,j) = abs(As_u(i,j)*u(i,j-1)+Aw_u(i,j)*u(i-1,j)+Ae_u(i,j)*u(i+1,j)+An_u(i,j)*u(i,j+1)+b_u(i,j)-Ap_u(i,j)*u(i,j))
        end if
      end if

    end do
  end do

  ! Compute the momentum residuals for v
  !call velocity_source2d("v")
  do i = 2,m-1
    do j = 2,n-1
      v_momentum_residual(i,j) = abs(As_v(i,j)*v(i,j-1)+Aw_v(i,j)*v(i-1,j)+Ae_v(i,j)*v(i+1,j)+An_v(i,j)*v(i,j+1)+b_v(i,j)-Ap_v(i,j)*v(i,j))
    end do
  end do

  R_u = sum(u_momentum_residual)
  R_v = sum(v_momentum_residual)

  return

end subroutine convergence2d
