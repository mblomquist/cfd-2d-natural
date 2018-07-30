! solver2d_bicg
!
! Written by Matt Blomquist
! LAxst Update: 2018-07-20 (YYYY-MM-DD)
!
! This program solves a three-dimensional finite volume discretization problem
! using the bi-conjugate gradients stabilized algorithm.
!
! Definition of input arguments
! Inputs:
!   Ab, As, Aw, Ap, Ae, An, At :: These arrays represent the coefficients for adjacent nodes
!   b :: This array represents the right-hand side of the equation Ax=b
!   phi :: This value represents the Axppropriate solution array (pressure, velocity, temperature)
!   m, n, l :: These values represent the number of nodes for i and j for the phi value
!   tol :: represents the solution tolerance
!   maxit :: represents the maximum number of iterations of the bicg Algorithm
!
! Outputs:
!   phi :: on exit, this value contains the updated solution
!   maxit :: on exit, this value contains the number of iterations of the bicg algorithm
!   tol :: on exit, this value represents the normalized residual

subroutine solver2d_bicg(As, Aw, Ap, Ae, An, b, phi, m, n, tol, maxit)

  ! Define implicit
  implicit none

  ! Include mkl functions
  include "mkl.fi"

  ! Define input variables
  integer, intent(in) :: m, n
  integer, intent(in) :: maxit
  real(8), intent(in) :: tol
  real(8), dimension(m,n), intent(in) :: As, Aw, Ap, Ae, An, b
  real(8), dimension(m,n), intent(inout) :: phi

  ! Define internal variables
  integer :: i, j, k, itr
  real(8), dimension(m*n,5) :: A_values
  integer, dimension(5) :: A_distance
  real(8), dimension(m*n) :: r, rt, u, ut, c, Axx, Atut, x, b_values
  real(8) :: rho, rho1, gamma, beta, alpha, r_norm

  A_distance = (/-m, -1, 0, 1, m/)

  ! Convert values into CDS Format
  do j = 1,n
    do i = 1,m

      ! Compress stiffness matrix values
      A_values(i+(j-1)*m,1) = -As(i,j)
      A_values(i+(j-1)*m,2) = -Aw(i,j)
      A_values(i+(j-1)*m,3) = Ap(i,j)
      A_values(i+(j-1)*m,4) = -Ae(i,j)
      A_values(i+(j-1)*m,5) = -An(i,j)

      ! Compress right-hand side values
      b_values(i+(j-1)*m) = b(i,j)+1.0e-8

      ! Compress preconditioning values
      x(i+(j-1)*m) = phi(i,j)

    end do
  end do

  ! ================================================================= !
  ! ====================== Start BiCG Algoritm ====================== !
  ! ================================================================= !

  ! Set x
  !x = 1.

  ! Compute r0
  call mkl_ddiagemv('N', m*n, A_values, m*n, A_distance, 5, x, Axx)
  r = b_values - Axx

  ! Set rt
  rt = r

  ! Set rho1
  rho1 = 1.

  ! Set u and ut
  u = 0.
  ut = 0.

  ! Start BiCG Loop
  do itr = 1,maxit

    rho = ddot(m*n, r, 1, rt, 1)

	  beta = -rho/rho1

	  u = r - beta * u

	  call mkl_ddiagemv('N', m*n, A_values, m*n, A_distance, 5, u, c)

	  ut = rt - beta * ut

	  gamma = ddot(m*n, c, 1, rt, 1)

	  alpha = rho / gamma

	  x = x + alpha * u

	  r = r - alpha * c

	  r_norm = dnrm2(m*n, r, 1)

    !print *, 'Iteration: ', itr
    !print *, 'Relative residual: ', r_norm

	  if (r_norm .le. tol) then

      !print *, 'BiCG Algorithm successfully converged!'
      !print *, 'Number of Iterations: ', itr
      !print *, 'Relative residual: ', r_norm

      ! Update phi with the solution
      do j = 1,n
        do i = 1,m
          phi(i,j) = x(i+(j-1)*m)
        end do
      end do

	    return

    elseif (itr .eq. maxit) then
      !print *, 'BiCG Algorithm did not converge!'
      !print *, 'Number of Iterations: ', itr
      !print *, 'Relative residual: ', r_norm
	  end if

	  call mkl_ddiagemv('T', m*n, A_values, m*n, A_distance, 5, ut, Atut)

	  rt = rt - alpha * Atut

	  rho1 = rho

  end do

  ! ================================================================= !
  ! ======================= End BiCG Algoritm ======================= !
  ! ================================================================= !

  return

end subroutine solver2d_bicg
