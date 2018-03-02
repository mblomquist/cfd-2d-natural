! solver2d_bicgstab
!
! Written by Matt Blomquist
! Last Update: 2018-02-02 (YYYY-MM-DD)
!
! This program solves a three-dimensional finite volume discretization problem
! using the bi-conjugate gradients stabilized algorithm.
!
! Definition of input arguments
! Inputs:
!   As, Aw, Ap, Ae, An :: These arrays represent the coefficients for adjacent nodes
!   b :: This array represents the right-hand side of the equation Ax=b
!   phi :: This value represents the appropriate solution array (pressure, velocity, temperature)
!   m, n :: These values represent the number of nodes for i and j for the phi value
!   tol :: represents the solution tolerance
!   maxit :: represents the maximum number of iterations of the BiCGStab Algorithm
!
! Outputs:
!   phi :: on exit, this value contains the updated solution
!   maxit :: on exit, this value contains the number of iterations of the BiCGStab algorithm
!   tol :: on exit, this value represents the normalized residual

subroutine solver2d_bicgstab(As_in, Aw, Ap_in, Ae, An, b, phi, m, n, tol, maxit)

  ! Define implicit
  implicit none

  ! Include mkl functions
  include "mkl.fi"

  ! Define input variables
  integer, intent(in) :: m, n
  integer, intent(in) :: maxit
  real(8), intent(in) :: tol
  real(8), dimension(m,n), intent(in) :: As_in, Aw, Ap_in, Ae, An, b
  real(8), dimension(m,n), intent(inout) :: phi

  ! Define internal variables
  integer :: i, j, itr
  real(8), dimension(m*n,5) :: A_values
  integer, dimension(5) :: A_distance
  real(8), dimension(m*n) :: r, r0, r1, x, p, Ap, As, s, Ax, b_values
  real(8) :: alpha, omega, beta, r_norm, s_norm


  A_distance = (/-m, -1, 0, 1, m/)

  ! Convert values into CDS Format
  do j = 1,n
    do i = 1,m

      ! Compress stiffness matrix values
      A_values(i+(j-1)*m,1) = -As_in(i,j)
      A_values(i+(j-1)*m,2) = -Aw(i,j)
      A_values(i+(j-1)*m,3) = Ap_in(i,j)
      A_values(i+(j-1)*m,4) = -Ae(i,j)
      A_values(i+(j-1)*m,5) = -An(i,j)

      ! Compress right-hand side values
      b_values(i+(j-1)*m) = b(i,j)

      ! Compress preconditioning values
      x(i+(j-1)*m) = phi(i,j)

    end do
  end do

  ! ================================================================= !
  ! ==================== Start BiCGStab Algoritm ==================== !
  ! ================================================================= !

  ! Set preconditioning array to 1
  if (sum(x) .eq. 0) then
    x = 1
  end if

  call mkl_ddiagemv('N', m*n, A_values, m*n, A_distance, 5, x, Ax)

  r = b_values - Ax

  r0 = 1

  if (m .eq. n) then
    r0(1) = 2
  end if

  p = r

  do itr = 1,maxit

    call mkl_ddiagemv('N', m*n, A_values, m*n, A_distance, 5, p, Ap)

    alpha = ddot(m*n,r,1,r0,1)/ddot(m*n,Ap,1,r0,1)

    s = r - alpha*Ap

    call mkl_ddiagemv('N', m*n, A_values, m*n, A_distance, 5, s, As)

    omega = ddot(m*n,As,1,s,1) / ddot(m*n,As,1,As,1)

    x = x + alpha*p + omega*s

    r1 = s - omega*As

    beta = ddot(m*n,r1,1,r0,1) / ddot(m*n,r,1,r0,1) * alpha/omega

    p = r1 + beta*(p-omega*Ap)

    r = r1

    ! Check convergence of BiCG
    call mkl_ddiagemv('N', m*n, A_values, m*n, A_distance, 5, x, Ax)
	  r_norm = abs(dnrm2(m*n, b_values-Ax, 1))

    if (r_norm < tol) then
        print *, 'BiCGStab Algorithm successfully converged!'
        print *, 'Number of Iterations: ', itr
        print *, 'Relative residual: ', r_norm
        exit
    end if

    if (itr .eq. maxit) then
      print *, '************************************'
      print *, '************************************'
      print *, 'BiCGStab Algorithm did not converge!'
      print *, 'Number of Iterations: ', itr
      print *, 'Relative residual: ', r_norm
      print *, '************************************'
      print *, '************************************'
    end if

  end do

  ! ================================================================= !
  ! ====================== End BiCG Algoritm ======================== !
  ! ================================================================= !

  ! Update phi with the solution
  do j = 1,n
    do i = 1,m
      phi(i,j) = x(i+(j-1)*m)
    end do
  end do

  return

end subroutine solver2d_bicgstab
