! solver3d_gmres
!
! Written by Matt Blomquist
! LAxst Update: 2018-07-23 (YYYY-MM-DD)
!
! This program solves a three-dimensional finite volume discretization problem
! using the generalized minimum residual (GMRES) algorithm.
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

subroutine solver2d_gmres(As, Aw, Ap, Ae, An, b, phi, m, n, tol, maxit, fault)

  ! Define implicit
  implicit none

  ! Include mkl functions
  include "mkl.fi"

  ! Define input variables
  integer, intent(in) :: m, n
  integer, intent(in) :: maxit
  integer, intent(inout) :: fault
  real(8), intent(in) :: tol
  real(8), dimension(m,n), intent(in) :: As, Aw, Ap, Ae, An, b
  real(8), dimension(m,n), intent(inout) :: phi

  ! Define internal variables :: Matrix Conversion
  integer :: i, j, k, ii, jj, itr_init
  real(8), dimension(m*n,5) :: A_values
  integer, dimension(5) :: A_distance
  real(8), dimension(m*n) :: b_values

  ! Define internal variables :: GMRES

  real(8) :: b_norm, err, r_norm, mult1, temp_tri_sol
  real(8), dimension(m*n) :: r, Axx, x, e1, Qy, beta
  real(8), dimension(maxit) :: sn, cs, y, cs2, sn2
  real(8), dimension(m*n, maxit+1) :: Q
  real(8), dimension(maxit+1, maxit) :: H



  ! ================================================================= !
  ! ==================== Start Matrix Conversion ==================== !
  ! ================================================================= !

  ! Set fault
  fault = 0

  ! Set input integer
  itr_init = 0

  do i = maxit-1, 3, -1
    itr_init = itr_init + i
  end do

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
  ! ==================== End Matrix Conversion ====================== !
  ! ================================================================= !

  ! ================================================================= !
  ! ====================== Start GMRES Algoritm ===================== !
  ! ================================================================= !

  ! Set x = 1
  !x = 1.

  ! Compute residual vector
  call mkl_ddiagemv('N', m*n, A_values, m*n, A_distance, 5, x, Axx)
  r = b_values - Axx

  ! Compute residual norm
  r_norm = dnrm2(m*n, r, 1)

  ! Compute norm of values vector
  b_norm = dnrm2(m*n, b_values, 1)

  ! Compute inital error
  err = r_norm / b_norm

  ! Initialize vectors
  sn = 0.
  cs = 0.
  e1 = 0.
  e1(1) = 1.

  ! Set Q to inital residual
  Q(:,1) = r / r_norm

  ! Set inital beta
  beta = r_norm * e1

  ! Start GMRES Loop
  do k = 1,maxit


    cs2 = cs
    sn2 = sn

    ! Run Arnoldi Loop
	  call arnoldi(A_values, A_distance, Q, H, m*n, k, maxit)

	  ! Eliminate the last element in H ith row and update rotation matrix
	  call givens_rotation(H, cs2, sn2, m*n, k, maxit)

    cs = cs2
    sn = sn2

	  ! Update the residual vector
	  beta(k+1) = -sn(k) * beta(k)
	  beta(k) = cs(k) * beta(k)
	  err = abs(beta(k+1)) / b_norm

	  if (err .le. tol) then

      y(1:k) = beta(1:k)

      call dtrsv('U', 'N', 'N', k, H(1:k,1:k), k, y(1:k), 1)

      Qy = 0.

      do i = 1,m*n
        do j = 1,k

          Qy(i) = Qy(i) + Q(i,j)*y(j)

        end do
      end do

      x = x + Qy


      do jj = 1,n
        do ii = 1,m
          phi(ii,jj) = x(ii+(jj-1)*m)
        end do
      end do


      !print *, 'GMRES Algorithm successfully converged!'
      !print *, 'Number of Iterations: ', k
      !print *, 'Relative residual: ', err
      fault = 1

      return

	  elseif (k .eq. maxit) then

      y(1:k) = beta(1:k)

      call dtrsv('U', 'N', 'N', k, H(1:k,1:k), k, y(1:k), 1)

      Qy = 0.

      do i = 1,m*n
        do j = 1,k

          Qy(i) = Qy(i) + Q(i,j)*y(j)

        end do
      end do

      x = x + Qy

      do jj = 1,n
        do ii = 1,m
          phi(ii,jj) = x(ii+(jj-1)*m)
        end do
      end do

      !print *, 'GMRES Algorithm did not converge!'
      !print *, 'Number of Iterations: ', k
      !print *, 'Relative residual: ', err
      fault = 0

      return

	  end if

  end do

  ! ================================================================= !
  ! ======================== End GMRES Algoritm ===================== !
  ! ================================================================= !



  return

end subroutine solver2d_gmres


subroutine arnoldi(A_values, A_distance, Q, H, m, k, maxit)

  ! Define implicit
  implicit none

  ! Include mkl functions
  include "mkl.fi"

  ! Define input variables
  integer, intent(in) :: m, k, maxit
  integer, dimension(5) :: A_distance
  real(8), dimension(m,5) :: A_values
  real(8), dimension(m,maxit+1) :: Q
  real(8), dimension(maxit+1,maxit) :: H

  ! Define internal variables
  integer :: i
  real(8), dimension(m) :: q_vec

  call mkl_ddiagemv('N', m, A_values, m, A_distance, 5, Q(:,k), q_vec)

  do i = 1,k

    H(i,k) =  ddot(m, q_vec, 1, Q(:,i), 1)
	  q_vec = q_vec - H(i,k) * Q(:,i)

  end do

  H(k+1,k) = dnrm2(m, q_vec, 1)
  q_vec = q_vec / H(k+1,k)
  Q(:,k+1) = q_vec

  return

end subroutine arnoldi

subroutine givens_rotation(H, cs, sn, m, k, maxit)

  ! Define implicit
  implicit none

  ! Define input variables
  integer, intent(in) :: m, k, maxit
  real(8), dimension(maxit) :: sn, cs
  real(8), dimension(maxit+1, maxit) :: H

  ! Define internal variables
  integer :: i
  real(8) :: temp

  do i = 1,k-1

    temp = cs(i) * H(i,k) + sn(i)*H(i+1,k)

	  H(i+1,k) = -sn(i)*H(i,k) + cs(i)*H(i+1,k)
	  H(i,k) = temp

  end do

  call calc_rotation(H, cs, sn, m, k, maxit)

  H(k,k) = cs(k)*H(k,k) + sn(k)*H(k+1,k)
  H(k+1,k) = 0.

  return

end subroutine givens_rotation


subroutine calc_rotation(H, cs, sn, m, k, maxit)

  ! Define implicit
  implicit none

  ! Define Input Variables
  integer, intent(in) :: m, k, maxit
  real(8), dimension(maxit) :: cs, sn
  real(8), dimension(maxit+1, maxit) :: H

  ! Define internal variables
  real(8) :: temp

  if (H(k,k) .eq. 0.) then

    cs(k) = 0.
	  sn(k) = 1.

  else

    temp = (H(k,k)**2.0+H(k+1,k)**2.0)**(0.5)
	  cs(k) = abs(H(k,k)) / temp
	  sn(k) = cs(k) * H(k+1,k) / H(k,k)

  end if

  return
end subroutine calc_rotation
