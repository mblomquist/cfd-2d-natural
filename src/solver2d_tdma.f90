! solver2d_tdma
!
! Written by Matt Blomquist
! Last Update: 2018-02-05 (YYYY-MM-DD)
!
! This program solves a two-dimensional discretization problem utilizing a line-by-line
! TDMA (tri-diagonal matrix algorithm). The algorithm sweeps from South to North then from
! West to East.

subroutine solver2d_tdma(As, Aw, Ap, Ae, An, b, phi, m, n, tol, maxit)

  !implicit none

  ! Include mkl header
  include "mkl.fi"

  ! Define input variables
  integer, intent(in) :: m, n
  integer, intent(in) :: maxit
  real(8), intent(in) :: tol
  real(8), dimension(m,n), intent(in) :: As, Aw, Ap, Ae, An, b
  real(8), dimension(m,n), intent(inout) :: phi

  ! Define internal variables
  integer :: i, j, k, itr
  integer, dimension(5) :: A_distance
  real(8) :: residual
  real(8), dimension(m) :: bx_line
  real(8), dimension(n) :: by_line
  real(8), dimension(m,n) :: x
  real(8), dimension(m*n) :: x_compressed, b_compressed, Ax
  real(8), dimension(m*n,5) :: A_compressed

  ! Set temporary x value
  x = phi

  ! Initilize b_line
  bx_line = 0
  by_line = 0

  ! Start TDMA Sweep Loop
  do itr = 1,maxit

    ! Iterate lines from South to North
    !do k = 1,itr

      do i = 1,m
	    do j = 1,n

	      if (i .eq. 1) then
		    by_line(j) = b(i,j)-An(i,j)*x(i+1,j)
		  else if (i .eq. m) then
		    by_line(j) = b(i,j)-As(i,j)*x(i-1,j)
		  else
		    by_line(j) = b(i,j)-As(i,j)*x(i-1,j)-An(i,j)*x(i+1,j)
		  end if

	    end do

	    call solver1d_tdma(Aw(i,:), Ap(i,:), Ae(i,:), by_line, x(i,:), m)

	  end do

    !end do

    ! Iterate lines from West to East
    !do k = 1,itr

      do j = 1,n
	    do i = 1,m

	      if (j .eq. 1) then
		    bx_line(i) = b(i,j)-Ae(i,j)*x(i,j+1)
		  else if (j .eq. n) then
		    bx_line(i) = b(i,j)-Aw(i,j)*x(i,j-1)
		  else
		    bx_line(i) = b(i,j)-Aw(i,j)*x(i,j-1)-Ae(i,j)*x(i,j+1)
		  end if

	  !end do

	  call solver1d_tdma(As(:,j), Ap(:,j), An(:,j), bx_line, x(:,j), m)

    end do

    ! Check Residual
    A_distance = (/-m, -1, 0, 1, m/)

    ! Convert values into CDS Format
    do j = 1,n
      do i = 1,m

        ! Compress stiffness matrix values
        A_compressed(i+(j-1)*m,1) = As(i,j)
        A_compressed(i+(j-1)*m,2) = Aw(i,j)
        A_compressed(i+(j-1)*m,3) = Ap(i,j)
        A_compressed(i+(j-1)*m,4) = Ae(i,j)
        A_compressed(i+(j-1)*m,5) = An(i,j)

        ! Compress right-hand side values
        b_compressed(i+(j-1)*m) = b(i,j)

        ! Compress preconditioning values
        x_compressed(i+(j-1)*m) = x(i,j)

      end do
	end do

	! Compute matrix-vector product of A_compressed and x_compressed
	call mkl_ddiagemv('N', m*n, A_compressed, m*n, A_distance, 5, x_compressed, Ax)

	! Compute norm of b-Ax
	residual = abs(dnrm2(m*n, b_compressed-Ax, 1))

	if (residual < tol) then
      print *, 'Completed!'
      print *, 'Number of Iterations: ', itr
      print *, 'Relative residual: ', residual
      exit
    elseif (k .eq. maxit) then
      print *, 'TDMA did not converge!'
      print *, 'Number of Iterations: ', itr
      print *, 'Relative residual: ', residual
    end if

  end do
  ! Update phi with the solution
  do j = 1,n
    do i = 1,m
      phi(i,j) = x(i,j)
    end do
  end do

  return

end subroutine solver2d_tdma
