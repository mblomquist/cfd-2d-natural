! solver3d_tdma
!
! Written by Matt Blomquist
! Last Update: 2018-07-18 (YYYY-MM-DD)
!
! This program solves a three-dimensional discretization problem utilizing a line-by-line
! TDMA (tri-diagonal matrix algorithm).
!
subroutine solver2d_tdma(As, Aw, Ap, Ae, An, b, phi, m, n, tol, maxit)

  integer, intent(in) :: m, n, maxit
  real(8), dimension(m,n), intent(in) :: Aw, Ae, As, An, Ap, b
  real(8), dimension(m,n), intent(inout) :: phi

  integer :: i, j, k
  real(8), dimension(m) :: awe, bwe, cwe, dwe, phiwe
  real(8), dimension(n) :: asn, bsn, csn, dsn, phisn
  real(8), dimension(m,n) :: r
  real(8) :: r_sum, tol

  do k = 1, maxit

    ! Start West - East Solve
    do j = 1, n, 1

      do i = 1, m, 1

	      awe(i) = Ap(i,j)
	      bwe(i) = -Ae(i,j)
	      cwe(i) = -Aw(i,j)

	      if (j .eq. 1) then
	        dwe(i) = b(i,j)+An(i,j)*phi(i,j+1)
	      elseif (j .eq. n) then
	        dwe(i) = b(i,j)+As(i,j)*phi(i,j-1)
	      else
	        dwe(i) = b(i,j)+As(i,j)*phi(i,j-1)+An(i,j)*phi(i,j+1)
	      end if

	    end do

	    call solver1d_tdma(awe, bwe, cwe, dwe, phiwe, m)

	    phi(:,j) = phiwe(:)

    end do

    ! Start South - North Solve
    do i = 1, m, 1

      do j = 1, n, 1

	    asn(j) = Ap(i,j)
	    bsn(j) = -An(i,j)
	    csn(j) = -As(i,j)

	    if (i .eq. 1) then
	      dsn(j) = b(i,j)+Ae(i,j)*phi(i+1,j)
	    elseif (i .eq. m) then
	      dsn(j) = b(i,j)+Aw(i,j)*phi(i-1,j)
	    else
	      dsn(j) = b(i,j)+Aw(i,j)*phi(i-1,j)+Ae(i,j)*phi(i+1,j)
	    end if

	  end do

	  call solver1d_tdma(asn, bsn, csn, dsn, phisn, n)

 	  phi(i,:) = phisn(:)

    end do

    ! Check Residual
    r = 0

    do j = 2, n-1
      do i = 2, m-1
        r(i,j) = Ap(i,j)*phi(i,j) - (Aw(i,j)*phi(i-1,j)+Ae(i,j)*phi(i+1,j)+As(i,j)*phi(i,j-1)+An(i,j)*phi(i,j+1)+b(i,j))
      end do
    end do

  	r_sum = 0.

  	do i = 2,m-1
  	  do j = 2,n-1
  		    r_sum = r_sum + (r(i,j))**2.0
  	  end do
  	end do

    r_sum = r_sum**(0.5)

  	if (r_sum .le. tol) then
  	  !print *, "TDMA Compelete."
      !print *, "r_sum:", r_sum
      !print *, "itrs:", k
        return
    else
      !print *, "r_sum:", r_sum
      !print *, "itrs:", k
  	end if

  end do

  return

end subroutine solver2d_tdma

subroutine solver1d_tdma(a, b, c, d, phi, n)

  ! Define input / output variables
  integer, intent(in) :: n
  real(8), intent(in) :: a(n), b(n), c(n), d(n)
  real(8), intent(out) :: phi(n)

  ! Define internal variables
  real(8), dimension(n) :: P, Q

  integer :: i

  ! Start forward-substitution
  P(1) = b(1)/a(1)
  Q(1) = d(1)/a(1)

  do i = 2, n, 1
    P(i) = b(i)/(a(i)-c(i)*P(i-1))
    Q(i) = (d(i)-c(i)*Q(i-1))/(a(i)-c(i)*P(i-1))

  end do

  phi(n) = Q(n)

  ! Start backward-substitution
  do i = n-1, 1, -1

    phi(i) = Q(i)-P(i)*phi(i+1)

  end do

  return

end subroutine solver1d_tdma
