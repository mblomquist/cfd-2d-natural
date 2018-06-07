! solver2d_tdma
!
! Written by Matt Blomquist
! Last Update: 2018-05-15 (YYYY-MM-DD)
!
! This program solves a two-dimensional discretization problem utilizing a line-by-line
! TDMA (tri-diagonal matrix algorithm).
!
subroutine solver2d_tdma(Aw, Ae, As, An, Ap, b, phi, m, n, tol, maxit)

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

    do j = 1, n

      if (j .eq. 1) then
        do i = 1, m

          if (i .eq. 1) then
            r(i,j) = Ap(i,j)*phi(i,j) - (Ae(i,j)*phi(i+1,j)+An(i,j)*phi(i,j+1)+b(i,j))
          elseif (i .eq. m) then
            r(i,j) = Ap(i,j)*phi(i,j) - (Aw(i,j)*phi(i-1,j)+An(i,j)*phi(i,j+1)+b(i,j))
          else
            r(i,j) = Ap(i,j)*phi(i,j) - (Aw(i,j)*phi(i-1,j)+Ae(i,j)*phi(i+1,j)+An(i,j)*phi(i,j+1)+b(i,j))
          end if

        end do
      elseif (j .eq. n) then
        do i = 1, m

          if (i .eq. 1) then
            r(i,j) = Ap(i,j)*phi(i,j) - (Ae(i,j)*phi(i+1,j)+As(i,j)*phi(i,j-1)+b(i,j))
          elseif (i .eq. m) then
            r(i,j) = Ap(i,j)*phi(i,j) - (Aw(i,j)*phi(i-1,j)+As(i,j)*phi(i,j-1)+b(i,j))
          else
            r(i,j) = Ap(i,j)*phi(i,j) - (Aw(i,j)*phi(i-1,j)+Ae(i,j)*phi(i+1,j)+As(i,j)*phi(i,j-1)+b(i,j))
          end if

        end do
      else
        do i = 1, m

          if (i .eq. 1) then
            r(i,j) = Ap(i,j)*phi(i,j) - (Ae(i,j)*phi(i+1,j)+As(i,j)*phi(i,j-1)+An(i,j)*phi(i,j+1)+b(i,j))
          elseif (i .eq. m) then
            r(i,j) = Ap(i,j)*phi(i,j) - (Aw(i,j)*phi(i-1,j)+As(i,j)*phi(i,j-1)+An(i,j)*phi(i,j+1)+b(i,j))
          else
            r(i,j) = Ap(i,j)*phi(i,j) - (Aw(i,j)*phi(i-1,j)+Ae(i,j)*phi(i+1,j)+As(i,j)*phi(i,j-1)+An(i,j)*phi(i,j+1)+b(i,j))
          end if

        end do
      end if

    end do

    r_sum = 0

    do j = 1, n
      do i = 1, m
        r_sum = r_sum + abs(r(i,j))
      end do
    end do

    !print *, "r_sum, itr:", r_sum, k

    if (r_sum < tol) then
      !print *, "TDMA Compelete."
      !print *, "r_sum:", r_sum
      !print *, "itrs:", k
      return
    end if

  end do

  return

end subroutine solver2d_tdma
