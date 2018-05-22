! solver1d_tdma
!
! Written by Matt Blomquist
! Last Update: 2018-05-15 (YYYY-MM-DD)
!
! This program solves a one-dimensional tri-diagonal matrix
! using the Thomas Algorithm.
!
! Inputs ::
!   a :: center diagonal
!   b :: upper diagonal
!   c :: lower diagonal
!   d :: right hand side
!   phi :: solution array
!   n :: vector length
!
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
