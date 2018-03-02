! solver1d_tdma
!
! Written by Matt Blomquist
! Last Update: 2018-01-31 (YYYY-MM-DD)
!
! This program solves a one-dimensional tri-diagonal matrix
! using the Thomas Algorithm.

subroutine solver1d_tdma(aa, bb, cc, dd, x, m)

  ! Define implicit
  implicit none

  ! Input variables
  integer, intent(in) :: m
  real(8), dimension(m), intent(in) :: aa, bb, cc, dd
  real(8), dimension(m), intent(out) :: x

  ! Internal variables
  integer :: i
  real(8) :: temp
  real(8), dimension(m) :: a, b, c, d

  ! Set x to 0
  x = 0

  a = aa
  b = bb
  c = cc
  d = dd

  ! Solve initial c and d
  c(1) = c(1)/b(1)
  d(1) = d(1)/b(1)

  ! Run forward sub loop
  do i = 2,m-1
    temp = b(i)-a(i)*c(i-1)
    c(i) = c(i)/temp
    d(i) = (d(i)-a(i)*d(i-1))/(b(i)-a(i)*c(i-1))
  end do

  x(m) = (d(m)-a(m)*d(m-1))/(b(m)-a(m)*c(m-1))

  do i = m, 1, -1
    x(i) = -c(i)*x(i+1)+d(i)
  end do

  return

end subroutine solver1d_tdma
