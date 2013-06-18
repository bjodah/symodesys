program main
use types, only: dp
use transform, only: perform
implicit none
integer, parameter :: n = 3
real(dp), dimension(n) :: x, y
real(dp), dimension(n,2) :: output
x = [1.0_dp, 2.0_dp, 3.0_dp]
y = [7.0_dp, 8.0_dp, 9.0_dp]
call perform(n, x, y, output)
print *, output
end program
