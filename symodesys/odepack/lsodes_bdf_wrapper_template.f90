module lsodes_bdf_wrapper

use iso_c_binding, only: c_double, c_int
use lsodes_bdf, only: integrate

implicit none

contains

subroutine c_integrate(y, t0, tend, atol, rtol, nt, h_init, h_max, nderiv, yres, tres) bind(c)
real(c_double), intent(in) :: y(${NY}+${NPARAM})
real(c_double), intent(in) :: t0
real(c_double), intent(in) :: tend
real(c_double), intent(in) :: atol
real(c_double), intent(in) :: rtol
integer(c_int), intent(in) :: nt
real(c_double), intent(in) :: h_init
real(c_double), intent(in) :: h_max
integer(c_int), intent(in) :: nderiv
real(c_double), intent(out) :: yres(${NY}, nderiv+1, nt)
real(c_double), intent(out) :: tres(nt)
call integrate(y, t0, tend, atol, rtol, nt, h_init, h_max, nderiv, yres, tres)
end subroutine

end module
