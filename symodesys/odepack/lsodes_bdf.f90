! Template for generating Fortran 90 code for calling lsodes from python (cython).
module lsodes_bdf
! ode: callbacks for calculating func and jac (generated from templates)
use ode, only: func, jac, d2ydt2, neq, nnz, nparams, ia, ja 
use iso_c_binding, only: c_double, c_int
implicit none
private
public lsodes_bdf_integrate

contains

subroutine lsodes_bdf_integrate(y, t0, tend, atol, rtol, nt, &
     h_init, h_max, nderiv, yres, tres) bind(c)
  ! This template is for BDF (5th order) with analytic Jac provided as callback

  ! mf: method flag, this means: MOSS=1, METH=2, MITER=1
  integer, parameter :: mf=121 !, maxord=5 
  ! MOSS=1: user supplied JAC (generated symbolically)
  ! METH=2: Backward Differentiation Formulas (stiff), this implies maxord 5
  ! MITER=1: chord iteration with a user-supplied sparse Jac, given by sub JAC

  ! Define iwork length
  integer, parameter :: liw=31+neq+nnz ! len(iwork)
  integer :: i, iwork(liw) ! MOSS=0 and MITER=1 or 2

  ! Define rwork length
  integer, parameter :: lenrat = 2 !hardcoded (x64): sizeof(c_double)/sizeof(integer)
  integer, parameter :: lwm = 2*nnz + 2*neq + (nnz+9*neq)/lenrat
  integer, parameter :: safety_factor = 2
  ! lrw: minimum len(rwork), formula for mf=121
  integer, parameter :: lrw=(20+9*neq+lwm)*safety_factor 
  real(c_double) :: rwork(lrw)

  ! Define input variables
  real(c_double), intent(in), dimension(neq+nparams) :: y ! params are stored in end of y
  real(c_double), intent(in) :: t0, tend, atol, rtol, h_init, h_max
  integer(c_int), intent(in) :: nt
  integer(c_int), intent(in) :: nderiv
  real(c_double), intent(out), dimension(neq, nderiv+1, nt) :: yres
  real(c_double), intent(out), dimension(nt) :: tres

  ! Define work variables
  real(c_double) :: tout, t
  integer :: itol, itask, istate, iopt
  real(c_double) :: ydot(neq), yddot(neq)
 
  itol = 1
  itask = 1
  istate = 1
  iopt = 1

  ! Initialize iwork, rwork to zero
  do i = 1,30
     iwork(i) = 0
  end do
  do i = 1,lrw
     rwork(i) = 0.0
  end do

  ! Set h_init and h_max
  rwork(5) = h_init
  rwork(6) = h_max

  ! Sparsity
  do i = 1,neq
     iwork(30+i) = ia(i)
  end do
  do i = 1,nnz
     iwork(30+neq+i) = ja(i)
  end do


  ! Save first data point
  i=0
  tres(i+1) = t0
  yres(:, 1, i+1) = y(1:neq)

  ! Save higher derivatives (for interpolation)
  if (nderiv > 0) then
     call func(neq, t, y, ydot)
     yres(:, 2, i+1) = ydot
  end if

  if (nderiv > 1) then
     call d2ydt2(neq, t, y, yddot, ia, ja)
     yres(:, 3, i+1) = yddot
  end if

  
  t = t0
  do i = 1,nt-1
     tout = t0 + (tend-t0)/real(nt-1, c_double)*real(i, c_double)

     call dlsodes (func, neq, y, t, tout, itol, rtol, atol, &
          itask, istate, iopt, rwork, lrw, iwork, liw, jac, mf)

     tres(i+1) = t
     yres(:, 1, i+1) = y(1:neq)

     ! Variables which could be of use if adding intelligence
     ! to this driving loop
     ! hu    = rwork(11) ! step sized used last
     ! hcur  = rwork(12) ! step size to be used on next step

     ! Save higher derivatives (for interpolation)
     if (nderiv > 0) then
        call func(neq, t, y, ydot)
        yres(:, 2, i+1) = ydot
     end if

     if (nderiv > 1) then
        call d2ydt2(neq, t, y, yddot, ia, ja)
        yres(:, 3, i+1) = yddot
     end if
  end do

end subroutine

end module
