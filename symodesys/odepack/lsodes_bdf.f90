! Template for generating Fortran 90 code to wrapped using Cython for calling lsodes from python.
! mako template variables: NY, NNZ, IA, JA
module lsodes_bdf
use types, only: dp ! dp === double precision
use ode, only: func, jac, d2ydt2, neq, nnz, nparams, ia, ja ! callbacks for calculating func and jac (generated from templates)
implicit none

private
public integrate

contains

subroutine integrate(y, t0, tend, atol, rtol, nt, nderiv, yres, tres)
  ! This template is for BDF (5th order) with analytic Jac provided as callback
  integer, parameter :: mf=121, maxord=5 ! method flag, this means: MOSS=1, METH=2, MITER=1
  ! MOSS=1: user supplied JAC (generated symbolically)
  ! METH=2: Backward Differentiation Formulas (stiff), this implies maxord 5
  ! MITER=1: chord iteration with a user-supplied sparse Jac, given by sub JAC

  ! Define iwork length
  integer, parameter :: liw=31+neq+nnz ! len(iwork)
  integer :: i, iwork(liw) ! MOSS=0 and MITER=1 or 2

  ! Define rwork length
  integer, parameter :: lenrat = 2 !hardcoded (x64): sizeof(dp)/sizeof(integer)
  integer, parameter :: lwm = 2*nnz + 2*neq + (nnz+9*neq)/lenrat
  integer, parameter :: safety_factor = 2
  integer, parameter :: lrw=(20+9*neq+lwm)*safety_factor ! minimum len(rwork), formula for mf=121
  real(dp) :: rwork(lrw)

  ! Define input variables
  real(dp), intent(in) :: y(neq+nparams) ! params are stored in end of y
  real(dp), intent(in) :: t0, tend, atol, rtol
  integer, intent(in) :: nt
  integer, intent(in) :: nderiv
  real(dp), intent(inout) :: yres(:, :, :)
  real(dp), intent(inout) :: tres(:)

  ! Define work variables
  real(dp) :: tout, t
  integer :: itol, itask, istate, iopt
  real(dp) :: ydot(neq), yddot(neq)
 
  itol = 1
  itask = 1
  istate = 1
  iopt = 1

  ! Initialize iwork, rwork to zero
  do i = 1,30
     iwork(i) = 0.0
  end do
  do i = 1,lrw
     rwork(i) = 0.0
  end do

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
  yres(:, 1, i+1) = y
  ! Save higher derivatives
  if (nderiv > 0) then
     call func(neq, t, y, ydot)
     yres(:, 2, i+1) = ydot
  end if

  if (nderiv > 1) then
     call d2ydt2(neq, t, y, yddot, ia, ja)
     yres(:, 3, i+1) = yddot
  end if

  
  t = t0
  do i = 1,nt
     tout = t0 + (tend-t0)/real(nt, dp)*real(i, dp)

     call dlsodes (func, neq, y, t, tout, itol, rtol, atol, &
          itask, istate, iopt, rwork, lrw, iwork, liw, jac, mf)

     ! TODO: determine if y corresponds to t (likely) or rwork(13)
     ! (rwork(13) is current value of t solver has actually reached)
     tres(i+1) = t ! rwork(13)
     yres(:, 1, i+1) = y
     ! Variables which could be of use if adding intelligence
     ! to this driving loop
     ! hu    = rwork(11) ! step sized used last
     ! hcur  = rwork(12) ! step size to be used on next step

     ! Save higher derivatives
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
