! Template for generating Fortran 90 code to wrapped using Cython for calling lsodes from python.
! mako template variables: NY, NNZ, IA, JA, NPARAM, f, cse_func, yale_jac_cse, yale_jac_expr, dfdt

module ode
use types, only: dp
implicit none
! Set problem specific values:
integer, parameter :: nnz=${NNZ}, neq=${NY}, nparams=${NPARAM}
!real(dp), save :: params(nparams)
! Sparsity structure
integer, parameter :: ia(neq) = ${IA} ! used in lsodes_bdf.f90
integer, parameter :: ja(nnz) = ${JA} ! used in lsodes_bdf.f90

public func, jac

contains


subroutine func(neq, t, y, ydot)
  integer, intent(in) :: neq
  real(dp), intent(in) :: t, y(${NY}+${NPARAM})
  real(dp), intent(inout) :: ydot(${NY})
% for cse_token, cse_expr in cse_func:
  real(dp) :: ${cse_token}
% endfor
% for cse_token, cse_expr in cse_func:
  ${cse_token} = ${cse_expr}
% endfor
% for i, expr in enumerate(f, 1):
  ydot(${i}) = ${expr}
% endfor
  return
end subroutine

subroutine jac(neq, t, y, j, ian, jan, pdj)
  integer, intent(in) :: neq, j, ian(${NY}), jan(${NNZ})
  real(dp), intent(in) :: t, y(${NY}+${NPARAM})
  real(dp), intent(inout) :: pdj(${NY})
% for i in range(1,NY+1):
% for cse_token, cse_expr in yale_jac_cse[i-1]:
  real(dp) :: ${cse_token}
% endfor
% endfor
  select case (j)
  ! e.g.: pdj(1) = -1.0_dp - 2.0_dp*y(1)*y(2)*y(nparams+1)
! i loops over column indices
% for i in range(1,NY+1):
  case (${i})
  % for cse_token, cse_expr in yale_jac_cse[i-1]:
     ${cse_token} = ${cse_expr}
  % endfor
  ! k loops over row indices
  % for k, expr in yale_jac_exprs[i-1]:
     pdj(${k}) = ${expr}
  % endfor
     return
% endfor
  end select
end subroutine

subroutine dfdt(t, y, pdt)
  real(dp), intent(in) :: t, y(${NY}+${NPARAM})
  real(dp), intent(inout) :: pdt(${NY})
  integer :: i
  ! used fo explicit calc of d2ydt2
% for k, expr in dfdt:
  pdt(${k}+1) = ${expr}
% endfor
end subroutine

subroutine d2ydt2(neq, t, y, yddot, ian, jan)
  integer, intent(in) :: neq, ian(${NY}), jan(${NNZ})
  real(dp), intent(in) :: t, y(${NY}+${NPARAM})
  real(dp), intent(inout) :: yddot(${NY})
  real(dp) :: pdj(${NY}), ydot(${NY}), pdt(${NY})
  integer :: i, j, k

  call func(neq, t, y, ydot)
  call dfdt(t, y, pdt)
  do i = 1,neq
     yddot(i) = pdt(i)
     do j = 1,neq
        pdj(j) = 0.0_dp
     end do
     call jac(neq, t, y, i, ian, jan, pdj)

     ! since ian doesn't store the upper boundary for 
     ! the last item, we need to handle this edge case
     if (i == neq) then
        k = ${NNZ}
     else
        k = iran(i+1)
     end if
     do j = ian(i),k
        yddot(i) = yddot(i) + pdj(jan(j))*ydot(jan(j))
     end do
  end do
end subroutine

end module
