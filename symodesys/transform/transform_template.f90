module types
    implicit none
    integer, parameter :: dp=kind(0.d0)  ! double precision
end module

module transform
    use types, only: dp
    implicit none

contains

    subroutine perform(n, ${ARGS_COMMA}, output)
        integer, intent(in) :: n
    %for arg_name in ARGS:
        real(dp), intent(in), dimension(n) :: ${arg_name}
    %endfor
        real(dp), intent(out), dimension(n, ${N_EXPRS}) :: output
    %for cse_token, cse_expr in CSES:
        real(dp), dimension(n) :: ${cse_token}
    %endfor
    %for cse_token, cse_expr in CSES:
        ${cse_token} = ${cse_expr}
    %endfor
    %for i, expr_in_cse in enumerate(EXPRS_IN_CSE):
        output(:, ${i}+1) = ${expr_in_cse}
    %endfor
    end subroutine
end module

module transform_wrapper
    use iso_c_binding, only: c_double, c_int
    use transform, only: perform
    implicit none

contains

    subroutine c_perform(n, ${ARGS_COMMA}, output) bind(c)
        integer(c_int), intent(in) :: n
    %for arg_name in ARGS:
        real(c_double), intent(in), dimension(n) :: ${arg_name}
    %endfor
        real(c_double), intent(out), dimension(n, ${N_EXPRS}) :: output
        call perform(n, ${ARGS_COMMA}, output)
    end subroutine

end module
