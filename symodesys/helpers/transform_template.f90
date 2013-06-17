module types
    implicit none
    integer, parameter :: dp=kind(0.d0)  ! double precision
end module

module transform
    use types, only: dp
    implicit none

contains

    subroutine perform(${ARGS_COMMA}, output)
    %for arg_name in ${ARGS}
        real(dp), intent(in) :: ${arg_name}(:)
    %endfor
        real(dp), intent(out) :: output(:,:)

    %for i, f90expr in enumerate(f90exprs)
        output(:, ${i}+1) = f90expr
    %endfor
end module

module transform_wrapper
    use iso_c_binding, only: c_double
    implicit none

contains

    subroutine c_perform(${ARGS_COMMA}) bind(c)
    %for arg_name in ${ARGS}
        real(c_double), intent(in) :: ${arg_name}(:)
    %endfor
        real(c_double), intent(out) :: output(:,:)
        call perform(${ARGS_COMMA})
    end subroutine

end module
