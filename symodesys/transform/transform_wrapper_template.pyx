from numpy cimport ndarray
from numpy import zeros, float64

cdef extern void c_perform(int * n,
%for arg_name in ARGS:
    double * ${arg_name},
%endfor
double * output)

def transform(
%for arg_name in ARGS:
    double [:] ${arg_name},
%endfor
    ):
    cdef ndarray[double, mode="fortran", ndim=2] output = zeros((len(${ARGS[0]}), ${N_EXPRS}), order='F', dtype=float64)
    cdef int n = len(${ARGS[0]})

    c_perform(&n,
%for arg_name in ARGS:
        &${arg_name}[0],
%endfor
        &output[0,0])
    return output