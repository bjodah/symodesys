from numpy cimport ndarray
from numpy import zeros, float64

cdef extern void c_perform(
% for arg_name in ARGS
    double * ${arg_name},
% endfor
double * output)

def transform(
% for arg_name in ARGS
    double [:] ${arg_name},
% endfor
    ):
    cdef ndarray[double] y = ascontiguousarray(concatenate((y0, params)), dtype=float64)
    cdef ndarray[double, mode="fortran", ndim=2] output = zeros((len(${ARGS[0]}), ${N_EXPRS}), order='F', dtype=float64)

    c_perform(
%for arg_name in ARGS:
        &${arg_name}[0],
%endfor
        &output[0,0])
    return output
