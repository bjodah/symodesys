from numpy cimport ndarray
from numpy import zeros, float64, asfortranarray

cdef extern void c_perform(int * n, int * m, double * input, double * output)

def transform(double [:,:] inp, int n_exprs):
    cdef int n = inp.shape[0]
    cdef int m = inp.shape[1]
    cdef ndarray[double, mode="fortran", ndim=2] inp_arr = \
        asfortranarray(inp)
    cdef ndarray[double, mode="fortran", ndim=2] output = \
        zeros((n, n_exprs), order='F', dtype=float64)

    c_perform(&n, &m, &inp_arr[0,0], &output[0,0])
    return output
