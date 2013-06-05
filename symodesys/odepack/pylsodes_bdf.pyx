from numpy cimport ndarray
from numpy import zeros, ascontiguousarray, float64, swapaxes, concatenate

cdef extern void c_integrate(double * y, double * t0, double * tend, double * atol, double * rtol, int * nt, int * nderiv, double * yres, double * tres)

def integrate(double [:] y0, double t0, double tend, double [:] params, double atol, double rtol, int nt, int nderiv):
    cdef ndarray[double] y = ascontiguousarray(concatenate((y0, params)), dtype=float64)
    cdef ndarray[double, mode="fortran", ndim=3] yres = zeros((len(y0), nderiv+1, nt+1), order='F', dtype=float64)
    cdef ndarray[double] tres = zeros(nt+1, dtype=float64)

    # integrate
    c_integrate(&y[0], &t0, &tend, &atol, &rtol, &nt, &nderiv, &yres[0,0,0], &tres[0])

    # yres axes are now: depv, nderiv, nt
    # swap to make them nt, depv, nderiv and change
    # ordering to C-order:
    Yres = swapaxes(yres, 0, 1)
    Yres = swapaxes(Yres, 0, 2)
    return tres, ascontiguousarray(Yres) #tres, yres
