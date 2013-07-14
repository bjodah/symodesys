from numpy cimport ndarray
from numpy import zeros, ascontiguousarray, float64, swapaxes, concatenate


# lsodes_bdf_integrate: lsodes_bdf.f90 (iso_c_binding module)
cdef extern void lsodes_bdf_integrate(
    double * y, double * t0, double * tend, double * atol,
    double * rtol, int * nt, double * h_init, double * h_max,
    int * nderiv, double * yres, double * tres)


def integrate_equidistant_output(
        double t0, double tend, double [:] y0, int nt,
        double h_init, double h_max, double atol,
        double rtol, double [:] params, int nderiv):
    """
    This routine allocates memory for output arrays and
    handles C/Fortran ordering issue of arrays and different
    structure of Yres during the integrationl.
    """
    cdef ndarray[double] y = ascontiguousarray(
        concatenate((y0, params)), dtype=float64)
    cdef ndarray[double, mode="fortran", ndim=3] yres = zeros(
        (len(y0), nderiv+1, nt), order='F', dtype=float64)
    cdef ndarray[double] tres = zeros(nt, dtype=float64)

    # integrate
    lsodes_bdf_integrate(&y[0], &t0, &tend, &atol, &rtol,
                         &nt, &h_init, &h_max, &nderiv,
                         &yres[0,0,0], &tres[0])

    # yres axes are now: depv, nderiv, nt
    # swap to make them nt, depv, nderiv and change
    # ordering to C-order:
    Yres = swapaxes(yres, 0, 1)
    Yres = swapaxes(Yres, 0, 2)
    print tres
    return tres, ascontiguousarray(Yres) #tres, yres
