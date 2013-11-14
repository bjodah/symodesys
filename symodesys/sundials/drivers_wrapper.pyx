import numpy as np
cimport numpy as cnp

cdef extern int integrate_fixed_step(
    double t, double t1, double * y0, int n_steps,
    double h_init, double h_max, double abstol,
    double reltol, void * params, size_t dim,
    int nderiv, double * tout, double * Yout,
    int step_type_idx, int mode, int mu, int ml)

step_types = ('adams','bdf')
modes = ('dense', 'band')

def integrate_equidistant_output(
        double t, double t1, double [::1] y0,
        int n_steps, double h_init,
        double h_max, double abstol,
        double reltol, double [::1] params,
        int nderiv, str step_type = 'bdf', str mode = 'dense'):
    """
    Runs the fixed step integration using SUNDIALS
    """
    # TODO: return statistics from solution, look at scipy.ode
    cdef int status
    cdef cnp.ndarray[cnp.float64_t, ndim=1] y0_arr = \
        np.array(y0, dtype=np.float64, copy=True)
    cdef cnp.ndarray[cnp.float64_t, ndim=1] tout = \
        np.ascontiguousarray(np.empty(n_steps, dtype=np.float64))
    cdef cnp.ndarray[cnp.float64_t, ndim=1] Yout = \
        np.ascontiguousarray(np.empty(
            n_steps * len(y0) * (nderiv + 1), dtype=np.float64))

    assert y0_arr.flags['C_CONTIGUOUS']
    status = integrate_fixed_step(
        t, t1, <double *>y0_arr.data, n_steps,
        h_init, h_max, abstol, reltol, &params[0],
        len(y0), nderiv, <double *>tout.data,
        <double *>Yout.data,
        step_types.index(step_type)+1,
        modes.index(mode),
        -1, # mu
        -1, # ml
        )

    if status != 0:
        raise RuntimeError(
            'Integration failed with SUNDIALS status {}'.format(status))

    return tout, Yout.reshape((n_steps, len(y0), nderiv + 1))
