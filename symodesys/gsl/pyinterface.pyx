import numpy as np
cimport numpy as np

cdef extern int integrate_ode_using_driver_fixed_step(double t, double t1, double y[], int n_steps,
			    double h_init, double h_max, double eps_abs,
			    double eps_rel, const void *params, size_t dim, int order, double tout[], double Yout[])

def integrate_odeiv2_driver(double t, double t1, np.ndarray[np.double_t, ndim=1] y0, int n_steps,
                            double h_init, double h_max, double eps_abs,
                            double eps_rel, np.ndarray[np.double_t, ndim=1] params, int order):
    cdef int status
    cdef np.ndarray[np.double_t, ndim=1] params_arr = np.ascontiguousarray(params)
    cdef np.ndarray[np.double_t, ndim=1] y0_arr = np.ascontiguousarray(y0)
    cdef np.ndarray[np.double_t, ndim=1] tout = np.ascontiguousarray(np.empty(n_steps, dtype=np.float64))
    cdef np.ndarray[np.double_t, ndim=1] Yout = np.ascontiguousarray(np.empty(n_steps * len(y0) * (order + 1), dtype=np.float64))
    status = integrate_ode_using_driver_fixed_step(
        t, t1, <double*>y0_arr.data, n_steps, h_init, h_max,
        eps_abs, eps_rel,<void*>params_arr.data,
        dim = len(y0), order = order,
        tout = <double*>tout.data, Yout = <double*>Yout.data
        )
    if status != 0:
        raise RuntimeError('Integration failed')
    print Yout
    print n_steps, y0, order + 1
    return tout, Yout.reshape((n_steps, len(y0), order + 1))

