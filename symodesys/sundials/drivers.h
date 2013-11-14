#ifndef _SYMODESYS_SUNDIALS_DRIVERS_H_
#define _SYMODESYS_SUNDIALS_DRIVERS_H_

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., CV_BDF, CV_ADAMS */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */

#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */


#ifdef USE_LAPACK
  #include <cvode/cvode_lapack.h>       /* prototype for CVDense */
  #define OUR_DENSE CVLapackDense
  #define OUR_BANDED CVLapackBand
  // Sundials 2.7.0 changed int -> long int, but BLAS uses int
  // We use SIZE_T defined here:
  #define SIZE_T int
#else
  #include <cvode/cvode_dense.h>       /* prototype for CVDense */
  #define OUR_DENSE CVDense
  #define OUR_BAND CVBand
  #define SIZE_T long int
#endif

int
func (realtype t, N_Vector y, N_Vector f, void * params);
int
dense_jac (long int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, void *params, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);


enum {
  DENSE_MODE = 0,
  BANDED_MODE = 1,
}


int
integrate_fixed_step (double t, double t1, double * y0, int n_steps,
			double h_init, double h_max, double abstol,
			double reltol, void * params, size_t dim,
			int nderiv, double * tout, double * Yout,
			int step_type_idx, int mode);

int
integrate_fixed_step_print(double t, double t1, double * y, int n_steps,
			     double h_init, double h_max, double abstol,
			     double reltol, void *params, size_t dim, 
			     int nderiv, int step_type_idx, int mode);


#endif // _SYMODESYS_SUNDIALS_DRIVERS_H_
