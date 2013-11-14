#ifndef _DRIVERS_H_
#define _DRIVERS_H_

#include <gsl/gsl_odeiv2.h>
#include "symodesys_util.h"

const gsl_odeiv2_step_type *
get_step_type(int index);

int
integrate_fixed_step (double t, double t1, double * y0, int n_steps,
		      double h_init, double h_max, double abstol,
		      double reltol, void * params, size_t dim,
		      int nderiv, double * tout, double * Yout, int step_type_idx);

int
integrate_fixed_step_print(double t, double t1, double * y, int n_steps,
			   double h_init, double h_max, double abstol,
			   double reltol, void *params, size_t dim, int nderiv,
			   int step_type_idx);

#endif
