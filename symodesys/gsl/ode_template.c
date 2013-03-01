#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_pow_int.h>

#include "ode.h"
#include "func.h"
#include "jac.h"

// Python Mako template of C file
// Variables: ${NY}

int
integrate_ode_using_driver_fixed_step (double t, double t1, double y[], int n_steps,
			    double h_init, double h_max, double eps_abs,
			    double eps_rel, void *params, int print_values)
{
  int i; // Counter in macro-step loop
  int j; // Counter in print loop
  int status;
  size_t dim = ${NY};
  double ti;
  double dt = (t1-t)/(double)n_steps;
  const gsl_odeiv2_step_type * T = gsl_odeiv2_step_msbdf;
  gsl_odeiv2_step * s = gsl_odeiv2_step_alloc (T, dim);
  //gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (eps_abs, eps_rel);
  //gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc (dim);

  gsl_odeiv2_system sys = {func, jac, dim, params};

  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_msbdf, h_init, eps_abs, eps_rel);
  gsl_odeiv2_step_set_driver(s, d);

  if (h_max > 0.0)
    {
      gsl_odeiv2_driver_set_hmax(d, h_max);
    }

  for (i = 0; i < n_steps; ++i)
    {
      // Macro-step loop
      ti = t + dt*(i+1);
      status = gsl_odeiv2_driver_apply (d, &t, ti, y);

      if (status != GSL_SUCCESS)
	{
	  printf ("error, return value=%d\n", status);
	  break;
	}
      if (print_values)
	{
	  printf(STRINGIFY(PRECISION), t);
	  for (j = 0; j < ${NY}; ++j)
	    {
	      printf(" " STRINGIFY(PRECISION), y[j]);
	    }
	  printf("\n");
	}

    }

  gsl_odeiv2_driver_free (d);
  gsl_odeiv2_step_free (s);

  return status;
}

