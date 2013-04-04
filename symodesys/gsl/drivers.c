#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_pow_int.h>

#include "drivers.h"
#include "func.h"
#include "jac.h"


int
integrate_ode_using_driver_fixed_step (double t, double t1, double y[], int n_steps,
                                       double h_init, double h_max, double eps_abs,
                                       double eps_rel, void * params, size_t dim, int order, double tout[], double Yout[])
{
  /* Order can be 0, 1 or 2 */
  size_t i; /* Counter in macro-step loop */
  size_t j; /* Counter in print loop */
  size_t k; /* Counter in dfdy loop */
  double temp;
  int status;
  double ti;
  double dt = (t1-t)/(double)n_steps;
  const gsl_odeiv2_step_type * T = gsl_odeiv2_step_msbdf;
  gsl_odeiv2_step * s = gsl_odeiv2_step_alloc (T, dim);
  /* gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (eps_abs, eps_rel); */
  /* gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc (dim); */

  gsl_odeiv2_system sys = {func, jac, dim, params};

  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_msbdf, h_init, eps_abs, eps_rel);
  gsl_odeiv2_step_set_driver(s, d);

  double *f = malloc(sizeof(double)*dim);
  double *dfdt = malloc(sizeof(double)*dim);
  gsl_matrix *dfdy = gsl_matrix_alloc(dim, dim);

  if (h_max > 0.0)
    {
      gsl_odeiv2_driver_set_hmax(d, h_max);
    }

  for (i = 0; i < n_steps; ++i)
    {
      /* Macro-step loop */
      ti = t + dt;//*(i+1);
      status = gsl_odeiv2_driver_apply (d, &t, ti, y);

      if (status != GSL_SUCCESS)
        {
          printf ("error, return value=%d\n", status);
          break;
        }
      tout[i] = t;
      if (order > 0)
        func(t, y, f, params);
      if (order > 1)
        jac(t, y, dfdy->data, dfdt, params);

      for (j = 0; j < dim; ++j)
        {
          Yout[i*dim*(order+1)+j*(order+1)+0] = y[j];
          if (order > 0)
            Yout[i*dim*(order+1)+j*(order+1)+1] = f[j];
          if (order > 1)
            {
              temp = 0;
              for (k=0; k<dim; ++k)
                {
                  temp += gsl_matrix_get(dfdy, j, k)*f[k];
                }
              Yout[i*dim*(order+1)+j*(order+1)+2] = dfdt[j]+temp;
            }
        }
    }

  /* Memory management */
  gsl_odeiv2_driver_free (d);
  gsl_odeiv2_step_free (s);
  free(f);
  free(dfdt);
  gsl_matrix_free(dfdy);

  return status;
}


int
integrate_ode_using_driver_fixed_step_print(double t, double t1, double y[], int n_steps,
                                            double h_init, double h_max, double eps_abs,
                                            double eps_rel, void *params, size_t dim, int order)
{
  int status;
  int i,j,k;
  double * tout;
  double * Yout;

  tout = malloc(sizeof(double)*n_steps);
  Yout = malloc(sizeof(double)*n_steps*dim*(order+1));
  status = integrate_ode_using_driver_fixed_step(t, t1, y, n_steps, h_init, h_max, eps_abs,
                                                 eps_rel, params, dim, order, tout, Yout);
  if (status != GSL_SUCCESS)
    {
      return status;
    }

  for (i = 0; i < n_steps; ++i)
    {
	  printf(STRINGIFY(PRECISION), tout[i]);
	  for (j = 0; j < dim; ++j)
	    {
          for (k = 0; k<=order; ++k)
            {
              printf(" " STRINGIFY(PRECISION), Yout[i*dim*(order+1)+j*(order+1)+k]);
            }
	    }
	  printf("\n");
    }

  free(tout);
  free(Yout);
  return status;
}

