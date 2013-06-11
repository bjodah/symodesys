#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_pow_int.h>
#include <gsl/gsl_block.h>

#include "drivers.h"
#include "func.h"
#include "jac.h"

const gsl_odeiv2_step_type * get_step_type(int index){
  switch(index){
  case 0:
    return gsl_odeiv2_step_rk2;
  case 1:
    return gsl_odeiv2_step_rk4;
  case 2:
    return gsl_odeiv2_step_rkf45;
  case 3:
    return gsl_odeiv2_step_rkck;
  case 4:
    return gsl_odeiv2_step_rk8pd;
  case 5:
    return gsl_odeiv2_step_rk2imp;
  case 6:
    return gsl_odeiv2_step_rk4imp;
  case 7:
    return gsl_odeiv2_step_bsimp;
  case 8:
    return gsl_odeiv2_step_rk1imp;
  case 9:
    return gsl_odeiv2_step_msadams;
  case 10:
    return gsl_odeiv2_step_msbdf;
  default:
    return gsl_odeiv2_step_bsimp;
  }
}

int
integrate_ode_using_driver_fixed_step (double t, double t1, double y[], int n_steps,
                                       double h_init, double h_max, double eps_abs,
                                       double eps_rel, void * params, size_t dim, int nderiv, double tout[], double Yout[], int step_type_idx)
{
  /* Nderiv can be 0, 1 or 2 */
  size_t i; /* Counter in macro-step loop */
  size_t j; /* Counter in print loop */
  size_t k; /* Counter in dfdy loop */
  double temp;
  int status;
  double ti;
  double dt = (t1-t)/(double)(n_steps-1);
  /* gsl_odeiv2_step_type * T = gsl_odeiv2_step_msbdf; */
  gsl_odeiv2_step * s = gsl_odeiv2_step_alloc (get_step_type(step_type_idx), dim);
  /* gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (eps_abs, eps_rel); */
  /* gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc (dim); */

  gsl_odeiv2_system sys = {func, jac, dim, params};

  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_msbdf, h_init, eps_abs, eps_rel);
  gsl_odeiv2_step_set_driver(s, d);

  gsl_block * f = gsl_block_calloc(dim);
  gsl_block * dfdt = gsl_block_calloc(dim);
  gsl_matrix *dfdy = gsl_matrix_calloc(dim, dim); // init to zero (calloc)

  if (h_max > 0.0)
    {
      gsl_odeiv2_driver_set_hmax(d, h_max);
    }

  for (i = 0; i < n_steps; ++i)
    {
      tout[i] = t;
      if (nderiv > 0)
        func(t, y, f->data, params);
      if (nderiv > 1)
        jac(t, y, dfdy->data, dfdt->data, params);

      for (j = 0; j < dim; ++j)
        {
          Yout[i*dim*(nderiv+1)+j*(nderiv+1)+0] = y[j];
          if (nderiv > 0)
            Yout[i*dim*(nderiv+1)+j*(nderiv+1)+1] = f->data[j];
          if (nderiv > 1)
            {
              temp = 0;
              for (k=0; k<dim; ++k)
                {
                  temp += gsl_matrix_get(dfdy, j, k)*f->data[k];
                }
              Yout[i*dim*(nderiv+1)+j*(nderiv+1)+2] = dfdt->data[j]+temp;
            }
        }
      /* Macro-step loop */
      ti = t + dt;//*(i+1);
      status = gsl_odeiv2_driver_apply (d, &t, ti, y);

      if (status != GSL_SUCCESS)
        {
          printf ("error, return value=%d\n", status);
          break;
        }

    }

  /* Memory management */
  gsl_odeiv2_driver_free (d);
  gsl_odeiv2_step_free (s);
  gsl_block_free(f);
  gsl_block_free(dfdt);
  gsl_matrix_free(dfdy);

  return status;
}


int
integrate_ode_using_driver_fixed_step_print(double t, double t1, double y[], int n_steps,
                                            double h_init, double h_max, double eps_abs,
                                            double eps_rel, void *params, size_t dim, int nderiv,
					    int step_type_idx)
{
  int status;
  int i,j,k;
  double * tout;
  double * Yout;

  tout = malloc(sizeof(double)*n_steps);
  Yout = malloc(sizeof(double)*n_steps*dim*(nderiv+1));
  status = integrate_ode_using_driver_fixed_step(t, t1, y, n_steps, h_init, h_max, eps_abs,
                                                 eps_rel, params, dim, nderiv, tout, Yout, 
						 step_type_idx);
  if (status != GSL_SUCCESS)
    {
      return status;
    }

  for (i = 0; i < n_steps; ++i)
    {
	  printf(STRINGIFY(PRECISION), tout[i]);
	  for (j = 0; j < dim; ++j)
	    {
          for (k = 0; k<=nderiv; ++k)
            {
              printf(" " STRINGIFY(PRECISION), Yout[i*dim*(nderiv+1)+j*(nderiv+1)+k]);
            }
	    }
	  printf("\n");
    }

  free(tout);
  free(Yout);
  return status;
}

