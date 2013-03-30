#include <stdio.h>
#include <gsl/gsl_errno.h>

#include "ode.h"

// Python Mako template of C file
// Variables: params, NY, Y0_COMMA_SEP_STR


int
main (void)
{
  int status;
  double	t	    = 0.0;
  double	t1	    = 10.0;
  double	h	    = 1e-6;
  double	eps_abs	    = 1e-6;
  double	eps_rel	    = 1e-6;
  int		print_values = 1;
  double        y[${NY}]    = {${Y0_COMMA_SEP_STR}};
  % for param_token, param_value in params:
  double ${param_token} = ${param_value};
  % endfor


  status = integrate_ode(t, t1, y, h, eps_abs, eps_rel, &params, print_values);

  if (status == GSL_SUCCESS)
    {
      return 0;
    }

  return 1;

}
