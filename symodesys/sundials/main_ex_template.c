#include <stdlib.h>
#include "drivers.h"

// Python Mako template of C file
// Variables: Y0_COMMA_SEP_STR, PARAM_VALS_COMMA_SEP_STR

// Note that analytic expressions are not evaluated here.


int
main (void)
{
  int status;
  int n = 10;
  size_t dim = ${NY};
  int nderiv = 2;
  double	t	    = 0.0;
  double	t1	    = 10.0;
  double	h	    = 1e-6;
  double	hmax    = 1e-2;
  double	eps_abs	= 1e-6;
  double	eps_rel	= 1e-6;
  double y[]    = {${y0_comma_sep_str}}; // ${y0_names}
  double params[] = {${param_vals_comma_sep_str}}; // ${param_names}
  int step_type_idx = 1;

  status = integrate_fixed_step_print(
    t, t1, y, n, h, hmax, eps_abs, eps_rel, params, dim, nderiv,
    step_type_idx, DENSE_MODE);

  if (status != CV_SUCCESS)
      return 1;
  return 0;

}
