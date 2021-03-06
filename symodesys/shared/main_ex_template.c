// ${_warning_in_the_generated_file_not_to_edit}
<%doc>
  // Python Mako template of C file
  // Variables: y0_comma_sep_str, param_vals_comma_sep_str

  // Note that analytic expressions are not evaluated here.
</%doc>
#include <stdlib.h>
#include <gsl/gsl_errno.h>

#include "drivers.h"



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
  int step_type_idx = 7;

  status = integrate_fixed_step_print(
    t, t1, y, n, h, hmax, eps_abs, eps_rel, params, dim, nderiv,
     step_type_idx);

  if (status == GSL_SUCCESS)
    {
      return 0;
    }

  return 1;

}
