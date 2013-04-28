/* #include <stdio.h> */
#include <gsl/gsl_errno.h>
#include <math.h>

// Python Mako template of C file
// Variables: f, cse_func
// CSE tokens: cse%d

int
func (double t, const double y[], double f[], void * params)
{
  /*
    Best is to name all parameters k[0] ... k[P]
   */
  /*int i;*/
  double *k = (double *) params;

  /*
    Define variables for common subexpressions
   */
% for cse_token, cse_expr in cse_func:
  double ${cse_token};
% endfor

  /*
    Calculate common subexpressions
   */
% for cse_token, cse_expr in cse_func:
  ${cse_token} = ${cse_expr};
% endfor

  /*
    Assign derivatives
   */
% for i, expr in enumerate(f):
  f[${i}] = ${expr};
% endfor

  return GSL_SUCCESS;
}
