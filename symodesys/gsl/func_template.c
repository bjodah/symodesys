/* #include <stdio.h> */
#include <gsl/gsl_errno.h>

// Python Mako template of C file
// Variables: ${f} ${cse_func}
// CSE tokens: cse_*


int
func (double t, const double y[], double f[], void * params)
{
  /*
    Best is to name all parameters k[0] ... k[P]
   */
  int i;
  double *k = (double *) params;

  /*
    Define variables for common subexpressions
   */
  % for cse_type, cse_token, cse_expr in cse_func
        ${cse_type} ${cse_token};
  % endfor


  /*
    Calculate common subexpressions
   */

  % for cse_type, cse_token, cse_expr in cse_func
        ${cse_token} = ${cse_expr};
  % endfor


  % for i expr in f
      f[${i}] = ${expr};
  % endfor

  return GSL_SUCCESS;
}
