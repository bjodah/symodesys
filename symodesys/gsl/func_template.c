/* #include <stdio.h> */
#include <gsl/gsl_errno.h>

// Python Mako template of C file
// Variables: ${f} ${cse_func} ${analytic_y}
// CSE tokens: cse_*

  /*
    Define analytic f
   */

/* % for y_token, cse_expr in analytic_y */
/* double */
/* ${y_token} (double t, const double y[], void * params) */
/* { */
/*   double *k = (double *) params; */
/*   return ${y_epxr} */
/* } */
/* % endfor */



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
  % for cse_token, cse_expr in cse_func
        double ${cse_token};
  % endfor

  /*
    Calculate common subexpressions
   */

  % for cse_token, cse_expr in cse_func
        ${cse_token} = ${cse_expr};
  % endfor

  /*
    Assign derivatives
   */


  % for i, expr in enumerate(f)
      f[${i}] = ${expr};
  % endfor

  return GSL_SUCCESS;
}
