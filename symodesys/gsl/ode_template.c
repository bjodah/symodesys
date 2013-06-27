#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <math.h>

// Python Mako template of C file
// Variables: f, cse_func
// Variables: jac, dfdt, NY, cse_jac
// CSE tokens: cse%d


int
func (double t, const double y[], double f[], void * params)
{
  /*
    Best is to name all parameters k[0] ... k[P]
   */
  /*int i;*/
  const double *k = (double *) params;

  /*
    Define variables for common subexpressions
   */
% for cse_token, cse_expr in cse_func:
  const double ${cse_token} = ${cse_expr};
% endfor

  /*
    Assign derivatives
   */
% for i, expr in enumerate(f):
  f[${i}] = ${expr};
% endfor

  return GSL_SUCCESS;
}


int
jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
  const double *k = (double *) params;
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, ${NY}, ${NY});
  gsl_matrix *m = &dfdy_mat.matrix;

  /*
    Define variables for common subexpressions
   */
% for cse_token, cse_expr in cse_jac:
  const double ${cse_token} = ${cse_expr};
% endfor

  /*
    Populate the NY times NY Jacobian matrix
   */
% for (i, j), expr in jac:
    gsl_matrix_set (m, ${i}, ${j}, ${expr});
% endfor


  /*
    Populate the array dfdt of length NY
   */
% for i, expr in dfdt:
  dfdt[${i}] = ${expr};
% endfor

  return GSL_SUCCESS;
}
