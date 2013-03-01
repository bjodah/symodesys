#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>

// Python Mako template of C file
// Variables: ${jac} ${dfdt} ${NY} ${cse_jac}

int
jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
  double *k = (double *) params;
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, ${NY}, ${NY});
  gsl_matrix *m = &dfdy_mat.matrix;

  /*
    Define variables for common subexpressions
   */
  % for cse_type, cse_token, cse_expr in cse_jac
        ${cse_type} ${cse_token};
  % endfor


  /*
    Calculate common subexpressions
   */

  % for cse_type, cse_token, cse_expr in cse_jac
        ${cse_token} = ${cse_expr};
  % endfor


  /*
    Populate the NY times NY Jacobian matrix
   */

  % for i, j, expr in jac
        gsl_matrix_set (m, ${i}, ${j}, ${expr});
  % endfor


  /*
    Populate the array dfdt of length NY
   */

  % for i expr in dfdt
      dfdt[${i}] = ${expr};
  % endfor

  return GSL_SUCCESS;
}
