#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_pow_int.h>

// Python Mako template of C file
// Variables: ${jac} ${dfdt} ${NY}

int
jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
  double *k = (double *) params;
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, ${NY}, ${NY});
  gsl_matrix *m = &dfdy_mat.matrix;

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
