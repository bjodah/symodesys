/* #include <stdio.h> */
#include <gsl/gsl_errno.h>
#include <gsl/gsl_pow_int.h>

// Python Mako template of C file
// Variables: ${f}


int
func (double t, const double y[], double f[], void * params)
{
  /*
    Best is to name all parameters k[0] ... k[P]
   */
  int i;
  double *k = (double *) params;

  % for i expr in f
      f[${i}] = ${expr};
  % endfor

  return GSL_SUCCESS;
}
