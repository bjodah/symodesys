#include <math.h>

#include "drivers.h" 

// Python Mako template of C file

int
band_jac (SIZE_T N, SIZE_T mu, SIZE_T ml,
	    realtype t, N_Vector u, N_Vector fu, 
	    DlsMat J, void *params,
	    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  const double *k = (double *) params;

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
  BAND_ELEM(J, ${i}, ${j}) = ${expr};
% endfor

  return 0;
}
