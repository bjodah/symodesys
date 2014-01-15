// ${_warning_in_the_generated_file_not_to_edit}
#include <math.h>
#include "drivers.h" 

int
dense_jac (DIM_T N, realtype t, N_Vector y, N_Vector fy,
	   DlsMat J, void *params, N_Vector tmp1,
	   N_Vector tmp2, N_Vector tmp3)
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
  DENSE_ELEM(J, ${i}, ${j}) = ${expr};
% endfor

  return 0;
}
