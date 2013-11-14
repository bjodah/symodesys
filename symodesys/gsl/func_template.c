#include <math.h>

#include "drivers.h" 

// Python Mako template of C file


int
func (realtype t, N_Vector y, N_Vector f, void * params)
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
  NV_Ith_S(f, ${i}) = ${expr};
% endfor

  return 0;
}
