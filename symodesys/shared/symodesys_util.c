#include "symodesys_util.h"

void print_state(double t, size_t dim, int nderiv, size_t idx, double * yout){
  printf(" " STRINGIFY(PRECISION), t);
  for (size_t j = 0; j < dim; ++j)
    {
      for (int k = 0; k<=nderiv; ++k)
	{
	  printf(" " STRINGIFY(PRECISION), yout[idx*dim*(nderiv+1)+j*(nderiv+1)+k]);
	}
    }
  printf("\n");  
}
