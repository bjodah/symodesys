#ifndef _SYMODESYS_UTIL_H_
#define _SYMODESYS_UTIL_H_

#include <stdio.h>

/* http://stackoverflow.com/questions/2740039/ \ */
/* using-c-preprocessor-to-construct-a-string-literal-for-scanf */
#define STR_EVALUATE(x)   #x
#define STRINGIFY(x)      STR_EVALUATE(x)

#define PRECISION %.9e

void print_state(double t, size_t dim, int nderiv, size_t idx, double * yout);

#endif // _SYMODESYS_SUNDIALS_DRIVERS_H_
