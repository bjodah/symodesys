#ifndef _ODE_H_
#define _ODE_H_

int func (double t, const double y[], double f[], void * params);
int jac (double t, const double y[], double *dfdy, double dfdt[], void *params);

#endif /* _FUNC_H_ */
