#ifndef _DRIVERS_H_
#define _DRIVERS_H_

/* http://stackoverflow.com/questions/2740039/ \ */
/* using-c-preprocessor-to-construct-a-string-literal-for-scanf */
#define STR_EVALUATE(x)   #x
#define STRINGIFY(x)      STR_EVALUATE(x)

#define PRECISION %.9e

gsl_odeiv2_step_type *
get_step_type(int index);

int
integrate_ode_using_driver_fixed_step (double t, double t1, double y[], int n_steps,
                                       double h_init, double h_max, double eps_abs,
                                       double eps_rel, void * params, size_t dim, int order, double tout[], double Yout[], int step_type_idx);

int
integrate_ode_using_driver_fixed_step_print(double t, double t1, double y[], int n_steps,
                                            double h_init, double h_max, double eps_abs,
                                            double eps_rel, void *params, size_t dim, int order,
					    int step_type_idx);


#endif
