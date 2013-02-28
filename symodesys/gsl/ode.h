#ifndef ODE_H_INCLUDED
#define ODE_H_INCLUDED

// http://stackoverflow.com/questions/2740039/ \
// using-c-preprocessor-to-construct-a-string-literal-for-scanf
#define STR_EVALUATE(x)   #x
#define STRINGIFY(x)      STR_EVALUATE(x)

#define PRECISION %.9e


int
integrate_ode_using_driver (double t, double t1, double y[], int n_steps,
			    double h_init, double h_max, double eps_abs,
			    double eps_rel, void *params, int print_values);

#endif
