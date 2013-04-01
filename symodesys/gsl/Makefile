
ode.so: drivers.c func.c jac.c
	gcc -fpic -O3 -c -Wall -std=c99 -DGSL_RANGE_CHECK_OFF drivers.c
	gcc -fpic -O3 -c -Wall -std=c99 -DGSL_RANGE_CHECK_OFF func.c
	gcc -fpic -O3 -c -Wall -std=c99 -DGSL_RANGE_CHECK_OFF jac.c
	ld -shared -soname ode.so -o ode.so -lgsl -lgslcblas -lm -lc $^

main_ex: main_ex.c ode.so
	gcc -DGSL_RANGE_CHECK_OFF -std=c99 $^ -O3 -o $@
