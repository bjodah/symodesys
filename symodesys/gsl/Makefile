ode.so: drivers.o func.o jac.o
	ld -shared -soname $@ -o $@ $^ -lgsl -lgslcblas -lm -lc

%.o: %.c
	gcc -fpic -O3 -c -Wall -std=c99 -DGSL_RANGE_CHECK_OFF $^

main_ex: main_ex.c ode.so
	gcc -DGSL_RANGE_CHECK_OFF -std=c99 $^ -O3 -o $@
