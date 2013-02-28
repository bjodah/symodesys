
ode.so: ode.c func.c jac.c
	gcc $^ -O3 -lgsl -lgslcblas -lm -o $@

main: main.c ode.so
	gcc $^ -O3 -lgsl -lgslcblas -lm -o $@
