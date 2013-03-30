
ode.so: ode.c func.c jac.c
	gcc -fpic -O3 -c -Wall ode.c
	gcc -fpic -O3 -c -Wall func.c
	gcc -fpic -O3 -c -Wall jac.c
	ld -shared -soname ode.so -o ode.so -lgsl -lgslcblas -lm -lc $^

main_ex: main_ex.c ode.so
	gcc $^ -O3 -o $@
