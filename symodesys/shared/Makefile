CC=gcc
CFLAGS=-O2 -funroll-loops -Wall -std=c99 -DGSL_RANGE_CHECK_OFF -DHAVE_INLINE
LIBS=-lgsl -lgslcblas -lm -lc

plot.png: plot.py data.dat
	python $< -o $@ data.dat

data.dat: main_ex
	./$< > $@

main_ex: main_ex.c ode.o drivers.o
	${CC} ${CFLAGS} -o $@ $^ ${LIBS}

%.o: %.c
	${CC} ${CFLAGS} -c $^

