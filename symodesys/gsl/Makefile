CFLAGS=-O2 -funroll-loops -Wall -std=c99 -DGSL_RANGE_CHECK_OFF
CC=gcc
LIBS=-lgsl -lgslcblas -lm -lc

main_ex: main_ex.c ode.o drivers.o
	${CC} ${CFLAGS} -o $@ $^ ${LIBS}

%.o: %.c
	${CC} ${CFLAGS} -c $^

