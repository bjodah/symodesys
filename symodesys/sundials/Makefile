CC=gcc
CFLAGS=-g -O0 -funroll-loops -Wall -std=c99
LIBS=-lm -lsundials_cvode -llapack -lsundials_nvecserial
INC=/usr/local/include


plot.png: plot.py data.dat
	python $< -o $@ data.dat

.PHONY: clean

clean:
	rm *.o
	rm main_ex

data.dat: main_ex
	./$< > $@

main_ex: dense_jac.o band_jac.o func.o drivers.o symodesys_util.o main_ex.c
	${CC} ${CFLAGS} -o $@ $^ ${LIBS}

%.o: %.c
	${CC} ${CFLAGS} -c $^

