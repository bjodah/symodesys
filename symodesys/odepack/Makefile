# GFortran
F90 = gfortran
F90FLAGS = -Wall -std=f2008 -Wextra -Wimplicit-interface -fPIC

# Debug flags:
#F90FLAGS += -g -fcheck=all
# Release flags:
#F90FLAGS += -O3 -march=native -ffast-math -funroll-loops
F77 = gfortran
F77FLAGS = -Wall -Wextra -Wimplicit-interface -fPIC
#F77FLAGS += -O3 -march=native -ffast-math -funroll-loops

# Intel ifort
#F90 = ifort-12.0.191
#F90FLAGS = -stand f95 -warn all
# Debug flags:
#F90FLAGS += -check all
# Release flags:
#F90FLAGS += -xHOST -O3 -no-prec-div -static

SPKG_LOCAL=/usr
#GCC
CC = gcc
#CFLAGS = -Wall -fPIC -I$(SPKG_LOCAL)/include/${PYTHONLIB} -I$(SPKG_LOCAL)/lib/${PYTHONLIB}/site-packages/numpy/core/include

# 'gcc -pthread',
# '-fno-strict-aliasing',
# '-DNDEBUG -g -fwrapv -O2 -Wall -Wstrict-prototypes',
# '-fPIC'

ODEPACKSRC = opkda1.f opkda2.f opkdmain.f 
ODEPACKOBJ = opkda1.o opkda2.o opkdmain.o
SYSOBJ=types.o ode.o lsodes_bdf.o lsodes_bdf_wrapper.o
PYTHONLIB=python2.7

CFLAGS = -pthread -fno-strict-aliasing -DNDEBUG -g -fwrapv -O2 -Wall -Wstrict-prototypes -fPIC -I$(SPKG_LOCAL)/include/${PYTHONLIB} -I$(SPKG_LOCAL)/lib/${PYTHONLIB}/site-packages/numpy/core/include



pylsodes_bdf.so: ${ODEPACKOBJ} ${SYSOBJ} pylsodes_bdf.o
	${F90} -shared -o $@ $^ -l${PYTHONLIB}

clean:
	rm *.o *.c *.so *.mod

.PHONY: clean


ode.f90: ode_template.f90 variables.pkl
	python enmako.py -p variables.pkl -o $@ $<

ode.o: ode.f90
	$(F90) $(F90FLAGS) -c $<

lsodes_bdf_wrapper.f90: lsodes_bdf_wrapper_template.f90 variables.pkl
	python enmako.py -p variables.pkl -o $@ $<

lsodes_bdf_wrapper.o: lsodes_bdf_wrapper.f90
	$(F90) $(F90FLAGS) -c $<

variables.pkl: mk_pickle.py
	python mk_pickle.py variables.pkl


pylsodes_bdf.c: pylsodes_bdf.pyx
	cython $<

%.o: %.f
	$(F77) $(F77FLAGS) -c $<

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

%.o: %.c
	$(CC) $(CFLAGS) -c $<
