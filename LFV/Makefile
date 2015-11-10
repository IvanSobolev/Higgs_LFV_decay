#!/bin/bash

#   The directory where the include files are.
#   By default it was installed into your home directory.
prefix = ${HOME}
name1 = LFV3_10F
name2 = LFV3_10F_new
#   include and lib directories
includedir = $(prefix)/include/cteq/
libdir = $(prefix)/lib

#   use you favorite fortran and C compilers 
#   (note, the package was tested with GNU compilers)
F77 = gfortran
CC = gcc

all: LFV3 LFV3_new

LFV3: Cteq6Pdf-2010.f LFV3.o nrutil.o jacobi.o spline.o splint.o
	$(F77) -O3 Cteq6Pdf-2010.f LFV3.o -L$(libdir) nrutil.o jacobi.o spline.o splint.o -o $(name1) -lm
LFV3_new: Cteq6Pdf-2010.f LFV3_new.o nrutil.o jacobi.o spline.o splint.o random.o
	$(F77) -O3 Cteq6Pdf-2010.f LFV3_new.o -L$(libdir) nrutil.o jacobi.o spline.o splint.o random.o -o $(name2) -lm
LFV3.o: LFV3.c
	gcc -c -O3 -I$(includedir) -I. LFV3.c
LFV3_new.o: LFV3_new.c
	gcc -c -O3 -I$(includedir) -I. LFV3_new.c
nrutil.o: nrutil.c
	gcc -c nrutil.c -I.
spline.o: spline.c
	gcc -c spline.c -I.
splint.o: splint.c
	gcc -c splint.c -I.
jacobi.o: jacobi.c
	gcc -c jacobi.c -I.
random.o: random.c
	gcc -c random.c
clean:
	rm -rf *.o $(name1) $(name2) rm *~
