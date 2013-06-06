.SUFFIXES: .o .c .f .cxx
options=-DLINUX -DBLAS -O3
CC=/usr/bin/g++
cc=/usr/bin/gcc
f77=gfortran
objects=ran1.o four1.o prod.o \
plgndr.o centroidstuff.o inter.o  sturmi.o inverr.o inverm.o isoev.o \
lancbis.o bisec.o genran.o scalar.o peckeris.o trivec.o hv.o eiganl.o
h2oh2: h2oh2.o $(objects) matvec.o BF.h
	$(CC)  -o h2oh2.x h2oh2.o matvec.o $(objects) \
-lblas -llapack -lgfortran
clean:
	rm *.o
.cxx.o: 
	$(CC)  -c  $(options) $<
.f.o:
	$(f77) -c  $<
.c.o:
	$(cc) -c $(options) $<
