F90=gfortran
FFLAGS=-g -O2 -fbounds-check -fcheck=all
LIBS=-llapack 

all: solver

solver: constants.o connectivity.o boundary.o leastSquares.o flow.o flowVars.o main.o   
	$(F90) $(FFLAGS) main.o flowVars.o flow.o leastSquares.o boundary.o connectivity.o constants.o -o solver $(LIBS)

constants.o: constants.f90
	$(F90) -c $(FFLAGS) constants.f90 

connectivity.o: connectivity.f90
	$(F90) -c $(FFLAGS) connectivity.f90 

boundary.o: boundary.f90
	$(F90) -c $(FFLAGS) boundary.f90 

leastSquares.o: leastSquares.f90
	$(F90) -c $(FFLAGS) leastSquares.f90 

flow.o: flow.f90
	$(F90) -c $(FFLAGS) flow.f90 

flowVars.o: flowVars.f90
	$(F90) -c $(FFLAGS) flowVars.f90 

main.o: main.f90
	$(F90) -c $(FFLAGS) main.f90

clean:
	rm -rf *.mod *.o
