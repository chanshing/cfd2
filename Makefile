OPT = 
# GFORTRAN
GWARN = -Wall -Wno-unused-variable -Wno-unused-dummy-argument
GOMP = $(OPT) -fopenmp 
GFLAGS = $(OPT) -O3 $(GWARN)
GDEBUG = $(OPT) -ffpe-trap=zero,invalid,overflow \
		 -fbacktrace -fbounds-check $(GWARN)
# IFORT
IFORTWARN = -warn -warn nodeclarations -warn nounused
IFORTOMP = $(OPT) -openmp 
IFORTFLAGS = $(OPT) -fast $(IFORTWARN)
IFORTDEBUG = $(OPT) -check all -traceback -fpe0 $(IFORTWARN)
IDBFLAGS = $(OPT) -debug -O0 -fpe0 $(IFORTWARN)

OBJECTS = util_mod.o\
		  mesh_mod.o\
		  general_mod.o\
		  bodies_mod.o\
		  IO_mod.o\
		  runge_kutta_mod.o\
		  iter_solver_mod.o\
		  solver_mod.o\
		  main.o

.PHONY: clean

default: FC = ifort
default: LFLAGS = 
default: CFLAGS = $(IFORTFLAGS)
default: ns

omp: FC = ifort
omp: LFLAGS = $(IFORTOMP)
omp: CFLAGS = $(IFORTFLAGS) $(IFORTOMP)
omp: ns

idb: FC = ifort
idb: LFLAGS =
idb: CFLAGS = $(IDBFLAGS)
idb: ns

debug: FC = ifort
debug: LFLAGS =
debug: CFLAGS = $(IFORTDEBUG)
debug: ns

gnu: FC = gfortran
gnu: LFLAGS =
gnu: CFLAGS = $(GFLAGS)
gnu: ns

gomp: FC = gfortran
gomp: LFLAGS = $(GOMP)
gomp: CFLAGS = $(GFLAGS) $(GOMP)
gomp: ns

gdebug: FC = gfortran
gdebug: LFLAGS =
gdebug: CFLAGS = $(GDEBUG)
gdebug: ns

ns: $(OBJECTS)
	$(FC) $(LFLAGS) $(OBJECTS) -o ns

%.o: %.f90
	$(FC) $(CFLAGS) -c $<

clean:
	rm -f $(OBJECTS) ns *.mod
