# icosio Makefile

# DEFS may include: -DSERIAL for no-MPI build
#                   -DNOGRIB for no-GRIB build

.PRECIOUS: %.f90

%.f90: %.F90
	$(CPP) $(CPP_FLAGS) $(DEFS) $< > $@

%.o: %.f90
	$(FC) -c $(FFLAGS) $<

all: icosio.o

icosio.o: stringlist.o

clean:
	$(RM) *.f90 *.o *.mod
