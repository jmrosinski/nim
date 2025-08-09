# dynamics Makefile (for NIM dynamics compiled for the CPU)

include ../macros.make
include NIM_DYNAMICS_OBJS

.SMS     = $(SMS)/bin/ppp modinfo
CNTLOBJS = $(addsuffix .o, $(basename $(CNTLSRCS)))
CNTLSRCS = $(wildcard $(CNTLDIR)/*.F90)
ifeq ($(ENABLE_AFFINITY),yes)
  COBJS  = runnable.o print_affinity.o set_affinity.o
endif
LIBS     = $(GPTLLIB)

.PHONY: clean objects1 objects2

%.o: %.f
	$(FCX) $(GPTLINCLUDE) -c $<

all: DEPENDENCIES objects1 objects2 $(NIMEXE)

CPUOBJS = $(ACCOBJS)

objects1:
	$(MAKE) $(GMAKEMINUSJ) -f Makefile.serial.cpu SFLAGS="$(SFLAGS)" $(SEROBJS)
	$(MAKE) $(GMAKEMINUSJ) -f Makefile.serial.cpu SFLAGS="" $(COMOBJS) $(CPUOBJS)

objects2:
ifeq ($(PAR),serial)
	$(MAKE) $(GMAKEMINUSJ) -f Makefile.serial.cpu SFLAGS="" $(SMSOBJS)
else
	$(MAKE) $(GMAKEMINUSJ) -f Makefile.sms $(SMSOBJS)
        LIBS += -L$(SMS)/lib -l$(SMSLIBNAME)
endif
OBJS = $(CNTLOBJS) $(SMSOBJS) $(SEROBJS) $(COMOBJS) $(CPUOBJS)

$(NIMEXE): $(OBJS) $(ICOSIOOBJS) $(UTILSOBJS) $(COBJS) $(GPTLOBJS) $(LIBshare) $(LIBcommon)
	$(FCX) -o $@ $(ACCFLAGS) $(MODDIRS) $(OMPFLAG) $(MICFLAGS) $(OBJS) \
                     $(ICOSIOOBJS) $(UTILSOBJS) $(COBJS) $(GPTLOBJS) \
                     $(LIBcommon) $(LIBshare) $(LIBS) $(LDFLAGS)

DEPENDENCIES:
	$(RM) Filepath Srcfiles
	echo "." > Filepath
	ls *.F90 *.f > Srcfiles
	$(MKDEPENDS) -m Filepath Srcfiles > $@

clean:
	$(RM) DEPENDENCIES .ACC* *.f90 *.m4 *.mod *.o *.sms
