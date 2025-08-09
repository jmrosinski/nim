# cntl Makefile (for NIM compiled for the CPU)

include ../macros.make

CNTLSRCS = $(wildcard *.F90)
CNTLOBJS = $(addsuffix .o,$(basename $(CNTLSRCS)))

all:
ifeq ($(PAR),serial)
	$(MAKE) -f Makefile.serial $(CNTLOBJS)
else
	$(MAKE) -f Makefile.sms    $(CNTLOBJS)
endif

DEPENDENCIES:
	$(RM) Filepath Srcfiles
	echo "." > Filepath
	ls *.F90 | grep -v sms__decomp > Srcfiles
	$(MKDEPENDS) -m -d sms__decomp.o Filepath Srcfiles > $@

-include DEPENDENCIES

clean:
	$(RM) DEPENDENCIES *.f90 *.mod *.o *.sms
