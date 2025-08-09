# NIM Makefile
include macros.make

all: info build

info:
	$(info Inside top level Makefile)
	$(info ARCH   = $(ARCH))
	$(info HW     = $(HW))
	$(info PAR    = $(PAR))
	$(info MAXNZ  = $(MAXNZ))
	$(info SRCDIR = $(SRCDIR))
	$(info DOUBLE_PRECISION = $(DOUBLE_PRECISION))

build:
	$(MAKE) $(GMAKEMINUSJ) $(ICOSIO_ARGS)               -C icosio
	$(MAKE) $(GMAKEMINUSJ) SRCDIR=$(SRCDIR)             -C utils
	$(MAKE) $(GMAKEMINUSJ)                              -C cntl
	$(MAKE) SRCDIR=$(SRCDIR)                            -C dynamics

clean:
	$(MAKE) -C icosio    clean
	$(MAKE) -C utils     clean
	$(MAKE) -C cntl      clean
	$(MAKE) -C dynamics  clean
	$(RM) bin/* lib/*
