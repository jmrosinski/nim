ROOTDIR   = $(PWD)
ARCH      = theiapgi
HW        = gpu
PAR       = serial
MAXNZ     = 96
OPENMP    = no
GPTLDIR   = /home/rosinski/install/gptl_cuda
GPTLINCLUDE = -I$(GPTLDIR)/include
DOUBLE_PRECISION = false
DEFINES = 
DEFINES += -DSERIAL
DEFINES += -DALIGN_OK

BUILD_DEFINES = -DHIDE_LONGLINES -DPGIOPENACC -DGPU_RESIDENT -D_OPENACC -DPAR_WRK=3 -DVEC_LEN=32 -DNZ=96
#JR Disable affinity stuff for GPU
BUILD_DEFINES += -UENABLE_AFFINITY

ICOSIO_FC   = nvfortran
ICOSIO_ARGS = CPP=cpp CPP_FLAGS="-C -P -nostdinc" DEFS="-DSERIAL -DNOGRIB" FC=$(ICOSIO_FC) FFLAGS="-g -O3"
SRCDIR      = $(ROOTDIR)

SMSLIBNAME = smsgpu
GMAKEMINUSJ = -j1
LDFLAGS     =

CC = nvcc
CFLAGS = -g -O3
CPP        = gfortran -E
CPP_OPTS   = -C -P -nostdinc
CPP_FLAGS  = $(CPP_OPTS) $(GPTLINCLUDE) $(DEFINES) $(BUILD_DEFINES)
FC         = nvfortran
FCserial   = $(FC)
FCfrontend = nvfortran
ACCFLAGS = -acc -Minfo=accel -gpu=ccnative,maxregcount:96 -cuda
FFLAGS        = $(DYNFLAGS)
UTILSDIR      = ../utils
ICOSIODIR     = ../icosio
CNTLDIR       = ../cntl
MODDIRS       = -I$(UTILSDIR) -I$(ICOSIODIR) -I$(CNTLDIR)

DYN_ACCFLAGS = $(ACCFLAGS)
DEBUG = no
ifeq ($(DEBUG),yes)
  PREPFLAGS  = -g -O0
  DYNFLAGS   = -g -O0 -mcmodel=medium -traceback -Meh_frame $(BUILD_DEFINES) $(DYN_ACCFLAGS) -Mchkfpstk -Mchkptr -Ktrap=fp -Mdclchk -Mbounds -Mchkstk
else
  PREPFLAGS  = -g -O3 -mcmodel=medium
  DYNFLAGS   = -g -O3 -mcmodel=medium -traceback -Meh_frame $(BUILD_DEFINES) $(DYN_ACCFLAGS)
endif
LDFLAGS =
LDFLAGS += $(DYNFLAGS) -lunwind
R8FLAG = -r8

GPTLLIB = -cuda -L$(GPTLDIR)/lib -lgptlf -lgptl -lgptl_cuda -lstdc++

NVCC = nvcc
NVCC_FLAGS = -pg -arch=sm_61 -ftz=true -fmad=false -w
ifeq ($(DEBUG),yes)
  NVCC_FLAGS += -g -G 
endif

DIVIDES_IN_SOLVEITHLS = yes
ifeq ($(DIVIDES_IN_SOLVEITHLS),no)
  BUILD_DEFINES += -DSOLVEITHLS_RECIPROCAL
endif

FCX    = $(FCserial)
SFLAGS =

ICOSIOOBJS     = $(shell ls $(ICOSIODIR)/*.o)
UTILSOBJS      = $(UTILSDIR)/taskinfo.o $(UTILSDIR)/namelistdata.o $(UTILSDIR)/wrappers.o
NIMEXE = nim

MKDEPENDS = $(SRCDIR)/tools/mkDepends
