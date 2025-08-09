# utils Makefile

include ../macros.make

BINDIR     = ../bin
EXES       = $(addprefix $(BINDIR)/,ConvData GetComputeTasks GetGLVL GetNZ GetPath GetPhysics GetQueueTime GetRestart)
CNTLF90S   = kinds.F90 ReadNamelist.F90 ToUpper.F90
COMMONOBJS = namelistdata.o kinds.o ReadNamelist.o headers.o ToUpper.o
#JR Need to filter out -xavx or -xAVX if present so don't have to be on a
#JR SNB to submit a SNB job
FFLAGS_LOC = $(filter-out -xavx -xAVX,$(shell echo "$(FFLAGS)"))

%.o: %.F90
	$(FCfrontend) -c $(DEFS) $(FFLAGS_LOC) $(MICFLAGS) $<

all: $(EXES) $(BINDIR)/get_task_info

$(EXES): $(BINDIR)/%: $(notdir %.F90) $(COMMONOBJS)
	$(FCfrontend) $(FFLAGS_LOC) $(MICFLAGS) $^ -o $@

$(BINDIR)/get_task_info: get_task_info.F90 $(COMMONOBJS) taskinfo.o taskinfo_comm.o wrappers.o
	$(FCfrontend) $(FFLAGS_LOC) $(MICFLAGS) $^ -o $@

wrappers.o: wrappers.F90
	$(CPP) $(CPP_FLAGS) $< > $*.f90
	$(FCfrontend) -c $(FFLAGS_LOC) $(MICFLAGS) $*.f90

taskinfo.o: taskinfo.F90 namelistdata.o
	$(CPP) $(CPP_FLAGS) $< > $*.f90
	$(FCfrontend) -c $(FFLAGS_LOC) $(MICFLAGS) $*.f90

kinds.o: kinds.F90
	$(CPP) $(CPP_FLAGS) $< > $*.f90
	$(FCfrontend) -c $(FFLAGS_LOC) $(MICFLAGS) $*.f90

taskinfo_comm.o: ../dynamics/taskinfo_comm.F90
	$(CPP) $(CPP_FLAGS) -DSERIAL $< > $*.f90
	$(FCfrontend) -c $(FFLAGS_LOC) $(MICFLAGS) $*.f90

$(CNTLF90S): %: ../cntl/%
	rsync -au $< .

DEPENDENCIES: $(wildcard *.F90) $(CNTLF90S)
	$(RM) Filepath Srcfiles
	echo "." > Filepath
	ls *.F90 > Srcfiles
	$(MKDEPENDS) -m Filepath Srcfiles > $@

-include DEPENDENCIES

clean:
	$(RM) *.o *.f90 *.mod DEPENDENCIES
