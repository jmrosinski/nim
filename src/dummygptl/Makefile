# dummygptl Makefile

include ../macros.make

SRCS = $(wildcard *.F90)
OBJS = $(addsuffix .o, $(basename $(SRCS)))

all: $(OBJS)

%.o: %.F90
	$(FCserial) -c $(MICFLAGS) $<

clean:
	$(RM) *.o
