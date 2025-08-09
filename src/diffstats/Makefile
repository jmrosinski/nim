FFLAGS = -g -O2
CFLAGS = -DFORTRANUNDERSCORE -g -O2

OBJS = check_header.o DiffStats.o jrdiff_vectors.o parse_hdr.o

diffstats: $(OBJS)
	$(FC) -o $@ $(OBJS) $(FFLAGS)

%.o: %.F90
	$(FC) -c $(FFLAGS) $<

%.o: %.c
	$(CC) -c $(CFLAGS) $<

clean:
	$(RM) $(OBJS) diffstats

check_header.o:
DiffStats.o:
jrdiff_vectors.o:
parse_hdr.o:
