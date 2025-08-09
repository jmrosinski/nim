#include <stdio.h>
#include <unistd.h>  // gethostname, getopt
#include <sched.h>
#ifdef _OPENMP
#include <omp.h>
#endif

extern void runnable (cpu_set_t *, int *, int *);

void print_affinity_ (int *rank)
{
  char hnbuf[64];
  int thread = 0;
  int lo;
  int hi;
  cpu_set_t coremask;

  (void) gethostname (hnbuf, sizeof (hnbuf));
#pragma omp parallel private (thread, coremask, lo, hi)
  {
#ifdef _OPENMP
    thread = omp_get_thread_num ();
#endif
    /*
    ** Disable this because nvcc 21.2 cannot find sched_getaffinity
    // Passing zero means use the calling process
    (void) sched_getaffinity (0, sizeof (coremask), &coremask);
    runnable (&coremask, &lo, &hi);
#pragma omp critical
    {
      printf ("Rank %d thread %d on %s. (Runnable range: lo=%d hi=%d)\n",
	      *rank, thread, hnbuf, lo, hi);
      fflush (stdout);
    }
    */
  }
}
