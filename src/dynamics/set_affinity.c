// Need to define _GNU_SOURCE sometimes to get value of CPU_SETSIZE
#define _GNU_SOURCE

#include <stdio.h>
#include <sched.h>

int set_affinity_ (int *rank, int *cpu_cores_per_node, int *max_compute_tasks_per_node, int *pin_to_single_core, int *root_on_socket1)
{
  int mincore;
  int maxcore;
  int cores_per_socket = *cpu_cores_per_node/2;
  int tasks_per_socket = *max_compute_tasks_per_node/2;
  int core;
  int retval = 0;
  cpu_set_t setmask;

  if(*pin_to_single_core) {
    // Pin to a sincle core
    mincore = *rank % (*max_compute_tasks_per_node);
    if (mincore >= tasks_per_socket){
      mincore = mincore + cores_per_socket - tasks_per_socket;
    }
    if (*root_on_socket1 ) {
      if (mincore < tasks_per_socket)
        mincore = mincore + cores_per_socket;
      else
        mincore = mincore - cores_per_socket;
    }
    maxcore = mincore;

  }else {
    // Pin to a socket
    if (*rank % *max_compute_tasks_per_node < tasks_per_socket){
      mincore = 0;
      if (*root_on_socket1 ) mincore = cores_per_socket;
    }else{
      mincore = cores_per_socket;
      if (*root_on_socket1 ) mincore = 0;
    }
    maxcore = mincore + cores_per_socket-1;
  }      

  CPU_ZERO (&setmask);
  for (core = mincore; core <= maxcore; ++core)
    CPU_SET (core, &setmask);
  if (sched_setaffinity (0, sizeof (setmask), &setmask) < 0) {
    printf ("setaff: bad return from sched_setaffinity\n");
    retval = -1;
  }
return retval;
}
