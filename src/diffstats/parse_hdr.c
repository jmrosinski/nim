/*
** Routines in this file parse the first 3 NIM header records.
** Since the Fortran "write" statement used i0 format, C was the only way I could
** manage to handle things like glvl= 4 (1 char) as well as glvl= 10 (2 chars),
** and also levels= 96 (2 chars) as well as levels= 192 (3 chars).
**
** Jim Rosinski November 15, 2013
*/

#include <string.h>
#include <stdio.h>

#if ( defined FORTRANUNDERSCORE )
#define parse_hdr_1 parse_hdr_1_
#define parse_hdr_2 parse_hdr_2_
#define parse_hdr_3 parse_hdr_3_
#elif ( defined FORTRANDOUBLEUNDERSCORE )
#define parse_hdr_1 parse_hdr_1__
#define parse_hdr_2 parse_hdr_2__
#define parse_hdr_3 parse_hdr_3__
#endif

void parse_hdr_1 (const char *hdr, char varname[4], char yyyymmddhhmm[12], 
		  int nc1, int nc2, int nc3)
{
  char vn[5];
  char yy[13];

  sscanf (hdr, "%*s %4s %*s %*s %*s %*s %12s",
	            vn,                 yy);
  if (nc2 != 4) {
    printf ("parse_hdr_1: expected length of 1st var=4 got %d\n", nc2);
    return;
  }
  strncpy (varname, vn, nc2);

  if (nc3 != 12) {
    printf ("parse_hdr_1: expected length of 2nd var=12 got %d\n", nc3);
    return;
  }
  strncpy (yyyymmddhhmm, yy, nc3);

  return;
}

void parse_hdr_2 (const char *hdr, int *nz, int *glvl, char timeunit[2], int *its, int *time,
		  int nc1, int nc2)
{
  char tu[2];

  sscanf (hdr, "%*s %d%*s %*s %d%*s %*s %*s %*s %2s%*s %*s %d%*s %d",
	  nz,       glvl,                       tu,        its,  time);
  if (nc2 != 2) {
    printf ("parse_hdr_2: expected length of 1st var=2 got %d\n", nc2);
    return;
  }
  strncpy (timeunit, tu, nc2);
  return;
}

void parse_hdr_3 (const char *hdr, int *dim1, int *nip, 
		  int nc1)
{
  sscanf (hdr, "%*s %d%*s %*s %d",
	            dim1,     nip);
  return;
}
