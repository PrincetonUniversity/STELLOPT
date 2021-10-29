/*==============================================================================
 * FILE: clean_hist.c
 *
 * PURPOSE: Goes through history file (or any ASCII data file with first column
 *          as time) and removes extraneous lines due to restart runs.
 *
 *       Code will output initial header, but will not output any subsequent
 *       comments.
 *
 * COMPILE USING: gcc -o clean_hist clean_hist.c
 *
 * USAGE: ./clean_hist <infile>
 *                
 *
 * WRITTEN BY: Denis St-Onge, August 2016
 *============================================================================*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

# define LENGTH 1000
# define MAX_LINES 500000

static void usage(const char *arg);

/* ========================================================================== */

int main(int argc, char* argv[])
{
  /* file variables */
  const char whitespace[] = " \f\n\r\t\v";
  char *in_name;
  char line[LENGTH];
  FILE *fidin;
  /* data variables */
  float t,box_size, thres;
  float x0,y0,z0,x1,y1,z1;
  float times[MAX_LINES], last_good_time;
  int i,expand = 0, first_char;
  long data_pos, nx=0,ny=0,nz=0,nlines=0;

  /* Read Arguments */

  if(argc != 2) usage(argv[0]);

  in_name = argv[1];

  fidin = fopen(in_name,"r");
  if (fidin == NULL){
    fprintf(stderr,"Failed to open input file %s.\n",in_name);
    exit(1);
  }

/* Grab header and write into output */
 
  fgets(line,LENGTH,fidin);
  first_char=strspn(line, whitespace);
  if(line[0] == '\0') exit(0);

  while(line[0] != '\0' && line[first_char] == '#'){
    fprintf(stdout,"%s",line);
    data_pos=ftell(fidin); // This is position in the file where the data starts
    fgets(line,LENGTH,fidin);
    first_char=strspn(line, whitespace);
  }

/* if the run has been restarted, clean the extraneous times */
  while(line[0] != '\0' && !feof(fidin)){
    if(line[first_char] != '#'){
      sscanf(line,"%f", &times[nlines++]);
    }
    fgets(line,LENGTH,fidin);
    first_char=strspn(line, whitespace);
  }

  fseek(fidin,data_pos,SEEK_SET);

  last_good_time=times[nlines-1];
  for(i = nlines-1; i >= 0; i--){
    if(times[i-1] < last_good_time)
      last_good_time=times[i-1];
    else
    times[i-1] = -1; 
  }  
 
  i=0;

  fgets(line,LENGTH,fidin);
  first_char=strspn(line, whitespace);
  while(line[0] != '\0' && !feof(fidin)){
    if(line[first_char] != '#'){
      if(times[i] != -1){
        fprintf(stdout,"%s",line);
      }
      i++;
    }
    fgets(line,LENGTH,fidin);
    first_char=strspn(line, whitespace);
  }
  return 0;
}


/* ========================================================================== */


static void usage(const char *arg)
{
  fprintf(stderr,"Required Usage: %s infile\n", arg);
  exit(0);
}
