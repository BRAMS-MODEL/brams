/*!############################# Change Log ##################################
! 2.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
*/

#include "utils_sub_names.h"

#include <stdio.h>
#include <malloc.h>
#include <math.h>
#ifdef CRAY
#include <stdlib.h>
#endif
/*#include <unistd.h>*/

#ifdef SGI
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#endif

/* DEBUG-ALF */
/* #include <time.h> */

/* Prototypes not needed except for C++
int vfscale(float *, long int, double *, double *);
void *malloc(int);
void *free(void *);
*/

/*********************************************************/

#define nlocrams 1000
void *iaddrams[nlocrams];

int iralloc(int *memtot,int *ia,int *ioff)
{
  extern void *iaddrams[nlocrams];
  int n, ifaddr, imaddr;
  void *iaddr;

/* Allocate the memory */

  iaddr = malloc((*memtot)*sizeof(float));

/* Compute the offset for Fortran */

#ifdef CRAY
  ifaddr = (int)ia;
  imaddr = (int)iaddr;
  *ioff=(imaddr-ifaddr);
#else
  ifaddr = (int)ia;
  imaddr = (int)iaddr;
  *ioff = (imaddr-ifaddr)/sizeof(float);
#endif

/* Find first empty location in address array */

  for(n=0; n < nlocrams ; n++) {
    if( iaddrams[n] == 0 ){
      iaddrams[n] = iaddr;
      break;
    }
  }

 /* printf("in ralloc- %i %i %i %i\n",n,*memtot,iaddr,*ioff);*/
  
/* Return this location to FORTRAN as an address "handle" */

  return(n);
}

/*********************************************************/

int irfree(int *nmem) 
{
  int ans=0;
  extern void *iaddrams[nlocrams];
  
  /* printf("in rfree- %i %i\n",*nmem,iaddrams[*nmem]); */
  free(iaddrams[*nmem]);                   
  iaddrams[*nmem]=0;
  return(ans);
  
}

/* ******************************************************* */

void irsleep(int *seconds)
{
   extern int sleep(int);

#if !defined (PC_NT1)
   sleep( *seconds );
#endif

   return;
}

/* ******************************************************* */

FILE *ramsfile;

#ifdef STARDENT 

int rams_c_open(filename,faccess)
     struct  { char *string; int len; } *filename,*faccess;
{
  extern FILE *ramsfile;
  
  /* printf(" C_open - %s %s \n",filename->string,faccess->string); */
  ramsfile=fopen(filename->string,faccess->string);
  return(0);
}
#else

int rams_c_open(char *filename,char *faccess)

{
  extern FILE *ramsfile;

  /* printf(" C_open-%s \n",filename); */
  ramsfile=fopen(filename,faccess);
 /* perror("rams_c_open"); */
  return(0);
}
#endif

/*********************************************************/

int rams_c_close()
{
  extern FILE *ramsfile;
  int istat;

  istat=fclose(ramsfile);
  return(istat);
}

/*********************************************************/

int rams_c_pos(long int *fbyte)
{ 
  int retcode;
  extern FILE *ramsfile;

  retcode=fseek(ramsfile,*fbyte,0);
  return(retcode);
}

/*********************************************************/

void rams_c_tell(int *pos)
{ 
  extern FILE *ramsfile;

  *pos=ftell(ramsfile);
}

/*********************************************************/

int rams_c_read(int *fbyte,int *numbytes,int *a)
{
  int retcode;
  extern FILE *ramsfile;

  retcode=fseek(ramsfile,*fbyte,0);
  fread(a,1,*numbytes,ramsfile);
  return(retcode);
}

/*********************************************************/
int rams_c_read_char(int *fbyte,int *numbytes,int *a)
 
{
  int retcode;
  extern FILE *ramsfile;

  retcode=fseek(ramsfile,*fbyte,0);
  fread(a,1,*numbytes,ramsfile);
  return(retcode);
}

/*********************************************************/

int rams_c_write(int *fbyte,int *numbytes,int *a)
{
  int retcode;
  extern FILE *ramsfile;

  retcode=fseek(ramsfile,*fbyte,0);
  fwrite(a,1,*numbytes,ramsfile);
  return(retcode);
}

/**********************************************************************/
/*   C versions of vfirecr and vforecr written by Peter Olsson, 1993 */

#include <ctype.h>
#define BITFIELDLENGTH 6
#define SMALL_OFFSET 1.e-20
#define VFORMMASK 63 

void vfirecr(int *unit, float *a, long int *n, char *type, char *b, 
	     long int *irec)    
{
  extern FILE *ramsfile;
  int i, j, nn, nbits, nchs;
  float bias, fact, inverse_fact;
  unsigned vnum, char_count;

/*   fseek( ramsfile, *irec, 0); */
  fseek( ramsfile, *irec, SEEK_SET);
  fread(b, 1, 80, ramsfile);
  sscanf(b,"%ld %d %f %f",&nn, &nbits, &bias, &fact);
  printf("Data found: nn=%ld, nbits=%d, bias=%f, fact=%f\n",
	 nn, nbits, bias, fact);

  inverse_fact = 1./fact;
  nchs=nbits/6;
  fread(b,1,*n * nchs,ramsfile);
  if( nn != *n )
    printf("Word count mismatch on vfirec record\n  Words on record - %d\n  Words expected  - %d\n ",nn,*n);

  for(i = 0, char_count = 0; i < *n; i++)
  {
    for(j = 0,vnum=0; j < nchs; j++, char_count++)
    {
      vnum = vnum << BITFIELDLENGTH;
      if( isdigit( b[char_count] ) ) 
               vnum = vnum | (unsigned)b[char_count] - 48;
      else if( isupper( b[char_count] ) ) 
               vnum = vnum |(unsigned)b[char_count] - 55;
      else 
               vnum = vnum | (unsigned)b[char_count] - 61;
    }
    a[i] = vnum*inverse_fact - bias;
  }
}

/*************************************************************************/

/* vforecr: write an ASCII coded (on nbits) portable version of array a[0:n-1] into */
/*          file ramsfile returning final record position irec.                     */
/*          Procedure arguments scr and cscr are scratch areas of size n.           */
/*          Procedure arguments unit and type are unused.                           */

char vc[65] = 
     "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz{|";

void vforecr(int *unit, float *a, long int *n, int *nbits, float *scr,
            char *cscr, char *type, long int *irec)
{
  extern FILE *ramsfile;
  extern char vc[];
  double amax, amin, bias, fact; 
  int i, j, char_count, nchs;
  float ftemp;
  extern int vfscale(float*, long int, double*, double* );

  vfscale( a, *n, &amin, &amax);

  bias = -amin + SMALL_OFFSET;
  fact = (pow( 2.0, (double)*nbits)-1 ) / ( bias + amax + SMALL_OFFSET);

  fprintf(ramsfile,"%8ld%8d%20.10e%20.10e                        "
	  ,*n,*nbits,bias,fact);

 /* assume for now that transformation is linear  */   

  for( i = 0; i < *n;  i++ )
    scr[i] = ( a[i] + bias ) * fact;

  nchs = *nbits / BITFIELDLENGTH;

  for( j = 0, char_count = 0 ; j < *n; j++ ) {
    ftemp=scr[j];
    for( i = nchs-1; i >= 0; i--, char_count++ )
      cscr[char_count]= 
          vc[( (unsigned)ftemp >> (i * BITFIELDLENGTH) ) & VFORMMASK];
  }

  fwrite(cscr , sizeof( char ), *n * nchs, ramsfile ); 
  *irec=ftell(ramsfile); 
}

/*************************************************************************/

int vfscale(float *a, long int n, double *min, double *max)
{ 
  int i;
  
  *min=  1.e20;
  *max = -*min;  

  for( i = 0; i < n; i++)
  {
    if( a[i] > *max ) *max = a[i];
    if( a[i] < *min ) *min = a[i];
  }
  return(0);
}


/*  ==========================================

    Rewrite basic i/o handling functions
    to allow multiple files simultaneously
    opened, controlled by fortran accessible
    integer "unit" (JP: unrelated to fortran io unit)

    ==========================================*/

#define MAXFILES 10

static FILE *filerams[MAXFILES];
static int  fileused[MAXFILES]={0,0,0,0,0,0,0,0,0,0};

void rams_c_open_u(char *filename, char *faccess, int *unit)

{
  int i;
  for (i=0; i<MAXFILES && fileused[i]; i++); /* search for unused entry */
  if (i < MAXFILES && (filerams[i]=fopen(filename,faccess)) != NULL) {
    fileused[i]=1;
    *unit=i;
  }
  else 
    *unit=-1;       /* return error */
}

void rams_c_close_u(int *unit)
{
  int istat;
  if (*unit >= 0 && *unit < MAXFILES && fileused[*unit]) {
    fileused[*unit]=0;
    *unit=fclose(filerams[*unit]);
  }
  else 
    *unit=-1;
}


void vfirecr_u(int *unit, float *a, long int *n, char *type, char *b, 
	     long int *irec)    
{
  int i, j, nn, nbits, nchs;
  float bias, fact, inverse_fact;
  unsigned vnum, char_count;

  fseek( filerams[*unit], *irec, SEEK_SET);
  fread(b, 1, 80, filerams[*unit]);
  sscanf(b,"%ld %d %f %f",&nn, &nbits, &bias, &fact);
  printf("Data found: nn=%ld, nbits=%d, bias=%f, fact=%f\n",
	 nn, nbits, bias, fact);

  inverse_fact = 1./fact;
  nchs=nbits/6;
  fread(b,1,*n * nchs,filerams[*unit]);
  if( nn != *n )
    printf("Word count mismatch on vfirec record\n  Words on record - %d\n  Words expected  - %d\n ",nn,*n);

  for(i = 0, char_count = 0; i < *n; i++)
  {
    for(j = 0,vnum=0; j < nchs; j++, char_count++)
    {
      vnum = vnum << BITFIELDLENGTH;
      if( isdigit( b[char_count] ) ) 
               vnum = vnum | (unsigned)b[char_count] - 48;
      else if( isupper( b[char_count] ) ) 
               vnum = vnum |(unsigned)b[char_count] - 55;
      else 
               vnum = vnum | (unsigned)b[char_count] - 61;
    }
    a[i] = vnum*inverse_fact - bias;
  }
}


void vforecr_u(int *unit, float *a, long int *n, int *nbits, float *scr,
            char *cscr, char *type, long int *irec)
{
  extern char vc[];
  double amax, amin, bias, fact; 
  int i, j, char_count, nchs;
  float ftemp;
  double den;
  extern int vfscale(float*, long int, double*, double* );

  vfscale( a, *n, &amin, &amax);

  bias = -amin + SMALL_OFFSET;
  den  = bias + amax;
  den += SMALL_OFFSET;
/*   printf("DEBUG-ALF:vforecr:bias=%g,amin=%g,SO=%g,amax=%g,den=%g,den2=%g\n", */
/* 	 bias,amin,SMALL_OFFSET,amax,(bias+amax+SMALL_OFFSET),den); */
/*   fact = (pow( 2.0, (double)*nbits)-1 ) / ( bias + amax + SMALL_OFFSET); */
  fact = (pow( 2.0, (double)*nbits)-1 )/den;

  fprintf(filerams[*unit],"%8ld%8d%20.10e%20.10e                        "
	  ,*n,*nbits,bias,fact);


  for( i = 0; i < *n;  i++ )
    scr[i] = ( a[i] + bias ) * fact;

  nchs = *nbits / BITFIELDLENGTH;

  for( j = 0, char_count = 0 ; j < *n; j++ ) {
    ftemp=scr[j];
    for( i = nchs-1; i >= 0; i--, char_count++ )
      cscr[char_count]= 
          vc[( (unsigned)ftemp >> (i * BITFIELDLENGTH) ) & VFORMMASK];
  }

  fwrite(cscr , sizeof( char ), *n * nchs, filerams[*unit] ); 
  *irec=ftell(filerams[*unit]); 
}
