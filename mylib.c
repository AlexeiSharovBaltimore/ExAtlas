/*************************************
* Programmer: Alexei Sharov   11/30/2000
* Department of Entomlogy, Virginia Tech
* (540)231-7316, e-mail sharov@vt.edu
* Home page: www.ento.vt.edu/~sharov
*
* The 'mylib.c' program is included into almost all other
* c-rpograms.
*
* This program is a part of the STS decision-support system. It
* has auxillary functions
***************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "mylib.h"

/*****************************************/
float             slip
/*****************************************
DESCRIPTION
   returns linear interpolation function, defined by a set of points with
   x-projection and y-projection
*/
(float x,             /* argument of a function */
float xy[][X_AND_Y],  /* coordinates of points */
int n_points)         /* number of points */
{
float y;           /* the interpolation */
float z;
int point;         /* number of the nearest point greater than argument */

   if (x <= xy [0][X_PROJECTION])
      y = xy [0][Y_PROJECTION];
   else if (x >= xy [n_points - 1][X_PROJECTION])
      y = xy [n_points - 1][Y_PROJECTION];
   else {   /*   F I N D I N G   'point': */
      for (point = 1; x > xy [point][X_PROJECTION]; ++point);
      z = xy [point][X_PROJECTION] - xy [point-1][X_PROJECTION];
      if (z <= 0.0)
	 printf("ERROR! in slip\n");
      else        /*   L I N E A R   INTERPOLATION: */
	 y = xy[point-1][Y_PROJECTION] + (x - xy[point-1][X_PROJECTION]) *
	   (xy[point][Y_PROJECTION]-xy[point-1][Y_PROJECTION]) / z;
      }
   return (y);
}
/*****************************************/
void        skip_comments      (FILE *fp)
/*****************************************
* DESCRIPTION: skips comments in the data file. Comments start with a single
* quotation mark ' and end either with an end-of-line sign or with a single
* quotation mark. */
{
int index;
int repeat = TRUE;

   while (repeat == TRUE){

      /* skip spaces and end-of-lines before comments start */
      index = ' ';
      while (index == ' ' || index == '\n')
	 index = fgetc (fp);

      /* if it is not a single quotation mark (not comments), return
       * the sign back to the file: */
      if (index != '\''){
         ungetc (index, fp);
	 repeat = FALSE;
	 }

      /* if comments are present, then find the end of comments: either
       * a single quotation mark, or the end-of-line */
      else{
	 index = ' ';
         while (index != '\'' && index != '\n')
	    index = fgetc (fp);
      }
   }
}

/**********************************************/
void    read_array_fl  (FILE *fp, float *x,int n)
/**********************************************
*DESCRIPTION: skips comments and then reads
*   n numbers (float format) to array x[].*/
{
float *px;

   for (px = x; px < x + n; ++px)
      *px = read_float (fp);
}

/**********************************************/
void   read_array_int  (FILE *fp, int *x, int n)
/**********************************************
DESCRIPTION: skips comments and then reads
   n numbers (integer format) to array x[].
*/
{
int *px;

	for (px = x; px < x + n; ++px)
		*px = read_int (fp);
}

/****************************/
void print_line (FILE *fp)
/***************************/
{
	char buffer[21];
	fgets(buffer, 20, fp);
	buffer[20] = '\0';
	printf("%s\n", buffer);
}

/**********************************************/
float      read_float       (FILE *fp)
/**********************************************
DESCRIPTION: skips comments and then reads
	a float number to the variable x.
*/
{
float x;

	skip_comments (fp);
	if (fscanf (fp,"%f",&x) != 1)
	{
		print_line(fp);
		error_message ("reading float");
	}
	return (x);
}
/**********************************************/
int            read_int             (FILE *fp)
/**********************************************
DESCRIPTION: skips comments and then reads
	an integer number to the variable x.
*/
{
int x, index;
char buffer[21];

	skip_comments (fp);
	if (fscanf (fp,"%d",&x) != 1)
	{
		print_line(fp);
		error_message ("reading int");
	}
	index = fgetc (fp);
	if (index == '.')   /* Float number is read instead of integer */
	{
		print_line(fp);
		error_message ("int is expected, not float");
	}
	else
		ungetc (index, fp);
	return (x);
}

/**********************************************/
void            clean_float      (float *x,int n)
/**********************************************
DESCRIPTION: cleans the array of floats; x - the pointer to array
* n - number of elements
*/
{
float *px;
   for (px=x; px < x+n; ++px)
      *px=0.0;
}
/**********************************************/
void            clean_alloc   (float **x, int repl, int n)
/**********************************************
DESCRIPTION: cleans allocated memory (float); x - the pointer to array
 of pointers to allocated memory, repl - size of array x, n - number of
 float elements allocated to each element of array x.
*/
{
int i,k;
   for (i=0; i<repl; ++i)
      for (k=0; k<n; ++k)
	 x[i][k] = 0.0;
}
/**********************************************/
void            clean_int         (int *x,int n)
/**********************************************
DESCRIPTION: cleans the array of integers; x - the pointer to array
* n - number of elements
*/
{
int *px;
   for (px=x; px < x+n; ++px)
      *px=0;
}

/**********************************************/
float         normal_distribution     (float x)
/**********************************************
DESCRIPTION:
   returns an integral of normal probability distribution
   (y changes from 0 to 1; y = 0.5 for x = 0)
   Function translated from FORTRAN code (Numerical Recipes, Press W.H.
   et al. 1986)
*/
{
   static double
      a1=-1.26551223,
      a2= 1.00002368,
		a3= 0.37409196,
      a4= 0.09678418,
      a5=-0.18628806,
      a6= 0.27886807,
      a7=-1.13520398,
      a8= 1.48851587,
      a9=-0.82215223,
      a10=0.17087277;
   double z, t;
   float y;

	z = fabs((double)x)/sqrt(2.);
   t = 1.0 / (1.0 + 0.5 * z);
   y = t*exp(-z * z + a1 + t * (a2 + t * (a3 + t * (a4 + t * (a5 + t *
     (a6 + t * (a7 + t * (a8 + t * (a9 + t * a10)))))))));
   if(x < 0.0) y = 2.0 - y;
   y = 1.0 - 0.5 * y;
   return(y);
}

/*****************************************/
float     divide       (float x, float y)
/*****************************************
DESCRIPTION: safe division of x by y. If y=0 then 0 is returned.
*/
{
float z;
   if (y != 0.0) z = x / y;
   else z = 0.0;
   return (z);
}

/********************************/
float    absfn    (float x)
/********************************/
{
float y;
   y = (float)fabs((double)x);
   return(y);
}

/********************************/
float    sqrtfn    (float x)
/********************************/
{
float y;
	if (x <= 0 && x > -0.0001)
		return(0.0);
	else if (x < 0)
	{
		error_message("sqrt < 0");
		return (0.0);
	}
	y = (float)sqrt((double)x);
	return(y);
}

/********************************/
float    expfn    (float x)
/********************************/
{
float y;
	if (x > FMAXEXP)
		x = FMAXEXP;
	else if (x < -FMAXEXP)
		return(0.0);
	y = (float)exp((double)x);
	return(y);
}

/********************************/
float    powfn    (float x, float p)
/********************************/
{
float t;

	t = p * log((double)x);
	return(expfn(t));
}

/********************************/
float    logfn    (float x)
/********************************/
{
float y;
	if (x > 0.0)
		y = (float)log((double)x);
	else
		error_message ("log < 0");
	return(y);
}

/********************************************/
void   find_string    (char *string, FILE *fp)
/********************************************/
{
int found = FALSE;
int c, l, n=0;

	l = strlen(string);
	while (found == FALSE){
		c = fgetc(fp);
	if (c == EOF){
		fclose(fp);
		error_message ("String not found");
	}
	if (c == *(string+n))
	 ++n;
      else
	 n = 0;
      if (n == l)
	 found = TRUE;
   }
}

/*****************************************/
char  *stringset   (char *x, int ch, int n)
/*****************************************
DESCRIPTION: creates a string of n characters 'ch'
*/
{
int i;
   if (n>0){
      for (i=0; i < n; ++i)
	 *(x+i) = ch;
      *(x+n) = '\0';
		}
   else
      x[0] = '\0';
   return (x);
}
/****************************************************/
void  get_mem(float **ptr, int repl, int elem)
/****************************************************/
{
int i;

	for(i=0; i<repl; ++i){
      check(ptr[i] = (float *)calloc(elem, sizeof(float)));
   }
}
/****************************************************/
void  free_mem(float **ptr, int repl)
/****************************************************/
{
int i;

   for(i=0; i<repl; ++i)
		free(ptr[i]);
}
/****************************************************/
void  get_mem_int(int **ptr, int repl, int elem)
/****************************************************/
{
int i;

   for(i=0; i<repl; ++i)
      check(ptr[i] = (int*)calloc(elem, sizeof(int)));
}
/****************************/
void   check     (void *ptr)
/****************************/
{
   if (ptr == NULL)
		error_message ("Out of memory");
}

/**********************************/
char    *truncate   (char buffer[])
/**********************************/
{
char *buf;

	buf = buffer + strlen(buffer) - 1;
	while ((*buf == ' ' || *buf == '\n') && buf > buffer)
		--buf;
	*(buf+1) = '\0';
	return (buffer);
}

/***************************************/
float  rand_norm  (float xm, float std)
/***************************************
DESCRIPTION: Generates normal distribution with mean xm and standard
deviation std.
*/
{
float  sum = 0.;
int i;

   for (i=0; i<12; ++i)
		sum += (float)rand()/RAND_MAX;
   return ((sum - 6.)*std + xm);
}

/***************************************/
float  rand_lognorm  (float xm, float cv)
/***************************************
DESCRIPTION: Generates log-normal distribution with mean xm and
coefficient of variation cv.
*/
{
float stdev_log,
   aver_log;

   stdev_log = sqrtfn (logfn (cv * cv + 1.0));
   aver_log = logfn(xm) - stdev_log * stdev_log / 2;
   return (expfn(rand_norm(aver_log, stdev_log)));
}

/********************************/
void   error_message   (char *message)
/********************************/
{
	printf("ERR: %s.\n", message);
	exit(1);
}

/*************************************/
char   *read_string     (char *string, int len, FILE *fp)
/*************************************/
{
	skip_comments(fp);
	return(truncate(fgets(string,len,fp)));
}
