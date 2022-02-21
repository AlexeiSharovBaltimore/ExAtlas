/*************************************
* Programmer: Alexei Sharov   11/30/2000
* Department of Entomlogy, Virginia Tech
* (540)231-7316, e-mail sharov@vt.edu
* Home page: www.ento.vt.edu/~sharov
*
* Header file for my library MYLIB.C
*
* This program is a part of the STS decision-support system.
***************************************/
#ifndef __MYLIB_H
#define __MYLIB_H
#define MAXINT	0x7FFF
#define FMAXEXP	38

#include <stdio.h>

/* ASCII control characters */
#define NO_INPUT       0
#define ASCNUL      (256+3)
#define ASCBEL         7
#define ASCBS          8
#define TAB            9
#define ASCLF        0xA
#define ASCFF        0xC
#define ENTER        0xD
#define ESC         0x1B
#define DEL         0x7F
#define SPACE       0x20

   /* special keys for IBM PC */
#define HOMEKEY     (256+71)
#define BACKTAB     (256+15)
#define UPARROW     (256+72)
#define LEFTARROW   (256+75)
#define RIGHTARROW  (256+77)
#define ENDKEY      (256+79)
#define DOWNARROW   (256+80)
#define PGUPKEY     (256+73)
#define PGDNKEY     (256+81)
#define INSERTKEY   (256+82)
#define DELETEKEY   (256+83)
#define CTLPRTSC    (256+114)
#define CTLLARROW   (256+115)
#define CTLRARROW   (256+116)
#define CTLEND      (256+117)
#define CTLPGDN     (256+118)
#define CTLHOME     (256+119)
#define CTLPGUP     (256+132)

  /* function key codes */
#define F1       (256+59)
#define F2       (256+60)
#define F3       (256+61)
#define F4       (256+62)
#define F5       (256+63)
#define F6       (256+64)
#define F7       (256+65)
#define F8       (256+66)
#define F9       (256+67)
#define F10      (256+68)

  /* alt-key + number key (top row) */
#define ALT1KEY     (256+120)
#define ALTSKEY     (256+31)

#define SQUARE(x)     ((x)*(x))
#define PI            3.1416
#define X_AND_Y       2
#define X_PROJECTION  0
#define Y_PROJECTION  1
#define TOO_BIG       200.0
#define FALSE         0
#define TRUE          1

/*   F U N C T I O N   P R O T O T Y P E S */

float slip               (float x,float xy[][X_AND_Y],int n);
void  read_array_fl      (FILE *fp,float *x,int n);
void  read_array_int     (FILE *fp,int *x,int n);
float read_float         (FILE *fp);
int   read_int           (FILE *fp);
void  clean_float        (float x[],int n);
void  clean_alloc   (float **x, int repl, int n);
void  clean_int          (int x[],int n);
float  normal_distribution  (float x);
float  divide (float x, float y);
float  expfn (float x);
float  powfn (float x, float p);
float  logfn (float x);
float  absfn (float x);
float  sqrtfn (float x);
void   skip_comments (FILE *fp);
void   clear_screen   (int textcol, int backgr);
void   find_string    (char *string, FILE *fp);
char  *stringset   (char *x, int ch, int n);
void  get_mem(float **ptr, int repl, int elem);
void  free_mem(float **ptr, int repl);
void  get_mem_int(int **ptr, int repl, int elem);
void  check      (void *ptr);
char  *truncate  (char buffer[]);
float rand_norm  (float xm, float std);
float rand_lognorm  (float xm, float cv);
void  remove_spaces (char *string);
void  beep   (int freq, int duration);
int   getkey (void);
void  initialize_mouse  (void);
void  mouse_off  (void);
void  mouse_on  (void);
int   get_mouse_click  (int *x, int *y);
void  error_message  (char *message);
float powint  (float x, int n);
char *truncate (char *buffer);
char *read_string  (char *string, int len, FILE *fp);

#endif
