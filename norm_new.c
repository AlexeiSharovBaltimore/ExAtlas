/********************************************
* Local normalization procedure. The distribution is
* adjusted by a piece-linear method so that quantiles
* coincide with quantiles of a given standard distribution
* It is possible to generate a standard distribution by
* averaging the quantiles from all available data.
* 
* The method was modified from the original method
* of Igor Sidorov et al. (NCI-NIH, Friderick)
*
* Author: Alexei A. Sharov 
*********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXCOL  5000
#define MISSING -9999

typedef struct data_st{
	long ncol;
	long nrow;
	float *data[MAXCOL];
	char **header_row;
	long missing[MAXCOL];
	char *comments_before[500];
	char *comments_after[500];
	int n_before;
	int n_after;
}DATA;

typedef struct param_st{
	long ncol;
	long nrow;
	long nquant;
	long ncol_target;
}PARAM;
float median1(float *a, long n);

/***********************************************/
void error_message (char *message)
/***********************************************/
{
printf("ERROR: %s\n", message);
exit(1);
}

/***********************************************/
void check (void *x)
/***********************************************/
{
if (!x){ error_message("ERROR: Out of memory\n"); }
return;
}

/*************************************************************************/
void sortem   (long ie, float *a, long *b)
/*************************************************************************/
{
long i, j, k, m, p, q, iring;
long lt[64], ut[64];
float ta, xa;
long tb, xb;
long i__1;

/* Function Body */
j = ie-1;
m = 1;
i = 0;
iring = 2;

/* If this segment has more than two elements  we split it */
L10:	 if ((i__1 = j - i - 1) < 0) {
	goto L100;
	 } else if (i__1 == 0) {
	goto L90;
	 } else {
	goto L15;
	 }

/* p is the position of an arbitrary element in the segment we choose the 
* middle element. Under certain circumstances it may be advantageous 
* to choose p at random. */

L15:
	p = (j + i) / 2;
	ta = a[p];
	a[p] = a[i];
L19:	if(b){
		tb = b[p];
		b[p] = b[i];
	}
L21: /* Start at the beginning of the segment, search for k such that a(k)>t */
    q = j;
    k = i;
L20:	++k;
	if (k > q) goto L60;
    if (a[k] <= ta) goto L20;
/* Such an element has now been found now search for a q such that a(q)<t 
* starting at the end of the segment. */
L30:  if (a[q] < ta) {
	goto L40;
	}
    --q;
    if (q > k) {
	goto L30;
    }
    goto L50;

/* a(q) has now been found. we interchange a(q) and a(k) */

L40: xa = a[k];
    a[k] = a[q];
    a[q] = xa;
L44:	if(b){
		xb = b[k];
		b[k] = b[q];
		b[q] = xb;
	}
L45:
/* Update q and search for another pair to interchange: */
    --q;
   goto L20;
L50:    q = k - 1;
L60:
/* The upwards search has now met the downwards search: */
    a[i] = a[q];
    a[q] = ta;
L64:    if(b){
		b[i] = b[q];
		b[q] = tb;
	}
L65:

/* The segment is now divided in three parts: (i,q-1),(q),(q+1,j) */
/* store the position of the largest segment in lt and ut */
    if (q << 1 <= i + j) {
	goto L70;
 }
 lt[m - 1] = i;
 ut[m - 1] = q - 1;
 i = q + 1;
 goto L80;
L70:	lt[m - 1] = q + 1;
	ut[m - 1] = j;
	j = q - 1;
/* Update m and split the new smaller segment */
L80:	++m;
	goto L10;

/* We arrive here if the segment has  two elements we test to see if */
/* the segment is properly ordered if not, we perform an interchange */
L90:
	if (a[i] <= a[j]) goto L100;
	xa = a[i];
	a[i] = a[j];
	a[j] = xa;
L94:    if(b){
		xb = b[i];
		b[i] = b[j];
		b[j] = xb;
	}
L95:

/* If lt and ut contain more segments to be sorted repeat process: */
L100:	 --m;
	if (m <= 0) goto L110;
	i = lt[m - 1];
	j = ut[m - 1];
	goto L10;
L110:	return;
} /* sortem_ */

/**************************/
float  median1  (float *a, long n)
/**************************/
{
float *b, median=MISSING,y;
long i, i1=0;

if(n<=0) return(median);
check(b = (float*)malloc(n*sizeof(float)));
memcpy(b,a,n*sizeof(float));
sortem(n,b,NULL);
while(i1<n && b[i1]<=MISSING) i1++;
n -= i1;
i = n/2;
if(n){
	median = b[i1+i];
	if(i*2==n){
		y = b[i1+i-1];
		median = (median+y)/2;
	}
}
free(b);
return(median);
}

/*****************************************/
long  count_missing (float *data, long n)
/*****************************************/
{
long i=0;
while(i<n && data[i]<=MISSING) i++;
return(i);
}

/*****************************************/
void  smooth_sorted (float *data, long n)
/*****************************************/
{
long i,i0,i1,d,nmax=50;
float b, x0, x1, x;

if(n<nmax) return;
i0 = 20;
i1 = nmax-1;
d = i1-i0;
x0 = data[i0];
x1 = data[i1];
if(x0<1.0e-6) x0=1.0e-6;
if(x1<1.0e-6) x1=1.0e-6;
x0 = log(x0);
x1 = log(x1);
b = (x1-x0)/d;
for(i=0; i<20; i++){
	x = exp(x0 + b*(i-i0));
	data[i] = x;
}
i0 = n-nmax;
i1 = n-20;
d = i1-i0;
x0 = data[i0];
x1 = data[i1];
if(x0<1.0e-6) x0=1.0e-6;
if(x1<1.0e-6) x1=1.0e-6;
x0 = log(x0);
x1 = log(x1);
b = (x1-x0)/d;
for(i=n-20; i<n; i++){
	x = exp(x1 + b*(i-i1));
	data[i] = x;
}
return;
}

/*****************************************/
void  normalize (DATA *d, PARAM *p)
/*****************************************/
{
float **data;
long **index, *numbers, *missingData;
float x, dq1, dq2;
float a, b;
long irow, icol, i, j;

check(data = (float**)malloc((p->ncol+1)*sizeof(float*)));
check(index = (long**)malloc((p->ncol+1)*sizeof(long*)));
check(missingData = (long*)malloc(p->ncol*sizeof(long)));
check(numbers = (long*)malloc(p->nrow*sizeof(long)));
check(data[d->ncol] = (float*)calloc(p->nrow,sizeof(float)));
for(irow=0; irow<d->nrow; irow++) numbers[irow]=irow;
for(icol=0; icol<d->ncol; ++icol){
	check(data[icol] = (float*)malloc(p->nrow*sizeof(float)));
	check(index[icol] = (long*)malloc(p->nrow*sizeof(long)));
	memcpy(data[icol],d->data[icol],d->nrow*sizeof(float));
	memcpy(index[icol],numbers,d->nrow*sizeof(long));
	sortem(d->nrow,data[icol],index[icol]);
	missingData[icol] = count_missing(data[icol],d->nrow);
	if(icol < p->ncol_target){
		for(irow=0; irow<d->nrow; ++irow){
			i = missingData[icol]+(irow+0.5)/d->nrow*(d->nrow-missingData[icol]);
			data[d->ncol][irow] += data[icol][i];
		}
	}
}
for(irow=0; irow<d->nrow; ++irow)
	data[d->ncol][irow] /= p->ncol_target;
smooth_sorted(data[d->ncol],d->nrow);
for(icol=0; icol<d->ncol; ++icol){
	for(irow=0; irow<d->nrow; ++irow){
		i = index[icol][irow];
		if(irow>=missingData[icol]){
			j = (irow-missingData[icol]+0.5)*d->nrow/(d->nrow-missingData[icol]);
			d->data[icol][i] = data[d->ncol][j];
		}
	}
}
for(icol=0; icol<=d->ncol; ++icol) free(data[icol]);
for(icol=0; icol <d->ncol; ++icol) free(index[icol]);
free(data);
free(missingData);
free(index);
free(numbers);
return;
}

/***********************************************/
void  read_to_the_end_of_line (FILE *fp)
/***********************************************/
{
char a;
while((a=fgetc(fp))!=EOF && a != '\n');
}

/***********************************************/
void  read_tabdelimited_item (FILE *fp, char *buffer)
/***********************************************/
{
char ch;
long i=0, j=0;
while((ch=fgetc(fp)) && ch != '\t' && ch != '\n'){
	if(j==0 && (ch=='\t' || ch=='\n')){
		j++;
		continue;
	}
	j++;
	if(i==0 && ch==' ') continue;
	if(i < 60-1) buffer[i++] = ch;	
}
buffer[i] = '\0';
if(ch == '\n' || ch != '\t') ungetc (ch, fp);
}

/***********************************************/
long  split_string (char *string, char *items[], long num)
/***********************************************/
{
char *ch;
long i=0;

ch = strchr(string,'\n');
if(ch) *ch = '\0';
ch = string;
while(1){
	items[i] = ch;
	ch = strchr(ch,'\t');
	if(ch) *ch = '\0';
	else break;
	if(i>=num-1) break;
	ch++;
	i++;
}
return(i+1);
}

/***********************************************/
char  *copy_string (char *x)
/***********************************************/
{
long len;
char *y;

len = strlen(x);
check(y = (char*)malloc((len+1)*sizeof(char)));
strcpy(y, x);
return(y);
}

/***********************************************/
DATA *read_data (char *filename, char *outputfile, PARAM *p) 
/***********************************************/
{
FILE *inputFile, *fp;
DATA *d;
long i, j, int_in_col1, nitems, count=0, norm=0, comments=0;
long buff_size = 100000;
char *buffer, **items;
float x, *data_row, median, sum, nn;

if(p->ncol >= MAXCOL-1) error_message("Increase MAXCOL!");
check(buffer = (char*)malloc(sizeof(char)*buff_size));
check(items = (char**)malloc(50000*sizeof(char*)));
check(d = (DATA*)malloc(sizeof(DATA)));
d->nrow = p->nrow;
d->ncol = p->ncol;
check(d->header_row = (char**)malloc(d->nrow*sizeof(char*)));
for(i=0; i<d->ncol; ++i){
	check(d->data[i] = (float*)malloc(d->nrow*sizeof(float)));
}
inputFile = fopen(filename, "r");
if(!inputFile) error_message("Input file not found\n");
while(fgets(buffer,buff_size-1,inputFile)){
	if(buffer[0] == '!') comments=1;
	if(!strncmp(buffer,"!Series_normalized",18)){
		if(norm==0){
			norm=1;
			strcpy(buffer,"!Series_normalized\t\"true\"\n");
		}else{
			continue;
		}
	}else if(norm==0 && comments && (strncmp(buffer,"!Series",7) || !strncmp(buffer,"!series_matrix_table_begin",26))){
		norm=1;
		d->comments_before[count++] = copy_string("!Series_normalized\t\"true\"\n");
	}
	d->comments_before[count++] = copy_string(buffer);
	if(count>490) error_message("Too many comment lines");
	if(buffer[0] != '!' && buffer[0] != '\t' && strlen(buffer)>2) break;
}
d->n_before = count;
for(i=0;i<d->nrow;++i){
	fgets(buffer,buff_size-1,inputFile);
	if(buffer[0] == '!') error_message("Unexpected comments\n");
	nitems = split_string (buffer,items,50000);
	d->header_row[i] = copy_string(items[0]);
	sum=0;
	nn=0;
	for(j=0; j<d->ncol; ++j){
		if(sscanf(items[j+1],"%f", &x) != 1){
			printf("Error reading: row %d, col %d\n", i+2, j);
			d->data[j][i] = MISSING;
		}else{
			if(x<=1.0e-8) x = MISSING;
			d->data[j][i] = x;
		}
	}
}
count=0;
while(fgets(buffer,buff_size-1,inputFile)){
	d->comments_after[count++] = copy_string(buffer);
	if(count>490) error_message("Too many comment lines after");
}
d->n_after = count;
fclose(inputFile);
free(buffer);
free(items);
return(d);
}

/***********************************************/
void print_output  (DATA *d, char *filename, PARAM *p)
/***********************************************/
{
FILE *outputFile;
long i,j;

outputFile= fopen(filename, "w");
if(!outputFile) error_message("ERROR: Cannot open ouput file");
for(j=0;j<d->n_before;++j){
	fprintf(outputFile,"%s",d->comments_before[j]);
}
for(j=0;j<d->nrow;++j){
	fprintf(outputFile,"%s",d->header_row[j]);
	for(i=0;i<d->ncol;++i){
		if(d->data[i][j]<=MISSING)
			fprintf(outputFile,"\t%.0f", d->data[i][j]);
		else if(d->data[i][j] > 1000)
			fprintf(outputFile,"\t%.1f", d->data[i][j]);
		else if(d->data[i][j] > 10)
			fprintf(outputFile,"\t%.3f", d->data[i][j]);
		else if(d->data[i][j] > 0.1)
			fprintf(outputFile,"\t%.5f", d->data[i][j]);
		else
			fprintf(outputFile,"\t%.7f", d->data[i][j]);
	}
	fprintf(outputFile,"\n");
}
for(j=0;j<d->n_after;++j){
	fprintf(outputFile,"%s",d->comments_after[j]);
}
fclose(outputFile);
}

/***********************************************/
long	count_fields(char *buffer, long len)
/***********************************************/
{
long i, n=0;
for(i=0; i<len && buffer[i] != '\n'; ++i){
	if(buffer[i] == '\t' && buffer[i+1] !='\n')
		++n;
}
return(n);
}

/***********************************************/
void  get_dimensions  (char *filename, long *nrow, long *ncol)
/***********************************************/
{
char *buffer,*ch;
FILE *fp;

fp = fopen(filename, "r");
if(!fp) error_message("Input file not found");
check(buffer = (char*)malloc(sizeof(char)*100000));
do{
	fgets(buffer,100000,fp);
	ch = buffer;
	while((*ch==' ' || *ch=='\t') && ch<buffer+10000) ch++;
}while(*ch=='!' || strlen(ch) < 2);
*ncol = count_fields(ch,100000);
*nrow = 0;
while(fgets(buffer,100000,fp) && strlen(buffer)>2){
	if(buffer[0]=='!') break;
	++(*nrow);
}
free(buffer);
fclose(fp);
}

/***********************************************/
int main (int argc, char **argv)
/***********************************************/
{
DATA *d;
PARAM *p;
long j;

if (argc < 2){
	error_message("norm_new inputFile outputFile [-n ncol_target]\n");
}
check(p = (PARAM*)calloc(1,sizeof(PARAM)));
get_dimensions(argv[1], &p->nrow, &p->ncol);
if(p->nrow<10 || p->ncol<1) error_message("Wrong table dimensions");
if(argc > 3 && !strcmp(argv[3],"-n")){
	p->ncol_target = atoi(argv[4]);
	if(p->ncol_target <= 0 || p->ncol_target > p->ncol)
		p->ncol_target=p->ncol;
}else{
	p->ncol_target = p->ncol;
}
d = read_data(argv[1], argv[2], p);
normalize(d, p);
print_output(d, argv[2], p);
return(0);
}

