/*************************************
* anova_oneway.c is a part of ExAtlas software
* http://lgsun.grc.nia.nih.gov/ANOVA
*
* Function: Analysis of variance (ANOVA) of microarray data; average error variance
* for genes with same intensity is used to adjust ANOVA for each gene (5 error models)
* FDR is estimated. This is a stand-alone version which accepts GEO-formatted matrix
* as input. It is different from the original version as follows:
* (1) Only one column per array is allowed. Thus, 2-dye arrays should be normalized
* by reference first, or considered as two arrays.
* (2) Parameters are entered at command line
* (3) Column headers are retrieved from GEO-formatted matrix
* (4) Permutation tests are not available
* 
* Syntax:
* anova_oneway -i inputFasta -o outputFile [-a annotationFile, -err errorModel,
* -z ZvalueOutliers, -FDR maxFDR]
* 
* Author: Alexei Sharov   10/25/2013
* Lab. of Genetics, National Institute on Aging (NIA/NIH)
* Phone: 410-558-8556   Email: sharoval@mail.nih.gov
* 
* If you use this program in your research, please cite:
* Sharov AA, Dudekula DB, Ko MS. 2005. A web-based tool for principal component 
* and significance analysis of microarray data. Bioinformatics. 21(10):2548-9. 
**************************************/

/*
anova_oneway -i matrix_NIAadult_cells.txt -o output.txt
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#define MAXIT 10000
#define EPS 3.0e-7
#define FPMIN 1.0e-30
#define M_SQRT2PI 2.50662827463100050242

#define MAXCOL  5000
#define MISSING -9999
#define STRINGLEN 50

#define ERR_ACTUAL		1
#define ERR_AVER		2
#define ERR_BAYESIAN		3
#define ERR_MAX_AVER_ACTUAL	4
#define ERR_MAX_AVER_BAYESIAN	5

#define DESIGN_ONEDYE	0
#define DESIGN_TWODYE	1

#define ASCENDING	0
#define DESCENDING	1

typedef struct data_st{
	long ncol;
	long nrow;
	float *data[MAXCOL];
	float *mean_ref;
	float *median_ref;
	long n_levels;		/* number of factor levels */
	long *n_reps;		/* number of replications for each level */
	long totrep;		/* total replications */
	char *header_col[MAXCOL];
	char *level_name[MAXCOL];
	char **header_row;
	char *header1;
}DATA;

typedef struct param_st{
	long window;		/* averaging window for smoothing error variance */
	float propErrVar;	/* proportion of highest error variances to be removed */
	long errorModel;		/* Error type */
	float thresh_z;		/* Threshold z-value for outlier detection */
	long bayesianDF;		/* Bayesian degrees of freedom */
	float FDR;		/* FDR threshold */
	float cutoff;		/* Threshold z-value for outlier detection */
	long use_probeID;
	char *inputFile;
	char *outputFile;
	char *annotationFile;
}PARAM;

typedef struct fstat_st{
	long nrow;
	float *F;
	float *var;
	float *errVar;
	float *smoothErrVar;
	float *errVarFinal;
	long *df1;
	long *df2;
	float *means[MAXCOL];
	long *n[MAXCOL];
	float *glob_mean;
	float *mean_ref;
	float *averIntensity;
	long *rank;
}FSTAT;

typedef struct annot_st{
	char *id;
	char *symbol;
	char *geneName;
	char *add1;
	char *add2;
	char *probeName;
}ANNOTATION;

typedef struct hash_st{
	long n;
	long n_alloc;
	long type;
	long increment;
	char **list;
	long *value1;
	float *value2;
	char **value3;
}HASH;

char  *copy_string (char *x);
double normal_inverse(double p);
double normal_distribution(double x);

/*************************************************************************/
void sortem   (long ie, float *a, long iperm, float *b,
		float *c, float *d, float *e, float *f, float *g, float *h)
/*************************************************************************/
{
	 long i, j, k, m, p, q, iring;
	 long lt[64], ut[64];
	 float ta, tb, tc, td, te, tf, tg, th, xa, xf, xg;
	 float xh, xe, xd, xc, xb;
	 long i__1;

	 /* Function Body */
	 j = ie-1;
	 m = 1;
	 i = 0;
	 iring = iperm + 1;
	 if (iperm > 7) {
		 iring = 1;
	 }

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
	 switch (iring) {
	case 1:  goto L21;
	case 2:  goto L19;
	case 3:  goto L18;
	case 4:  goto L17;
	case 5:  goto L16;
	case 6:  goto L161;
	case 7:  goto L162;
	case 8:  goto L163;
	 }
L163:	 th = h[p];
	 h[p] = h[i];
L162:	 tg = g[p];
    g[p] = g[i];
L161:	 tf = f[p];
    f[p] = f[i];
L16:	 te = e[p];
	 e[p] = e[i];
L17:	 td = d[p];
    d[p] = d[i];
L18:	 tc = c[p];
    c[p] = c[i];
L19:	 tb = b[p];
	 b[p] = b[i];
L21: /* Start at the beginning of the segment, search for k such that a(k)>t */
    q = j;
    k = i;
L20:	 ++k;
	 if (k > q) {
	goto L60;
	 }
    if (a[k] <= ta) {
	goto L20;
    }
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
    switch (iring) {
	case 1:  goto L45;
	case 2:  goto L44;
	case 3:  goto L43;
	case 4:  goto L42;
	case 5:  goto L41;
	case 6:  goto L411;
	case 7:  goto L412;
	case 8:  goto L413;
    }
L413:     xh = h[k];
    h[k] = h[q];
	 h[q] = xh;
L412:	 xg = g[k];
    g[k] = g[q];
	 g[q] = xg;
L411:    xf = f[k];
	 f[k] = f[q];
    f[q] = xf;
L41:	 xe = e[k];
	 e[k] = e[q];
    e[q] = xe;
L42:	 xd = d[k];
    d[k] = d[q];
    d[q] = xd;
L43:    xc = c[k];
	 c[k] = c[q];
	 c[q] = xc;
L44:	 xb = b[k];
    b[k] = b[q];
	 b[q] = xb;
L45:
/* Update q and search for another pair to interchange: */
    --q;
   goto L20;
L50:    q = k - 1;
L60:
/* The upwards search has now met the downwards search: */
    a[i] = a[q];
    a[q] = ta;
   switch (iring) {
	case 1:  goto L65;
	case 2:  goto L64;
	case 3:  goto L63;
	case 4:  goto L62;
	case 5:  goto L61;
	case 6:  goto L611;
	case 7:  goto L612;
	case 8:  goto L613;
    }
L613:	 h[i] = h[q];
	 h[q] = th;
L612:    g[i] = g[q];
	 g[q] = tg;
L611:    f[i] = f[q];
    f[q] = tf;
L61:    e[i] = e[q];
	 e[q] = te;
L62:	 d[i] = d[q];
    d[q] = td;
L63:    c[i] = c[q];
    c[q] = tc;
L64:    b[i] = b[q];
    b[q] = tb;
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
L70:	 lt[m - 1] = q + 1;
	 ut[m - 1] = j;
	 j = q - 1;
/* Update m and split the new smaller segment */
L80:	 ++m;
	 goto L10;

/* We arrive here if the segment has  two elements we test to see if */
/* the segment is properly ordered if not, we perform an interchange */
L90:
    if (a[i] <= a[j]) {
	goto L100;
    }
	 xa = a[i];
	 a[i] = a[j];
	 a[j] = xa;
	 switch (iring) {
	case 1:  goto L95;
	case 2:  goto L94;
	case 3:  goto L93;
	case 4:  goto L92;
	case 5:  goto L91;
	case 6:  goto L911;
	case 7:  goto L912;
	case 8:  goto L913;
    }
L913:	 xh = h[i];
	 h[i] = h[j];
	 h[j] = xh;
L912:    xg = g[i];
    g[i] = g[j];
	 g[j] = xg;
L911:	 xf = f[i];
    f[i] = f[j];
	 f[j] = xf;
L91:	 xe = e[i];
    e[i] = e[j];
    e[j] = xe;
L92:	 xd = d[i];
	 d[i] = d[j];
    d[j] = xd;
L93:	 xc = c[i];
	 c[i] = c[j];
	 c[j] = xc;
L94:    xb = b[i];
    b[i] = b[j];
	 b[j] = xb;
L95:

/* If lt and ut contain more segments to be sorted repeat process: */
L100:	 --m;
	 if (m <= 0) {
	goto L110;
	 }
	 i = lt[m - 1];
	 j = ut[m - 1];
	 goto L10;
L110:	 return;
} /* sortem_ */

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
if(!x){ error_message("Out of memory\n"); }
return;
}

/***********************************************/
long hash_word_pos (char *input, HASH *h) 
/***********************************************/
{
long n, j, i, cmp;
char *value;

if(!input) error_message("null in hash_word_pos");
if(!h) error_message("null hash in hash_word_pos");
n = h->n;
if(!n) return(-1);
j = -1;
i = n/2;
while(1){
	value = h->list[i];
	cmp = strcmp(value,input);
	if(cmp==0){
		break;
	}else if(cmp<0){
		if(n-i<=1){ break; }
		j = i;
		i = (n+i)/2;
	}else{
		if(i-j<=1){ --i; break; }
		n = i;
		i = (j+n)/2;
	}
}
return(i);
}

/***********************************************/
HASH *new_hash (long type, long increment) 
/***********************************************/
{
HASH *h;
check(h = (HASH*)calloc(1,sizeof(HASH)));
h->n_alloc = increment;
h->type = type;
h->increment = increment;
check(h->list = (char**)calloc(increment,sizeof(char*)));
if(type==1) check(h->value1 = (long*)calloc(increment,sizeof(long)));
else if(type==2) check(h->value2 = (float*)calloc(increment,sizeof(float)));
else if(type==3) check(h->value3 = (char**)calloc(increment,sizeof(char*)));
else if(type!=0) error_message("wrong hash type");
return(h);
}

/***********************************************/
void destroy_hash (HASH *h) 
/***********************************************/
{
long i;
if(!h) return;
if(!h->n_alloc){ free(h); return; }
for(i=0; i<h->n; i++){
	free(h->list[i]);
	if(h->type==3) free(h->value3[i]);
}
free(h->list);
if(h->type==1) free(h->value1);
else if(h->type==2) free(h->value2);
else if(h->type==3) free(h->value3);
else if(h->type!=0) error_message("wrong hash type");
free(h);
return;
}

/***********************************************/
void  add_hash_long (char *word, long value, long add, HASH *h)
/***********************************************/
{
long j;
if(!h) error_message("null hash in add_hash_long");
if(h->type>1) error_message("incorrect hash type long");
if(!word) error_message("null in add_hash_long");
if(h->n == h->n_alloc-1){
	h->n_alloc += h->increment;
	check(h->list = (char**)realloc(h->list,h->n_alloc*sizeof(char*)));
	if(h->type==1) check(h->value1 = (long*)realloc(h->value1,h->n_alloc*sizeof(long)));
}
j = hash_word_pos(word, h);
//printf("A %d %d %s\n",h->n,j,word);
if(j>=0 && !strcmp(word,h->list[j])){
	if(h->type==1){
		if(add) h->value1[j] += value;
		else h->value1[j] = value;
	}
	return;
}
if(j<h->n-1) memmove(h->list+j+2,h->list+j+1,(h->n-j-1)*sizeof(char*));
h->list[j+1] = copy_string(word);
if(h->type==1){
	if(j<h->n-1) memmove(h->value1+j+2,h->value1+j+1,(h->n-j-1)*sizeof(long));
	h->value1[j+1] = value;
}
h->n++;
return;
}

/***********************************************/
long hash_value_long (char *input, HASH *h) 
/***********************************************/
{
long pos;
if(!h) error_message("null hash in hash_value_long");
if(h->type>1) error_message("incorrect hash type in hash_value_long");
pos = hash_word_pos (input, h);
if(pos>=0 && !strcmp(input,h->list[pos])){
	if(h->type==1) return(h->value1[pos]);
	return(1);
}
if(h->type==1) return(-1);
return(0);
}

/***********************************************/
void  reverse (long *index, long n)
/***********************************************/
{
long i, swap;
for(i=0; i<n/2; ++i){
	swap = index[i];
	index[i] = index[n-1-i];
	index[n-1-i] = swap;
}
}

/***********************************************/
long  accountNum (char *word)
/***********************************************/
{
long i, dot=0, len, count=0;
len = strlen(word);
for(i=0; i<len; ++i){
	if(isdigit(word[i])) count++;
	if(word[i]=='.' && !dot && i>0) dot=i;
}
if(count<3) dot=0;
return dot;
}

/***********************************************/
FSTAT *init_fstat (DATA *d)
/***********************************************/
{
FSTAT *fs;
long i;

fs = (FSTAT*)malloc(sizeof(FSTAT));
fs->nrow = d->nrow;
fs->F = (float*)calloc(d->nrow,sizeof(float));
fs->var = (float*)calloc(d->nrow,sizeof(float));
fs->df1 = (long*)calloc(d->nrow,sizeof(long));
fs->df2 = (long*)calloc(d->nrow,sizeof(long));
fs->errVarFinal = (float*)calloc(d->nrow,sizeof(float));
fs->errVar      = (float*)calloc(d->nrow,sizeof(float));
fs->smoothErrVar= (float*)calloc(d->nrow,sizeof(float));
for(i=0; i<d->n_levels; ++i){
	fs->means[i] = (float*)calloc(d->nrow,sizeof(float));
	fs->n[i] = (long*)calloc(d->nrow,sizeof(long));
}
fs->glob_mean = (float*)calloc(d->nrow,sizeof(float));
fs->mean_ref = (float*)calloc(d->nrow,sizeof(float));
fs->averIntensity = (float*)calloc(d->nrow,sizeof(float));
fs->rank = (long*)calloc(d->nrow,sizeof(long));

return(fs);
}

/***********************************************/
void	smooth_error_variance(FSTAT *fs, long *index, long window, float propErrVar)
/***********************************************/
{
long ii, i, irow, imax, imin, old_irow=-1, filled=0;
float sum;
long n;
float *errVar, errMax, var, var1, maxErrVar=0;

errVar = (float*)calloc(fs->nrow,sizeof(float));
for(ii=0; ii<fs->nrow; ++ii)
	errVar[ii] = fs->errVar[ii];
sortem(fs->nrow, errVar, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

imax = (fs->nrow-1)*(1.0-propErrVar);
if(imax < 0) imax = 0;
if(imax >= fs->nrow) imax = fs->nrow-1;
errMax = errVar[imax];

sum = 0;
n = 0;
for(ii=0; ii<window/2 && ii<fs->nrow; ++ii){
	var = fs->errVar[index[ii]];
	if (var > MISSING && var <= errMax){
		sum += var;
		n++;
	}
}
for(ii=0; ii<fs->nrow; ++ii){
	irow = index[ii];
	imax = ii + (window+0.5)/2;
	imin = floor(ii - (window-0.5)/2);
	if (imax < fs->nrow){
		var = fs->errVar[index[imax]];
		if (var > MISSING && var <= errMax){
			sum += var;
			n++;
		}
	}
	if (imin >= 0){
		var1 = fs->errVar[index[imin]];
		if (var1 > MISSING && var1 <= errMax){
			sum -= var1;
			n--;
		}
	}
	if(n){
		fs->smoothErrVar[irow] = sum/n;
		if(!filled){
			i = 0;
			while(i<ii && fs->smoothErrVar[index[i]]==0){
				fs->smoothErrVar[index[i++]] = fs->smoothErrVar[irow];
			}
			filled=1;
		}
		old_irow = irow;
	}
	else fs->smoothErrVar[irow] = fs->smoothErrVar[old_irow];
}
for(ii=fs->nrow/2; ii<fs->nrow; ++ii){
	irow = index[ii];
	if(maxErrVar < fs->smoothErrVar[irow]) maxErrVar = fs->smoothErrVar[irow];
	else fs->smoothErrVar[irow] = maxErrVar;
}
free(errVar);
}

/***********************************************/
long	*get_index	(float *x, long n)
/***********************************************/
{
float *index1, *x1;
long *index, i;

index1 = (float*)malloc(n*sizeof(float));
x1 = (float*)malloc(n*sizeof(float));
index = (long*)malloc(n*sizeof(long));

for(i=0; i<n; ++i){
	index1[i]=i;
	x1[i] = x[i];
}
sortem(n, x1, 1, index1, NULL, NULL, NULL, NULL, NULL, NULL);
for(i=0; i<n; ++i)
	index[i]=index1[i];
free(index1);
free(x1);
return(index);
}

/***********************************************/
long	*get_rank	(float *x, long n, long direction, float *index3)
/***********************************************/
{
float *index1, *index2, *x1;
long *index, i;

index1 = (float*)malloc(n*sizeof(float));
index2 = (float*)malloc(n*sizeof(float));
x1 = (float*)malloc(n*sizeof(float));
index = (long*)malloc(n*sizeof(long));

for(i=0; i<n; ++i){
	index1[i]=i;
	index2[i]=i;
	x1[i] = x[i];
}
sortem(n, x1, 1, index1, NULL, NULL, NULL, NULL, NULL, NULL);
for(i=0; i<n; ++i)
	index[i]=index1[i];
if(direction==DESCENDING)
	reverse(index, n);
for(i=0; i<n; ++i){
	index3[i]=index[i];
	index1[i]=index[i];
}
sortem(n, index1, 1, index2, NULL, NULL, NULL, NULL, NULL, NULL);
for(i=0; i<n; ++i)
	index[i]=index2[i];
free(index1);
free(index2);
free(x1);
return(index);
}

/***********************************************/
float median (float *xr, long n)
/***********************************************/
{
float median;
long start=0, j;

if(!n){ return(0); }
sortem(n, xr, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

while(start<n && xr[start]==MISSING){ start++; }
if(start==n){ return(0); }
j = (n-start)/2;
median = xr[start+j];
if(n-start > 2*j){
	median = (median + xr[start+j+1])/2;
}
return(median);
}

/***********************************************/
void	error_variance_noreps (FSTAT *fs, DATA *d, long *index, long window)
/***********************************************/
{
long i, k, iw, irow, icol, imax, imin;
float *deviate, x, M, errVar, maxErrVar;

maxErrVar=MISSING;
//printf("Window %d\n",window);
check(deviate = (float*)malloc(d->ncol*window*sizeof(float)));
for(iw=0; iw < d->nrow/(window/2)+1; ++iw){
	imin = iw*window/2-window/4;
	imax = imin+window;
	if(imin <0) imin=0;
	if(imax >= d->nrow) imax = d->nrow-1;
	if(imax <= imin) break;
	k = 0;
	for(i=imin; i < imax; ++i){
		irow = index[i];
		M = d->median_ref[irow];
		for(icol=0; icol < d->ncol; ++icol){
			x = d->data[icol][irow];
			if(x<=MISSING) continue;
			x = fabs(x - M);
			if(x > 2) continue;
			deviate[k++] = x;
		}
	}
	errVar = MISSING;
	if(d->ncol>1 && k>=1){
		errVar = median(deviate,k)/0.675;
		errVar = errVar*errVar*d->ncol/(d->ncol-1);
		if(iw >= d->nrow/(window/2)/2){
			if(maxErrVar < errVar) maxErrVar = errVar;
			else errVar = maxErrVar;
		}
	}
	for(i=iw*window/2; i < (iw+1)*window/2 && i<d->nrow; ++i){
		fs->errVar[index[i]] = errVar;
		fs->smoothErrVar[index[i]] = errVar;
	}
}
free(deviate);
}

/***********************************************/
void  anova (long irow, DATA *d, FSTAT *fs, PARAM *p)
/***********************************************/
{
float x, y, sumIntensity;
long icol, ilevel, n, naver;

/* Estimate means and df*/
fs->mean_ref[irow] = d->mean_ref[irow];
for(ilevel=0; ilevel<d->n_levels; ++ilevel){
	fs->means[ilevel][irow] = 0;
	fs->n[ilevel][irow] = 0;
}
n = 0;
ilevel = 0;
fs->df1[irow] = 0;
fs->df2[irow] = 0;
fs->glob_mean[irow] = 0;
for(icol=0;icol<d->ncol;++icol){
	x = d->data[icol][irow];
	if(x > MISSING){
		fs->means[ilevel][irow] += x;
		fs->n[ilevel][irow]++;
		fs->df2[irow]++;
	}
	++n;
	if(n == d->n_reps[ilevel]){
		if(fs->n[ilevel][irow]>0) { 
			fs->glob_mean[irow] += fs->means[ilevel][irow];
			fs->means[ilevel][irow] /= fs->n[ilevel][irow];
			fs->df1[irow]++;
		}
		else fs->means[ilevel][irow] = MISSING;	
		n = 0;
		++ilevel;
	}
}
if(fs->df2[irow])
	fs->glob_mean[irow] /= fs->df2[irow];
else
	fs->glob_mean[irow] = MISSING;
fs->df2[irow] -= fs->df1[irow];
fs->df1[irow]--;

/* Estimate variances */
fs->var[irow] = 0;
fs->errVar[irow] = 0;
n = 0;
ilevel = 0;
for(icol=0;icol<d->ncol;++icol){
	x = d->data[icol][irow];
	if(x > MISSING){
		x -= fs->means[ilevel][irow];
		fs->errVar[irow] += x*x;
		y = fs->means[ilevel][irow] - fs->glob_mean[irow];
		fs->var[irow] += y*y;
	}
	++n;
	if(n == d->n_reps[ilevel]){
		n = 0;
		++ilevel;
	}
}
if(fs->df2[irow]){
	fs->errVar[irow] /= fs->df2[irow];
}else{
	fs->errVar[irow] = MISSING;
}
if(fs->df1[irow])
	fs->var[irow] /= fs->df1[irow];
else
	fs->var[irow] = 0;
sumIntensity = 0;
naver = 0;
for(ilevel=0; ilevel<d->n_levels; ++ilevel){
	if(fs->means[ilevel][irow] > MISSING){ 
		sumIntensity += fs->means[ilevel][irow];
		naver++;
	}
}
if(naver>0) fs->averIntensity[irow] = sumIntensity/naver;
else fs->averIntensity[irow] = MISSING;
}

/***********************************************/
void   estimate_z (long irow, float *z, DATA *d, FSTAT *fs)
/***********************************************/
{
float x;
long icol, ilevel, n;

n = 0;
ilevel = 0;
for(icol=0;icol<d->ncol;++icol){
	z[icol] = 0;
	x = d->data[icol][irow];
	if(x > MISSING && fs->smoothErrVar[irow] > 0){
		x -= fs->means[ilevel][irow];
		z[icol] = x/sqrt(fs->smoothErrVar[irow]);
		if(z[icol] < 0) z[icol] *= -1;
	}
	++n;
	if(n==d->n_reps[ilevel]){
		n = 0;
		++ilevel;
	}
}
}

/***********************************************/
long   remove_outliers (long irow, float *z, DATA *d, PARAM *p)
/***********************************************/
{
long i, j, icol, n=0, nremoved = 0;
long iremove;
float M, x, dmax, dist;

M = d->mean_ref[irow];
icol=0;
for(i=0; i<d->n_levels; i++){
	if(d->n_reps[i]==2 && z[icol] > p->thresh_z){
		if(d->data[icol][irow]==MISSING || d->data[icol+1][irow]==MISSING){ icol+=2; continue; }
		if(fabs(d->data[icol][irow]-M) > fabs(d->data[icol+1][irow]-M)){ d->data[icol][irow] = MISSING; }
		else{ d->data[icol+1][irow] = MISSING; }
		icol+=2;
		++nremoved;
		continue;
	}
	dmax=0;
	iremove=-1;
	for(j=0; j<d->n_reps[i]; j++){
		if(z[icol] > p->thresh_z && d->data[icol][irow] != MISSING){
			dist = fabs(d->data[icol][irow]-M);
			if(dmax<dist){ dmax=dist; iremove=icol; }
		}
		icol++;
	}
	if(iremove>=0){
		d->data[iremove][irow] = MISSING;
		++nremoved;
	}
}
if(nremoved){
	d->mean_ref[irow]=0;
	for(icol=0; icol<d->ncol; icol++){
		x = d->data[icol][irow];
		if(x > MISSING){
			d->mean_ref[irow] += x;
			n++;
		}
	}
	d->mean_ref[irow] /= n;
}
return(nremoved);
}

/***********************************************/
void  recursive_estimation (DATA *d, FSTAT *fs, PARAM *p)
/***********************************************/
{
long icol, irow;
float *z;
long ncycles;
long noutliers=0, nmissing=0, noutliers_old=-1;
long *index;

for(irow=0; irow<d->nrow; ++irow)
	for(icol=0; icol<d->ncol; ++icol)
		if(d->data[icol][irow] <= MISSING) ++nmissing;

z = (float*)malloc(d->ncol*sizeof(float));
ncycles = 0;
do{
	for(irow=0; irow<d->nrow; ++irow)
		anova(irow, d, fs, p);
	index = get_index(fs->averIntensity, d->nrow);
	reverse(index, d->nrow);
	if(d->ncol > d->n_levels){
		smooth_error_variance(fs, index, p->window, p->propErrVar);
	}else if (d->ncol > 1){
		error_variance_noreps(fs, d, index, p->window);
	}
	free(index);
	//printf("C\n");
	if(noutliers==noutliers_old || ++ncycles >= 5) 
		break;
	for(irow=0; irow<d->nrow; ++irow){
		estimate_z(irow, z, d, fs);
		noutliers += remove_outliers (irow, z, d, p);
	}
	noutliers_old = noutliers;
}while(1);
if(nmissing) printf("Number of missing values = %d\n", nmissing);
if(noutliers) printf("Number of outliers removed = %d\n", noutliers);
free(z);
}

/***********************************************/
void  fisher_statistics (FSTAT *fs, PARAM *p)
/***********************************************/
{
long irow;
float errVar;

for(irow =0; irow<fs->nrow; ++irow){
	if(fs->var[irow] <= MISSING){
		fs->F[irow] = MISSING;
		fs->errVarFinal[irow] = MISSING;
		continue;
	}
	switch(p->errorModel){
		case ERR_ACTUAL: errVar = fs->errVar[irow]; break;
		case ERR_AVER: errVar = fs->smoothErrVar[irow]; break;
		case ERR_MAX_AVER_ACTUAL: errVar = fs->smoothErrVar[irow];
			if(errVar < fs->errVar[irow])
				errVar = fs->errVar[irow];
			break;
		case ERR_BAYESIAN:
		case ERR_MAX_AVER_BAYESIAN: 
			if(fs->df2[irow] >= p->bayesianDF)
				errVar = fs->errVar[irow];
			else{
				if(fs->errVar[irow] <= MISSING)
					errVar = fs->smoothErrVar[irow];
				else
					errVar = (fs->df2[irow]*fs->errVar[irow] 
					 +(p->bayesianDF-fs->df2[irow])*fs->smoothErrVar[irow])/p->bayesianDF;
			}
			if(p->errorModel==ERR_MAX_AVER_BAYESIAN && errVar < fs->smoothErrVar[irow])
				errVar = fs->smoothErrVar[irow];
			break;
	}
	fs->errVarFinal[irow] = errVar;
	if(errVar <= MISSING){
		if(fs->var[irow]==0) fs->F[irow] = 0;
		else fs->F[irow] = MISSING;
	}else{
		if(errVar <= 0){
			errVar = 0.00001;
			fs->errVarFinal[irow] = errVar;
		}
		fs->F[irow] = fs->var[irow]/errVar;
	}
}
}

/***********************************************/
void  read_to_the_end_of_line (FILE *fp)
/***********************************************/
{
char a;
while((a=fgetc(fp))!=EOF && a != '\n');
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
char *buffer;
FILE *fp;
long n;

fp = fopen(filename, "r");
if(!fp){
	error_message("Input file not found\n");
}
check(buffer = (char*)malloc(sizeof(char)*100000));
fgets(buffer,100000,fp);		/* read headers */
*ncol = count_fields(buffer, 100000);
*nrow = 0;
while(fgets(buffer,100000,fp)){
	n = count_fields(buffer, 100000);
	if(n == *ncol && strlen(buffer) - n > 2)
		++(*nrow);
	else
		break;
}
free(buffer);
fclose(fp);
}

/***********************************************/
void  read_tabdelimited_item (FILE *fp, char *buffer)
/***********************************************/
{
char ch;
long i=0, j=0;
while((ch=fgetc(fp)) && (j==0 || ch != '\t' && ch != '\n') && ch != EOF){
	if(j==0 && (ch=='\t' || ch=='\n')){
		j++;
		continue;
	}
	j++;
	if(i==0 && ch==' ') continue;
	if(i < STRINGLEN-1) buffer[i++] = ch;	
}
buffer[i] = '\0';
if(ch == '\n' || ch == '\t') ungetc (ch, fp);
}

/***********************************************/
int  split_string (char *string, char *items[], int num)
/***********************************************/
{
char *ch;
int i=0;

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
void  lower_case (char *s, char *low)
/***********************************************/
{
long i;
for(i=0;i<=strlen(s);i++){
	if(s[i]>=65 && s[i]<=90) s[i] = s[i]+32;
}
return;
}

/***********************************************/
char  *string_replace (char *s, char *find, char *replace)
/***********************************************/
{
long i, len_s, len_f, len_r, j, len_new;
char *new_string;

len_s = strlen(s);
len_f = strlen(find);
len_r = strlen(replace);
if(!len_s || !len_f){ return(s); }
len_new = len_s + 1;
if(len_r > len_f) len_new = (len_r-len_f)*(len_s/len_f+1)+1;
if(len_new > len_s+100000){ len_new = len_s+100000; }

check(new_string = (char*)calloc(len_new,sizeof(char)));
j=0;
for(i=0; i<len_s-len_f+1; ++i){
	if(!strncmp(find,s+i,len_f)){
		strcpy(new_string+j,replace);
		j+=len_r;
		i+=len_f-1;
	}else{
		new_string[j++]=s[i];
	}
}
new_string[j] = '\0';
new_string = realloc(new_string,(j+1)*sizeof(char));
free(s);
return(new_string);
}

/***********************************************/
DATA *read_data (PARAM *p) 
/***********************************************/
{
FILE *inputFile;
DATA *d;
long i, j, n, k, nitems, ilevel, pos, len, reference=0, *reorder, *used, count;
float x, xmin=1E10, xmax=-1E10, *values, logBase, adjustment, lowvalue;
char *buffer, **items, **headers=NULL, number[10], *ch;
long buff_size = 100000;

check(buffer = (char*)malloc(buff_size*sizeof(char)));
check(items = (char**)malloc(50000*sizeof(char*)));
check(d = (DATA*)calloc(1,sizeof(DATA)));

inputFile = fopen(p->inputFile, "r");
if(!inputFile){
	error_message("Input file not found\n");
}
while(fgets(buffer,buff_size-1,inputFile)){
	len = strlen(buffer);
	if(buffer[len-1]=='\n'){ buffer[len-1]='\0'; --len; }
	if(buffer[len-1]=='\r'){ buffer[len-1]='\0'; --len; }
	if(!len){ continue; }
	if(strncmp(buffer,"!",1)) break;
	if(!strncmp(buffer,"!Sample_title",13) || !strncmp(buffer,"!sample_title",13)){
		nitems = split_string (buffer,items,50000);
		d->ncol = nitems-1;
		check(headers = (char**)malloc(d->ncol*sizeof(char*)));
		for(i=0; i<d->ncol; ++i){
			headers[i] = copy_string(items[i+1]);
			headers[i] = string_replace(headers[i],"\"","");
			//ch = strstr(headers[i],"_rep");
			//if(ch && ch > headers[i]+strlen(headers[i])-15) *ch = '\0';
			//ch = strstr(headers[i]," rep");
			//if(ch && ch > headers[i]+strlen(headers[i])-15) *ch = '\0';
		}
	}
	if(!strncmp(buffer,"!Sample_data_row_count",22) || !strncmp(buffer,"!sample_data_row_count",22)){
		split_string (buffer,items,50000);
		if(items[1][0]=='\"') items[1]++;
		if(sscanf(items[1],"%d",&d->nrow)==0 || d->nrow <10){
			printf("N rows = %d is not correct",d->nrow); exit(1);
		}
	}
}
if(!headers){
	nitems = split_string (buffer,items,50000);
	d->ncol = nitems-1;
	check(headers = (char**)malloc(d->ncol*sizeof(char*)));
	for(i=0; i<d->ncol; ++i){
		headers[i] = copy_string(items[i+1]);
		headers[i] = string_replace(headers[i],"\"","");
		//ch = strstr(headers[i],"_rep");
		//if(ch && ch > headers[i]+strlen(headers[i])-15) *ch = '\0';
		//ch = strstr(headers[i]," rep");
		//if(ch && ch > headers[i]+strlen(headers[i])-15) *ch = '\0';
	}
}
if(!d->nrow){
	while(fgets(buffer,buff_size-1,inputFile)){
		if(!strncmp(buffer,"!",1)) break;
		d->nrow++;
	}
	rewind(inputFile);
	while(fgets(buffer,buff_size-1,inputFile)){
		if(!strncmp(buffer,"\n",1)) continue;
		if(strncmp(buffer,"!",1)) break;
	}
}
if(p->window > d->nrow){ p->window = d->nrow; }
	
//printf("A %d %d\n",d->nrow,d->ncol);

/* Determine the order of samples */
check(reorder = (long*)malloc(d->ncol*sizeof(long)));
check(used = (long*)calloc(d->ncol,sizeof(long)));
check(d->n_reps = (long*)calloc(d->ncol,sizeof(long)));
ilevel=0;
pos=0;
for(i=0;i<d->ncol;++i){
	if(used[i]) continue;
	d->level_name[ilevel] = copy_string(headers[i]);
	d->n_reps[ilevel]++;
	reorder[i] = pos++;
	for(k=i+1; k<d->ncol; ++k){
		if(!strcmp(headers[k],headers[i])){
			reorder[k] = pos++;
			used[k] = 1;
			d->n_reps[ilevel]++;
		}
	}
	ilevel++;
}
d->n_levels = ilevel;
d->totrep = d->ncol;
j = 0;
for(i=0;i<d->n_levels;++i){
	for(k=0; k<d->n_reps[i]; ++k){
		check(d->header_col[j] = (char*)calloc((strlen(d->level_name[i])+10),sizeof(char)));
		strcpy(d->header_col[j],d->level_name[i]);
		sprintf(number,"_rep%d",k+1);
		strcat(d->header_col[j],number);
	}
}
for(i=0;i<d->ncol;++i){
	check(d->data[i] = (float*)malloc(d->nrow*sizeof(float)));
}

check(d->mean_ref = (float*)calloc(d->nrow,sizeof(float)));
check(d->median_ref = (float*)calloc(d->nrow,sizeof(float)));
check(d->header_row = (char**)malloc(d->nrow*sizeof(char*)));
check(values = (float*)malloc(d->ncol*sizeof(float)));
d->header1 = copy_string("SampleID");

k=0;
for(i=0;i<d->nrow;++i){
	if(!fgets(buffer,buff_size-1,inputFile)){
		d->nrow = i;
		break;
	}
	if(!strncmp(buffer,"!",1) || strlen(buffer)<3){
		d->nrow = i;
		break;
	}
	split_string (buffer,items,50000);
	d->header_row[i] = copy_string(items[0]);
	d->header_row[i] = string_replace(d->header_row[i],"\"","");
	n = 0;
	xmax = -1.0e10;
	for(j=0;j<d->ncol;++j){
		if(sscanf(items[j+1],"%f", &x) != 1){
			printf("Error reading: row %d, col %d\n", i+2, j);
			d->data[reorder[j]][i] = MISSING;
		}else if(x==MISSING){
			d->data[reorder[j]][i] = MISSING;
		}else{
			if(p->cutoff>=0.000001 && x < p->cutoff){  x = p->cutoff*(0.3+0.7*(float)rand()/RAND_MAX); }
			if(xmax<x) xmax=x;
			x = log10(x);
			d->data[reorder[j]][i] = x;
			d->mean_ref[i] += x;
			values[n++] = x;
		}
	}
	if(n>0) d->median_ref[i] = median(values,n);
	if(n<1 || xmax<p->cutoff){
		d->mean_ref[i]=0;
		free(d->header_row[i]);
		i--;
		d->nrow--;
	}
	if(n>0)	d->mean_ref[i] /= n;
}
fclose(inputFile);
//printf("Data loaded\n");

for(i=0;i<d->ncol;++i){ free(headers[i]); }
free(headers);
free(buffer);
free(items);
free(reorder);
free(used);
free(values);
return(d);
}

/***********************************************/
char *uc (char *a) 
/***********************************************/
{
long i=0;
if(!a) return(a);
while(a[i]){ a[i]=toupper(a[i]); i++; }
return(a);
}

/***********************************************/
ANNOTATION *read_annotation (char *filename, DATA *d, long *nAnnot, HASH **hashAnnot) 
/***********************************************/
{
FILE *fp;
long i, n, nrows=0, nmeta=0, len, i1, i2;
long probeNames=0, N1=0, N4=0;
char *buffer, **items;
long buff_size = 10000;
ANNOTATION *annot;
HASH *hashID;
HASH *hashName;
HASH *hashRowID;

check(buffer = (char*)malloc(buff_size*sizeof(char)));
check(items = (char**)malloc(50*sizeof(char*)));
hashID = new_hash(1,10000);
hashName = new_hash(1,10000);
hashRowID = new_hash(1,10000);
for(i=0; i<d->nrow; ++i){
	if(d->header_row[i]) add_hash_long(d->header_row[i],1,0,hashRowID);
}
fp = fopen(filename, "r");
if(!fp){
	error_message("Input file not found\n");
}
while(fgets(buffer,buff_size-1,fp)){
	if(buffer[0]=='!' && nrows==0) nmeta++;
	else nrows++;
}
check(annot = (ANNOTATION*)calloc(nrows,sizeof(ANNOTATION)));
rewind(fp);
for(i=0; i<nmeta; ++i){
	fgets(buffer,buff_size-1,fp);
}
for(i=0; i<nrows; ++i){
	fgets(buffer,buff_size-1,fp);
	len = strlen(buffer);
	if(buffer[len-1]=='\n'){ buffer[len-1]='\0'; --len; }
	if(!len){ continue; }
	n = split_string (buffer,items,50);
	if(i==0 && n>3 && !strcmp(uc(items[3]),"PROBENAME")){
		probeNames = 1;
	}
	if(n>0) annot[i].id = copy_string(items[0]);
	if(n>1) annot[i].symbol = copy_string(items[1]);
	if(n>2) annot[i].geneName = copy_string(items[2]);
	if(annot[i].id && hash_value_long(annot[i].id,hashRowID)>=0) add_hash_long(annot[i].id,i,0,hashID);
	if(probeNames){
		if(n>3) annot[i].probeName = copy_string(items[3]);
		if(n>4) annot[i].add1 = copy_string(items[4]);
		if(n>5) annot[i].add2 = copy_string(items[5]);
		if(annot[i].probeName && hash_value_long(annot[i].probeName,hashRowID)>=0) add_hash_long(annot[i].probeName,i,0,hashName);
	}else{	
		if(n>3) annot[i].add1 = copy_string(items[3]);
		if(n>4) annot[i].add2 = copy_string(items[4]);
	}
}
fclose(fp);
N1 = hashID->n;
N4 = hashName->n;
if(probeNames && N4 > N1){
	for(i=0; i<nrows; ++i){
		if(annot[i].id) free(annot[i].id);
		annot[i].id = annot[i].probeName;
	}
	*hashAnnot = hashName;
	destroy_hash(hashID);
}else{
	*hashAnnot = hashID;
	destroy_hash(hashName);
}
*nAnnot = nrows;
free(buffer);
free(items);
destroy_hash(hashRowID);
return(annot);
}

/***********************************************/
float gammln(float xx)
/***********************************************
* Returns the value ln[Ã(xx)] for xx > 0. */
{
double x,y,tmp,ser;
static double cof[6]={76.18009172947146,-86.50532032941677,
24.01409824083091,-1.231739572450155,
0.1208650973866179e-2,-0.5395239384953e-5};
long j;

y=x=xx;
tmp=x+5.5;
tmp -= (x+0.5)*log(tmp);
ser=1.000000000190015;
for (j=0;j<=5;j++) 
	ser += cof[j]/++y;
return -tmp+log(2.5066282746310005*ser/x);
}

/***********************************************/
float betacf(float a, float b, float x)
/***********************************************
* Evaluates continued fraction for incomplete beta function */
{
long m,m2;
float aa,c,d,del,h,qab,qam,qap;

qab=a+b;
qap=a+1.0;
qam=a-1.0;
c=1.0;
d=1.0-qab*x/qap;
if (fabs(d) < FPMIN) d=FPMIN;
d=1.0/d;
h=d;
for (m=1;m<=MAXIT;m++) {
	m2=2*m;
	aa=m*(b-m)*x/((qam+m2)*(a+m2));
	d=1.0+aa*d;
	if (fabs(d) < FPMIN) d=FPMIN;
	c=1.0+aa/c;
	if (fabs(c) < FPMIN) c=FPMIN;
	d=1.0/d;
	h *= d*c;
	aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
	d=1.0+aa*d;
	if (fabs(d) < FPMIN) d=FPMIN;
	c=1.0+aa/c;
	if (fabs(c) < FPMIN) c=FPMIN;
	d=1.0/d;
	del=d*c;
	h *= del;
	if (fabs(del-1.0) < EPS || m > MAXIT) break;
}
//if (m > MAXIT) printf("No convergence in betacf\n");
return h;
}

/***********************************************/
float betai(float a, float b, float x)
/***********************************************
* Returns the incomplete beta function Ix(a, b). */
{
float bt;

if (x < 0.0 || x > 1.0){
	printf("Bad x in routine betai\n");
	return(1);
}
if (x == 0.0 || x == 1.0) bt=0.0;
else bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
if (x < (a+1.0)/(a+b+2.0))
	return bt*betacf(a,b,x)/a;
return 1.0-bt*betacf(b,a,1.0-x)/b;
}

/**********************************************/
double     normal_distribution     (double x)
/**********************************************/
{
static double a1=-1.26551223, a2= 1.00002368, a3= 0.37409196, a4= 0.09678418,
      a5=-0.18628806, a6= 0.27886807, a7=-1.13520398, a8= 1.48851587, a9=-0.82215223,
      a10=0.17087277;
double z, t, y;

z = fabs((double)x)/sqrt(2.);
   t = 1.0 / (1.0 + 0.5 * z);
   y = t*exp(-z * z + a1 + t * (a2 + t * (a3 + t * (a4 + t * (a5 + t *
     (a6 + t * (a7 + t * (a8 + t * (a9 + t * a10)))))))));
   if(x < 0.0) y = 2.0 - y;
   y = 1.0 - 0.5 * y;
   return(y);
}

/**********************************************/
double normal_inverse(double p)
/**********************************************/
{
const double a[6] = {  -3.969683028665376e+01,  2.209460984245205e+02, -2.759285104469687e+02,  1.383577518672690e+02,
  -3.066479806614716e+01,  2.506628277459239e+00};
const double b[5] = {  -5.447609879822406e+01,  1.615858368580409e+02, -1.556989798598866e+02,  6.680131188771972e+01,
  -1.328068155288572e+01 };
const double c[6] = { -7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00,
   4.374664141464968e+00,  2.938163982698783e+00 };
const double d[4] = { 7.784695709041462e-03,  3.224671290700398e-01, 2.445134137142996e+00,  3.754408661907416e+00 };

register double q, t, u;

if(p >= 1 || p <= 0){ error_message("p-value out of range\n"); }
q = p; if(q>1-p) q=1-p;
if (q > 0.02425) {
  /* Rational approximation for central region. */
  u = q-0.5;
  t = u*u;
  u = u*(((((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4])*t+a[5])
    /(((((b[0]*t+b[1])*t+b[2])*t+b[3])*t+b[4])*t+1);
} else {
  /* Rational approximation for tail region. */
  t = sqrt(-2*log(q));
  u = (((((c[0]*t+c[1])*t+c[2])*t+c[3])*t+c[4])*t+c[5])
   /((((d[0]*t+d[1])*t+d[2])*t+d[3])*t+1);
}
/* The relative error of the approximation has absolute value less
    than 1.15e-9.  One iteration of Halley's rational method (third
    order) gives full machine precision... */
t = normal_distribution(u)-q;    /* error */
t = t*M_SQRT2PI*exp(u*u/2);   /* f(u)/df(u) */
u = u-t/(1+u*t/2);     /* Halley's method */

if(p > 0.5) return(-u);
return(u);
}

/***********************************************/
void print_data (DATA *d, FILE *outputFile)
/***********************************************/
{
long i, j;

fprintf(outputFile,"%s\t",d->header1);
for(i=0;i<d->ncol;++i)
	fprintf(outputFile,"%s\t", d->header_col[i]);
fprintf(outputFile,"\n");
for(j=0;j<d->nrow;++j){
	fprintf(outputFile,"%s\t",d->header_row[j]);
	for(i=0;i<d->ncol;++i)
		fprintf(outputFile,"%.4f\t", d->data[i][j]);
	fprintf(outputFile,"\n");
}
return;
}

/***********************************************/
void print_output  (DATA *d, FSTAT *fs, PARAM *p, FILE *outputFile, ANNOTATION *annot, long nAnnot, HASH *hashAnnot)
/************************************************/
{
long irow, i, nn=0, nsignif=0, iii, icol, nAnnotated=0;
long rank, ilevel, n, DF2, iannot=1, iannot1, point;
float x, y, prob, FDR1=1;
float *F1, *pvalue, *FDR, *index;
char buffer[100];

F1 = (float*)calloc(fs->nrow,sizeof(float));
pvalue = (float*)calloc(fs->nrow,sizeof(float));
FDR = (float*)calloc(fs->nrow,sizeof(float));
index = (float*)calloc(fs->nrow,sizeof(float));
n = d->n_levels*(d->n_levels-1)/2;

for(irow =0; irow<fs->nrow; ++irow){
	F1[irow] = fs->F[irow];
	x = fs->averIntensity[irow];
	++nn;
}
//printf("Number of columns = %d\n", d->totrep);
//printf("Number of factor levels = %d\n", d->n_levels);
//printf("Number of rows = %d\n", nn);
fs->rank = get_rank(F1, fs->nrow, DESCENDING, index);
DF2 = (d->ncol - d->n_levels)*p->window;
if(!DF2) DF2 = d->ncol*p->window/2;
if(p->errorModel == ERR_BAYESIAN) DF2 = p->bayesianDF;
for(i=d->nrow-1; i>=0; --i){
	irow = index[i];
	rank = fs->rank[irow]+1;
	if(p->errorModel == ERR_ACTUAL){
		DF2 = fs->df2[irow];
		if(!DF2) DF2 = d->ncol/2;
	}
	if(F1[irow] <= MISSING) prob = MISSING;
	else if (F1[irow] <1.05) prob = 1;
	else{
		//printf("%f %d %d\n",fs->F[irow],fs->df1[irow],DF2);
		prob = betai(0.5*DF2, 0.5*fs->df1[irow], DF2/(DF2+fs->df1[irow]*fs->F[irow]));
	}
	pvalue[irow] = prob;
	if(prob <= MISSING){
		FDR[irow] = MISSING;
	}
	else{
		FDR[irow] = prob*nn/rank;
		if(FDR[irow] > FDR1) FDR[irow] = FDR1;
		else FDR1 = FDR[irow];
		if(FDR[irow] <= p->FDR) ++nsignif;
	}
}

/* Print ANOVA results */
fprintf(outputFile,"%s\tAverIntensity",d->header1);
for(i=0;i<d->n_levels;++i){
	fprintf(outputFile,"\t%s (%d)", d->level_name[i],d->n_reps[i]);
}
fprintf(outputFile,"\tVar(Factor)");
fprintf(outputFile,"\tVar(Err)\tSmoothVar(Err)\tFinalMSE");
fprintf(outputFile,"\tF");
fprintf(outputFile,"\tP\tFDR\trank");
if(nAnnot){
	if(annot[0].symbol) fprintf(outputFile,"\t%s", annot[0].symbol);
	if(annot[0].geneName) fprintf(outputFile,"\t%s", annot[0].geneName);
	if(annot[0].add1) fprintf(outputFile,"\t%s", annot[0].add1);
	if(annot[0].add2) fprintf(outputFile,"\t%s", annot[0].add2);
}else if(p->use_probeID){
	fprintf(outputFile,"\tSymbol");
}
fprintf(outputFile,"\n");
for(irow=0; irow<d->nrow; ++irow){
	fprintf(outputFile,"%s",d->header_row[irow]);
	if(fs->averIntensity[irow]>MISSING) fprintf(outputFile,"\t%.4f", fs->averIntensity[irow]);
	else fprintf(outputFile,"\t%d", MISSING);

	/* Print the means */
	for(ilevel=0; ilevel<d->n_levels; ++ilevel){
		if(fs->means[ilevel][irow] <= MISSING) fprintf(outputFile,"\t%d", MISSING);
		else fprintf(outputFile,"\t%.4f", fs->means[ilevel][irow]);
	}
	fprintf(outputFile,"\t%.5f", fs->var[irow]);
	fprintf(outputFile,"\t%.5f\t%.5f\t%.5f", fs->errVar[irow], fs->smoothErrVar[irow], fs->errVarFinal[irow]);
	if(fs->F[irow] <= MISSING)
		fprintf(outputFile,"\t%d", MISSING);
	else 
		fprintf(outputFile,"\t%.3f", fs->F[irow]);
	rank = fs->rank[irow]+1;
	fprintf(outputFile,"\t%.5f\t%.5f\t%d", pvalue[irow], FDR[irow], rank);
	if(nAnnot){
		iannot = hash_value_long(d->header_row[irow],hashAnnot);
		//printf("%s %d \n",d->header_row[irow],iannot);
		if(iannot<0){
			point = accountNum(d->header_row[irow]);
			if(point && point<100){
				strcpy(buffer,d->header_row[irow]);
				buffer[point]='\0';
				iannot = hash_value_long(buffer,hashAnnot);
			}
		}
		if(p->use_probeID==2){
			fprintf(outputFile,"\t%s",d->header_row[irow]);
			nAnnotated++;
			if(iannot>=0 && annot[iannot].id && !strcmp(annot[iannot].symbol, d->header_row[irow])){
				if(annot[iannot].geneName) fprintf(outputFile,"\t%s", annot[iannot].geneName);
				if(annot[iannot].add1) fprintf(outputFile,"\t%s", annot[iannot].add1);
				if(annot[iannot].add2) fprintf(outputFile,"\t%s", annot[iannot].add2);
			}
		}else if(iannot>=0 && annot[iannot].id){
			nAnnotated++;
			if(annot[iannot].symbol && strlen(annot[iannot].symbol)>0) fprintf(outputFile,"\t%s", annot[iannot].symbol);
			if(annot[iannot].geneName) fprintf(outputFile,"\t%s", annot[iannot].geneName);
			if(annot[iannot].add1) fprintf(outputFile,"\t%s", annot[iannot].add1);
			if(annot[iannot].add2) fprintf(outputFile,"\t%s", annot[iannot].add2);
		}else if(p->use_probeID==1){
			fprintf(outputFile,"\t%s",d->header_row[irow]);
		}
	}else if(p->use_probeID){
		fprintf(outputFile,"\t%s",d->header_row[irow]);
	}
	fprintf(outputFile,"\n");
}
//printf("%d significant genes found\n", nsignif);
free(F1);
free(FDR);
free(pvalue);
free(index);
return;
}

/****************************/
void print_line (FILE *fp)
/***************************/
{
char buffer[51];
fgets(buffer, 50, fp);
buffer[50] = '\0';
printf("%s\n", buffer);
}

/**********************************************/
float      read_float   (FILE *fp, char *description)
/**********************************************/
{
float x;

if (fscanf (fp,"%f",&x) != 1)
{
	print_line(fp);
	error_message("Reading parameter\n");
}
//printf("%f\t%s\n", x, description);
return (x);
}

/**********************************************/
long   read_int  (FILE *fp, char *description)
/**********************************************/
{
long x, index;
int nnn;

if (fscanf (fp,"%d",&nnn) != 1)
{
	print_line(fp);
	error_message("Error reading parameter\n");
}
x = nnn;
index = fgetc (fp);
if (index == '.')   /* Float number is read instead of integer */
{
	print_line(fp);
	printf("Integer is expected, in parameter %s", description);
	exit(1);
}
else
	ungetc (index, fp);
//printf("%d\t%s\n", x, description);
return (x);
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
PARAM *read_parameters (int nargs, char **argv)
/***********************************************/
{
PARAM *p;
int iarg=1, score=0, i;

if(nargs < 5){ error_message("anova_oneway -i inputFile -o outputFile [-a annotationFile, -err errorModel, -z ZvalueOutliers, -FDR maxFDR, -prop propErrVar, -win window, -bayes bayesianDF\n"); }
p = (PARAM*)calloc(1,sizeof(PARAM));
p->window = 500;
p->propErrVar = 0.01;
p->errorModel=4;
p->thresh_z=8;
p->bayesianDF=10;
p->FDR=0.05;
p->cutoff=0;
p->use_probeID=0;

while(iarg < nargs){
	if(!strcmp(argv[iarg],"-i") && iarg < nargs-1) p->inputFile = copy_string(argv[++iarg]);
	else if(!strcmp(argv[iarg],"-o") && iarg < nargs-1) p->outputFile = copy_string(argv[++iarg]);
	else if(!strcmp(argv[iarg],"-a") && iarg < nargs-1) p->annotationFile = copy_string(argv[++iarg]);
	else if(!strcmp(argv[iarg],"-err") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->errorModel);
	else if(!strcmp(argv[iarg],"-z") && iarg < nargs-1) sscanf(argv[++iarg],"%f",&p->thresh_z);
	else if(!strcmp(argv[iarg],"-cutoff") && iarg < nargs-1){
		sscanf(argv[++iarg],"%f",&p->cutoff);
		if(p->cutoff<1.0e-8) p->cutoff=1.0e-8;
	}
	else if(!strcmp(argv[iarg],"-FDR") && iarg < nargs-1) sscanf(argv[++iarg],"%f",&p->FDR);
	else if(!strcmp(argv[iarg],"-bayes") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->bayesianDF);
	else if(!strcmp(argv[iarg],"-prop") && iarg < nargs-1) sscanf(argv[++iarg],"%f",&p->propErrVar);
	else if(!strcmp(argv[iarg],"-win") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->window);
	else if(!strcmp(argv[iarg],"-probe") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->use_probeID);
	else{
		printf("ERROR: Wrong option %s\n", argv[iarg]);
		exit(1); 
	}
	++iarg;
}
if(p->errorModel>5 || p->errorModel<0){ p->errorModel=4; printf("Wrong err. model (0-5 are acceptable). Restore default = 4\n"); }
if(p->thresh_z<=2){ p->thresh_z=8; printf("Wrong z-threshold (2 and up are acceptable). Restore default = 8\n"); }
if(p->FDR<0 || p->FDR>1){ p->FDR=0.05; printf("Wrong FDR-threshold (0-1 are acceptable). Restore default = 0.05\n"); }
if(p->bayesianDF<1){ p->bayesianDF=10; printf("Wrong Bayesian DF (>0 are acceptable). Restore default = 10\n"); }
if(p->propErrVar<0 || p->propErrVar>0.5){ p->propErrVar=0.01; printf("Wrong proportion of highest err. variances to be removed (0 - 0.5 are acceptable). Restore default = 0.01\n"); }
if(p->window<2 || p->window>3000){ p->window=500; printf("Wrong err. variance averaging window (2-3000 are acceptable). Restore default = 500\n"); }

/*
printf("%f\tThreshold z-value to remove outliers\n",p->thresh_z);
printf("%f\tProportion of highest err. variances to be removed\n",p->propErrVar);
printf("%f\tFDR threshold\n",p->FDR);
printf("%d\tErr. variance averaging window\n",p->window);
printf("%d\tErr. Model\n",p->errorModel);
if(p->cutoff>-0.5e10) printf("%f\tCutoff value\n",p->cutoff);
if(p->errorModel==3 || p->errorModel==5) printf("%d\tBayesian degrees of freedom\n",p->bayesianDF);
*/
return(p);
}

/***********************************************/
int main (long argc, char **argv) 
/***********************************************/
{
DATA *d, *d1;
PARAM *p;
FSTAT *fs;
FILE *fpout;
float *dist;
long n, i, nAnnot=0;
ANNOTATION *annot=NULL;
HASH *hashAnnot=NULL;

p = read_parameters(argc, argv);
d = read_data(p);
if(p->annotationFile){
	annot = read_annotation(p->annotationFile, d, &nAnnot, &hashAnnot);
}

/* ANOVA */
fpout = fopen(p->outputFile, "w");
if(!fpout){
	printf("Output file %s not opened!\n", p->outputFile);
	exit(1);
}
fs = init_fstat(d);
recursive_estimation (d, fs, p);
fisher_statistics (fs, p);
print_output(d, fs, p, fpout, annot, nAnnot, hashAnnot);
fclose(fpout);
return(0);
}

