/********************************************
* pairwise comparison of microarray averages based on ANOVA.
* All cell types (or tissues) are compared to one selected
* cell type/tissue or median of all cell types/tissues
* Programmer: Alexei Sharov (sharoval@mail.nih.gov)
* National Institute on Aging, Genetics Lab
pairwise anova-KeioTFs_Mar2015_combined.txt output.txt 0 -symbols
*********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define MAXCOL  5000
#define MISSING -9999

typedef struct data_st{
	long ncol;
	long nrow;
	float *data[MAXCOL];
	float *F;
	float *FDR;
	float *MSE;
	long *nrepl;
	char **header_col;
	char **header_row;
}DATA;
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
}

/*************************************************************************/
void sort_char   (long ie, char **a, long *b)
/*************************************************************************/
{
long i, j, k, m, p, q, iring, cmp;
long lt[64], ut[64];
char *ta, *xa;
long tb, xb;
long i1;

/* Function Body */
j = ie-1;
m = 1;
i = 0;
iring = 2;

/* If this segment has more than two elements  we split it */
L10:	 if ((i1 = j - i - 1) < 0) {
	goto L100;
	 } else if (i1 == 0) {
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
    if (strcmp(a[k],ta)<=0) goto L20;
/* Such an element has now been found now search for a q such that a(q)<t 
* starting at the end of the segment. */
L30:  if (strcmp(a[q],ta)<0) {
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
	if (strcmp(a[i],a[j])<=0) goto L100;
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
}

/***********************************************/
long hash_word_pos (char *input, HASH *h) 
/***********************************************/
{
long n, j, i, cmp, value_len=0;
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
long i,j,k;
if(!h) error_message("null hash in add_hash_long");
if(h->type>1) error_message("incorrect hash type long");
if(!word) error_message("null in add_hash_long");
if(h->n == h->n_alloc-1){
	h->n_alloc += h->increment;
	check(h->list = (char**)realloc(h->list,h->n_alloc*sizeof(char*)));
	if(h->type==1) check(h->value1 = (long*)realloc(h->value1,h->n_alloc*sizeof(long)));
}
j = hash_word_pos(word, h);
if(j>=0 && !strcmp(word,h->list[j])){
	if(h->type==1){
		if(add) h->value1[j] += value;
		else h->value1[j] = value;
	}
	return;
}
//if(j>0 && j<h->n-1) printf("F %d %d %s %s %s\n",j,h->n,word,h->list[j],h->list[j+1]);
//else if(j>0) printf("F %d %d %s %s\n",j,h->n,word,h->list[j]);
//else printf("F %d %s\n",j,word);

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
void  get_dimensions  (char *filename, long *nrow, long *ncol)
/***********************************************/
{
char *buffer, **items;
long buffer_size = 100000, i, nitems;
FILE *fp;

fp = fopen(filename, "r");
if(!fp) error_message("Input file not found");
check(buffer = (char*)malloc(sizeof(char)*buffer_size));
check(items = (char**)malloc(50000*sizeof(char*)));
fgets(buffer,buffer_size-1,fp);
nitems = split_string (buffer,items,50000);
for(i=2; i<nitems; i++){
	if(!strncmp("Var(",items[i],4)){ *ncol = i-2; break; }
}
*nrow = 0;
while(fgets(buffer,buffer_size,fp) && strlen(buffer)>2){
	if(buffer[0]=='!' || strlen(buffer)<2+(*ncol)*2) break;
	++(*nrow);
}
free(buffer);
free(items);
fclose(fp);
return;
}

/***********************************************/
DATA *read_data (char *filename, long ncol, long nrow, long backgr, long useSymbols, 
                 long redundant, float *minMSE, float exprThresh) 
/***********************************************/
{
FILE *inputFile, *fp;
DATA *d;
HASH *symbol2row;
long i, j, nitems, count=0, norm=0, comments=0, nRep=0;
long buff_size = 100000, len, irow, irow1;
char *buffer, **items, *ch, *name, *symbol=NULL;
float x, *data_row, median, F, MSE, baseline;

if(ncol >= MAXCOL-1) error_message("Increase MAXCOL!");
check(buffer = (char*)malloc(sizeof(char)*buff_size));
check(items = (char**)malloc(50000*sizeof(char*)));
check(data_row = (float*)malloc(ncol*sizeof(float)));
check(d = (DATA*)malloc(sizeof(DATA)));
d->nrow = nrow;
d->ncol = ncol;
check(d->header_row = (char**)calloc(nrow,sizeof(char*)));
check(d->header_col = (char**)calloc(ncol,sizeof(char*)));
for(i=0; i<ncol; ++i){
	check(d->data[i] = (float*)malloc(d->nrow*sizeof(float)));
}
check(d->nrepl = (long*)calloc(ncol,sizeof(long)));
check(d->F = (float*)malloc(nrow*sizeof(float)));
check(d->FDR = (float*)malloc(nrow*sizeof(float)));
check(d->MSE = (float*)malloc(nrow*sizeof(float)));
symbol2row = new_hash(1,100);

inputFile = fopen(filename, "r");
if(!inputFile) error_message("Input file not found\n");
fgets(buffer,buff_size-1,inputFile);
nitems = split_string (buffer,items,50000);
for(i=0;i<ncol;++i){
	ch = items[i+2];
	if(!strncmp("Mean(",ch,5)){
		ch+=5;
		len = strlen(ch);
		ch[len-1] = '\0';
		d->header_col[i] = copy_string(ch);
		d->nrepl[i] = 1;
	}else{
		len = strlen(ch);
		j = len-1;
		if(ch[j] == ')') ch[j--] = '\0';
		while(j>0 && isdigit(ch[j])) j--;
		if(j<len-2) d->nrepl[i] = atoi(ch+j+1);
		if(!d->nrepl[i]) d->nrepl[i]=1;
		if(ch[j] == '(') ch[j--] = '\0';
		if(ch[j] == ' ') ch[j]='\0';
		d->header_col[i] = copy_string(ch);
	}
	nRep += d->nrepl[i];
}
if(!redundant){ nrow=0; }
for(i=0;i<d->nrow;++i){
	fgets(buffer,buff_size-1,inputFile);
	nitems = split_string (buffer,items,50000);
	irow = i;
	F = atof(items[ncol+4+2]);
	if(nRep<=d->ncol){
		F = atof(items[1]);
	}
	MSE = atof(items[ncol+3+2]);
	if(MSE>0 && *minMSE>MSE) *minMSE=MSE;
	name = items[0];
	if(nitems > ncol+8+2){
		symbol = items[ncol+8+2];
	}
	if(!redundant){
		if(nitems <= ncol+8+2){
			if(useSymbols) continue;
			symbol = name;
		}
		if(strlen(symbol)==0) continue;
		irow1 = hash_value_long(symbol,symbol2row);
		if(irow1>=0){
			if(d->F[irow1] > F) continue;
			irow = irow1;
		}else{
			irow = nrow;
			add_hash_long(symbol,nrow++,0,symbol2row);		
		}
		if(d->header_row[irow]) free(d->header_row[irow]);
		d->header_row[irow] = copy_string(symbol);
		//printf("D %d %d\n",i,irow);
	}else{
		if(!symbol) d->header_row[irow] = copy_string(name);
		else{
			strcpy(buffer,name);
			strcat(buffer," (");
			strcat(buffer,symbol);
			strcat(buffer,")");
			d->header_row[irow] = copy_string(buffer);
		}
	}
	d->F[irow] = F;
	d->MSE[irow] = atof(items[ncol+3+2]);
	d->FDR[irow] = atof(items[ncol+6+2]);
	for(j=0; j<ncol; ++j){
		if(sscanf(items[j+2],"%f", &x) != 1){
			printf("Error reading: row %d, col %d\n", i+1, j+3);
			x = MISSING;
		}
		data_row[j] = x;
	}
	if(!backgr){
		baseline = median1(data_row,ncol);
	}else{
		baseline = data_row[backgr-1];
	}
	for(j=0; j<ncol; ++j){
		if(data_row[j]>MISSING && baseline>MISSING && (data_row[j]>=exprThresh || baseline>=exprThresh)){
			d->data[j][irow] = data_row[j]-baseline;
		}else{
			d->data[j][irow] = MISSING;
		}
	}
}
//printf("D1 %d\n",nrow);
if(!redundant){ d->nrow=nrow; }
fclose(inputFile);
destroy_hash(symbol2row);
free(data_row);
free(buffer);
free(items);
return(d);
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
float *get_distance_matrix (DATA *d)
/***********************************************/
{
float *dist, sum, x, y;
int i,j,irow,nn;

check(dist = (float*)calloc(d->ncol*d->ncol,sizeof(float)));
for(i=0; i<d->ncol; i++){
	for(j=i+1; j<d->ncol; j++){
		sum=0;
		nn=0;
		for(irow=0; irow<d->nrow; irow++){
			if(d->FDR[irow] > 0.05) continue;
			x = d->data[i][irow];
			y = d->data[j][irow];
			if(x>MISSING && y>MISSING){
				sum += (x-y)*(x-y);
				nn++;
			}
		}
		dist[i*d->ncol+j] = sqrt(sum/nn);
		dist[j*d->ncol+i] = dist[i*d->ncol+j];
	}
}
return(dist);
}

/***********************************************/
void  fill_output_table (DATA *d, long icol, long backgr, float **zvalue, float **FDR, 
	float **pvalue, float minMSE)
/***********************************************/
{
float x, y, z, x1, MSE;
float *zcopy, coeff, FDR1;
long *index, nn=0, n;
long i,j,irow;

check(zvalue[icol] = (float*)calloc(d->nrow,sizeof(float)));
check(FDR[icol] = (float*)calloc(d->nrow,sizeof(float)));
check(pvalue[icol] = (float*)calloc(d->nrow,sizeof(float)));
check(zcopy = (float*)calloc(d->nrow,sizeof(float)));
check(index = (long*)calloc(d->nrow,sizeof(long)));
for(irow=0; irow<d->nrow; irow++){
	index[irow]=irow;
	x = d->data[icol][irow];
	if(x==MISSING) continue;
	nn++;
	coeff = 1.0/d->nrepl[icol];
	if(backgr>0) coeff += 1.0/d->nrepl[backgr-1];
	else coeff += 1.0/d->ncol;
	MSE = d->MSE[irow];
	if(MSE <= 0) MSE=minMSE;
	zvalue[icol][irow] = fabs(x)/sqrt(MSE*coeff);
}
memcpy(zcopy,zvalue[icol],d->nrow*sizeof(float));
sortem(d->nrow,zcopy,index);
FDR1=1;
for(i=d->nrow-nn; i<d->nrow; i++){
	irow = index[i];
	z = zvalue[icol][irow];
	if(z==0){
		pvalue[icol][irow] = 1;
		FDR[icol][irow] = 1;
	}else{
		pvalue[icol][irow] = 2*(1 - normal_distribution(z));
		FDR[icol][irow] = pvalue[icol][irow]*(float)nn/(float)(d->nrow-i);
		if(FDR[icol][irow] > FDR1){ FDR[icol][irow] = FDR1; }
		else{ FDR1 = FDR[icol][irow]; }
	}
}
free(zcopy);
free(index);
}

/***********************************************/
void print_output_table (char *outputFile, DATA *d, long testcol, long backgr, float **zvalue, 
	float **FDR, float **pvalue, char *sourceFile)
/***********************************************/
{
FILE *fp;
long i, j;
float x, log10;
char buffer[12];

log10 = log(10.0);
fp = fopen(outputFile,"w");
fprintf(fp, "!Output_name\tPairwise comparison of gene expression\n");
if(sourceFile){
	fprintf(fp, "!Output_data_set1\t%s\n",sourceFile);
	fprintf(fp, "!Output_data_type1\tmatrix\n");
}
fprintf(fp, "!Output_method\tComparioson of expression profiles, ANOVA\n");
if(backgr<=0) fprintf(fp, "!Output_baseline\tmedian\n");
else          fprintf(fp, "!Output_baseline\t%s\n",d->header_col[backgr]);
if(testcol) fprintf(fp, "!Output_testcolumn\t%s\n",d->header_col[testcol-1]);
fprintf(fp, "!Output_N_columns\t%d\n",d->ncol);
fprintf(fp, "!Output_N_rows\t%d\n",d->nrow);
fprintf(fp, "!Output_N_matrixes\t5\n");
fprintf(fp, "!Output_matrix1_name\tz-value\n");
fprintf(fp, "!Output_matrix2_name\tlog-ratio\n");
fprintf(fp, "!Output_matrix3_name\tfold_change\n");
fprintf(fp, "!Output_matrix4_name\tFDR\n");
fprintf(fp, "!Output_matrix5_name\tp-value\n");
fprintf(fp, "!Output_matrix1_type\tnumbers\n");
fprintf(fp, "!Output_matrix2_type\tnumbers\n");
fprintf(fp, "!Output_matrix3_type\tnumbers\n");
fprintf(fp, "!Output_matrix4_type\tnumbers\n");
fprintf(fp, "!Output_matrix5_type\tnumbers\n");
fprintf(fp,"!Matrix1_start\n");
fprintf(fp,"Headers");
for(i=0; i<d->ncol; ++i){
	fprintf(fp,"\t%s",d->header_col[i]);
}
fprintf(fp,"\n");
for(i=0; i<d->nrow; ++i){
	fprintf(fp,"%s",d->header_row[i]);
	for(j=0; j<d->ncol; ++j){
		x = zvalue[j][i];
		if(d->data[j][i] != MISSING && d->data[j][i]<0) x = -x;
		fprintf(fp,"\t%.4f",x);
	}
	fprintf(fp,"\n");
}
fprintf(fp,"!Matrix1_end\n");
fprintf(fp,"!Matrix2_start\n");
fprintf(fp,"Headers");
for(i=0; i<d->ncol; ++i){
	fprintf(fp,"\t%s",d->header_col[i]);
}
fprintf(fp,"\n");
for(i=0; i<d->nrow; ++i){
	fprintf(fp,"%s",d->header_row[i]);
	for(j=0; j<d->ncol; ++j){
		x = d->data[j][i];
		if(x==MISSING) x=0;
		fprintf(fp,"\t%.4f",x);
	}
	fprintf(fp,"\n");
}
fprintf(fp,"!Matrix2_end\n");
fprintf(fp,"!Matrix3_start\n");
fprintf(fp,"Headers");
for(i=0; i<d->ncol; ++i){
	fprintf(fp,"\t%s",d->header_col[i]);
}
fprintf(fp,"\n");
for(i=0; i<d->nrow; ++i){
	fprintf(fp,"%s",d->header_row[i]);
	for(j=0; j<d->ncol; ++j){
		x = d->data[j][i];
		if(x==MISSING) x=0;
		x = exp(x*log10);
		if(x<0.01) fprintf(fp,"\t%.6f",x);
		else fprintf(fp,"\t%.4f",x);
	}
	fprintf(fp,"\n");
}
fprintf(fp,"!Matrix3_end\n");
fprintf(fp,"!Matrix4_start\n");
fprintf(fp,"Headers");
for(i=0; i<d->ncol; ++i){
	fprintf(fp,"\t%s",d->header_col[i]);
}
fprintf(fp,"\n");
for(i=0; i<d->nrow; ++i){
	fprintf(fp,"%s",d->header_row[i]);
	for(j=0; j<d->ncol; ++j){
		if(FDR[j][i]<1.0e-15) fprintf(fp,"\t0");
		else if(FDR[j][i]<1.0e-5){
			sprintf(buffer,"%.3e",FDR[j][i]);
			if(strlen(buffer)==10){
				buffer[7]=buffer[8]; buffer[8]=buffer[9]; buffer[9]=buffer[10]; 
			}
			fprintf(fp,"\t%s",buffer);
		}
		else if(FDR[j][i]<0.01) fprintf(fp,"\t%.6f",FDR[j][i]);
		else fprintf(fp,"\t%.4f",FDR[j][i]);
	}
	fprintf(fp,"\n");
}
fprintf(fp,"!Matrix4_end\n");
fprintf(fp,"!Matrix5_start\n");
fprintf(fp,"Headers");
for(i=0; i<d->ncol; ++i){
	fprintf(fp,"\t%s",d->header_col[i]);
}
fprintf(fp,"\n");
for(i=0; i<d->nrow; ++i){
	fprintf(fp,"%s",d->header_row[i]);
	for(j=0; j<d->ncol; ++j){
		if(pvalue[j][i]<1.0e-15) fprintf(fp,"\t0");
		else if(pvalue[j][i]<1.0e-5){
			sprintf(buffer,"%.3e",pvalue[j][i]);
			if(strlen(buffer)==10){
				buffer[7]=buffer[8]; buffer[8]=buffer[9]; buffer[9]=buffer[10]; 
			}
			fprintf(fp,"\t%s",buffer);
		}
		else if(pvalue[j][i]<0.01) fprintf(fp,"\t%.6f",pvalue[j][i]);
		else fprintf(fp,"\t%.4f",pvalue[j][i]);
	}
	fprintf(fp,"\n");
}
fprintf(fp,"!Matrix5_end\n");
fclose(fp);
return;
}

/***********************************************/
long  get_significant (DATA *d, float *dist, long icol, long backgr, float FDRthresh, float logratio_thresh, 
	float specific, FILE *fout, float minMSE)
/***********************************************/
{
float maxdist = 0, x, y, z, x1, MSE, sqrtMSE, m, sd;
float *zvalue, *zcopy, *FDRlist, FDR, FDR1=1;
float *zspecif, p, logratio, coeff;
long *index, *set_up, *set_dwn, nn=0, n, nset_up=0, nset_dwn=0;
long i,j,irow,n_geneset=0;

check(zvalue = (float*)calloc(d->nrow,sizeof(float)));
check(zcopy = (float*)calloc(d->nrow,sizeof(float)));
check(zspecif = (float*)calloc(d->nrow,sizeof(float)));
check(FDRlist = (float*)calloc(d->nrow,sizeof(float)));
check(index = (long*)calloc(d->nrow,sizeof(long)));
check(set_up = (long*)calloc(d->nrow,sizeof(long)));
check(set_dwn = (long*)calloc(d->nrow,sizeof(long)));
for(i=0; i<d->ncol; i++){
	for(j=i+1; j<d->ncol; j++){
		if(maxdist<dist[i*d->ncol+j]) maxdist=dist[i*d->ncol+j];
	}
}
for(irow=0; irow<d->nrow; irow++){
	index[irow]=irow;
	MSE = d->MSE[irow];
	sqrtMSE = sqrt(MSE);
	x = d->data[icol][irow];
	if(x==MISSING) continue;
	nn++;
	coeff = 1.0/d->nrepl[icol];
	if(backgr>0) coeff += 1.0/d->nrepl[backgr-1];
	else coeff += 1.0/d->ncol;
	if(MSE <= 0) MSE=minMSE;
	zvalue[irow] = fabs(x)/sqrt(MSE*coeff);
	m=0; sd=0; n=0;
	for(i=0; i<d->ncol; i++){
		if(i==icol || dist[i*d->ncol+icol] < maxdist/3) continue;
		x1 = d->data[i][irow];
		if(x1==MISSING) continue;
		if(x1<=MISSING) continue;
		m += x1;
		sd += x1*x1;
		n++;
	}
	if(n>1 && x*(x*n-m) > 0){
		sd = sqrt(fabs(sd-m*m/n)/(n-1));
		if(sd < sqrtMSE) sd =sqrtMSE;
		m /= n;
		if(sd>0) zspecif[irow] = floor(1000*fabs(x-m)/sd+0.5)/1000;
	}
}
memcpy(zcopy,zvalue,d->nrow*sizeof(float));
sortem(d->nrow,zcopy,index);
for(i=d->nrow-nn; i<d->nrow; i++){
	irow = index[i];
	z = zvalue[irow];
	if(!z) continue;
	p = 2*(1 - normal_distribution(z));
	FDR = p*(float)nn/(d->nrow-i);
	if(FDR > FDR1){ FDR = FDR1; }
	else{ FDR1 = FDR; }
	FDRlist[irow] = FDR;
	logratio = d->data[icol][irow];
	if(FDR <= FDRthresh && fabs(logratio) >= logratio_thresh && (d->ncol<3 || zspecif[irow] >= specific)){
		if(logratio > 0){
			set_up[nset_up++] = irow;
		}else{
			set_dwn[nset_dwn++] = irow;
		}
	}
}
if(nset_up){
	n_geneset++;
	fprintf(fout,"%s_up\t",d->header_col[icol]);
	for(i=nset_up-1; i>=0; i--)
		fprintf(fout,"\t%s",d->header_row[set_up[i]]);
	fprintf(fout,"\n\tFDR");
	for(i=nset_up-1; i>=0; i--)
		fprintf(fout,"\t%.4f",FDRlist[set_up[i]]);
	fprintf(fout,"\n\tlogratio");
	for(i=nset_up-1; i>=0; i--)
		fprintf(fout,"\t%.4f",d->data[icol][set_up[i]]);
	fprintf(fout,"\n\tspecific");
	for(i=nset_up-1; i>=0; i--)
		fprintf(fout,"\t%.4f",zspecif[set_up[i]]);
	fprintf(fout,"\n");
}
if(nset_dwn){
	n_geneset++;
	fprintf(fout,"%s_down\t",d->header_col[icol]);
	for(i=nset_dwn-1; i>=0; i--)
		fprintf(fout,"\t%s",d->header_row[set_dwn[i]]);
	fprintf(fout,"\n\tFDR");
	for(i=nset_dwn-1; i>=0; i--)
		fprintf(fout,"\t%.4f",FDRlist[set_dwn[i]]);
	fprintf(fout,"\n\tlogratio");
	for(i=nset_dwn-1; i>=0; i--)
		fprintf(fout,"\t%.4f",d->data[icol][set_dwn[i]]);
	fprintf(fout,"\n\tspecific");
	for(i=nset_dwn-1; i>=0; i--)
		fprintf(fout,"\t%.4f",zspecif[set_dwn[i]]);
	fprintf(fout,"\n");
}
free(zvalue);
free(zcopy);
free(zspecif);
free(FDRlist);
free(index);
free(set_up);
free(set_dwn);
return(n_geneset);
}

/***********************************************/
int main (int argc, char **argv)
/***********************************************/
{
DATA *d;
FILE *fout;
long ncol, nrow, backgr, testcol=0, i,j, useSymbols=0,redundant=0,delta=1;
float FDRthresh=0.05, foldchange=2, specific=0, minMSE=1.0e15, exprThresh=MISSING;
float *DIST, **FDR, **pvalue, **zvalue;
long n_geneset=0;
char *filemode = "w", *logFile=NULL, *outputFile=NULL, *command, *sourceFile=NULL;
if (argc < 4){
	error_message("pairwise inputFile outputFile backgr_column [-col testcolumn]\n");
}

backgr = atoi(argv[3]);
for(i=4; i<argc; i++){
	if(!strcmp(argv[i],"-col") && i<argc-1) testcol = atoi(argv[++i]);
	else if(!strcmp(argv[i],"-fdr") && i<argc-1) sscanf(argv[++i],"%f",&FDRthresh);
	else if(!strcmp(argv[i],"-fold") && i<argc-1) sscanf(argv[++i],"%f",&foldchange);
	else if(!strcmp(argv[i],"-spec") && i<argc-1) sscanf(argv[++i],"%f",&specific);
	else if(!strcmp(argv[i],"-expr") && i<argc-1) sscanf(argv[++i],"%f",&exprThresh);
	else if(!strcmp(argv[i],"-log") && i<argc-1) logFile = copy_string(argv[++i]);
	else if(!strcmp(argv[i],"-out") && i<argc-1) outputFile = copy_string(argv[++i]);
	else if(!strcmp(argv[i],"-source") && i<argc-1) sourceFile = copy_string(argv[++i]);
	else if(!strcmp(argv[i],"-symbol")) useSymbols=1;
	else if(!strcmp(argv[i],"-add")) filemode = "a";
	else if(!strcmp(argv[i],"-red")) redundant = 1;
	else error_message("Wrong option");
}
if(foldchange<1.0e-15) error_message("Fold change is zero or negative");
foldchange = log(foldchange)/log(10);
if(exprThresh>1.0e-15){ exprThresh = log(exprThresh)/log(10); }

get_dimensions(argv[1], &nrow, &ncol);
if(!nrow || !ncol){ error_message("Data undefined\n"); }
if(ncol>200) delta=20;
else if(ncol>100) delta=10;
else if(ncol>30) delta=5;
else if(ncol>15) delta=2;
check(FDR = (float**)malloc(ncol*sizeof(float*)));
check(pvalue = (float**)malloc(ncol*sizeof(float*)));
check(zvalue = (float**)malloc(ncol*sizeof(float*)));
check(command = (char*)malloc(sizeof(char)*10000));
if(nrow<10 || ncol<1) error_message("Wrong table dimensions");
d = read_data(argv[1],ncol,nrow,backgr,useSymbols,redundant,&minMSE,exprThresh);
if(minMSE==1.0e15) minMSE = 0.00025;
else minMSE *= 0.25;
DIST = get_distance_matrix(d);
fout = fopen(argv[2],filemode);
for(i=0; i<ncol; i++){
	if(!outputFile && (i==backgr-1 || testcol>0 && i!=testcol-1)) continue;
	if(logFile && i%delta==0){
		sprintf(command,"echo \"Column %d (N=%d)\" >> %s",i,ncol,logFile);
		system(command);
	}
	if(i != backgr-1 && (testcol==0 || i==testcol-1)){
		n_geneset += get_significant(d,DIST,i,backgr,FDRthresh,foldchange,specific,fout,minMSE);
	}
	if(outputFile){
		fill_output_table(d,i,backgr,zvalue,FDR,pvalue,minMSE);
	}
}
free(DIST);
if(outputFile){
	print_output_table(outputFile,d,testcol,backgr,zvalue,FDR,pvalue,sourceFile);
}
printf("%d",n_geneset);
return(0);
}


