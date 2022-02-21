/* correlation_exatlas input.txt output.txt */

//./page_exatlas_test ../output/5097.txt output.txt 79 23608 32 -cut 0.15 -attrib -epfp 0.5 -fold 1.5 -genes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define EPS 3.0e-7
#define FPMIN 1.0e-30
#define MAXCOL  50000
#define MISSING -9999
#define STRINGLEN 50

static long *blank;

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
L19:	tb = b[p];
	b[p] = b[i];
L21: /* Start at the beginning of the segment, search for k such that a(k)>t */
    q = j;
    k = i;
L20:	++k;
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
L44:	xb = b[k];
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
L70:	lt[m - 1] = q + 1;
	ut[m - 1] = j;
	j = q - 1;
/* Update m and split the new smaller segment */
L80:	++m;
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
L110:	return;
} /* sortem_ */

/***********************************************/
void error_message (char *message)
/***********************************************/
{
printf("Error: %s\n", message);
exit(0);
}

/***********************************************/
void check (void *x)
/***********************************************/
{
if (!x){ error_message("Out of memory\n"); }
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

/**  Allocation of float matrix storage  *****************************/
float **new_matrix(long n, long m)
/**********************************************************/
{
    long i;
    float **mat;

    check(mat = (float **) malloc((unsigned)n*sizeof(float*)));
    for (i = 0; i < n; i++){
        check(mat[i] = (float *) calloc((unsigned)m,sizeof(float)));
    }
    return mat;
}

/**  Allocation of long matrix storage  *****************************/
long **new_matrix_long(long n, long m)
/**********************************************************/
{
    long i;
    long **mat;

    check(mat = (long **) malloc(n*sizeof(long*)));
    for (i = 0; i < n; i++){
        check(mat[i] = (long *) calloc(m,sizeof(long)));
    }
    return mat;
}

/*******************************************/
float **read_data (FILE *fp, int ncol, long nrows, char **headers, char **symbols)
/*******************************************/
{
char *buffer, **items, *pchar;
long buff_size=500000, i, j, len, nitems, iline=0;
float x, **a;

a = new_matrix(ncol,nrows);
check(buffer = (char*)malloc(buff_size*sizeof(char)));
check(items = (char**)malloc(MAXCOL*sizeof(char*)));
fgets(buffer,buff_size-1,fp);
nitems = split_string(buffer,items,MAXCOL);
for(i=0; i<ncol; ++i){
	if(i>=nitems-1) break;
	headers[i] = copy_string(items[i+1]);
}
while(fgets(buffer,buff_size-1,fp)){
	len = strlen(buffer);
	if(buffer[len-1]=='\n'){ buffer[len-1]='\0'; --len; }
	if(!len){ break; }
	nitems = split_string (buffer,items,buff_size-1);
	symbols[iline] = copy_string(items[0]);
	for(i=0; i<ncol; ++i){
		if(i>=nitems-1) break;
		if(sscanf(items[i+1],"%f",&x)){
			a[i][iline] = x;
		}else{
			a[i][iline] = MISSING;
		}
	}
	iline++;
}
if(iline != nrows) error_message("Wrong number of lines");
free(items);
free(buffer);
return(a);
}

/*******************************************/
long **sort_matrix (float **a, int ncol, long nrows)
/*******************************************/
{
long **a_sorted, *index, *index1, i, j;
char *genelist;
float zzz, *response;

check(response = (float*)malloc(sizeof(float)*nrows));
check(index = (long*)malloc(sizeof(long)*nrows));
check(index1 = (long*)malloc(sizeof(long)*nrows));
for(j=0; j<nrows; ++j){
	index1[j] = j;
}
a_sorted = new_matrix_long(ncol,nrows);
for(i=0; i<ncol; ++i){
	memcpy(index,index1,nrows*sizeof(long));
	for(j=0; j<nrows; ++j){
		response[j] = a[i][j];
	}
	sortem(nrows, response, index);
	memcpy(a_sorted[i],index,nrows*sizeof(long));
}
return(a_sorted);
}

/*******************************************/
void  read_geneset (FILE *fp, long *b, int iGS, long nrows, char **geneset_headers, long *Ngenes)
/*******************************************/
{
char *buffer, **items, *pchar;
long buff_size=1000000, i, j, len, nitems, iline=0;
long x;

check(buffer = (char*)malloc(buff_size*sizeof(char)));
check(items = (char**)malloc(MAXCOL*sizeof(char*)));
fgets(buffer,buff_size-1,fp);
nitems = split_string(buffer,items,MAXCOL);
geneset_headers[iGS] = copy_string(items[0]);
memcpy(b,blank,nrows*sizeof(long));
for(i=1; i<nitems; ++i){
	x = atoi(items[i]);
	if(x >= nrows || x<0) error_message("Invalid gene number");
	b[x] = i;
}
*Ngenes = nitems-1;
free(items);
free(buffer);
return;
}

/*******************************************************/
char  *get_coregulated_genes (char **symbols, float *x, long *xsorted, long *member, 
	float propControl, long nrows, long nSet, long idir, float fold_thresh, float EPFP_thresh, float cutoff)
/*******************************************************/
{
float EPFP, EPFP1;
float *EPFPlist, resp;
long i, j, n, k, dir;
long *genelist, alloc, len1, len2, nTargets, nsubset;
char *output=NULL, *buffer, number[20];

if(nSet < 3) return(output);
check(EPFPlist = (float*)malloc(sizeof(float)*nSet));
check(genelist = (long*)malloc(sizeof(long)*nSet));
nsubset = nrows*cutoff;
dir = idir*2-1;
nTargets=0;
n=0;
for(i=0; i<nsubset; ++i){
	j = i;
	if(idir){ j=nrows-1-i; }
	k = xsorted[j];
	resp = x[k];
	if(resp==MISSING){ continue; }
	if(resp*dir<0.0792 || resp*dir<fold_thresh){ break; }
	n++;
	if(member[k]){
		//printf("T %d %d %f %f\n",j,k,resp,propControl);
		EPFPlist[nTargets] = propControl/(nTargets+0.5)*n;
		//printf("T %d %f\n",j,EPFPlist[nTargets]);
		genelist[nTargets++] = k;
		if(nTargets >= nSet) break;
	}
}
EPFP1=1;
n = 0;
for(i=nTargets-1; i>=0; i--){
	if(EPFPlist[i] > EPFP1){
		EPFPlist[i] = EPFP1;
	}else{
		EPFP1 = EPFPlist[i];
	}
	if(n==0 && EPFPlist[i] <= EPFP_thresh){ n=i+1; }
}
nTargets = n;
alloc = 100000;
if(nTargets>0){
	check(output = (char*)malloc(sizeof(char)*alloc));
	check(buffer = (char*)malloc(sizeof(char)*alloc));
	strcpy(output,"genes\t");
	strcpy(buffer,"\nEPFP\t");
	len1 = strlen(output);
	len2 = strlen(buffer);
	for(i=0; i<nTargets; i++){
		strcat(output,symbols[genelist[i]]);
		strcat(output,",");
		sprintf(number,"%.4f,",EPFPlist[i]);
		strcat(buffer,number);
		len1 = strlen(output);
		len2 = strlen(buffer);
		if(len1>alloc-100 || len2>alloc-100 ){
			alloc += 100000;
			check(output = (char*)realloc(output,sizeof(char)*alloc));
			check(buffer = (char*)realloc(buffer,sizeof(char)*alloc));
		}
	}
	output[len1-1]='\0';
	buffer[len2-1]='\0';
	if(len1+len2 > alloc){
		check(output = (char*)realloc(output,sizeof(char)*(len1+len2+10)));
	}
	strcat(output,buffer);
	free(buffer);
}
free(EPFPlist);
free(genelist);
return (output);
}

/*******************************************/
long  page_analysis (float **a, long **asorted, long *b, float **z_up, float **z_dn, char **database, 
	long **id_up, long **id_dn, long ndata, int ncol, int nrows, long iGS,
	char **symbols, int option_genes, float EPFP_thresh,float fold_thresh, float cutoff)
/*******************************************/
{
long i, j1, j, k, n,nSet, nsubset, idir, dir, nControl, nSetControl;
char *genelist;
float zzz, *propControl, resp;
float mean,sd,meanset;

nsubset = nrows*cutoff;
check(propControl = (float*)malloc(sizeof(float)*ncol));
for(i=0; i<ncol && option_genes; ++i){
	nControl=0;
	nSetControl =0;
	for(j=0; j<nrows; ++j){
		k = asorted[i][j];
		resp = a[i][k];
		if(resp < -0.0792){ continue; }
		if(resp > 0.0792){ break; }
		nControl++;
		if(b[k]){
			nSetControl++;
		}
	}
	propControl[i] = (float)nSetControl/nControl;
}
for(idir=0; idir<2; ++idir){
	dir = idir*2-1;
	for(i=0; i<ncol; ++i){
		n=0; nSet=0; mean=0; sd=0; meanset=0;
		for(j=0; j<nsubset; ++j){
			j1 = j;
			if(idir){ j1=nrows-1-j; }
			k = asorted[i][j1];
			resp = a[i][k];
			if(resp==MISSING){ continue; }
			if(resp*dir<=0){ break; }
			n++;
			mean += resp;
			sd += resp*resp;
			if(b[k]){
				//printf("%d %d %d %d %s\n",iGS,idir,j,k,symbols[k]);
				nSet++;
				meanset += resp;
			}
		}
		zzz = 0;
		if(nSet >= 5){
			sd = sd-mean*mean/n;
			if(sd<EPS) sd = EPS;
			sd = sqrt(sd/(n-1));
			mean /= n;
			meanset /= nSet;
			if(idir==0){
				mean *= -1;
				meanset *= -1;
			}
			zzz = (meanset - mean)*sqrt((float)nSet)/sd;
		}
		if(idir==0){
			z_dn[i][iGS] = zzz;
			if(zzz >=2 && option_genes){
				genelist = get_coregulated_genes(symbols,a[i],asorted[i],b,propControl[i],nrows,nSet,idir,fold_thresh,EPFP_thresh,cutoff);
				if(genelist){
					id_dn[i][iGS] = ndata+1;
					database[ndata++] = genelist;
				}
			}
		}else{
			z_up[i][iGS] = zzz;
			if(zzz >=2 && option_genes){
				genelist = get_coregulated_genes(symbols,a[i],asorted[i],b,propControl[i],nrows,nSet,idir,fold_thresh,EPFP_thresh,cutoff);
				if(genelist){
					id_up[i][iGS] = ndata+1;
					database[ndata++] = genelist;
				}
			}
		}
	}
}
free(propControl);
return(ndata);
}

/*******************************************************/
char  *get_coregulated_genes_attrib (char **symbols, float *x, long *xsorted, long *member, 
	float *propControl, long nrows, long idir, float fold_thresh, float EPFP_thresh, float cutoff, long nGenes)
/*******************************************************/
{
float EPFP, EPFP1;
float *EPFPlist, resp;
long i, j, n, k, k1, i1, dir;
long *genelist, alloc, len1, len2, nTargets, nTargetsAttr, nsubset, *stack;
char *output=NULL, *buffer, number[20];

check(EPFPlist = (float*)malloc(sizeof(float)*nGenes));
check(genelist = (long*)malloc(sizeof(long)*nGenes));
check(stack = (long*)malloc(sizeof(long)*nGenes));
nTargets=0;
nsubset = nrows*cutoff;
dir = idir*2-1;
n=0;
for(i=0; i<nsubset; ++i){
	j = i;
	if(idir){ j=nrows-1-i; }
	k = xsorted[j];
	resp = x[k];
	if(resp==MISSING){ continue; }
	if(resp*dir<0.0792 || resp*dir<fold_thresh){ break; }
	n++;
	k1 = member[k];
	if(k1){
		stack[nTargets] = k1;
		nTargetsAttr=1;
		for(i1=0; i1<nTargets; ++i1){
			if(stack[i1]<k1) nTargetsAttr++;
		}
		EPFPlist[nTargets] = propControl[k1-1]/(nTargetsAttr+0.5)*n;
		//printf("%d %d %d %d %f %f\n",i,k1,nTargetsAttr,n,propControl[k1-1],EPFPlist[nTargets]);
		for(i1=0; i1<nTargets; ++i1){
			if(stack[i1]<k1 && EPFPlist[i1]>EPFPlist[nTargets]){
				EPFPlist[i1]=EPFPlist[nTargets];
			}
		}
		genelist[nTargets++] = k;
	}
}
n = 0;
alloc = 100000;
if(nTargets>0){
	check(output = (char*)malloc(sizeof(char)*alloc));
	check(buffer = (char*)malloc(sizeof(char)*alloc));
	strcpy(output,"genes\t");
	strcpy(buffer,"\nEPFP\t");
	len1 = strlen(output);
	len2 = strlen(buffer);
	for(i=0; i<nTargets; i++){
		if(EPFPlist[i] > EPFP_thresh) continue;
		n++;
		strcat(output,symbols[genelist[i]]);
		strcat(output,",");
		sprintf(number,"%.4f,",EPFPlist[i]);
		strcat(buffer,number);
		len1 = strlen(output);
		len2 = strlen(buffer);
		//printf("A1 %d %d %d %f\n",i, len1,len2,EPFPlist[i]);
		if(len1>alloc-100 || len2>alloc-100 ){
			alloc += 100000;
			check(output = (char*)realloc(output,sizeof(char)*alloc));
			check(buffer = (char*)realloc(buffer,sizeof(char)*alloc));
		}
	}
	output[len1-1]='\0';
	buffer[len2-1]='\0';
	if(len1+len2 > alloc){
		check(output = (char*)realloc(output,sizeof(char)*(len1+len2+10)));
	}
	strcat(output,buffer);
	free(buffer);
	if(n==0){
		free(output);
		output = NULL;
	}
}
free(EPFPlist);
free(genelist);
free(stack);
return (output);
}

/*******************************************/
long  page_analysis_attribute (float **a, long **asorted, long *b, float **z_up, float **z_dn, char **database, 
	long **id_up, long **id_dn, long ndata, int ncol, int nrows, long iGS, long Ngenes,
	char **symbols, int option_genes, float EPFP_thresh,float fold_thresh, float cutoff)
/*******************************************/
{
long i, j1, j, k, n, *nSet, nsubset, idir, dir, nControl, *nSetControl, sum, imax, jmax, nSetMax;
long *blanc;
char *genelist;
float zzz, **propControl, resp, sumfloat, xmax, ymax, x;
float mean,sd,*meanset,zmax;

nsubset = nrows*cutoff;
check(propControl = (float**)malloc(sizeof(float*)*ncol));
check(nSetControl = (long*)calloc(Ngenes,sizeof(long)));
check(blanc = (long*)calloc(Ngenes,sizeof(long)));
for(i=0; i<ncol && option_genes; ++i){
	memcpy(nSetControl,blanc,Ngenes*sizeof(long));
	check(propControl[i] = (float*)malloc(sizeof(float)*Ngenes));
	nControl=0;
	for(j=0; j<nrows; ++j){
		k = asorted[i][j];
		resp = a[i][k];
		if(resp < -0.0792){ continue; }
		if(resp > 0.0792){ break; }
		nControl++;
		if(b[k]){
			nSetControl[b[k]-1]++;
		}
	}
	sum=0;
	for(j=0; j<Ngenes; ++j){
		sum += nSetControl[j];
		if(nControl) propControl[i][j] = (float)sum/nControl;
	}
}
for(idir=0; idir<2; ++idir){
	dir = idir*2-1;
	for(i=0; i<ncol; ++i){
		n=0; mean=0; sd=0; xmax=0; ymax=0; jmax=0;
		check(meanset = (float*)calloc(Ngenes,sizeof(float)));
		check(nSet = (long*)calloc(Ngenes,sizeof(long)));
		for(j=0; j<nsubset; ++j){
			j1 = j;
			if(idir) j1=nrows-1-j;
			k = asorted[i][j1];
			resp = a[i][k];
			if(resp==MISSING) continue;
			if(resp*dir<=0) break;
			n++;
			mean += resp;
			sd += resp*resp;
			if(b[k]){
				j1 = b[k]-1;
				nSet[j1]++;
				meanset[j1] += resp;
				x = resp*dir/b[k];
				if(xmax<x){ xmax=x; ymax=resp; jmax=j1; }
			}
		}
		sum=0;
		sumfloat=0;
		sd = sd-mean*mean/n;
		if(sd<EPS) sd = EPS;
		sd = sqrt(sd/(n-1));
		mean /= n;
		zmax = 0;
		imax = -1;
		nSetMax = 0;
		for(j=0; j<Ngenes; ++j){
			sum += nSet[j];
			sumfloat += meanset[j];
			zzz = 0;
			if(sum >= 5){
				if(j>=jmax){
					meanset[j] = (sumfloat-ymax)/(sum-1);
				}else{
					meanset[j] = sumfloat/sum;
				}
				zzz = dir*(meanset[j] - mean)*sqrt((float)sum)/sd;
			}
			if(zmax < zzz){ zmax=zzz; imax=j; nSetMax=sum; }
		}
		if(idir==0){
			z_dn[i][iGS] = zmax;
			if(zmax >=2 && option_genes){
				genelist = get_coregulated_genes_attrib(symbols,a[i],asorted[i],b,propControl[i],nrows,idir,fold_thresh,EPFP_thresh,cutoff,Ngenes);
				if(genelist){
					id_dn[i][iGS] = ndata+1;
					database[ndata++] = genelist;
				}
			}
		}else{
			//if(i<3) printf("A1 %d %f\n",i,zmax);
			z_up[i][iGS] = zmax;
			if(zmax >=2 && option_genes){
				genelist = get_coregulated_genes_attrib(symbols,a[i],asorted[i],b,propControl[i],nrows,idir,fold_thresh,EPFP_thresh,cutoff,Ngenes);
				if(genelist){
					id_up[i][iGS] = ndata+1;
					database[ndata++] = genelist;
				}
			}
		}
		free(meanset);
		free(nSet);
	}
}
for(i=0; i<ncol && option_genes; ++i){
	free(propControl[i]);
}
free(propControl);
free(nSetControl);
free(blanc);
return(ndata);
}

/***********************************************/
int main (long argc, char **argv) 
/***********************************************/
{
int ncol, nGeneSet, option_genes=0, use_attribute=0, geneset_delta=50;
long nrows, **id_up, **id_dn, **a_sorted, ndata=0, *b, i, j, Ngenes;
char option, *logfile=NULL, *input_file, *output_file;
char **headers, **symbols, **geneset_names, **database, *command=NULL;
float **a, **z_up, **z_dn, fold_thresh=2, cutoff=0.25, EPFP_thresh=0.5, x, y;
FILE *fpout, *fp;

if (argc < 6){
	printf("Syntax: page_exatlas input_file output_file ncol nrows nGeneSet [-genes, -attrib, -l logfile, -epfp EPFP, -cut cutoff]\n");
	printf(" Input file should be tab-delimited\n");
	printf(" -genes: identify associated genes\n");
	printf(" -attrib: genes in geneset are sorted by some attribute (e.g., score)\n");
	printf(" -epfp:  Expected Proportion of False Positives threshold (default 0.5)\n");
	printf(" -fold:  Fold threshold for associated genes (default 2)\n");
	printf(" -cut:   Cutoff proportion of genes used for analysis (default 0.25)\n");
	exit(1);
}
input_file = copy_string(argv[1]);
output_file = copy_string(argv[2]);
ncol = atoi(argv[3]);
nrows = atoi(argv[4]);
nGeneSet = atoi(argv[5]);
if(ncol<1 || ncol>MAXCOL) error_message("Number of columns 1 - wrong value");
if(nrows<10) error_message("Not enough rows < 10");
for(i=6; i<argc; i++){
	if(!strcmp(argv[i],"-genes")) option_genes = 1;
	else if(!strcmp(argv[i],"-attrib")) use_attribute = 1;
	else if(!strcmp(argv[i],"-fold")){
		if(i+1==argc) error_message("fold change not specified");
		if(!sscanf(argv[++i],"%f",&fold_thresh)) error_message("fold change non-readable");
		if(fold_thresh<1 || fold_thresh>20) error_message("fold change should be from 1 to 20");
		fold_thresh = log(fold_thresh)/log(10.0);
	}
	else if(!strcmp(argv[i],"-epfp")){
		if(i+1==argc) error_message("EPFP threshold not specified");
		if(!sscanf(argv[++i],"%f",&EPFP_thresh)) error_message("EPFP threshold non-readable");
		if(EPFP_thresh<=0.0001 || EPFP_thresh>1) error_message("EPFP threshold should be from 0.0001 to 1");
	}
	else if(!strcmp(argv[i],"-cut")){
		if(i+1==argc) error_message("cutoff threshold not specified");
		if(!sscanf(argv[++i],"%f",&cutoff)) error_message("EPFP threshold non-readable");
		if(cutoff<=0.001 || cutoff>0.5) error_message("Cutoff should be from 0.001 to 0.5");
	}
	else if(!strcmp(argv[i],"-l")){
		if(i+1==argc) error_message("logfile name not specified");
		logfile = copy_string(argv[++i]);
	}
	else{
		error_message("unrecognized command-line parameter");
	}
}
a = new_matrix(ncol,nrows);
check(blank = (long*)calloc(nrows,sizeof(long)));
check(b = (long*)malloc(2*nrows*sizeof(long)));
check(database = (char**)malloc(sizeof(char*)*ncol*nGeneSet*2));
check(headers = (char**)calloc(ncol,sizeof(char*)));
check(geneset_names = (char**)calloc(nGeneSet,sizeof(char*)));
check(symbols = (char**)calloc(nrows,sizeof(char*)));
fp = fopen(input_file, "r");
a = read_data(fp,ncol,nrows,headers,symbols);
a_sorted = sort_matrix(a,ncol,nrows);
z_up = new_matrix(ncol,nGeneSet);
z_dn = new_matrix(ncol,nGeneSet);
id_up = new_matrix_long(ncol,nGeneSet);
id_dn = new_matrix_long(ncol,nGeneSet);
if(nGeneSet<10) geneset_delta=1;
else if (nGeneSet<100) geneset_delta=5;
else if (nGeneSet<500) geneset_delta=20;
for(i=0; i<nGeneSet; i++){
	if(logfile && i%geneset_delta==0){
		if(!command)
			check(command = malloc(sizeof(char)*(strlen(logfile)+200)));
		sprintf(command,"echo \"Geneset %d (N=%d)\" >> ",i+1,nGeneSet);
		strcat(command,logfile);
		system(command);
	}
	read_geneset(fp,b,i,nrows,geneset_names,&Ngenes);
	//printf("A2 %d %s\n",i,geneset_names[i]);
	if(use_attribute)
		ndata = page_analysis_attribute(a,a_sorted,b,z_up,z_dn,database,id_up,id_dn,ndata,ncol,nrows,i,Ngenes,symbols,option_genes,EPFP_thresh,fold_thresh,cutoff);
	else
		ndata = page_analysis(a,a_sorted,b,z_up,z_dn,database,id_up,id_dn,ndata,ncol,nrows,i,symbols,option_genes,EPFP_thresh,fold_thresh,cutoff);
}
fclose(fp);

fpout = fopen(output_file, "a");
if(!fpout) error_message("Output file not opened!");
fprintf(fpout,"!Matrix1_start\n");
fprintf(fpout,"Headers");
for(i=0; i<ncol; ++i){
	fprintf(fpout,"\t%s",headers[i]);
}
fprintf(fpout,"\n");
for(i=0; i<nGeneSet; ++i){
	fprintf(fpout,"%s",geneset_names[i]);
	for(j=0; j<ncol; ++j){
		fprintf(fpout,"\t%.4f",z_up[j][i]);
	}
	fprintf(fpout,"\n");
}
fprintf(fpout,"!Matrix1_end\n");
fprintf(fpout,"!Matrix2_start\n");
fprintf(fpout,"Headers");
for(i=0; i<ncol; ++i){
	fprintf(fpout,"\t%s",headers[i]);
}
fprintf(fpout,"\n");
for(i=0; i<nGeneSet; ++i){
	fprintf(fpout,"%s",geneset_names[i]);
	for(j=0; j<ncol; ++j){
		fprintf(fpout,"\t%.4f",z_dn[j][i]);
	}
	fprintf(fpout,"\n");
}
fprintf(fpout,"!Matrix2_end\n");
fprintf(fpout,"!Matrix3_start\n");
fprintf(fpout,"Headers");
for(i=0; i<ncol; ++i){
	fprintf(fpout,"\t%s",headers[i]);
}
fprintf(fpout,"\n");
for(i=0; i<nGeneSet; ++i){
	fprintf(fpout,"%s",geneset_names[i]);
	for(j=0; j<ncol; ++j){
		x = z_up[j][i];
		if(x<0) x=0;
		y = z_dn[j][i];
		if(y<0) y=0;
		fprintf(fpout,"\t%.4f",x-y);
	}
	fprintf(fpout,"\n");
}
fprintf(fpout,"!Matrix3_end\n");
if(!option_genes){ fclose(fpout); exit(0); }
fprintf(fpout,"!Matrix4_start\n");
fprintf(fpout,"Headers");
for(i=0; i<ncol; ++i){
	fprintf(fpout,"\t%s",headers[i]);
}
fprintf(fpout,"\n");
for(i=0; i<nGeneSet; ++i){
	fprintf(fpout,"%s",geneset_names[i]);
	for(j=0; j<ncol; ++j){
		fprintf(fpout,"\t%d",id_up[j][i]);
	}
	fprintf(fpout,"\n");
}
fprintf(fpout,"!Matrix4_end\n");
fprintf(fpout,"!Matrix5_start\n");
fprintf(fpout,"Headers");
for(i=0; i<ncol; ++i){
	fprintf(fpout,"\t%s",headers[i]);
}
fprintf(fpout,"\n");
for(i=0; i<nGeneSet; ++i){
	fprintf(fpout,"%s",geneset_names[i]);
	for(j=0; j<ncol; ++j){
		fprintf(fpout,"\t%d",id_dn[j][i]);
	}
	fprintf(fpout,"\n");
}
fprintf(fpout,"!Matrix5_end\n");
fprintf(fpout,"!Database_start\n");
for(i=0; i<ndata; i++){
	fprintf(fpout,">%d\n%s\n",i+1,database[i]);
}
fprintf(fpout,"!Database_end\n");
fclose(fpout);
return(0);
}

