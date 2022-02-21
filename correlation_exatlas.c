/* correlation_exatlas 2384.txt output.txt 62 91 14883 P -genes -epfp 0.5 -angle 0.5*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define EPS 3.0e-7
#define FPMIN 1.0e-30
#define MAXCOL  5000
#define MISSING -9999
#define STRINGLEN 50

void error_message (char *);
static char *logfile=NULL;

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
void check (void *x)
/***********************************************/
{
if (!x){
	error_message("Out of memory");
	exit(0);
}
}

/***********************************************/
void error_message (char *message)
/***********************************************/
{
char *command=NULL;
if(logfile){
	check(command = malloc(sizeof(char)*(strlen(logfile)+200)));
	sprintf(command,"echo \"Error: %s\" >> ",message);
	strcat(command,logfile);
	system(command);
}else{
	printf("Error: %s\n", message);
}
exit(0);
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

/**************************/
float   pearson_correlation (float *x1, float *y1, long n1, long *n, float *cov)
/**************************/
{
double sx=0, sxx=0, sxy=0, sy=0, syy=0, x, y, vx, vy, vxy;
float r=0;
long i;

*n = 0;
if(!x1 || !y1 || !n1){ return(0); }
for(i=0; i<n1; ++i){
	x = x1[i];
	y = y1[i];
	//if(i<50) printf("%f %f\n",x,y);
	if(x<=MISSING || y<=MISSING){ continue; }
	sx += x;
	sxx += x*x;
	sxy += x*y;
	sy += y;
	syy += y*y;
	++(*n);
}
if(*n<3){ return(0); }
vx = sxx-sx*sx/(*n);
vy = syy-sy*sy/(*n);
*cov = sxy-sx*sy/(*n);
vxy = vx*vy;
if(vxy<0.0000001) vxy=0.0000001;
r = *cov/sqrt(vxy);
return (r);
}

/**************************/
float   spearman_correlation (float *x1, float *y1, long n1, long *n)
/**************************/
{
double sx=0, sxx=0, sxy=0, sy=0, syy=0, x, y, vx, vy, vxy;
float r=0, *x2, *y2, *x3, *y3;
long i, j=0, *index, *index1;

if(!x || !y){ return(0); }
check(x2 = (float*)malloc(sizeof(float)*n1));
check(y2 = (float*)malloc(sizeof(float)*n1));
check(index = (long*)malloc(sizeof(long)*n1));
check(index1 = (long*)malloc(sizeof(long)*n1));
for(i=0; i<n1; ++i){
	if(x1[i] > MISSING && y1[i] > MISSING){
		x2[j] = x1[i];
		y2[j] = y1[i];
		index[j] = j;
		j++;
	}
}
*n = j;
memcpy(index1,index,j*sizeof(long));
sortem(j, x2, index);
sortem(j, y2, index1);
for(i=0; i<*n; ++i){
	x2[index[i]] = i;
	y2[index1[i]] = i;
}
free(index);
free(index1);
for(i=0; i<*n; ++i){
	x = x2[i];
	y = y2[i];
	sx += x;
	sxx += x*x;
	sxy += x*y;
	sy += y;
	syy += y*y;
}
free(x2);
free(y2);
if(*n<3){ return(0); }
vx = sxx-sx*sx/(*n);
vy = syy-sy*sy/(*n);
vxy = vx*vy;
if(vxy<0.0000001) vxy=0.0000001;
r = (sxy-sx*sy/(*n))/sqrt(vxy);
return (r);
}

/*******************************************/
void  read_data  (char *input_file, float **a, int ncol1, int ncol2, long nrows, 
		char **headers, char **row_names, char **symbols)
/*******************************************/
{
FILE *fp;
char *buffer, **items;
long buff_size=500000, maxcol, i, len, nitems, iline=0;
float x;

maxcol = ncol1+ncol2+50;
check(buffer = (char*)malloc(buff_size*sizeof(char)));
check(items = (char**)calloc(maxcol,sizeof(char*)));
fp = fopen(input_file, "r");
if(!fp) error_message("Input file not opened!");
fgets(buffer,buff_size-1,fp);
len = strlen(buffer);
if(buffer[len-1]=='\n'){ buffer[len-1]='\0'; --len; }
nitems = split_string(buffer,items,maxcol);
for(i=0; i<ncol1; ++i){
	if(i>=nitems-1) break;
	headers[i] = copy_string(items[i+1]);
	/* printf("%s\n",headers[i]); */
}
while(fgets(buffer,buff_size-1,fp)){
	len = strlen(buffer);
	if(buffer[len-1]=='\n'){ buffer[len-1]='\0'; --len; }
	if(!len){ break; }
	nitems = split_string (buffer,items,50000);
	symbols[iline] = copy_string(items[0]);
	for(i=0; i<ncol1; ++i){
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
fgets(buffer,buff_size-1,fp);
len = strlen(buffer);
if(buffer[len-1]=='\n'){ buffer[len-1]='\0'; --len; }
nitems = split_string (buffer,items,maxcol);
for(i=0; i<ncol2; ++i){
	if(i>=nitems) break;
	row_names[i] = copy_string(items[i+1]);
}
iline=0;
while(fgets(buffer,buff_size-1,fp)){
	len = strlen(buffer);
	if(buffer[len-1]=='\n'){ buffer[len-1]='\0'; --len; }
	if(!len){ break; }
	nitems = split_string (buffer,items,50000);
	if(strcmp(symbols[iline],items[0])) error_message("symbols do not match");
	for(i=0; i<ncol2; ++i){
		if(sscanf(items[i+1],"%f",&x)){
			a[i+ncol1][iline] = x;
		}else{
			a[i+ncol1][iline] = MISSING;
		}
	}
	/* printf("%d %s %f\n",iline,symbols[iline],a[ncol1][iline]); */
	iline++;
}
fclose(fp);
free(items);
free(buffer);
return;
}

/*******************************************************/
char  *get_coregulated_genes (char **symbols, float *x, float *y, int nrows, int sign,
	float angle, float EPFP_threshold)
/*******************************************************/
{
float x2, y2, x1, y1, z, wedgeArea, contArea;
float *xvalue, *yvalue, *xrank, *yrank, *score, *EPFPlist, *EPFP, *EPFP_profile;
long i, j, i1, j1, ix, iy, ix1, iy1, ix_last, n1=0, n2=0, count, *isymbol;
long *index, *index1, *genelist, *xlink, *ylink, *freqy, alloc, len1, len2;
char *output=NULL, *buffer, number[20];

check(xvalue = (float*)malloc(sizeof(float)*nrows));
check(yvalue = (float*)malloc(sizeof(float)*nrows));
check(index = (long*)malloc(sizeof(long)*nrows));
check(isymbol = (long*)malloc(sizeof(long)*nrows));
for(i=0; i<nrows; ++i){
	x1 = x[i];
	y1 = y[i];
	if(y1<=MISSING){ continue; }
	if(sign<0){ x1 = -x1; y1 = -y1; }
	if(x1<=0 || y1<=0){ continue; }
	xvalue[n1] = x1;
	yvalue[n1] = y1;
	isymbol[n1] = i;
	index[n1]=n1;
	n1++;
}
if(n1<5){
	free(xvalue);  free(yvalue); free(index); free(isymbol);
	return (NULL);
}
//printf("n1=%d %f\n",n1,angle);
check(index1 = (long*)malloc(sizeof(long)*n1));
check(xrank = (float*)malloc(sizeof(float)*n1));
check(yrank = (float*)malloc(sizeof(float)*n1));
check(xlink = (long*)malloc(sizeof(long)*n1));
check(ylink = (long*)malloc(sizeof(long)*n1));
check(score = (float*)malloc(sizeof(float)*n1));
check(EPFP = (float*)malloc(sizeof(float)*n1));
check(EPFPlist = (float*)malloc(sizeof(float)*n1));
check(freqy = (long*)calloc(n1,sizeof(long)));
check(genelist = (long*)malloc(sizeof(long)*n1));
memcpy(xrank,xvalue,n1*sizeof(float));
memcpy(yrank,yvalue,n1*sizeof(float));

/* Estimate rank */
memcpy(xlink,index,n1*sizeof(long));
sortem(n1, xrank, xlink);
for(i=0; i<n1; ++i){ xrank[xlink[i]]= i; }
memcpy(ylink,index,n1*sizeof(long));
sortem(n1, yrank, ylink);
for(i=0; i<n1; ++i){ yrank[ylink[i]]= i; }
if(angle>0){
	/* score = distance (vert. or horiz) from wedge border that goes through origin */
	for(i=0; i<n1; ++i){
		x1 = xrank[i]+0.5;
		y1 = yrank[i]+0.5;
		if(x1>y1){ score[i] = y1-(1.0-angle)*x1; freqy[(long)x1]++; }
		else     { score[i] = x1-(1.0-angle)*y1; freqy[(long)y1]++; }
	}
	count=0;
	for(i=0; i<n1; i++){
		count += freqy[i];
		freqy[i] = count;
	}
	sortem(n1,score,index);
	x2 = 1;
	for(i=n1/20; i<n1; ++i){
		j = index[i];
		if(x2 > EPFP_threshold+0.2 && (xrank[j]>yrank[j]*1.2 || yrank[j]>xrank[j]*1.2)) continue;
		if(x2 <= EPFP_threshold && (xrank[j]>yrank[j]*1.05 || yrank[j]>xrank[j]*1.05)){
			xlink[n2] = isymbol[j];
			EPFP[n2++] = x2;
			continue;
		}

		/* (n1-i) = actual num. of points in the wedge; z=hieght of wedge */
		if(score[i]>0){
			z = n1-score[i]/angle;
			wedgeArea = angle*z*z;
			contArea = 2*n1*z-wedgeArea-z*z;
			x1 = (i-freqy[n1-(long)z]+1)*wedgeArea/contArea/(n1-i);
		}else{
			z = n1-n1*angle+score[i];
			contArea = z*z/angle;
			wedgeArea = n1*n1-contArea;
			x1 = i*wedgeArea/contArea/(n1-i);
		}
		if(x1 > x2) x1 = x2;
		else x2 = x1;
		if(x1 <= EPFP_threshold){
			xlink[n2] = isymbol[j];
			EPFP[n2++] = x1;
		}
		if(x2 > EPFP_threshold+0.2){
			for(i1=i+30; i1<n1; ++i1){
				j = index[i1];
				if(xrank[j]<yrank[j]*1.2 && yrank[j]<xrank[j]*1.2){
					i = i1-1;
					break;
				}
			}
		}
	}
	for(i=0; i<n2; ++i){
		genelist[i] = xlink[n2-i-1];
		EPFPlist[i] = EPFP[n2-i-1];
	}
}else{
	check(EPFP_profile = (float*)malloc(sizeof(float)*n1));
	for(iy=0; iy<n1; iy++){
		EPFP_profile[iy] = 1;
		EPFP[iy] = 1;
	}
	for(ix=n1-1; ix>0; --ix){
		j = xlink[ix];
		iy = yrank[j];
		freqy[iy]++;
		count=0;
		for(i=n1-1; i>=iy; --i){
			count += freqy[i];
		}
		x1 = (float)(n1-iy-count)*(n1-ix)/ix/count;
		x2 = (float)(n1-ix-count)*(n1-iy)/iy/count;
		if(x1 < x2 && (n1-ix-count>=10||ix<iy+3)){ x1=x2; }
		EPFP[j] = x1;
		//printf("S %d %d %f\n",ix,iy,x1);
		if(iy > 0.9*n1 && count > 20 && x1 > EPFP_threshold*1.2){ break; }
		ix_last = ix;
	}
	for(ix=ix_last; ix<n1; ix++){
		j = xlink[ix];
		x1 = EPFP[j];
		iy = yrank[j];
		//if(j>=n1 || iy>=n1){ printf("Err\n"); exit(0); }
		if(x1 > EPFP_threshold && EPFP_profile[iy]>EPFP_threshold) continue;
		//printf("S %d %d %f %f\n",ix,iy,x1,EPFP_profile[iy]);
		if(EPFP_profile[iy] > x1){
			EPFP_profile[iy] = x1;
			for(iy1=iy+1; iy1<n1; iy1++){
				if(EPFP_profile[iy1] > EPFP_profile[iy]){ 
					EPFP_profile[iy1] = EPFP_profile[iy];
					//printf("B %d %d %d %f\n",ix,iy,iy1,EPFP_profile[iy]);
				}
				else{ break; }
			}
		}else if(EPFP_profile[iy] < x1){
			EPFP[j] = EPFP_profile[iy];
			//printf("T %d %f\n",iy,EPFP_profile[iy]);
		}
	}
	memcpy(EPFP_profile,EPFP,n1*sizeof(float));
	sortem(n1, EPFP_profile, index);
	for(i=0; i<n1; i++){
		j = index[i];
		if(EPFP[j] > EPFP_threshold) break;
		genelist[n2] = isymbol[j];
		EPFPlist[n2++] = EPFP[j];
	}
	free(EPFP_profile);
}
alloc = 100000;
if(n2>0){
	check(output = (char*)malloc(sizeof(char)*alloc));
	check(buffer = (char*)malloc(sizeof(char)*alloc));
	strcpy(output,"genes\t");
	strcpy(buffer,"\nEPFP\t");
	for(i=0; i<n2; i++){
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
free(xvalue);
free(yvalue);
free(index);
free(index1);
free(isymbol);
free(xrank);
free(yrank);
free(xlink);
free(ylink);
free(EPFP);
free(EPFPlist);
free(freqy);
free(score);
free(genelist);
return (output);
}

/*******************************************/
long  correlation_analysis (float **a, float **r, float **z, char **database, 
	long **id_up, long **id_dn, int ncol1, int ncol2, long nrows, char option, char **symbols, 
	int autocorrelation, int option_genes, float angle, float EPFP_threshold)
/*******************************************/
{
long data_id = 0, i, j, k, n1;
char *command=NULL, *genelist;
float cor, zzz, cov;

for(i=0; i<ncol1; ++i){
	if(logfile){
		if(!command)
			check(command = malloc(sizeof(char)*(strlen(logfile)+200)));
		sprintf(command,"echo \"Column %d (N=%d)\" >> ",i+1,ncol1);
		strcat(command,logfile);
		system(command);
	}
	for(j=0; j<ncol2; ++j){
		if(autocorrelation && j<=i){
			if(i==j){
				r[i][j] = 1;
				r[j][i] = 1;
				z[i][j] = 1000;
				z[j][i] = 1000;
			}else{
				r[i][j] = r[j][i];
				z[i][j] = z[j][i];
				id_up[i][j] = id_up[j][i];
				id_dn[i][j] = id_dn[j][i];
			}
			continue;
		}
		cor=0;
		cov=0;
		zzz=0;
		if(option == 'P' || option == 'C'){
			cor = pearson_correlation(a[i],a[ncol1+j],nrows,&n1,&cov);
		}else if(option == 'S'){
			cor = spearman_correlation(a[i],a[ncol1+j],nrows,&n1);
		}else{
			error_message("Unknown correlation method");
		}
		//printf("COR %d %d %d %f\n",i,j,n1,cor);
		if(n1 < 15){
			cor = 0;
		}
		if(option == 'P' || option == 'S'){
			r[i][j] = floor(10000*cor+0.5)/10000;
		}else{
			r[i][j] = floor(10000*cov+0.5)/10000;
		}
		if(n1>=15){
			if(cor>0.9999){ cor=0.9999; }
			if(cor<-0.9999){ cor=-0.9999; }
			zzz = cor/sqrt((1-cor*cor)/(n1-2));
		}
		z[i][j] = floor(1000*zzz+0.5)/1000;
		if(zzz >=2 && option_genes){
			genelist = get_coregulated_genes(symbols,a[i],a[ncol1+j],nrows,1,angle,EPFP_threshold);
			if(genelist){
				id_up[i][j] = data_id+1;
				database[data_id++] = genelist;
			}
			genelist = get_coregulated_genes(symbols,a[i],a[ncol1+j],nrows,-1,angle,EPFP_threshold);
			if(genelist){
				id_dn[i][j] = data_id+1;
				database[data_id++] = genelist;
			}
		}
	}
}
if(command) free(command);
return(data_id);
}

/***********************************************/
int main (long argc, char **argv) 
/***********************************************/
{
int ncol1, ncol2, autocorrelation=0, option_genes=0;
long nrows, **id_up, **id_dn, ndata, i, j;
char option, *input_file, *output_file, *error_text;
char **headers, **row_names, **symbols, **database;
float **a, **r, **z, angle=0, EPFP_thresh=0.3;
FILE *fpout;

if (argc < 7){
	check(error_text=(char*)malloc(sizeof(char*)*500));
	sprintf(error_text, "Error: insufficient number of arguments.\n"
	" Syntax: correlation_exatlas input_file output_file ncol1 ncol2 nrows option [-genes, -auto, -l logfile, -angle angle, -epfp EPFP]\n"
	" Input file should be tab-delimited\n"
	" option: P = Pearson correlation, S = Spearman correlation, C = covariance\n"
	" -genes: identify coregulated genes\n"
	" -auto:  autocorrelation\n"
	" -epfp:  Expected Proportion of False Positives threshold (default 0.3)\n"
	" -angle: angle for estimating EPFP\n");
	error_message(error_text);
}
input_file = copy_string(argv[1]);
output_file = copy_string(argv[2]);
ncol1 = atoi(argv[3]);
ncol2 = atoi(argv[4]);
nrows = atoi(argv[5]);
option = argv[6][0];
if(ncol1<1 || ncol1>MAXCOL) error_message("Number of columns 1 - wrong value");
if(ncol2<1 || ncol2>MAXCOL) error_message("Number of columns 2 - wrong value");
if(nrows<10) error_message("Not enough rows < 10");
if(option != 'P' && option != 'S' && option != 'C') error_message("Option should be P or S or C");
for(i=7; i<argc; i++){
	if(!strcmp(argv[i],"-auto")) autocorrelation = 1;
	else if(!strcmp(argv[i],"-genes")) option_genes = 1;
	else if(!strcmp(argv[i],"-angle")){
		if(i+1==argc) error_message("angle not specified");
		if(!sscanf(argv[++i],"%f",&angle)) error_message("angle non-readable");
		if(angle<0 || angle>1) error_message("angle should be from 0 to 1");
	}
	else if(!strcmp(argv[i],"-epfp")){
		if(i+1==argc) error_message("EPFP threshold not specified");
		if(!sscanf(argv[++i],"%f",&EPFP_thresh)) error_message("EPFP threshold non-readable");
		if(EPFP_thresh<=0.0001 || EPFP_thresh>1) error_message("EPFP threshold should be from 0.0001 to 1");
	}
	else if(!strcmp(argv[i],"-l")){
		if(i+1==argc) error_message("logfile name not specified");
		logfile = copy_string(argv[++i]);
	}
	else{
		error_message("unrecognized command-line parameter");
	}
}
a = new_matrix(ncol1+ncol2,nrows);
check(database = (char**)malloc(sizeof(char*)*ncol1*ncol2*2));
check(headers = (char**)calloc(ncol1,sizeof(char*)));
check(row_names = (char**)calloc(ncol2,sizeof(char*)));
check(symbols = (char**)calloc(nrows,sizeof(char*)));
read_data(input_file,a,ncol1,ncol2,nrows,headers,row_names,symbols);
r = new_matrix(ncol1,ncol2);
z = new_matrix(ncol1,ncol2);
id_up = new_matrix_long(ncol1,ncol2);
id_dn = new_matrix_long(ncol1,ncol2);
ndata = correlation_analysis(a,r,z,database,id_up,id_dn,ncol1,ncol2,nrows,option,symbols,autocorrelation,option_genes,angle,EPFP_thresh);
fpout = fopen(output_file, "a");
if(!fpout) error_message("Output file not opened!");
fprintf(fpout,"!Matrix1_start\n");
fprintf(fpout,"Headers");
for(i=0; i<ncol1; ++i){
	fprintf(fpout,"\t%s",headers[i]);
}
fprintf(fpout,"\n");
for(i=0; i<ncol2; ++i){
	fprintf(fpout,"%s",row_names[i]);
	for(j=0; j<ncol1; ++j){
		fprintf(fpout,"\t%.4f",z[j][i]);
	}
	fprintf(fpout,"\n");
}
fprintf(fpout,"!Matrix1_end\n");
fprintf(fpout,"!Matrix2_start\n");
fprintf(fpout,"Headers");
for(i=0; i<ncol1; ++i){
	fprintf(fpout,"\t%s",headers[i]);
}
fprintf(fpout,"\n");
for(i=0; i<ncol2; ++i){
	fprintf(fpout,"%s",row_names[i]);
	for(j=0; j<ncol1; ++j){
		fprintf(fpout,"\t%.4f",r[j][i]);
	}
	fprintf(fpout,"\n");
}
fprintf(fpout,"!Matrix2_end\n");
if(!option_genes){ fclose(fpout); exit(0); }
fprintf(fpout,"!Matrix3_start\n");
fprintf(fpout,"Headers");
for(i=0; i<ncol1; ++i){
	fprintf(fpout,"\t%s",headers[i]);
}
fprintf(fpout,"\n");
for(i=0; i<ncol2; ++i){
	fprintf(fpout,"%s",row_names[i]);
	for(j=0; j<ncol1; ++j){
		fprintf(fpout,"\t%d",id_up[j][i]);
	}
	fprintf(fpout,"\n");
}
fprintf(fpout,"!Matrix3_end\n");
fprintf(fpout,"!Matrix4_start\n");
fprintf(fpout,"Headers");
for(i=0; i<ncol1; ++i){
	fprintf(fpout,"\t%s",headers[i]);
}
fprintf(fpout,"\n");
for(i=0; i<ncol2; ++i){
	fprintf(fpout,"%s",row_names[i]);
	for(j=0; j<ncol1; ++j){
		fprintf(fpout,"\t%d",id_dn[j][i]);
	}
	fprintf(fpout,"\n");
}
fprintf(fpout,"!Matrix4_end\n");
fprintf(fpout,"!Database_start\n");
for(i=0; i<ndata; i++){
	fprintf(fpout,">%d\n%s\n",i+1,database[i]);
}
fprintf(fpout,"!Database_end\n");
fclose(fpout);
return(0);
}


