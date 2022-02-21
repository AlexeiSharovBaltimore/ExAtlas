/*********************************/
/* Principal Components Analysis */
/*********************************/

//pca 859283918885497.txt output.txt 10 0 V 1

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#define SIGN(a, b) ( (b) < 0 ? -fabs(a) : fabs(a) )

typedef struct cluster_st{
	long c1;
	long c2;
	long sign;
	long weight;
	float distance;
	float *v;
}CLUSTER;

void tqli(float *d, float *e, long n, float **z);
void tred2(float **a, long n, float *d, float *e);
void free_matrix(float **mat,long n, long m);
void free_vector(float *v, long n);
void scpcol(float **data, long n, long m, float **symmat);
void covcol(float **data, long n, long m, float **symmat, float *mean);
void corcol(float **data, long n, long m, float **symmat, float *mean);
void read_tabdelimited_item (FILE *fp, char *buffer);
long count_fields(char *buffer, long len);
float **new_matrix (long n, long m);
float *new_vector (long n);
void erhand (char *err_msg);
void free_index(long *v, long n);
void shellSort(float eigv[], const long array_size);
void resort (float var[], const long arr_size);
float **get_distance_matrix (float **symmat, float *evals, long m, long npc);
void combine_clusters (float **a, CLUSTER *c, long ic, long c1, long c2, long n1, long n);
CLUSTER *hierarchical_clustering (float **a, int n);
void print_clusters(CLUSTER *c, int n);
void check (void *x);

long* sortIndex;
FILE *OUTPUT;

/******************************************/
int main(int argc, char **argv)
/******************************************/
{
   FILE *stream;
   long  n, m,  i, j, k, k2, pc_number, *nn, nnn, center=0, nrow=0;
   float **data, **matrix(), **symmat, **symmat1, **symmat2, *vector(), *evals, *interm;
   float *evals_modif, *evals_modif1, *mean, *meanrow, **dist;
   float in_value, alpha, sum;
   char option;
   char **rowID, **colID, *buffer;
   CLUSTER *c1, *c2;

   if (argc <  5){
      printf("Syntax: pca.exe inputfile outputfile #PC alpha option [center nrows]\n\n");
      printf(" File should be tab-delimited in Stanford format\n");
      printf(" option   -- R for correlation analysis,\n");
      printf("             V for variance/covariance analysis\n");
      printf("             S for SSCP analysis.)\n");
      printf(" center   -- 0 non-centered data\n");
      printf("             1 centered rows\n");
      exit(1);
   }

   pc_number = atoi(argv[3]);
   sscanf(argv[4],"%f",&alpha);    /* alpha=1 if row metric, alpha=0 if column metric preserved */
   strncpy(&option,argv[5],1);     /* Analysis option */
   if(argc >= 7) center = atoi(argv[6]);
   if(argc >= 8) nrow = atoi(argv[7]);
   if(nrow < 0) nrow = 0;

   if((stream = fopen(argv[1],"r")) == NULL){ exit(1); }
   if((OUTPUT = fopen(argv[2],"w")) == NULL){ exit(1); }
   if(pc_number<2){ exit(1); }
   buffer = (char*)malloc(10000*sizeof(char));
   fgets(buffer,10000,stream);
   m = count_fields(buffer,strlen(buffer))-1;    /* number of columns */
   n = 0;
   while(fgets(buffer,10000,stream)) ++n;     /* count number of rows */
   if(nrow && n>nrow) n = nrow;
   if(pc_number > m || pc_number < 1) pc_number=m;
   rewind(stream);

   /* Now read in data. */

   data = new_matrix(n, m);  /* Storage allocation for input data */
   mean = new_vector(m);
   meanrow = new_vector(n);
   nn = (long*)calloc((m+1),sizeof(long*));
   rowID = (char**)malloc((n+1)*sizeof(char*));
   colID = (char**)malloc((m+1)*sizeof(char*));
   read_tabdelimited_item (stream, buffer);
   for (j = 1; j <= m; j++){
	read_tabdelimited_item (stream, buffer);
	colID[j] = (char*)malloc((strlen(buffer)+1)*sizeof(char*));
	strcpy(colID[j], buffer);
   }
   fgets(buffer,10000,stream);
   for (i = 1; i <= n; i++){
	read_tabdelimited_item (stream, buffer);
	if(!strlen(buffer)){
		n = i;
		break;
	}
        rowID[i] = (char*)malloc((strlen(buffer)+1)*sizeof(char*));
        strcpy(rowID[i], buffer);
	nnn = 0;
	meanrow[i] = 0;
        for (j = 1; j <= m; j++){
            fscanf(stream, "%f", &in_value);
            if(in_value <= -9999) in_value = -9999;
            else{
               meanrow[i] += in_value;
               nn[j]++;
               nnn++;
            }
            data[i][j] = in_value;
        }
        meanrow[i] /= nnn;
   }
   free(buffer);

   /* Substitute missing values, and center data. */
   for (i = 1; i <= n; i++){
      for (j = 1; j <= m; j++){
         if(data[i][j] == -9999) data[i][j] = meanrow[i];
         if(center)
            data[i][j] -= meanrow[i];
         mean[j] += data[i][j];
      }
   }
   if(center) fprintf(OUTPUT,"Data centered. ");

   /* Estimate means */
   for (j = 1; j <= m; j++)
      mean[j] /= n;

   /* Check on (part of) input data.
   for (i = 1; i <= 18; i++){
      for (j = 1; j <= 8; j++){
        printf("%.4f\t", data[i][j]);
      }
      printf("\n");
   }
   */
   symmat = new_matrix(m, m);  /* Allocation of correlation (etc.) matrix */

   /* Look at analysis option; branch in accordance with this. */

   switch(option){
          case 'R':
          case 'r':
              fprintf(OUTPUT,"Analysis of correlations chosen.\n");
              corcol(data, n, m, symmat, mean);

                   /* Output correlation matrix.
                   for (i = 1; i <= m; i++) {
                      for (j = 1; j <= 8; j++) {
                         printf("%7.4f", symmat[i][j]);
                      }
                      printf("\n");  }
                   */
              break;
          case 'V':
          case 'v':
              fprintf(OUTPUT,"Analysis of variances-covariances chosen.\n");
              covcol(data, n, m, symmat, mean);

                          /* Output variance-covariance matrix.
                          for (i = 1; i <= m; i++) {
                          for (j = 1; j <= m; j++)  {
                            printf("%.4f\t", symmat[i][j]);  }
                            printf("\n");  }
                          */
              break;
          case 'S':
          case 's':
              fprintf(OUTPUT,"Analysis of sums-of-squares-cross-products matrix chosen.\n");
              scpcol(data, n, m, symmat);

                         /* Output SSCP matrix.
                         for (i = 1; i <= m; i++) {
                          for (j = 1; j <= 8; j++)  {
                            printf("%.4f\t", symmat[i][j]);  }
                            printf("\n");  }
                         */
              break;
          default:
              printf("Option: %c\n", option);
              printf("For option, please type R, V, or S\n");
              printf("(upper or lower case).\n");
              printf("Exiting to system.\n");
              exit(1);
              break;
   }

   /* Allocate storage for dummy and new vectors. */
   evals = new_vector(m);     /* Storage alloc. for vector of eigenvalues */
   evals_modif = new_vector(m);     /* Modified eigenvalues*/
   evals_modif1 = new_vector(m);     /* Modified eigenvalues*/
   interm = new_vector(m);    /* Storage alloc. for 'intermediate' vector */
   symmat1 = new_matrix(m, m);  /* needed for eigenvectors */
   symmat2 = new_matrix(m, m);  /* Duplicate of correlation (etc.) matrix */
   for (i = 1; i <= m; i++) {
     for (j = 1; j <= m; j++) {
      symmat2[i][j] = symmat[i][j]; /* Needed below for col. projections */
     }
   }
   tred2(symmat, m, evals, interm);  /* Triangular decomposition */
   tqli(evals, interm, m, symmat);   /* Reduction of sym. trid. matrix */

   /* evals now contains the eigenvalues,
       columns of symmat now contain the associated eigenvectors. */

   sortIndex = (long*)calloc(m+1,sizeof(long));  //initializes the sort index to hold the relative order of values
   shellSort(evals, m); //sorts the eigenvalues and retains relative values in sortIndex

   fprintf(OUTPUT,"Eigenvalues:\n");
   for (j = m; j >= 1; j--) {
       fprintf(OUTPUT,"%.5f\n", evals[j]); 
       evals_modif[m+1-j] = evals[j];
       evals_modif1[m+1-j] = evals[j];
   }

   fprintf(OUTPUT,"\nEigenvectors: (in terms of original vbes.)\n");
   for (j = 1; j <= m; j++) {
       fprintf(OUTPUT,"%s\t", colID[j]);
       resort(symmat[j], m); // sorted based on sortIndex
       for (i = 1; i <= pc_number; i++){
            symmat1[j][i] = symmat[j][m-i+1];
          fprintf(OUTPUT,"%.5f", symmat[j][m-i+1]);
          if(i<pc_number) fprintf(OUTPUT,"\t"); 
       }
       fprintf(OUTPUT,"\n");  
   }
   dist = get_distance_matrix (symmat1,evals_modif,m,pc_number);
   c1 = hierarchical_clustering (dist,m);
   free_matrix(dist,m+1,m+1);

   /* Form projections of row-points on first pc_number prin. components. */
   /* Store in 'data', overwriting original data. */
   for (i = 1; i <= n; i++) {
      for (j = 1; j <= m; j++) {
        interm[j] = data[i][j] - mean[j];   /* data[i][j] will be overwritten */
      }
      for (k = 1; k <= pc_number; k++) {
          data[i][k] = 0.0;
          for (k2 = 1; k2 <= m; k2++){
            data[i][k] += interm[k2] * symmat[k2][m-k+1]; 
          }
      }
   }

   /* Modify eigenvalues */
   for (j = 1; j <= m; j++) {
     if(evals_modif[j] > 0.000000001)
       evals_modif[j] = exp(log(evals_modif[j])*(alpha-1.0)/2);
     else if (alpha == 1)
       evals_modif[j] = 1;
     else
       evals_modif[j] = 0.000000001;
   }

   fprintf(OUTPUT,"\nRow means (2nd column) and projections of row-points on PC (alpha=%.1f):\n", alpha);
   for (i = 1; i <= n; i++) {
       fprintf(OUTPUT,"%s\t%.4f", rowID[i], meanrow[i]);
       for (j = 1; j <= pc_number; j++)  {
          data[i][j] *= evals_modif[j];
          fprintf(OUTPUT,"\t%.5f", data[i][j]);  
       }
       fprintf(OUTPUT,"\n");  
   }
   fprintf(OUTPUT,"\nSorted columns by hierarchical clustering\n");
   print_clusters(c1,m);
   free(c1);
   if(n <= 2000){
      dist = get_distance_matrix (data,evals_modif1,n,pc_number);
      c2 = hierarchical_clustering (dist,n);
      free_matrix(dist,n+1,n+1);
      fprintf(OUTPUT,"\nSorted rows by hierarchical clustering\n");
      print_clusters(c2,n);
      free(c2);
   }
   fclose(OUTPUT);
   free_matrix(data, n, m);
   free_matrix(symmat, m, m);
   free_matrix(symmat1, m, m);
   free_matrix(symmat2, m, m);
   free_vector(evals, m);
   free_vector(evals_modif, m);
   free_vector(interm, m);
   return(0);
}

/***********************************************/
long	count_fields(char *buffer, long len)
/***********************************************/
{
	long i, n=1;
	for(i=0; i<len && buffer[i] != '\n'; ++i){
		if(buffer[i] == '\t') ++n;
	}
	return(n);
}

/***********************************************/
void  read_tabdelimited_item (FILE *fp, char *buffer)
/***********************************************/
{
	char ch;
	long i=0, j=0;

	while(ch=fgetc(fp)){
		if(j>0 && (ch=='\t' || ch=='\n')){
			ungetc (ch, fp);
			break;
		}
		if(j==0 && (ch=='\t' || ch=='\n')){
			j++;
			continue;
		}
		j++;
		if(i==0 && ch==' ') continue;
		if(i < 80-1) buffer[i++] = ch;	
	}
	buffer[i] = '\0';
}

/******* Correlation matrix: creation  ********************/
void corcol(float **data, long n, long m, float **symmat, float *mean)
/**********************************************************/
/* Create m * m correlation matrix from given n * m data matrix. */
{
float eps = 0.005;
float x, *stddev, *vector();
long i, j, j1, j2;

/* Allocate storage for std. dev. vectors */
stddev = new_vector(m);

/* Determine standard deviations of column vectors of data matrix. */

for (j = 1; j <= m; j++){
    stddev[j] = 0.0;
    for (i = 1; i <= n; i++){
        stddev[j] += (data[i][j]-mean[j]) * (data[i][j]-mean[j]);
    }
    stddev[j] /= (float)n;
    stddev[j] = sqrt(stddev[j]);
    /* The following in an inelegant but usual way to handle
    near-zero std. dev. values, which below would cause a zero-
    divide. */
    if (stddev[j] <= eps) stddev[j] = eps;
}

/* fprintf(OUTPUT,"\nStandard deviations of columns:\n");
   for (j = 1; j <= m; j++) { fprintf(OUTPUT,"%.4f %.4f\t", stddev[j], mean[j]); }
   fprintf(OUTPUT,"\n");
*/

/* Center and reduce the column vectors. */

for (i = 1; i <= n; i++) {
    for (j = 1; j <= m; j++){
        x = sqrt((float)n);
        x *= stddev[j];
        data[i][j] -= mean[j];
        data[i][j] /= x;
    }
}

/* Calculate the m * m correlation matrix. */
for (j1 = 1; j1 <= m-1; j1++){
    symmat[j1][j1] = 1.0;
    for (j2 = j1+1; j2 <= m; j2++){
        symmat[j1][j2] = 0.0;
        for (i = 1; i <= n; i++) {
            symmat[j1][j2] += ( data[i][j1] * data[i][j2]);
        }
        symmat[j2][j1] = symmat[j1][j2];
    }
}
symmat[m][m] = 1.0;
return;
}

/*************  Variance-covariance matrix: creation  *****************/
void covcol(float **data, long n, long m, float **symmat, float *mean)
/**********************************************************/
/* Create m * m covariance matrix from given n * m data matrix. */
{
float *vector();
long i, j, j1, j2;

/* Calculate the m * m covariance matrix. */
for (j1 = 1; j1 <= m; j1++){
    for (j2 = j1; j2 <= m; j2++){
        symmat[j1][j2] = 0.0;
        for (i = 1; i <= n; i++){
            symmat[j1][j2] += (data[i][j1]-mean[j1]) * (data[i][j2]-mean[j2]);
        }
        symmat[j1][j2] /= (n-1);
        symmat[j2][j1] = symmat[j1][j2];
    }
}
return;
}

/************  Sums-of-squares-and-cross-products matrix: creation  ************/
void scpcol(float **data, long n, long m, float **symmat)
/**********************************************************/
/* Create m * m sums-of-cross-products matrix from n * m data matrix. */
{
long i, j1, j2;

/* Calculate the m * m sums-of-squares-and-cross-products matrix. */

for (j1 = 1; j1 <= m; j1++)
    {
    for (j2 = j1; j2 <= m; j2++)
        {
        symmat[j1][j2] = 0.0;
        for (i = 1; i <= n; i++)
            {
            symmat[j1][j2] += data[i][j1] * data[i][j2];
            }
        symmat[j2][j1] = symmat[j1][j2];
        }
    }
return;
}

/******************************************/
float **get_distance_matrix (float **symmat, float *evals, long m, long npc)
/******************************************/
{
float **a;
float x, y, sxy;
long i, j, k, n;
long ncol, nrow;

a = new_matrix(m+1,m+1);
/* for(i=1; i<=19; ++i){
	fprintf(OUTPUT,"%d\t%f\n",i,evals[i]);
}
*/
for(i=1; i<=m; ++i){
	for(j=i+1; j<=m; ++j){
		sxy = 0;
		n = 0;
		for(k=1; k<=npc; ++k){
			x = symmat[i][k];
			y = symmat[j][k];
			sxy += (x-y)*(x-y)*evals[k];
			++n;
		}
		a[i][j] = sqrt(sxy);
		a[j][i] = a[i][j];
		if(a[i][j] > 10000000000){
			fprintf(OUTPUT,"%d %d\t%f\n",i,j,a[i][j]);
			for(k=1; k<=npc; ++k){
				x = symmat[i][k];
				y = symmat[j][k];
				fprintf(OUTPUT,"%d\t%f\t%f\n",k,x,y);
			}
		}
	}
}
return(a);
}

/******************************************/
void combine_clusters  (float **a, CLUSTER *c, long ic, long c1, long c2, long n1, long n)
/******************************************/
{
long i, j, i1, i2, ic1, ic2;
float w1=1, w2=1, x;

ic1 = c[ic].c1;
if(ic1 > n){
	ic1 -= n;
	w1 = c[ic1].weight;
	if(c[ic1].v[c2] < 0){ c[ic1].sign = -1; }
	else{ c[ic1].sign = 1; }
}
ic2 = c[ic].c2;
if(ic2 > n){
	ic2 -= n;
	w2 = c[ic2].weight;
	if(c[ic2].v[c1] < 0){ c[ic2].sign = -1; }
	else{ c[ic2].sign = 1; }
}
c[ic].weight = w1+w2;
c[ic].sign = 1;
for(i=1; i<=n1; i++){
	c[ic].v[i] = a[c1][i]-a[c2][i];
	a[i][n1+1] = (w1*a[c1][i] + w2*a[c2][i])/(w1+w2);
	a[n1+1][i] = a[i][n1+1];
}
for(i=1; i<=n1+1; ++i){
	if(c2-c1>1) memmove(&(a[i][c1]),&(a[i][c1+1]),(c2-c1-1)*sizeof(float));
	memmove(&(a[i][c2-1]),&(a[i][c2+1]),(n1-c2+1)*sizeof(float));
}
for(i=c1; i<=c2-2; ++i){
	memcpy(a[i],a[i+1],(n1+2)*sizeof(float));
}
for(i=c2-1; i<=n1-1; ++i){
	memcpy(a[i],a[i+2],(n1+2)*sizeof(float));
}
for(i=1; i<=n-n1+1; i++){
	x = (w1*c[i].v[c1] + w2*c[i].v[c2])/(w1+w2);
	if(c2-c1>1) memmove(&(c[i].v[c1]),&(c[i].v[c1+1]),(c2-c1-1)*sizeof(float));
	memmove(&(c[i].v[c2-1]),&(c[i].v[c2+1]),(n1-c2)*sizeof(float));
	c[i].v[n1-1] = x;
}
return;
}

/******************************************/
CLUSTER  *hierarchical_clustering (float **a, int n)
/******************************************/
{
long i, j, n1, i1, c1, c2, is=1, k=1, ic, ic1, ic2, sign;
float distance, s=1, maxx;
long *icluster;
CLUSTER *c;

check(icluster = (long*)calloc((2*n+1),sizeof(long)));
check(c = (CLUSTER*)calloc(n+1,sizeof(CLUSTER)));  // This is an array of clusters

/* Initialize clusters */
for(i=1;i<=n;++i){
	icluster[i] = i;
	icluster[i+n] = i+n;
	c[i].v = new_vector(n+1);
}
n1 = n;
while(n1 >= 2){
	/* Find minimal distance */
	distance = 10000000000;
	maxx = 0;
	for(i=1;i<=n1;++i){
		for(j=i+1;j<=n1;++j){
			if(distance > a[i][j]){
				distance = a[i][j];
				c1 = i;
				c2 = j;
			}
			//if(maxx < a[i][j]) maxx = a[i][j];
			//if(a[i][j] > 10000000000) fprintf(OUTPUT,"%d %d\t%f\n",i,j,a[i][j]);
		}
	}
	if(distance==10000000000){ fprintf(OUTPUT,"Dist=10000000000\n"); exit(0); }
	ic = n-n1+1;

	c[ic].c1 = icluster[c1]; // Mark memebers
	c[ic].c2 = icluster[c2];
	c[ic].distance = distance; // Distance between members

	/* if(n1 == n){
	for(j=1;j<=n1;++j){
		fprintf(OUTPUT,"\t%d",icluster[j]);
	}
	fprintf(OUTPUT,"\n");
	for(i=1; i<=n1; ++i){
		fprintf(OUTPUT,"%d",icluster[i]);
		for(j=1;j<=n1;++j){
			fprintf(OUTPUT,"\t%f",a[i][j]);
		}
		fprintf(OUTPUT,"\n");
	}
	fprintf(OUTPUT,"\n");
	}*/
	
	//fprintf(OUTPUT,"%d\t%d\t%d\t%f\t%f\n",ic+n,icluster[c1],icluster[c2],distance,maxx);
	combine_clusters(a, c, ic, c1, c2, n1, n);  // shrink matrix a to n1-1

	for(i=c1+1;i<=2*n;++i){
		i1 = i-1;
		if(i>c2) i1--;
		icluster[i1] = icluster[i];
	}
	--n1;
}
free (icluster);
for(i=1;i<=n;++i) free(c[i].v);
return(c);
}

/****************************************************/
void print_clusters(CLUSTER *c, int n)
/**********************************************************/
{
long *stackCl, *stackSign, is=1, ic, ic1, ic2, sign, nClusters;
float dist_tresh;

nClusters = (long)sqrt((float)n)*3;
if(nClusters > n-1) nClusters = n-1;
check(stackCl = (long*)calloc(3*n+1,sizeof(long)));
check(stackSign = (long*)calloc(3*n+1,sizeof(long)));

/*for(ic=1;ic<=n;++ic)
	printf("D: %d %d %d %f\n",ic,c[ic].c1,c[ic].c2,c[ic].distance);
*/
stackCl[is] = 2*n-1;
stackSign[is++] = 1;
while(is>1){
	ic = stackCl[--is];
	sign = stackSign[is];
	if(ic <= n){
		fprintf(OUTPUT,"%d,",ic);
		continue;
	}
	ic -= n;
	ic1 = c[ic].c1;
	ic2 = c[ic].c2;

//printf("%d - %d %d %d\n",is,ic,ic1,ic2);

	if(c[ic].sign*sign > 0){
		stackCl[is] = ic2;
		stackSign[is++] = 1;
		stackCl[is] = ic1;
		stackSign[is++] = -1;
	}else{
		stackCl[is] = ic1;
		stackSign[is++] = -1;
		stackCl[is] = ic2;
		stackSign[is++] = 1;
	}
}
fprintf(OUTPUT,"\n");

/* for(ic=2*n-1; ic>n; ic--){
	is = ic-n;
	printf("%d\t%d\t%d\t%d\t%f\n",ic,c[is].c1,c[is].c2,c[is].sign,c[is].distance);
}*/

free (stackCl);
free (stackSign);
return;
}

/**  Error handler  **************************************************/
void erhand(char *err_msg)
/**********************************************************/
{
    fprintf(stderr,"ERROR %s\n", err_msg);
    exit(1);
}

/**  Allocation of vector storage  ***********************************/
float *new_vector(long n)
/**********************************************************/
/* Allocates a float vector with range [1..n]. */
{
    float *v;
    v = (float *)calloc ((unsigned)(n+1),sizeof(float));
    if (!v) erhand("Allocation failure in vector().");
    return v;
}

/**  Allocation of float matrix storage  *****************************/
float **new_matrix(long n, long m)
/**********************************************************/
/* Allocate a float matrix with range [1..n][1..m]. */
{
    long i;
    float **mat;

    /* Allocate pointers to rows. */
    mat = (float **) malloc((unsigned)(n+1)*sizeof(float*));
    if (!mat) erhand("Allocation failure 1 in matrix().");
    /* Allocate rows and set pointers to them. */
    for (i = 1; i <= n; i++)
        {
        mat[i] = (float *) calloc((unsigned)(m+1),sizeof(float));
        if (!mat[i]) erhand("Allocation failure 2 in matrix().");
        }
     /* Return pointer to array of pointers to rows. */
     return mat;
}

/**  Deallocate vector storage  *********************************/
void free_vector(float *v, long n)
/**********************************************************/
/* Free a float vector allocated by vector(). */
{
   free(v);
}

/**  Deallocate float matrix storage  ***************************/
void free_matrix(float **mat,long n,long m)
/**********************************************************/
/* Free a float matrix allocated by matrix(). */
{
   long i;
   for (i = n; i >= 1; i--){
       free(mat[i]);
   }
   free(mat);
}

/**  Reduce a real, symmetric matrix to a symmetric, tridiag. matrix. */
void tred2(float **a, long n, float *d, float *e)
/**********************************************************/
/* Householder reduction of matrix a to tridiagonal form.
   Algorithm: Martin et al., Num. Math. 11, 181-195, 1968.
   Ref: Smith et al., Matrix Eigensystem Routines -- EISPACK Guide
        Springer-Verlag, 1976, pp. 489-494.
        W H Press et al., Numerical Recipes in C, Cambridge U P,
        1988, pp. 373-374.  */
{
long l, k, j, i;
float scale, hh, h, g, f;

for (i = n; i >= 2; i--)
    {
    l = i - 1;
    h = scale = 0.0;
    if (l > 1)
       {
       for (k = 1; k <= l; k++)
           scale += fabs(a[i][k]);
       if (scale == 0.0)
          e[i] = a[i][l];
       else
          {
          for (k = 1; k <= l; k++)
              {
              a[i][k] /= scale;
              h += a[i][k] * a[i][k];
              }
          f = a[i][l];
          g = f>0 ? -sqrt(h) : sqrt(h);
          e[i] = scale * g;
          h -= f * g;
          a[i][l] = f - g;
          f = 0.0;
          for (j = 1; j <= l; j++)
              {
              a[j][i] = a[i][j]/h;
              g = 0.0;
              for (k = 1; k <= j; k++)
                  g += a[j][k] * a[i][k];
              for (k = j+1; k <= l; k++)
                  g += a[k][j] * a[i][k];
              e[j] = g / h;
              f += e[j] * a[i][j];
              }
          hh = f / (h + h);
          for (j = 1; j <= l; j++)
              {
              f = a[i][j];
              e[j] = g = e[j] - hh * f;
              for (k = 1; k <= j; k++)
                  a[j][k] -= (f * e[k] + g * a[i][k]);
              }
         }
    }
    else
        e[i] = a[i][l];
    d[i] = h;
    }
d[1] = 0.0;
e[1] = 0.0;
for (i = 1; i <= n; i++)
    {
    l = i - 1;
    if (d[i])
       {
       for (j = 1; j <= l; j++)
           {
           g = 0.0;
           for (k = 1; k <= l; k++)
               g += a[i][k] * a[k][j];
           for (k = 1; k <= l; k++)
               a[k][j] -= g * a[k][i];
           }
       }
       d[i] = a[i][i];
       a[i][i] = 1.0;
       for (j = 1; j <= l; j++)
           a[j][i] = a[i][j] = 0.0;
    }
}

/**  Tridiagonal QL algorithm -- Implicit  **********************/
void tqli(float *d, float *e, long n, float **z)
/**********************************************************/
{
long m, l, iter, i, k;
float s, r, p, g, f, dd, c, b;

for (i = 2; i <= n; i++)
    e[i-1] = e[i];
e[n] = 0.0;
for (l = 1; l <= n; l++)
    {
    iter = 0;
    do
      {
      for (m = l; m <= n-1; m++)
          {
          dd = fabs(d[m]) + fabs(d[m+1]);
          if (fabs(e[m]) + dd == dd || iter>100) break;
          }
          if (m != l)
             {
             if (iter++ == 60) // erhand("No convergence in TLQI.");
             g = (d[l+1] - d[l]) / (2.0 * e[l]);
             r = sqrt((g * g) + 1.0);
             g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
             s = c = 1.0;
             p = 0.0;
             for (i = m-1; i >= l; i--)
                 {
                 f = s * e[i];
                 b = c * e[i];
                 if (fabs(f) >= fabs(g))
                    {
                    c = g / f;
                    r = sqrt((c * c) + 1.0);
                    e[i+1] = f * r;
                    c *= (s = 1.0/r);
                    }
                 else
                    {
                    s = f / g;
                    r = sqrt((s * s) + 1.0);
                    e[i+1] = g * r;
                    s *= (c = 1.0/r);
                    }
                 g = d[i+1] - p;
                 r = (d[i] - g) * s + 2.0 * c * b;
                 p = s * r;
                 d[i+1] = g + p;
                 g = c * r - b;
                 for (k = 1; k <= n; k++)
                     {
                     f = z[k][i+1];
                     z[k][i+1] = s * z[k][i] + c * f;
                     z[k][i] = c * z[k][i] - s * f;
                     }
                 }
                 d[l] = d[l] - p;
                 e[l] = g;
                 e[m] = 0.0;
             }
          }  while (m != l);
      }
 }

void shellSort(float eigv[], const long array_size)
{	//sorts the eigenvalues in evals
 //will make the index for resorting "sortIndex"
  long i, j, increment,temp2;
  float temp;

  for (i=1; i <= array_size; i++){
    sortIndex[i] = i;
  }
  increment = 3;
  while (increment > 0)
  {
    for (i=1; i <= array_size; i++)
    {
      j = i;
      temp = eigv[i];
      temp2 = sortIndex[i];
      while ((j > increment) && (eigv[j-increment] > temp))
      {
        eigv[j] = eigv[j - increment];
        sortIndex[j] = sortIndex[j - increment];

        j = j - increment;
      }
      eigv[j] = temp;
      sortIndex[j] = temp2;
    }
    if (increment/2 != 0)
      increment = increment/2;
    else if (increment == 1)
      increment = 0;
    else
      increment = 1;
  }
}

void resort (float var[], const long arr_size)
{//resort a vector using the index of relative order in sortIndex
   float *temp1 = new_vector(arr_size);
   long i;
   for (i=1; i<= arr_size; i++){
       temp1[i] = var[sortIndex[i]];
   }
   for (i=1; i<= arr_size; i++){
      var[i] =  temp1[i];
   }
   free(temp1);
}

/***********************************************/
void check (void *x)
/***********************************************/
{
if (!x){
	fprintf(OUTPUT,"Out of memory\n");
	exit(1);
}
}


