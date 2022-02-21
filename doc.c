/*************************************
* Programmer: Alexei Sharov   11/03/2000
* Department of Entomlogy, Virginia Tech
* (540)231-7316, e-mail sharov@vt.edu
* Home page: www.ento.vt.edu/~sharov
* 
* This file is a library of functions for handling raster data
*
* It is a part of the Slow-the-Spread
* Decision-support system.
**************************************/

#define  NDATATYPES  3
#define  BYTE    0
#define  INTEGER 1
#define  REAL    2

#define  NFILETYPES  2
#define  BINARY  0
#define  ASCII   1

#define  NITEMS  14	/* items in documentation file */
#define  MAXCOL 256	/* maximum number of colors in a map */

/* Structure that holds documentation for a file with raster data
* in IDRISI format */

typedef struct doc_st{
	char title[80];
	int data_type;
	int file_type;
	int nx, ny;		/* number of columns and rows */
	float xmin, ymin;	/* coordinates of left lower corner */
	float xmax, ymax;	/* coordinates of right upper corner */
	float thresh[MAXCOL];	/* upper threshold for a particular color in the map */
	int col[MAXCOL];	/* color number that corresponds to each threshold */
	int ncol;		/* number of colors */
}DOC;

static char *data_types[NDATATYPES] = {"byte","integer","real"};
static char *file_types[NFILETYPES] = {"binary","ascii"};
static char *items[NITEMS] =
{"image title", "data type", "file type", "rows", "columns",
 "min. X", "max. X", "min. Y", "max. Y", "legend",
 "cell x", "cell y", "colors"};

/*************************************/
void  modify_name  (char *nameold, char *namenew, char *modif)
/*************************************
* Replaces extention in the file name. It is used to generate
* ".doc" file names from ".img" file names, and vice versa.
*/
{
	char *ps;

	strcpy(namenew, nameold);
	ps = strrchr(namenew,'.');
	if (ps) *ps = '\0';
	strcat(namenew, modif);
}

/*********************************/
void	print_doc (DOC *doc, char *docfile)
/*********************************
* Prints a documentation file
*/
{
	FILE *fp;
	int i;

	fp = fopen(docfile,"w");
	fprintf(fp,"%-12s: %s\n", items[0], doc->title);
	fprintf(fp,"%-12s: %s\n", items[1], data_types[doc->data_type]);
	fprintf(fp,"%-12s: %s\n", items[2], file_types[doc->file_type]);
	fprintf(fp,"%-12s: %d\n", items[3], doc->ny);
	fprintf(fp,"%-12s: %d\n", items[4], doc->nx);
	fprintf(fp,"%-12s: %f\n", items[5], doc->xmin);
	fprintf(fp,"%-12s: %f\n", items[6], doc->xmax);
	fprintf(fp,"%-12s: %f\n", items[7], doc->ymin);
	fprintf(fp,"%-12s: %f\n", items[8], doc->ymax);
	fprintf(fp,"%-12s: %d\n", items[9], 0);
	if (doc->ncol > 0 && doc->ncol < MAXCOL)
	{
		fprintf(fp,"%12s: %d\n", items[13], doc->ncol);
		for (i=0; i<doc->ncol; ++i)
			fprintf(fp,"\t%.4f\t%d\n", doc->thresh[i], doc->col[i]);
	}
	fclose (fp);
}

/*********************************/
int get_type (char *value, char **list, int n)
/*********************************/
{
	int i;
	for(i=0;i<n;++i)
		if (strcmp(value, list[i]) == 0)
			return (i);
	return(-1);
}

/*********************************/
int  get_item (char *buffer, char *attr, char *value)
/*********************************/
{
	char *c;
	int len;

	c = strchr(buffer,':');
	if (!c)
		return(0);
	len = c-buffer;
	strncpy(attr,buffer,len);
	while(len > 0 && attr[len-1] == ' ')
		--len;
	if (len == 0)
		return(0);
	attr[len] = '\0';
	len = c-buffer;
	while (len < (int)strlen(buffer) && buffer[len+1] == ' ')
		++len;
	if (len >= (int)strlen(buffer)-1)
		return(0);
	strcpy(value, buffer+len+1);
	return(1);
}

/*********************************/
void  read_doc (DOC *doc, char *filename)
/*********************************
* Reads raster documentation from file
*/
{
	char buffer[80], attr[20], value[20];
	FILE *fp;
	int it, k, i;
	float cellx=0, celly=0;

	doc->xmin = -1.0e30;
	doc->xmax = -1.0e30;
	doc->ymin = -1.0e30;
	doc->ymax = -1.0e30;
	doc->nx = 0;
	doc->ny = 0;
	doc->ncol = 0;
	fp = fopen(filename,"r");
	if (!fp)
		error_message("File '*.doc' not found!");
	while (fgets(buffer, 79, fp))
	{
		i = strlen(buffer);
		while (i > 0 && (buffer[i-1]=='\n' || buffer[i-1]=='_' ||
			buffer[i-1]=='\t'))
				--i;
		buffer[i] = '\0';
		if (!get_item(buffer, attr, value))
			continue;
		it = get_type(attr,items, NITEMS);
		switch(it)
		{
			case 0: strcpy(doc->title, value); break;
			case 1: doc->data_type = get_type(value,data_types,NDATATYPES); break;
			case 2: doc->file_type = get_type(value,file_types,NFILETYPES); break;
			case 3: doc->ny = atoi(value); break;
			case 4: doc->nx = atoi(value); break;
			case 5: sscanf(value,"%f", &doc->xmin); break;
			case 6: sscanf(value,"%f", &doc->xmax); break;
			case 7: sscanf(value,"%f", &doc->ymin); break;
			case 8: sscanf(value,"%f", &doc->ymax); break;
			case 10: sscanf(value,"%f", &cellx); break;
			case 11: sscanf(value,"%f", &celly); break;
			case 12: sscanf(value,"%f", &doc->ncol);
				if (doc->ncol > 0 && doc->ncol < MAXCOL)
				{
					for (i=0; i<doc->ncol; ++i)
						fscanf(fp,"%f%d", &doc->thresh[i], &doc->col[i]);
				}
/*				else
					printf("Bad number of colors!\n"); */
				break;
		}
	}
	fclose(fp);
	if (doc->nx <=0 || doc->ny <= 0)
		error_message("nx or ny not defined!");
	if (doc->xmin < -1.0e29)
		doc->xmin = 0;
	if (doc->xmax < -1.0e29)
	{
		if(cellx <= 0)
			cellx = 1;
		doc->xmax = doc->xmin + doc->nx*cellx;
	}
	if (doc->ymin < -1.0e29)
		doc->ymin = 0;
	if (doc->ymax < -1.0e29)
	{
		if(celly <= 0)
			celly = 1;
		doc->ymax = doc->ymin + doc->ny*celly;
	}
}

/*********************************/
void compare_docs  (DOC *d1, DOC *d2)
/*********************************
* Checks if documentation files for 2 raster maps are compatible
*/
{
	if (d1->nx != d2->nx || d1->ny != d2->ny ||
	 d1->xmin != d2->xmin || d1->ymin != d2->ymin ||
	 d1->xmax != d2->xmax || d1->ymax != d2->ymax)
		error_message("Docs are different\n");
}

