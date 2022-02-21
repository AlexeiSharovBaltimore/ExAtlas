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

#include "doc.c"

typedef struct tagBOARD{
	float **statef; /* pointer to a grid (raster) with float values */
	int  **statei;	/* pointer to a grid with integer values */
	char **statec;	/* pointer to a grid with 'char' values */
	DOC doc;
}BOARD;

FILE *open_board_file (char *filename, DOC *d, char *attr);
void   read_board (BOARD *b, char *filename);

/*********************************/
BOARD *get_board (void)
/*********************************
* Allocates memory for structure BOARD
*/
{
	BOARD *b;

	check(b = (BOARD*)calloc(1, sizeof(BOARD)));
	return(b);
}

/***************************************/
DOC *get_doc (BOARD *b)
/***************************************
* Returns a pointer tothe DOC structure
*/
{
	return(&(b->doc));
}

/***************************************/
int correct_entry (int x, int y, DOC *d)
/***************************************/
{
	if (x >=0 && x < d->nx && y >=0 && y < d->ny)
		return(1);
	else
		return(0);
}

/*********************************/
void unshrink_board (BOARD *b)
/*********************************
* Allocates memory for raster data
*/
{
	int i;
	DOC *d;

	d = get_doc(b);
	switch(d->data_type)
	{
		case BYTE:
			check(b->statec = (char**)malloc(d->nx*sizeof(char*)));
			for (i=0; i<d->nx; ++i)
				check(b->statec[i] = (char*)calloc(d->ny,sizeof(char)));
			break;
		case INTEGER:
			check(b->statei = (int**)malloc(d->nx*sizeof(int*)));
			for (i=0; i<d->nx; ++i)
				check(b->statei[i] = (int*)calloc(d->ny,sizeof(int)));
			break;
		case REAL:
			check(b->statef = (float**)malloc(d->nx*sizeof(float*)));
			for (i=0; i<d->nx; ++i)
				check(b->statef[i] = (float*)calloc(d->ny,sizeof(float)));
			break;
	}
}

/*********************************/
BOARD *new_board (char *filename)
/*********************************
* Loads raster data from file 'filename'. Returns a pointer
* to the new BOARD structure
*/
{
	BOARD *b;
	char docfile[80];

	b = get_board();
	modify_name(filename, docfile, ".doc");
	read_doc(&b->doc, docfile);
	unshrink_board(b);
	read_board(b, filename);
	return(b);
}

/*********************************/
void shrink_board (BOARD *b)
/*********************************
* Free memory
*/
{
	int i;
	DOC *d;

	d = get_doc(b);
	switch(d->data_type)
	{
		case BYTE:
			for (i=0; i<d->nx; ++i)
				free(b->statec[i]);
			free(b->statec);
			break;
		case INTEGER:
			for (i=0; i<d->nx; ++i)
				free(b->statei[i]);
			free(b->statei);
			break;
		case REAL:
			for (i=0; i<d->nx; ++i)
				free(b->statef[i]);
			free(b->statef);
			break;
	}
}

/*********************************/
void dispose_board (BOARD *b)
/*********************************/
{
	shrink_board(b);
	free(b);
}

/*********************************/
float get_cell (int ix, int iy, BOARD *b)
/*********************************
* Returns cell value in column ix and row iy.
*/
{
	DOC *d;
	float val=0;

	d = get_doc(b);
	switch(d->data_type)
	{
		case BYTE:    val = b->statec[ix][iy]; break;
		case INTEGER: val = b->statei[ix][iy]; break;
		case REAL:    val = b->statef[ix][iy]; break;
	}
	return(val);
}

/*********************************/
void put_cell (float val, int ix, int iy, BOARD *b)
/*********************************
* Sets the value of a cell in column ix and row iy to val
*/
{
	DOC *d;
	float valint;

	d = get_doc(b);
	valint = val;
	if (valint >= 0)
		valint = (int)(valint + 0.5);
	else
		valint = (int)(valint - 0.5);
	switch(d->data_type)
	{
		case BYTE:    b->statec[ix][iy] = valint; break;
		case INTEGER: b->statei[ix][iy] = valint; break;
		case REAL:    b->statef[ix][iy] = val; break;
	}
}

/*********************************/
void   print_board (BOARD *b, char *filename)
/*********************************
* Prints raster information to a file in IDRISI format
*/
{
	int jx, jy;
	FILE *fp;
	float val;
	DOC *d;
	char docfile[80];

	d = get_doc(b);
	fp = open_board_file (filename, d, "w");
	for (jy = d->ny-1; jy >= 0; --jy)
	{
		for (jx = 0; jx < d->nx; ++jx)
		{
			switch(d->file_type)
			{
				case BINARY:
					switch(d->data_type)
					{
						case BYTE:    fwrite(&b->statec[jx][jy], sizeof(char), 1, fp); break;
						case INTEGER: fwrite(&b->statei[jx][jy], sizeof(int), 1, fp); break;
						case REAL:    fwrite(&b->statef[jx][jy], sizeof(float), 1, fp); break;
					}
					break;
				case ASCII:
					switch(d->data_type)
					{
						case BYTE:    fprintf(fp,"%d\n", (int)b->statec[jx][jy]); break;
						case INTEGER: fprintf(fp,"%d\n", b->statei[jx][jy]); break;
						case REAL:    
							val = b->statef[jx][jy];
							if (val == (long)val)
								fprintf(fp, "%.0f\n", val);
							else if (val*10 == (long)(val*10))
								fprintf(fp, "%.1f\n", val);
							else if (val*100 == (long)(val*100))
									fprintf(fp, "%.2f\n", val);
							else
								fprintf(fp, "%.4f\n", val);
							break;
					}
			}
		}
	}
	fclose(fp);
	modify_name(filename, docfile, ".doc");
	print_doc(d, docfile);
}

/*********************************/
void   modify_board (BOARD *b, float oldcolor, float newcolor)
/*********************************
* Change cells with value 'oldcolor' to 'newcolor'. It is used
* for example, to recolor a map.
*/
{
	int jx, jy;

	for (jy = 0; jy < b->doc.ny; ++jy)
		for (jx = 0; jx < b->doc.nx; ++jx)
			if (get_cell(jx,jy,b) == oldcolor)
				put_cell(newcolor,jx,jy,b);
}

/*********************************/
void   clear_board (BOARD *b)
/*********************************/
{
	int jx, jy;

	for (jy = 0; jy < b->doc.ny; ++jy)
		for (jx = 0; jx < b->doc.nx; ++jx)
			put_cell(0.0,jx,jy,b);
}

/*********************************/
void   read_board (BOARD *b, char *filename)
/*********************************
* Read raster information from a file
*/
{
	int jx, jy;
	FILE *fp;
	DOC *d;
	int val;

	d = get_doc(b);
	fp = open_board_file (filename, d, "r");
	for (jy = d->ny-1; jy >= 0; --jy)
	{
		for (jx = 0; jx < d->nx; ++jx)
		{
			switch(d->file_type)
			{
				case BINARY:
					switch(d->data_type)
					{
						case BYTE:    fread(&b->statec[jx][jy], sizeof(char), 1, fp); break;
						case INTEGER: fread(&b->statei[jx][jy], sizeof(int), 1, fp); break;
						case REAL:    fread(&b->statef[jx][jy], sizeof(float), 1, fp); break;
					}
					break;
				case ASCII:
					switch(d->data_type)
					{
						case BYTE:    fscanf(fp,"%d", &val); b->statec[jx][jy]=val; break;
						case INTEGER: fscanf(fp,"%d", &b->statei[jx][jy]); break;
						case REAL:    fscanf(fp,"%f", &b->statef[jx][jy]); break;
					}
			}
		}
	}
	fclose(fp);
}

/*********************************/
FILE *open_board_file (char *filename, DOC *d, char *attr)
/*********************************
* Open raster file for sequential reading (it is needed to
* save memory, because the values of all cells are not stored)
*/
{
	char binattr[5];
	FILE *fp;

	strcpy(binattr,attr);
	strcat(binattr,"b");
	switch(d->file_type)
	{
		case BINARY: fp = fopen(filename,binattr); break;
		case ASCII:  fp = fopen(filename,attr); break;
	}
   return(fp);
}

/*********************************/
float  read_value (DOC *d, FILE *fp)
/*********************************
* Reads raster values sequentially
*/
{
	float valf;
	int vali;
	char valc;
	int n;

	switch(d->data_type)
	{
		case BYTE:
			switch(d->file_type)
			{
				case BINARY: fread(&valc, sizeof(char), 1, fp);
					valf = (float)valc; break;
				case ASCII:  fscanf(fp,"%d", &vali);
					valf = (float)vali; break;
			}
			break;
		case INTEGER:
			switch(d->file_type)
			{
				case BINARY: fread(&vali, sizeof(int), 1, fp); break;
				case ASCII:  fscanf(fp,"%d", &vali); break;
			}
			valf = (float)vali;
			break;
		case REAL:
			switch(d->file_type)
			{
				case BINARY: fread(&valf, sizeof(float), 1, fp); break;
				case ASCII:  fscanf(fp,"%f", &valf); break;
			}
			break;
	}
	return(valf);
}

/*********************************/
void  write_value (float valf, DOC *d, FILE *fp)
/*********************************
* Writes raster values to a file sequentially
*/
{
	int vali;
	char valc;

	if (valf >= 0)
		vali = valf + 0.5;
	else
		vali = valf - 0.5;
	valc = vali;
	switch(d->data_type)
	{
		case BYTE:
			switch(d->file_type)
			{
				case BINARY: fwrite(&valc, sizeof(char), 1, fp); break;
				case ASCII:  fprintf(fp,"%d\n", vali); break;
			}
			break;
		case INTEGER:
			switch(d->file_type)
			{
				case BINARY: fwrite(&vali, sizeof(int), 1, fp); break;
				case ASCII:  fprintf(fp,"%d\n", vali); break;
			}
			break;
		case REAL:
			switch(d->file_type)
			{
				case BINARY: fwrite(&valf, sizeof(float), 1, fp); break;
				case ASCII:  fprintf(fp,"%.4f\n", valf); break;
			}
			break;
	}
}

