/*************************************
* Programmer: Alexei Sharov   11/03/2000
* Department of Entomlogy, Virginia Tech
* (540)231-7316, e-mail sharov@vt.edu
* Home page: www.ento.vt.edu/~sharov
* 
* The 'togif.c' program is called from 'stsdecis', 'treat',
* and 'detailmp' PERL scripts
*
* Program GIFLIB.C is a part of the Slow-the-Spread
* Decision-support system. It is used to generate GIF maps
* from a script that specifies map elements.
**************************************/

#include <stdlib.h>
#include "mylib.c"
#include "board.c"
#include "giflib.c"

/* Map objects */
#define OBJ_POINT    1	/* individual point */
#define OBJ_LINE     2	/* individual line */
#define OBJ_POINTF   3	/* file with multiple points */
#define OBJ_LINEF    4	/* file with multiple lines (ARC/INFO ungenerate format) */
#define OBJ_POINTFV  5	/* file with points that are colored according to their values */
#define OBJ_BITMAP   6	/* bitmap object */
#define OBJ_TEXTF    7	/* text objects in a file */
#define OBJ_FILL     8	/* fill flood the map at a particular point */
#define OBJ_SCALE    9	/* draw a space scale */
#define OBJ_TEXT    10	/* one line of text */
#define OBJ_TRANSPAR 11	/* set transparent color */
#define OBJ_POINTFV_TXT 12	/* points with values printed as text */
#define OBJ_GIFIMAGE 13	/* load a gif image */
#define OBJ_BOX      14	/* draw a box */
#define RASTER       1	/* raster input file */
#define VECTOR       2	/* drawing script as input file */
#define RIGHT  1	/* usedfor flood fill */
#define LEFT  -1	/* usedfor flood fill */
#define IDRISI	0	/* IDRISI-formated raster file */
#define ASCIIGRID	1	/* ASCIIGRID format for ARC/INFO */

typedef struct elem_st{	/* Element of stack for fill flood */
	int a, b, c, d;
}ELEM;

typedef struct stack_st{	/* Stack for fill flood */
	ELEM *elem;
	long num;	/* number of elements in the stack */
	long allc;	/* number of memory elements allocated */
}STACK;

typedef struct param_st{
	char rastfile[80];
	int pix, n_chng;
	int drawtype;
	int backgr;
	float thresh[MAXCOL];
	int col[MAXCOL];
	int ncol;
	float xcent, ycent, dx, dy;
	float tmin;
	DOC doc;
}PARAM;

static gdImagePtr image;
static int *palette;
static int XSCREEN =  600;
static int YSCREEN =  500;
static gdFontPtr gdFontGiant;
static gdFontPtr gdFontSmall;
static gdFontPtr gdFontTiny;
static char *gdFontGiantData;
static char *gdFontSmallData;
static char *gdFontTinyData;
static STACK stack;
static int fill_color;

void gdImageChar(gdImagePtr im, gdFontPtr f, int x, int y,
	int c, int color);
void gdImageString(gdImagePtr im, gdFontPtr f,
	int x, int y, unsigned char *s, int color);
void  scale (float x, float y, int *ix, int *iy, PARAM *p);
int gdImageColorClosest(gdImagePtr im, int r, int g, int b);
gdImagePtr gdImageCreateFromGif(FILE *fd);
int gdImageColorExact(gdImagePtr im, int r, int g, int b);
void gdImageStringVertical(gdImagePtr im, gdFontPtr f, int x, int y, unsigned char *s, int color);

/*********************************************/
void  stack_push  (int a, int b, int c, int d)
/*********************************************
* Used for fill flood
*/
{
	stack.elem[stack.num].a = a;
	stack.elem[stack.num].b = b;
	stack.elem[stack.num].c = c;
	stack.elem[stack.num].d = d;
	stack.num++;
	if (stack.num >= stack.allc)
		error_message("Stack too large");
}

/*********************************************/
void  stack_pop  (int *a, int *b, int *c, int *d)
/*********************************************
* Used for fill flood
*/
{
	--stack.num;
	if (stack.num < 0)
		error_message("Empty stack\n");
	*a = stack.elem[stack.num].a;
	*b = stack.elem[stack.num].b;
	*c = stack.elem[stack.num].c;
	*d = stack.elem[stack.num].d;
}

/*********************************************/
int  getbound (int x, int y, int dir, int color)
/*********************************************
* Used for fill flood
*/
{
	int xbeg;

	if (x < 0 || x >= XSCREEN || y < 0 || y >= YSCREEN) {
		printf("%d %d\n", x, y);
		error_message ("getbound");
	}
	xbeg = x;
	while (x >= 0 && x < XSCREEN && gdImageGetPixel(image,x,y) == fill_color)
	{
		if (dir > 0 || x != xbeg)
			gdImageSetPixel(image, x, y, color);
		x += dir;
	}
	if (x != xbeg)
		x -= dir;
	return (x);
}

/*********************************************/
void tilt  (int color)
/*********************************************
* Used for fill flood
*/
{
	int updwn, left, right, row, negupdwn;
	int newl, newr, center, left1, right1;

	stack_pop (&updwn, &left, &right, &row);
	negupdwn = -updwn;
	do{
		if (row < 0 || row >= YSCREEN)
			break;
		center = left;
		while (center <= right && gdImageGetPixel(image,center,row) != fill_color)
			++center;
		if (center > right)
			break;
		newl = getbound (center, row, LEFT, color);
		newr = getbound (center, row, RIGHT, color);

		if (newl < left-1)
			stack_push(negupdwn, newl, left, row-updwn);
		else if (newl > left+1)
			stack_push(updwn, left, newl, row);
		left = newl;

		if (newr > right+1)
			stack_push(negupdwn, right, newr, row-updwn);
		else if (newr < right-1)
			stack_push(updwn, newr, right, row);
		right = newr;

		row += updwn;
	}while (newr - newl >= 0);
}

/*********************************************/
void  fill_flood  (float x0, float y0, int color, PARAM *p)
/*********************************************
* Used for fill flood
*/
{
	int left, right, iter = 0, x, y;

	scale (x0, y0, &x, &y, p);
	if (x < 0 || x >= XSCREEN || y < 0 || y >= YSCREEN)
		return;
	fill_color = gdImageGetPixel(image,x,y);
	if (color == fill_color)
		return;
	left = getbound (x, y, LEFT, color);
	right = getbound (x, y, RIGHT, color);
	if (x >=0 && x < XSCREEN)
	{
		if (y+1 >=0 && y+1<YSCREEN)
			stack_push(1, left, right, y+1);
		if (y-1 >=0 && y-1<YSCREEN)
			stack_push(-1, left, right, y-1);
	}
	while(stack.num > 0 && iter < 1000)
	{
		tilt(color);
		++iter;
	}
}

/********************************/
void circle (int xc, int yc, int r, int color, int fill)
/********************************/
{
	int x, y, dec, iy;
	y = r;
	dec = 3-2*r;
	for (x = 0; x <= y; x++)
	{
		for (iy = -y+1; iy<y && fill >=0; ++iy)
		{
			gdImageSetPixel(image, xc+x, yc+iy, fill);
			gdImageSetPixel(image, xc+iy, yc+x, fill);
			gdImageSetPixel(image, xc-x, yc+iy, fill);
			gdImageSetPixel(image, xc+iy, yc-x, fill);
		}
		gdImageSetPixel(image, xc+x, yc+y, color);
		gdImageSetPixel(image, xc+x, yc-y, color);
		gdImageSetPixel(image, xc-x, yc+y, color);
		gdImageSetPixel(image, xc-x, yc-y, color);
		gdImageSetPixel(image, xc+y, yc+x, color);
		gdImageSetPixel(image, xc+y, yc-x, color);
		gdImageSetPixel(image, xc-y, yc+x, color);
		gdImageSetPixel(image, xc-y, yc-x, color);
		if (r >= 20) {
			gdImageSetPixel(image, xc+x, yc+y-1, color);
			gdImageSetPixel(image, xc+x, yc-y+1, color);
			gdImageSetPixel(image, xc-x, yc+y-1, color);
			gdImageSetPixel(image, xc-x, yc-y+1, color);
			gdImageSetPixel(image, xc+y-1, yc+x, color);
			gdImageSetPixel(image, xc+y-1, yc-x, color);
			gdImageSetPixel(image, xc-y+1, yc+x, color);
			gdImageSetPixel(image, xc-y+1, yc-x, color);
		}
		if ( dec >= 0 )
			dec += -4*(y--)+4;
		dec += 4*x+2;
	}
}

/********************************/
void  scale (float x, float y, int *ix, int *iy, PARAM *p)
/********************************/
{
	*ix = (float)XSCREEN/2 + (x - p->xcent)/p->dx;
	*iy = (float)YSCREEN/2 - (y - p->ycent)/p->dy;
}

/********************************/
void  scale_back (int ix, int iy, float *x, float *y, PARAM *p)
/********************************/
{
	*x = (ix - (float)XSCREEN/2 +0.5)*p->dx + p->xcent;
	*y = ((float)YSCREEN/2 - iy - 0.5)*p->dy + p->ycent;
}

/********************************/
void  draw_point (float x, float y, int color, PARAM *p)
/********************************/
{
	int ix, iy, dx, dy, d0;

	if (color < 0)
		return;
	scale (x, y, &ix, &iy, p);
	ix -= p->pix/2;
	iy -= p->pix/2;
	for(dx=0;dx<p->pix;++dx)
		for(dy=0;dy<p->pix;++dy)
			gdImageSetPixel(image, ix+dx, iy+dy, palette[color]);
}

/********************************/
void  draw_circle (float x, float y, int radius, int color, int fill, PARAM *p)
/********************************/
{
	int ix, iy, dx, dy, d0;

	if (color < 0)
		return;
	scale (x, y, &ix, &iy, p);
	if (radius == 0)
		gdImageSetPixel(image, ix, iy, color);
	else
		circle(ix, iy, radius,  color, fill);
}

/********************************/
void  draw_line (float x, float y, float x1, float y1, PARAM *p,
			int color, int width, int style)
/********************************/
{
	int ix, iy, ix1, iy1;

	scale (x, y, &ix, &iy, p);
	scale (x1, y1, &ix1, &iy1, p);
	gdImageLine(image, ix, iy, ix1, iy1, color, width, style);
}

/***************************************/
void copy_colors (PARAM *p)
/***************************************/
{
	int i;
	p->doc.ncol = p->ncol;
	for(i=0; i<p->ncol; ++i)
	{
		p->doc.col[i] = p->col[i];
		p->doc.thresh[i] = p->thresh[i];
	}
}

/***************************************/
void skip_line (DOC *d, FILE *fp, int raster_type)
/***************************************/
{
	int i, n;
	char buffer[256];

	if (raster_type == ASCIIGRID){
		for(i=0; i<d->nx; ++i)
			fscanf(fp,"%s",buffer);
		return;
	}
	switch(d->file_type)
	{
		case BINARY:
			switch(d->data_type)
			{
				case BYTE: n = d->nx*sizeof(char); break;
				case INTEGER: n = d->nx*sizeof(int); break;
				case REAL: n = d->nx*sizeof(float); break;
			}
			while(n > 0)
			{
				i = 256;
				if (n < i) i = n;
				fread(buffer, 1, i, fp);
				n -= i;
			}
			break;
		case ASCII:
			for(i=0; i<d->nx; ++i)
				fscanf(fp,"%s",buffer);
			break;
	}
}

/***************************************/
void  plot_raster  (PARAM *p)
/***************************************/
{
	FILE *fp;
	int i, color, ix, iy, dupx, dupy, key, mult;
	float x, y, value, dx, dy, ymin, ymax, xmin;
	char docname[100], buffer[100], *c;
	int raster_type;

	if (sscanf(p->rastfile,"%s%d",buffer,&mult)==1)
		mult = 1;
	if (strstr(buffer,".img")){
		raster_type = IDRISI;
		modify_name(buffer, docname, ".doc");
		fp = fopen(buffer, "r");
		if (!fp){
			printf("Image file not found\n");
			return;
		}
		fclose(fp);
		read_doc(&p->doc, docname);
	}
	else if (strstr(buffer,".asc")){
		raster_type = ASCIIGRID;
		fp = fopen(buffer, "r");
		if (!fp){
			printf("ASCIIGRID file not found\n");
			return;
		}
		fscanf(fp, "%s%d", buffer, &p->doc.nx);
		fscanf(fp, "%s%d", buffer, &p->doc.ny);
		fscanf(fp, "%s%f", buffer, &p->doc.xmin);
		fscanf(fp, "%s%f", buffer, &p->doc.ymin);
		fscanf(fp, "%s%f", buffer, &dx);
		fscanf(fp, "%s%f", buffer, &p->tmin);
		p->doc.xmax = p->doc.xmin + dx*p->doc.nx;
		p->doc.ymax = p->doc.ymin + dx*p->doc.ny;
		p->doc.ncol = 0;
	}
	else
		error_message("Unknown raster type");
	if (p->dx <= 0)
	{
		p->xcent = (p->doc.xmax+p->doc.xmin)/2;
		p->ycent = (p->doc.ymax+p->doc.ymin)/2;
		p->dx = (p->doc.xmax - p->doc.xmin)/p->doc.nx;
		p->dy = (p->doc.ymax - p->doc.ymin)/p->doc.ny;
	}
	if (p->doc.ncol <= 0 && p->ncol > 0)
		copy_colors(p);
	dx = (p->doc.xmax-p->doc.xmin)/p->doc.nx;
	dy = (p->doc.ymax-p->doc.ymin)/p->doc.ny;

	printf("dx=%.0f dy=%.0f xmin=%.0f xmax=%.0f ymin=%.0f ymax=%.0f",
	 dx,dy,p->doc.xmin,p->doc.xmax,p->doc.ymin,p->doc.ymax);

	p->pix = dx/(p->dx + 0.001)+1;
	if (dy/p->dy > dx/p->dx)
		p->pix = dy/(p->dy + 0.001)+1;
	dupx = p->dx/dx;
	dupy = p->dy/dy;
	scale_back (0, YSCREEN, &xmin, &ymin, p);
	scale_back (0, 0, &xmin, &ymax, p);
	if (raster_type == IDRISI) {
		fp = open_board_file(buffer, &p->doc, "r");
		if (!fp)
			error_message("img file not found");
	}
	for (iy = p->doc.ny-1; iy >= 0; --iy)
	{
		y = p->doc.ymin + (iy+0.5)*dy;
		while (iy >=0 && (y < ymin || y > ymax || (dupy && iy%dupy)))
		{
			skip_line(&p->doc, fp, raster_type);
			--iy;
			y = p->doc.ymin + (iy+0.5)*dy;
		}
		for (ix = 0; ix < p->doc.nx; ++ix)
		{
			x = p->doc.xmin + (ix+0.5)*dx;
			if (raster_type == IDRISI)
				value = read_value(&p->doc, fp);
			else
				fscanf(fp,"%f",&value);
			if (value >= -1000 && mult != 1)
				value *=mult;
			if (dupx && ix%dupx)
				continue;
			i = 0;
			while (i < p->doc.ncol-1 && value >= p->doc.thresh[i])
				++i;
			draw_point(x, y, p->doc.col[i], p);
		}
	}
	fclose(fp);
}

/********************************/
void  draw_box (float x, float y, float x1, float y1, PARAM *p,
			int color, int fill)
/********************************/
{
	int ix, iy, ix1, iy1, i, j;

	scale (x, y, &ix, &iy, p);
	scale (x1, y1, &ix1, &iy1, p);
	if(color >=0){
		color=palette[color];
		gdImageLine(image, ix, iy, ix, iy1, color, 1, 0);
		gdImageLine(image, ix1, iy, ix1, iy1, color, 1, 0);
		gdImageLine(image, ix, iy, ix1, iy, color, 1, 0);
		gdImageLine(image, ix, iy1, ix1, iy1, color, 1, 0);
	}
	if(fill >= 0){
		fill=palette[fill];
		if(ix1 < ix){
			int swap=ix; ix=ix1; ix1=swap;
		}
		if(iy1 < iy){
			int swap=iy; iy=iy1; iy1=swap;
		}
		for(i=ix+1; i<ix1; ++i){
			for(j=iy+1; j<iy1; ++j){
				gdImageSetPixel(image, i, j, fill);
			}
		}
	}
}


/***************************************/
void  draw_linef(int color, char *objfile, int linewid, int linestyle, PARAM *p)
/***************************************/
{
	float x, y, x1, y1;
	int i, ix, iy;
	FILE *fp;
	long num=0, done = FALSE;
	char buffer[100];

	fp = fopen(objfile,"r");
	if (!fp) {
		printf("file %s not found\n", objfile);
		return;
	}
	while (!done)
	{
		num = 0;
		fgets(buffer, 80, fp);
		while (fscanf(fp,"%f", &x) == 1)
		{
			fscanf(fp,"%f", &y);
			if (num)
				draw_line (x, y, x1, y1,p,palette[color], linewid, linestyle);
			x1 = x;
			y1 = y;
			++num;
		}
		if (!fgets(buffer, 80, fp))
			done = TRUE;
	}
	fclose(fp);
}

/***************************************/
void  draw_pointfv (int radius, char *objfile, PARAM *p)
/***************************************/
{
	float x, y, value;
	int i, ix, iy, color1;
	FILE *fp;
	char buffer[100];

	fp = fopen(objfile,"r");
	if (!fp) {
		printf("file %s not found\n", objfile);
		return;
	}
	while (fscanf(fp,"%f", &x) == 1)
	{
		fscanf(fp,"%f", &y);
		fscanf(fp,"%f", &value);
		fgets(buffer, 80, fp);
		scale (x, y, &ix, &iy, p);
		if (ix < 0 || ix >= XSCREEN || iy < 0 || iy >= YSCREEN)
			continue;

		i = 0;
		while (i < p->doc.ncol-1 && value >= p->doc.thresh[i])
			++i;
		color1 = p->doc.col[i];
		if (color1 < 0)
			continue;
		if (radius >= 0)
			draw_circle(x,y,radius,palette[color1], palette[color1], p);
	}
	fclose(fp);
}

/***************************************/
void  draw_pointfv_txt (PARAM *p, FILE *fp1)
/***************************************/
{
	float x, y, value;
	int i, ix, iy, len;
	FILE *fp;
	char buffer[100];
	float threshold[10];
	int color[10];
	int fontsize[10];
	int n_thresholds;
	gdFontPtr font;

	fscanf(fp1,"%s", buffer);
	fp = fopen(buffer,"r");
	if (!fp) {
		printf("file %s not found\n", buffer);
		return;
	}
	n_thresholds = read_int(fp1);
	if (n_thresholds > 10)
		error_message("Too many thresholds in draw_pointfv_txt");
        for(i=0; i<n_thresholds; ++i){
		threshold[i] = read_float(fp1);
		color[i] = read_int(fp1);
		fontsize[i] = read_int(fp1);
	}
	while (fscanf(fp,"%f", &x) == 1)
	{
		fscanf(fp,"%f", &y);
		fscanf(fp,"%f", &value);
		fgets(buffer, 80, fp);
		scale (x, y, &ix, &iy, p);
		if (ix < 0 || ix >= XSCREEN || iy < 0 || iy >= YSCREEN)
			continue;

		i = 0;
		while (i < n_thresholds-1 && value >= threshold[i])
			++i;
		if (color[i] < 0 || fontsize[i] < 0)
			continue;
		if (fontsize == 0)
			draw_circle(x,y,1,palette[color[i]],-1, p);
		else {
			sprintf(buffer, "%d", (int)value);
			if (fontsize[i] == 1) font = gdFontSmall;
			else font = gdFontGiant;
			len = strlen(buffer);
			gdImageString(image, font, ix-len*font->w/2, iy-font->h/2, buffer, color[i]);
		}
	}
	fclose(fp);
}

/***************************************/
void  draw_pointf(int color, int fill, int radius,char *objfile, PARAM *p)
/***************************************/
{
	float x, y;
	int i, ix, iy, color1, len;
	FILE *fp;
	char buffer[100];

	fp = fopen(objfile,"r");
	if (!fp) {
		printf("file %s not found\n", objfile);
		return;
	}
	while (fscanf(fp,"%f", &x) == 1)
	{
		fscanf(fp,"%f", &y);
		fgets(buffer, 80, fp);
		if (fill >=0)
			fill = palette[fill];
		if (color >= 0)
			draw_circle(x,y,radius,palette[color],fill, p);
	}
	fclose(fp);
}

/***************************************/
void  draw_textf (int color, int size, char *objfile, PARAM *p)
/***************************************/
{
	float x, y;
	int i, ix, iy, color1, len, ix1, iy1, xmin, xmax, ymin, ymax;
	FILE *fp;
	char buffer[100], *ch;
	gdFontPtr f;

	fp = fopen(objfile,"r");
	if (!fp) {
		printf("file %s not found\n", objfile);
		return;
	}
	while (fscanf(fp,"%f", &x) == 1)
	{
		fscanf(fp,"%f", &y);
		fgets(buffer, 99, fp);
		scale (x, y, &ix, &iy, p);
		i = strlen(buffer)-1;
		while (i >= 0 && (buffer[i]=='\n'||buffer[i]==' '||buffer[i]=='\t'))
			--i;
		buffer[i+1] = '\0';
		ch = buffer;
		while (ch-buffer < i && (ch[0]==' '||ch[0]=='\t'))
			++ch;
		len = strlen(ch);
		if (len == 0 || ix < 0 || ix >= XSCREEN || iy < 0 || iy >= YSCREEN)
			continue;
		if (size%3 == 1) f = gdFontSmall;
		else if (size%3 == 2) f = gdFontGiant;
		else f = gdFontTiny;
		if (size >= 3 && size < 6) {
			xmin = ix-len*f->w/2-2;
			xmax = ix+len*f->w/2+1;
			ymin = iy-f->h/2-1;
			ymax = iy+f->h/2+1;
			for (ix1=xmin; ix1<=xmax; ++ix1) {
				for (iy1=ymin; iy1<=ymax; ++iy1) {
					if (ix1==xmin || ix1==xmax || iy1==ymin || iy1==ymax)
						gdImageSetPixel(image, ix1, iy1, color);
					else
						gdImageSetPixel(image, ix1, iy1, 0);
				}
			}
		}
		if(size < 6){
			gdImageString(image, f, ix-len*f->w/2, iy-f->h/2, ch, color);
		}else{
			gdImageStringVertical(image, f, ix-f->h/2, iy, ch, color);
		}
	}
	fclose(fp);
}

/***************************************/
void  draw_text (int color, int size, float x, float y, char *text, PARAM *p)
/***************************************/
{
	int i, ix, iy, color1, len, ix1, iy1, xmin, xmax, ymin, ymax;
	gdFontPtr f;

	scale (x, y, &ix, &iy, p);
	len = strlen(text);
	if (len == 0 || ix < 0 || ix >= XSCREEN || iy < 0 || iy >= YSCREEN)
		return;
	if (size%3 == 1) f = gdFontSmall;
	else if (size%3 == 2) f = gdFontGiant;
	else f = gdFontTiny;
	if (size >= 3 && size < 6) {
		xmin = ix-len*f->w/2-4;
		xmax = ix+len*f->w/2+1;
		ymin = iy-f->h/2-1;
		ymax = iy+f->h/2+1;
		for (ix1=xmin; ix1<=xmax; ++ix1)
			for (iy1=ymin; iy1<=ymax; ++iy1)
				gdImageSetPixel(image, ix1, iy1, 0);
		for (ix1=xmin; ix1<=xmax; ++ix1){
			gdImageSetPixel(image, ix1, ymin, color);
			gdImageSetPixel(image, ix1, ymax, color);
		}
		for (iy1=ymin; iy1<=ymax; ++iy1){
			gdImageSetPixel(image, xmin, iy1, color);
			gdImageSetPixel(image, xmax, iy1, color);
		}
	}
	if(size < 6){
		gdImageString(image, f, ix-len*f->w/2, iy-f->h/2, text, color);
	}else{
		gdImageStringVertical(image, f, ix-f->h/2, iy, text, color);
	}
}

/***************************************/
void  draw_simple (int obj_type, int color, int r,
		 int linewid, int linestyle, FILE *fp, PARAM *p)
/***************************************/
{
	float x, y, x1, y1, dx, dy;
	int fill;

	if (obj_type == OBJ_POINT || obj_type == OBJ_BOX)
	{
		fill = read_int(fp);
		if (fill >=0)
			fill = palette[fill];
	}
	x = read_float(fp);
	y = read_float(fp);
	if (obj_type == OBJ_BOX || obj_type == OBJ_LINE){
		x1 = read_float(fp);
		y1 = read_float(fp);
	}
	if (obj_type == OBJ_POINT)
		draw_circle (x, y, r, color, fill, p);
	else if (obj_type == OBJ_LINE)
		draw_line (x, y, x1, y1,p,palette[color], linewid, linestyle);
	else if (obj_type == OBJ_BOX)
		draw_box (x, y, x1, y1,p,color,fill);
}

/***************************************/
void   draw_gifimage(PARAM *p, float x, float y, float resolution,
	 char *filename)
/***************************************
* x,y = coords of the center of image, resolution=m/pixel
*/
{
	FILE *fp;
	gdImagePtr imin;
	int ix, iy, xmin, xmax, ymin, ymax;
	int ix1, iy1, i, c, nc;
	float dx, dy, x1, y1;

	int colorMap[gdMaxColors];
	for (i=0; (i<gdMaxColors); i++) {
		colorMap[i] = (-1);
	}

	fp = fopen(filename,"rb");
	if (!fp){
		printf("Gif file %s not found!\n", filename);
		return;
	}
	imin = gdImageCreateFromGif(fp);
	fclose(fp);
	dx = imin->sx*resolution;
	dy = imin->sy*resolution;
	scale (x-dx/2, y+dy/2, &xmin, &ymin, p);
	scale (x+dx/2, y-dy/2, &xmax, &ymax, p);
	if (xmin<0) xmin = 0;
	if (ymin<0) ymin = 0;
	if (xmax>XSCREEN-1) xmax = XSCREEN-1;
	if (ymax>YSCREEN-1) ymax = YSCREEN-1;
	for(ix=xmin; ix<xmax; ++ix){ 
		for(iy=ymin; iy<ymax; ++iy){ 
			scale_back(ix, iy, &x1, &y1, p);
			ix1 = (float)imin->sx/2 + (x1 - x)/resolution;
			iy1 = (float)imin->sy/2 - (y1 - y)/resolution;
			if (!gdImageBoundsSafe(imin, ix1, iy1))
				continue;
			c = gdImageGetPixel(imin,ix1,iy1);

			/* Have we established a mapping for this color? */
			if (colorMap[c] == (-1)) {
				/* First look for an exact match */
				nc = gdImageColorExact(image,imin->red[c], 
					imin->green[c],	imin->blue[c]);
				if (nc == (-1)) {
					/* No, so try to allocate it */
					nc = gdImageColorAllocate(image, imin->red[c],
					 imin->green[c], imin->blue[c]);
					/* If we're out of colors, go for the closest color */
					if (nc == (-1)) {
						nc = gdImageColorClosest(image,imin->red[c],
						 imin->green[c],imin->blue[c]);
					}
				}
				colorMap[c] = nc;
			}
			gdImageSetPixel(image, ix, iy, colorMap[c]);
		}
	}
	gdImageDestroy(imin);
}

/***************************************/
void  drawing  (char *filename, PARAM *p)
/***************************************/
{
	int nobj, obj_type, nx, ny, radius;
	char objfile[81], buffer[90];
	int color, fill, linewid, linestyle, i;
	float xcent, ycent, dx, dy, x, y, resolution;
	FILE *fp;

	fp = fopen(filename,"r");
	if (!fp)
	{
		printf("Drawing file not found!\n");
		return;
	}
	dx = read_float(fp);
	fgets(buffer,80,fp);
	if (sscanf(buffer,"%f",&dy) != 1)
		dy = dx;
	xcent = read_float(fp);
	ycent = read_float(fp);
	if (p->dx <=0)
	{
		p->dx = dx;
		p->dy = dy;
		p->xcent = xcent;
		p->ycent = ycent;
		p->doc.xmin = xcent-dx*((float)XSCREEN/2);
		p->doc.ymin = ycent-dy*((float)YSCREEN/2);
		p->doc.xmax = p->doc.xmin + dx*XSCREEN;
		p->doc.ymax = p->doc.ymin + dy*YSCREEN;
	}
	nobj = read_int(fp);
	for (i=0; i<nobj; ++i)
	{
		skip_comments(fp);
		if (fscanf(fp,"%d", &obj_type) != 1)
			break;
		if (obj_type == OBJ_LINE || obj_type == OBJ_LINEF || obj_type == OBJ_FILL
		 || obj_type == OBJ_TEXTF || obj_type == OBJ_POINT || obj_type == OBJ_POINTF
		 || obj_type == OBJ_TEXT || obj_type == OBJ_TRANSPAR || obj_type == OBJ_BOX)
			color = read_int(fp);
		else
			color = 0;

		if (obj_type == OBJ_POINTFV && p->doc.ncol <= 0 && p->ncol > 0)
			copy_colors(p);

		if (obj_type == OBJ_LINE || obj_type == OBJ_LINEF)
		{
			linewid = read_int(fp);
			linestyle = read_int(fp);
		}

		if (obj_type == OBJ_POINT || obj_type == OBJ_POINTF || obj_type == OBJ_POINTFV 
		 || obj_type == OBJ_TEXTF || obj_type == OBJ_TEXT)
		{
			radius = read_int(fp); /* or size of font */
		}

		if (obj_type == OBJ_TEXT || obj_type == OBJ_GIFIMAGE){
			x = read_float(fp);
			y = read_float(fp);
		}
		if (obj_type == OBJ_GIFIMAGE){
			resolution = read_float(fp);
		}

		if (obj_type == OBJ_LINEF || obj_type == OBJ_POINTFV || obj_type == OBJ_GIFIMAGE 
		  || obj_type == OBJ_BITMAP || obj_type == OBJ_TEXTF  || obj_type == OBJ_TEXT)
			read_string(objfile, 80, fp);

		switch(obj_type)
		{
			case OBJ_TRANSPAR:
				image->transparent = color;
				break;
			case OBJ_POINTF:
				fill = read_int(fp);
				read_string(objfile, 80, fp);
				draw_pointf(color, fill, radius, objfile, p);
				break;
			case OBJ_LINEF:
				draw_linef(color, objfile, linewid, linestyle, p);
				break;
			case OBJ_POINTFV:
				draw_pointfv(radius, objfile, p);
				break;
			case OBJ_TEXTF:
				draw_textf(color, radius, objfile, p);
				break;
			case OBJ_TEXT:
				draw_text(color, radius, x, y, objfile, p);
				break;
			case OBJ_BITMAP:
				strcpy(p->rastfile, objfile);
				plot_raster(p);
				break;
			case OBJ_FILL:
				xcent = read_float(fp);
				ycent = read_float(fp);
				check(stack.elem = (ELEM*)malloc(5000*sizeof(ELEM)));
				stack.allc=5000;
				fill_flood(xcent, ycent, color, p);
				free(stack.elem);
				break;
			case OBJ_POINTFV_TXT:
				draw_pointfv_txt(p, fp);
				break;
			case OBJ_GIFIMAGE:
				draw_gifimage(p, x, y, resolution, objfile);
				break;
			default:
				draw_simple(obj_type, color, radius, linewid, linestyle, fp,p);
				break;
		}
	}
	fclose(fp);
}

/***************************************/
int  get_drawtype (char *filename)
/***************************************/
{
	char buffer[80], *ps;

	strcpy(buffer, filename);
	ps = strrchr(buffer,'.');
	ps++;
	if (ps == NULL) exit(0);
	if (strcmp(ps,"img")==0)
		return(RASTER);
	else if (!strcmp(ps,"dra")||!strcmp(ps,"txt"))
		return(VECTOR);
	else
		exit(0);
	return(-1);
}

/***************************************/
void  set_palette  (char *palfile)
/***************************************/
{
	FILE *fp;
	int r, g, b, i, n;

	check(palette = (int*)calloc(256,sizeof(int)));
	if (!palfile)
	{
		palette[0] = gdImageColorAllocate(image, 0, 0, 0);
		palette[1] = gdImageColorAllocate(image, 0, 0, 127);
		palette[2] = gdImageColorAllocate(image, 0, 127, 0);
		palette[3] = gdImageColorAllocate(image, 0, 127, 127);
		palette[4] = gdImageColorAllocate(image, 127, 0, 0);
		palette[5] = gdImageColorAllocate(image, 127, 0, 127);
		palette[6] = gdImageColorAllocate(image, 127, 63, 0);
		palette[7] = gdImageColorAllocate(image, 127, 127, 127);
		palette[8] = gdImageColorAllocate(image, 63, 63, 63);
		palette[9] = gdImageColorAllocate(image, 63, 63, 255);
		palette[10] = gdImageColorAllocate(image, 63, 255, 63);
		palette[11] = gdImageColorAllocate(image, 63, 255, 255);
		palette[12] = gdImageColorAllocate(image, 255, 63, 63);
		palette[13] = gdImageColorAllocate(image, 255, 63, 255);
		palette[14] = gdImageColorAllocate(image, 255, 255, 63);
		palette[15] = gdImageColorAllocate(image, 255, 255, 255);
	}
	else
	{
		fp = fopen(palfile, "r");
		if (!fp)
			error_message("palette not found");
		n = read_int(fp);
		if (n > 256)
			error_message("Maximum 256 colors allowed");
		for(i=0; i<n; ++i)
		{
			fscanf(fp,"%d%d%d",&r,&g,&b);
			palette[i] = gdImageColorAllocate(image, r, g, b);
		}
		fclose(fp);
	}
}

/*********************************/
PARAM *read_parameters (char *filename)
/*********************************/
{
	PARAM *p;
	FILE *fp;
	int i;

	check(p=(PARAM*)calloc(1,sizeof(PARAM)));
	fp = fopen(filename,"r");
	if (!fp)
	{
		printf("Param. file not found\n");
		return(p);
	}
	p->pix = read_int(fp);
	p->backgr = read_int(fp);
	p->ncol = read_int(fp);
	if (p->ncol <= 0)
		p->ncol = 0;
	if (p->ncol >= MAXCOL)
		p->ncol = MAXCOL-1;
	for (i=0; i<p->ncol; ++i)
	{
		p->thresh[i] = read_float(fp);
		p->col[i] = read_int(fp);
	}
	fclose(fp);
	p->doc.nx = XSCREEN;
	p->doc.ny = YSCREEN;
	p->doc.data_type = INTEGER;
	p->doc.file_type = ASCII;
	return(p);
}

/***************************************/
void initialize_fonts (char *filename)
/***************************************/
{
	FILE *fp;
	char fontfile[80], *ptr;

	check(gdFontTinyData=(char*)malloc(sizeof(char)*256*8*5));
	check(gdFontSmallData=(char*)malloc(sizeof(char)*19968));
	check(gdFontGiantData=(char*)malloc(sizeof(char)*34560L));
	check(gdFontTiny=(gdFontPtr)malloc(sizeof(gdFont)));
	check(gdFontSmall=(gdFontPtr)malloc(sizeof(gdFont)));
	check(gdFontGiant=(gdFontPtr)malloc(sizeof(gdFont)));

	/* Look for the file with fonts 'fonts.bin' in the same */
	/* directory where is the program 'togif' */
	strcpy(fontfile, filename);
	if (!(ptr = strstr(fontfile, "TOGIF")))
		ptr = strstr(fontfile, "togif");
	if (!ptr)
		error_message("program name");
	strcpy(ptr,"fonts.bin");
	fp=fopen(fontfile,"rb");
	if (!fp)
		error_message("font file not found");
	fread(gdFontTinyData, sizeof(char), 256*8*5, fp);
	fread(gdFontSmallData, sizeof(char), 19968, fp);
	fread(gdFontGiantData, sizeof(char), 34560L, fp);
	fclose (fp);
	gdFontTiny->nchars=256;
	gdFontTiny->offset=0;
	gdFontTiny->w=5;
	gdFontTiny->h=8;
	gdFontTiny->data=gdFontTinyData;
	gdFontSmall->nchars=256;
	gdFontSmall->offset=0;
	gdFontSmall->w=6;
	gdFontSmall->h=13;
	gdFontSmall->data=gdFontSmallData;
	gdFontGiant->nchars=256;
	gdFontGiant->offset=0;
	gdFontGiant->w=9;
	gdFontGiant->h=15;
	gdFontGiant->data=gdFontGiantData;
}

/***************************************/
void gdImageChar(gdImagePtr im, gdFontPtr f, int x, int y, 
	int c, int color)
/***************************************/
{
	int cx, cy;
	int px, py;
	int fline;
	cx = 0;
	cy = 0;
	if ((c < f->offset) || (c >= (f->offset + f->nchars))) {
		return;
	}
	fline = (c - f->offset) * f->h * f->w;
	for (py = y; (py < (y + f->h)); py++) {
		for (px = x; (px < (x + f->w)); px++) {
			if (f->data[fline + cy * f->w + cx]) {
				gdImageSetPixel(im, px, py, color);	
			}
			cx++;
		}
		cx = 0;
		cy++;
	}
}

/***************************************/
void gdImageString(gdImagePtr im, gdFontPtr f, 
	int x, int y, unsigned char *s, int color)
/***************************************/
{
	int i;
	int l;

	l = strlen(s);

	for (i=0; (i<l); i++) {
		gdImageChar(im, f, x, y, s[i], color);
		x += f->w;
	}
}

/***************************************/
void gdImageCharVertical(gdImagePtr im, gdFontPtr f, int x, int y, 
	int c, int color)
/***************************************/
{
	int cx, cy;
	int px, py;
	int fline;
	cx = 0;
	cy = 0;
	if ((c < f->offset) || (c >= (f->offset + f->nchars))) {
		return;
	}
	fline = (c - f->offset) * f->h * f->w;
	for (py = y-f->w; py < y; py++) {
		for (px = x; (px < (x + f->h)); px++) {
			if (f->data[fline + cx * f->w + (f->w-cy-1)]) {
				gdImageSetPixel(im, px, py, color);	
			}
			cx++;
		}
		cx = 0;
		cy++;
	}
}

/***************************************/
void gdImageStringVertical(gdImagePtr im, gdFontPtr f, 
	int x, int y, unsigned char *s, int color)
/***************************************/
{
	int i;
	int l;

	l = strlen(s);

	for (i=0; (i<l); i++) {
		gdImageCharVertical(im, f, x, y, s[i], color);
		y -= f->w;
	}
}

/***************************************/
void   writeImage(gdImagePtr im, char *filename, PARAM *p)
/***************************************/
{
	FILE *fp;
	char docfile[80];
        int ix, iy;

        fp=fopen(filename,"w");
	modify_name(filename, docfile, ".doc");
	print_doc(&p->doc, docfile);
	for(iy=0; iy<YSCREEN; ++iy)
        	for(ix=0;ix<XSCREEN; ++ix)
			fprintf(fp,"%d\n", gdImageGetPixel(im, ix, iy));
        fclose(fp);
}

static int
ReadColorMap(FILE *fd, int number, unsigned char (*buffer)[256])
{
       int             i;
       unsigned char   rgb[3];


       for (i = 0; i < number; ++i) {
               if (! ReadOK(fd, rgb, sizeof(rgb))) {
                       return TRUE;
               }
               buffer[CM_RED][i] = rgb[0] ;
               buffer[CM_GREEN][i] = rgb[1] ;
               buffer[CM_BLUE][i] = rgb[2] ;
       }


       return FALSE;
}

static void
ReadImage(gdImagePtr im, FILE *fd, int len, int height, unsigned char (*cmap)[256], int interlace, int ignore)
{
       unsigned char   c;      
       int             v;
       int             xpos = 0, ypos = 0, pass = 0;
       int i;
       /* Stash the color map into the image */
       for (i=0; (i<gdMaxColors); i++) {
               im->red[i] = cmap[CM_RED][i];	
               im->green[i] = cmap[CM_GREEN][i];	
               im->blue[i] = cmap[CM_BLUE][i];	
               im->open[i] = 1;
       }
       /* Many (perhaps most) of these colors will remain marked open. */
       im->colorsTotal = gdMaxColors;
       /*
       **  Initialize the Compression routines
       */
       if (! ReadOK(fd,&c,1)) {
               return; 
       }
       if (LWZReadByte(fd, TRUE, c) < 0) {
               return;
       }

       /*
       **  If this is an "uninteresting picture" ignore it.
       */
       if (ignore) {
               while (LWZReadByte(fd, FALSE, c) >= 0)
                       ;
               return;
       }

       while ((v = LWZReadByte(fd,FALSE,c)) >= 0 ) {
               /* This how we recognize which colors are actually used. */
               if (im->open[v]) {
                       im->open[v] = 0;
               }
               gdImageSetPixel(im, xpos, ypos, v);
               ++xpos;
               if (xpos == len) {
                       xpos = 0;
                       if (interlace) {
                               switch (pass) {
                               case 0:
                               case 1:
                                       ypos += 8; break;
                               case 2:
                                       ypos += 4; break;
                               case 3:
                                       ypos += 2; break;
                               }

                               if (ypos >= height) {
                                       ++pass;
                                       switch (pass) {
                                       case 1:
                                               ypos = 4; break;
                                       case 2:
                                               ypos = 2; break;
                                       case 3:
                                               ypos = 1; break;
                                       default:
                                               goto fini;
                                       }
                               }
                       } else {
                               ++ypos;
                       }
               }
               if (ypos >= height)
                       break;
       }

fini:
       if (LWZReadByte(fd,FALSE,c)>=0) {
               /* Ignore extra */
       }
}

gdImagePtr
gdImageCreateFromGif(FILE *fd)
{
       int imageNumber;
       int BitPixel;
       int ColorResolution;
       int Background;
       int AspectRatio;
       int Transparent = (-1);
       unsigned char   buf[16];
       unsigned char   c;
       unsigned char   ColorMap[3][MAXCOLORMAPSIZE];
       unsigned char   localColorMap[3][MAXCOLORMAPSIZE];
       int             imw, imh;
       int             useGlobalColormap;
       int             bitPixel;
       int             imageCount = 0;
       char            version[4];
       gdImagePtr im = 0;
       ZeroDataBlock = FALSE;

       imageNumber = 1;
       if (! ReadOK(fd,buf,6)) {
		return 0;
	}
       if (strncmp((char *)buf,"GIF",3) != 0) {
		return 0;
	}
       strncpy(version, (char *)buf + 3, 3);
       version[3] = '\0';

       if ((strcmp(version, "87a") != 0) && (strcmp(version, "89a") != 0)) {
		return 0;
	}
       if (! ReadOK(fd,buf,7)) {
		return 0;
	}
       BitPixel        = 2<<(buf[4]&0x07);
//       ColorResolution = (int) (((buf[4]&0x70)>>3)+1);
//       Background      = buf[5];
//       AspectRatio     = buf[6];

       if (BitSet(buf[4], LOCALCOLORMAP)) {    /* Global Colormap */
               if (ReadColorMap(fd, BitPixel, ColorMap)) {
			return 0;
		}
       }
       for (;;) {
               if (! ReadOK(fd,&c,1)) {
                       return 0;
               }
               if (c == ';') {         /* GIF terminator */
                       int i;
                       if (imageCount < imageNumber) {
                               return 0;
                       }
                       /* Terminator before any image was declared! */
                       if (!im) {
                              return 0;
                       }
		       /* Check for open colors at the end, so
                          we can reduce colorsTotal and ultimately
                          BitsPerPixel */
                       for (i=((im->colorsTotal-1)); (i>=0); i--) {
                               if (im->open[i]) {
                                       im->colorsTotal--;
                               } else {
                                       break;
                               }
                       } 
                       return im;
               }

               if (c == '!') {         /* Extension */
                       if (! ReadOK(fd,&c,1)) {
                               return 0;
                       }
                       DoExtension(fd, c, &Transparent);
                       continue;
               }

               if (c != ',') {         /* Not a valid start character */
                       continue;
               }

               ++imageCount;

               if (! ReadOK(fd,buf,9)) {
	               return 0;
               }

               useGlobalColormap = ! BitSet(buf[8], LOCALCOLORMAP);

               bitPixel = 1<<((buf[8]&0x07)+1);

               imw = LM_to_uint(buf[4],buf[5]);
               imh = LM_to_uint(buf[6],buf[7]);
	       if (!(im = gdImageCreate(imw, imh))) {
			 return 0;
	       }
               im->interlace = BitSet(buf[8], INTERLACE);
               if (! useGlobalColormap) {
                       if (ReadColorMap(fd, bitPixel, localColorMap)) { 
                                 return 0;
                       }
                       ReadImage(im, fd, imw, imh, localColorMap, 
                                 BitSet(buf[8], INTERLACE), 
                                 imageCount != imageNumber);
               } else {
                       ReadImage(im, fd, imw, imh,
                                 ColorMap, 
                                 BitSet(buf[8], INTERLACE), 
                                 imageCount != imageNumber);
               }
               if (Transparent != (-1)) {
                       gdImageColorTransparent(im, Transparent);
               }	   
       }
}

int gdImageColorClosest(gdImagePtr im, int r, int g, int b)
{
	int i;
	long rd, gd, bd;
	int ct = (-1);
	long mindist = 0;
	for (i=0; (i<(im->colorsTotal)); i++) {
		long dist;
		if (im->open[i]) {
			continue;
		}
		rd = (im->red[i] - r);	
		gd = (im->green[i] - g);
		bd = (im->blue[i] - b);
		dist = rd * rd + gd * gd + bd * bd;
		if ((i == 0) || (dist < mindist)) {
			mindist = dist;	
			ct = i;
		}
	}
	return ct;
}

int gdImageColorExact(gdImagePtr im, int r, int g, int b)
{
	int i;
	for (i=0; (i<(im->colorsTotal)); i++) {
		if (im->open[i]) {
			continue;
		}
		if ((im->red[i] == r) && 
			(im->green[i] == g) &&
			(im->blue[i] == b)) {
			return i;
		}
	}
	return -1;
}

/***************************************/
int  main   (int argc, char *argv[])
/***************************************/
{
	int ix, iy;
	FILE *fout;
	PARAM *p;
	DOC *d;

	if (argc < 4)
		error_message("togif infile outfile paramfile [palette [nx ny]]");
	if (argc >= 7)
	{
		XSCREEN = atoi(argv[5]);
		YSCREEN = atoi(argv[6]);
		if (XSCREEN <= 0 || XSCREEN > 5000 || YSCREEN <= 0 || YSCREEN > 5000)
			error_message ("Wrong image size");
	}
	p = read_parameters(argv[3]);
	p->drawtype = get_drawtype(argv[1]);
	image = gdImageCreate(XSCREEN, YSCREEN);

	initialize_fonts(argv[0]);
	if (argc >= 5)
		set_palette(argv[4]);
	else
		set_palette(NULL);

	/* draw background */
	for(ix=0; ix < XSCREEN; ++ix)
		for(iy=0; iy < YSCREEN; ++iy)
			gdImageSetPixel(image, ix, iy, palette[p->backgr]);

	/* draw map elements */
	switch (p->drawtype)
	{
		case RASTER:
			strcpy(p->rastfile, argv[1]);
			plot_raster (p);
			break;
		default:
			drawing (argv[1], p); /* read the script and draw */
			break;
	}
	/* write output */
	if (strstr(argv[2],".gif")){
		fout = fopen(argv[2], "wb");
		if(!fout) error_message("Cannot write to gif file");
		gdImageGif(image, fout);
		fclose(fout);
	}
	else if (strstr(argv[2],".img")) {
		writeImage(image, argv[2], p);
	}
	else
		error_message("Wrong extension in output file name");
	gdImageDestroy(image);
	return(0);
}


