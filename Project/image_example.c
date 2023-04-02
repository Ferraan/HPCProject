/*
	READ and WRITE images easily from C/C++ to 
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define AC_IMG(img,i,j,c) img[(int)(ch*height*(i)+ch*(j)+(c))]


double *readImage(const char *fname, int *width, int *height, int *ch) {
/*
	This subroutine reads an image stored as a double pointer
	into a binary file. The image can be either multi channel (RGB)
	or single channel (B&W). It also needs to read an info file with the
	dimensions of the image. This is a normal text file.

	Inputs:
		> fname:  filename as a C string
		> with:   with dimensions
		> height: height dimensions
		> ch:     number of channels (1: B&W, 3: RGB)

	Outputs:
		> img:    image arrays already allocated by this subroutine
*/
	FILE *myfile;
	char aux[256];
	// Read info file
	sprintf(aux,"%s.info",fname); myfile = fopen(aux,"r");
	char *line = NULL;
	size_t len = 0;
	// Read header line and discard it
	getline(&line,&len,myfile); 
	// Read image metadata
	getline(&line,&len,myfile); *ch     = atoi(line); // Number of channels 
	getline(&line,&len,myfile); *width  = atoi(line); // Width
	getline(&line,&len,myfile); *height = atoi(line); // Height 
	fclose(myfile);
	// Read image data
	size_t size = (*ch)*(*width)*(*height);
	sprintf(aux,"%s.bin",fname); myfile = fopen(aux,"rb");
	// Read binary data
	double *img;
	img = (double*)malloc(size*sizeof(double));
	fread(img,sizeof(double),size,myfile);
	fclose(myfile);
	// Return image
	return img;
}


void writeImage(const char *fname, const double *img, 
	const int width, const int height, const int ch) {
/*
	This subroutine writes an image stored as a double pointer
	into a binary file. The image can be either multi channel (RGB)
	or single channel (B&W). It also saves an info file with the
	dimensions of the image. This is a normal text file.

	Inputs:
		> fname:  filename as a C string
		> img:    pointer to the image
		> with:   with dimensions
		> height: height dimensions
		> ch:     number of channels (1: B&W, 3: RGB)

	Outputs:
		writes fname.info and fname.bin to the disk
		they can be visualized using MATLAB or Python
*/
	FILE *myfile;
	char aux[256];
	// Write info file
	sprintf(aux,"%s.info",fname); myfile = fopen(aux,"w");
	fprintf(myfile,"# Image information - channels, width, height\n"); // Header
	fprintf(myfile,"%d\n",ch);     // Number of channels
	fprintf(myfile,"%d\n",width);  // Width
	fprintf(myfile,"%d\n",height); // Height
	fclose(myfile);
	// Write binary file
	sprintf(aux,"%s.bin",fname); myfile = fopen(aux,"wb");
	fwrite(img,sizeof(double),ch*height*width,myfile);
	fclose(myfile);
}


// Simple example how to generate an image
int main(int argc, char *argv[]) {

	const double width=800, height=600, ch=3; // A 800x600 monochrome image
	double *img; 

	// Allocate memory
	img = (double*)malloc(ch*width*height*sizeof(double));

	// Fill the image with random values
	for (int ii=0;ii<ch*width*height;++ii)
		img[ii] = drand48();

	// Store on disk
	writeImage("test",img,width,height,ch);

	return 0;
}