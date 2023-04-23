#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#define ac_mat(X,i,j,n) X[(i)*(n)+(j)]
#define AC_IMG(img,i,j,c) img[(int)(ch*height*(i)+ch*(j)+(c))]

double *readImage(const char *fname, int *width, int *height, int *ch);
void writeImage(const char *fname, const double *img, const int width, const int height, const int ch);

int main(int argc, char *argv[]){
    int n, iter;
    double WindowX[2], WindowY[2];
    if (argc == 3) {
        n = atoi(argv[1]);
        iter = atoi(argv[2]);
        WindowX[0] = -2;
        WindowX[1] = 0.47;
        WindowY[0] = -1.12;
        WindowY[1] = 1.12;
    } else if (argc == 7) {
        n = atoi(argv[1]);
        iter = atoi(argv[2]);
        WindowX[0] = atof(argv[3]);
        WindowX[1] = atof(argv[4]);
        WindowY[0] = atof(argv[5]);
        WindowY[1] = atof(argv[6]);
    } else {
        printf("Add arguments\n");
        exit(0);
    }
    
    clock_t time1, time2;
    double dub_time;

    const int ch=1;

    const double deltax=(WindowX[1]-WindowX[0])/(double)(n-1);
    const double deltay=(WindowY[1]-WindowY[0])/(double)(n-1);
    double* k=(double*)calloc(n*n,sizeof(double)); //Pot ser un static array
    double* k2=(double*)calloc(n*n,sizeof(double)); 
    time1=clock();
    //Naive code

    //Compute the mandelbrot set
    for (int px = 0; px < n; ++px)
    {   
        double cx=WindowX[0]+px*deltax; //Posar dintre altre for per OPEN*
        for (int py = 0; py < n; ++py)
        {
            
            double cy=WindowY[0]+py*deltay;
            double zx=0.;
            double zy=0.;
            double zz=0.;
            //printf("X:%d Y:%d ",px,py);
            for (int i = 1; i < iter+1; ++i) //To be the same as the matlab code
            {
                double zx2=cx+zx*zx-zy*zy;
                double zy2=cy+zx*zy*2;
                zx=zx2;
                zy=zy2;
                if (zx*zx+zy*zy>4)
                {
                    zz=iter-i;
                    break;

                }
                
            }
            ac_mat(k,py,px,n) = zz;
            
        }
    }
    time2=clock();
    dub_time = (time2 - time1)/(double) CLOCKS_PER_SEC;
    printf("Time for n=%d and iter=%d -----> %lf \n", n,iter,dub_time);

    //Optimized code
    time1=clock();
    //Compute the mandelbrot set
    for (int px = 0; px < n; ++px)
    {   
        double x0=WindowX[0]+px*deltax; //Posar dintre altre for per OPEN*
        for (int py = 0; py < n; ++py)
        {
            
            double y0=WindowY[0]+py*deltay;
            double x2=0;
            double y2=0;
            double x=0;
            double y=0;
            //printf("X:%d Y:%d ",px,py);
            int iteration=0;
            while (x2+y2<=4 && iteration<iter)
            {
                y=(x+x)*y+y0;
                x=x2-y2+x0;
                x2=x*x;
                y2=y*y;
                iteration++;
            }
            
            ac_mat(k2,py,px,n) = iter-iteration; //Transposed, explained in report
            
        }
    }
    time2=clock();
    dub_time = (time2 - time1)/(double) CLOCKS_PER_SEC;
    printf("Time for n=%d and iter=%d optimized -----> %lf \n", n,iter,dub_time);

    

    #ifdef DEBUG
        for (int px = 0; px < n; px++)
        {
            for (int py = 0; py < n; py++)
            {
                printf("%f ", ac_mat(k,px,py,n));
            }
            printf("\n");
        }
        printf("Mat2\n");
        for (int px = 0; px < n; px++)
        {
            for (int py = 0; py < n; py++)
            {
                printf("%f ", ac_mat(k2,px,py,n));
            }
            printf("\n");
        }
    #endif
    writeImage("MandelbrotImageNaive",k,n,n,ch);
    writeImage("MandelbrotImageOptimized",k2,n,n,ch);

    free(k);
    free(k2);
}

////////
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