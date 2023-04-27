#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <openacc.h>


#define AC_IMG(img,i,j,c) img[(int)(ch*height*(i)+ch*(j)+(c))]
#define ac_mat(X,i,j,n) X[(i)*(n)+(j)]

double *readImage(const char *fname, int *width, int *height, int *ch);
void writeImage(const char *fname, const double *img, const int width, const int height, const int ch);

__global__ void MandelbrotGPU(int n, double deltax, double deltay,double WindowX0, double WindowY0, int iter,double *k);

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
    
    double time1, time2;
    double dub_time;
    const int ch=1;
    const double deltax=(WindowX[1]-WindowX[0])/(double)(n-1);
    const double deltay=(WindowY[1]-WindowY[0])/(double)(n-1);
    double* k=(double*)calloc(n*n,sizeof(double)); //Pot ser un static array? No, memoria estatica petita

    //Mandelcuda
    double *kGPU;
    cudaMalloc((void **)&kGPU,n*n*sizeof(double));
    int ntx = 16;
    int nty = 16;
    int nBlocksX = (n+ntx-1)/ntx;
    int nBlocksY = (n+nty-1)/nty;
    time1=omp_get_wtime();
    dim3 block2(ntx, nty);
    dim3 grid2(nBlocksX, nBlocksY);
    MandelbrotGPU<<<grid2,block2>>>(n,deltax,deltay,WindowX[0],WindowY[0],iter,kGPU);
    time2=omp_get_wtime();
    dub_time = time2 - time1 ;
    printf("Time for n=%d and iter=%d CUDA -----> %lf \n", n,iter,dub_time);
    cudaMemcpy(k,kGPU,n*n*sizeof(double),cudaMemcpyDeviceToHost);
    //Compute the mandelbrot set 

    #ifdef DEBUG
        printf("GPU:\n");
        for (int px = 0; px < n; px++)
        {
            for (int py = 0; py < n; py++)
            {
                printf("%f ", ac_mat(k,px,py,n));
            }
            printf("\n");
        }
    #endif
    #ifdef WRITE
        writeImage("MandelbrotImage",k,n,n,ch);
    #endif
    free(k);
    cudaFree(kGPU);
}


__global__ void MandelbrotGPU(int n, double deltax, double deltay,double WindowX0, double WindowY0, int iter,double *k){
    int px = blockIdx.x * blockDim.x + threadIdx.x;
    int py = blockIdx.y * blockDim.y + threadIdx.y;
    if (px < n && py < n) {
            double x0=WindowX0+px*deltax;
            double y0=WindowY0+py*deltay;
            double x2=0;
            double y2=0;
            double x=0;
            double y=0;
            int iteration=0;
            while (x2+y2<=4 && iteration<iter)
            {
                y=(x+x)*y+y0;
                x=x2-y2+x0;
                x2=x*x;
                y2=y*y;
                iteration++;
            }
            
            ac_mat(k,py,px,n) = iter-iteration;
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
