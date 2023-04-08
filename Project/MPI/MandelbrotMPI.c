#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "mpi.h"
#include "Image.h"

#define AC_IMG(img,i,j,c) img[(int)(ch*height*(i)+ch*(j)+(c))]
#define ac_mat(X,i,j,n) X[(i)*(n)+(j)]

double *readImage(const char *fname, int *width, int *height, int *ch);
void writeImage(const char *fname, const double *img, const int width, const int height, const int ch);
void worksplit(int *mystart, int *myend, int proc, int nproc, int start, int end);
void checkr(int r,char *txt);

int main(int argc, char *argv[]){
    //MPI
    int myrank; /* from 0 to P-1 */
    int P; /* number of processes */
    char name[MPI_MAX_PROCESSOR_NAME+1]; /* processor name */
    int r,l; /* detalls */
    r=MPI_Init(&argc,&argv); checkr(r,"init");
    r=MPI_Comm_rank(MPI_COMM_WORLD,&myrank); checkr(r,"rank");
    r=MPI_Comm_size(MPI_COMM_WORLD,&P); checkr(r,"size");
    r=MPI_Get_processor_name(name,&l); checkr(r,"name");
    
    const int n=10000;  //Up to n=40000 using 16Gb RAM and a little bit of swap
    const int iter=1000;
    const int ch=1;
    const double WindowX[2]={-2,0.47};
    const double WindowY[2]={-1.12, 1.12};
    const double deltax=(WindowX[1]-WindowX[0])/(double)(n-1);
    const double deltay=(WindowY[1]-WindowY[0])/(double)(n-1);
    
    double starttime, endtime;
    int *mystart, *myend;
    mystart = (int*)malloc(sizeof(int));
	myend = (int*)malloc(sizeof(int));
    double* coordX=(double*)calloc(n,sizeof(double)); //Can be arrays, not changing potser no, el heap? es petit
    double* coordY=(double*)calloc(n,sizeof(double));
    starttime = MPI_Wtime();
    
    worksplit(mystart,myend,myrank,P,0,n*n-1);
    int ntasks=*myend-*mystart+1;
    #ifdef DEBUG
        printf("ntasks=%d, *mystart=%d, *myend=%d\n",ntasks,*mystart,*myend);
    #endif
    double* kLOC=(double*)calloc(ntasks,sizeof(double)); //Pot ser un static array
    int posk=0;
    for(int i = *mystart; i <= *myend; i++){
		//printf("rank: %d, task: %d\n",myrank, i);
        int posx = i/n;
        int posy = i % n;
        double x0=WindowX[0]+posx*deltax;
        double y0=WindowY[0]+posy*deltay;        
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
            
        kLOC[posk] = iter-iteration; //Transposed, explained in report
        //printf("cx:%f  cy:%f  rank:%d \n",x0,y0,myrank);
        //printf("k[%d,%d]=%f  rank:%d\n",posy,posx,ac_mat(k,posy,posx,n),myrank);
        posk++;
	}
    endtime   = MPI_Wtime();
    printf("Time for rank:%d before gather ----> %f seconds\n",myrank,endtime-starttime);

    double* kGLO=NULL; //Initialized for all variables
    int* tasksGLO=(int*)malloc(P*sizeof(int));
    int* displ=(int*)malloc(P*sizeof(int));
    
    r=MPI_Gather(&ntasks,1,MPI_INT,tasksGLO,1,MPI_INT,0,MPI_COMM_WORLD); checkr(r,"gather"); //Faster to communicate, or to calculate?
    if (myrank==0)
    {
        kGLO=(double*)realloc(kGLO,n*n*sizeof(double));
        displ[0]=0;
        for (int i = 1; i < P; i++)
        {
            #ifdef DEBUG
                printf("displ[i]=%d displ[i-1]=%d tasksGLO[i-1]=%d\n",displ[i],displ[i-1],tasksGLO[i-1]);
            #endif
            displ[i]=displ[i-1]+tasksGLO[i-1];
            
        }
    }

    r=MPI_Gatherv(kLOC,ntasks,MPI_DOUBLE,kGLO,tasksGLO,displ,MPI_DOUBLE,0,MPI_COMM_WORLD); checkr(r,"gather");
    endtime   = MPI_Wtime();
    printf("Time for rank:%d ----> %f seconds\n",myrank,endtime-starttime);
    
    if (myrank==0)
    {   
        #ifdef DEBUG
            for (int i = 0; i < P; i++)
            {
                printf("recvcounts=%d  displs=%d\n",tasksGLO[i],displ[i]);
            }
            for (int px = 0; px < n; px++)
            {
                for (int py = 0; py < n; py++)
                {
                    printf("%f ", ac_mat(kGLO,px,py,n));
                }
                printf("\n");
            }
        #endif
        writeImage("test",kGLO,n,n,ch);
    }
    
    free(mystart);
	free(myend);
    free(kLOC);
    free(kGLO);
    free(displ);
    free(tasksGLO);
    MPI_Finalize();
}

void worksplit(int *mystart, int *myend, int proc, int nproc, int start, int end){
	// Comprovar que totes les variables tinguin valor permesos
	int nTasks = end - start + 1;
	int q = nTasks/nproc;
	int r = nTasks - nproc*q;
	//printf("nTasks = %d\n", nTasks);
	//printf("q = %d\n", q);
	//printf("r = %d\n", r);

	if(proc <= r-1){
		*mystart = start + proc*(q + 1);
		*myend = *mystart + q;
	} else {
		int S = start + (r - 1)*(q + 1) + q;
		*mystart = S + (proc - r)*q + 1;
		*myend = *mystart + q - 1;
	}
}

void checkr(int r,char *txt) {
if (r!=MPI_SUCCESS) {
fprintf(stderr,"Error: %s\n",txt);
exit(-1);
}
}

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