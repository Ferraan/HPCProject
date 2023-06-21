#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "mpi.h"


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
    int r,l; 
    r=MPI_Init(&argc,&argv); checkr(r,"init");
    r=MPI_Comm_rank(MPI_COMM_WORLD,&myrank); checkr(r,"rank");
    r=MPI_Comm_size(MPI_COMM_WORLD,&P); checkr(r,"size");
    r=MPI_Get_processor_name(name,&l); checkr(r,"name");
    
    //Arguments
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

    
    
    const int nx=3,ny=2;
    if ((nx*ny)!=P)
    {
        printf("Different number of divisions than ranks\n");
        exit(1);
    }
    
    
    const double deltax=(WindowX[1]-WindowX[0])/(double)(n-1);
    const double deltay=(WindowY[1]-WindowY[0])/(double)(n-1);
    
    double starttime, endtime;
    int mystartx, myendx,mystarty, myendy;
    double* coordX=(double*)calloc(n,sizeof(double)); 
    double* coordY=(double*)calloc(n,sizeof(double));
        
    int count=0;
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {               
            worksplit(&mystartx,&myendx,i,nx,0,n-1);
            worksplit(&mystarty,&myendy,j,ny,0,n-1);
            if (count==myrank)
            {
               goto endLoops;
            }
            count++;
        }
    }
    endLoops:;

    int ntasksx=myendx-mystartx+1;
    int ntasksy=myendy-mystarty+1;
    int ntasks=0;
    ntasks=ntasksx*ntasksy;   

    //printf("Myrank:%d Tasksx:%d Tasksy:%d\n",myrank,ntasksx,ntasksy);
    //printf("Myrank:%d Mystartx:%d Myendx:%d Mystarty:%d Myendy:%d Ntasks:%d\n",myrank,mystartx,myendx,mystarty,myendy,ntasks);
    
    int* kLOC=(int*)calloc(ntasks,sizeof(int)); 
    starttime = MPI_Wtime();
    for(int i = mystartx, posk = 0; i <= myendx; i++){
		for (int j = mystarty; j <= myendy; j++,posk++)
        {
            //printf("Myrank:%d i:%d j:%d \n",myrank,i,j);
            double x0=WindowX[0]+j*deltax;
            double y0=WindowY[0]+i*deltay;        
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
            kLOC[posk]= iter-iteration;
            
            //printf("Myrank:%d X:%f Y:%f i:%d j:%d value:%d \n",myrank,x0,y0,i,j,kLOC[pos]);
            
        }
	}
    checkr(MPI_Barrier(MPI_COMM_WORLD),"Barrier");
    endtime   = MPI_Wtime();
    printf("Time for rank:%d before gather ----> %f seconds\n",myrank,endtime-starttime,kLOC[0]); fflush(stdout);
    //printf("myrank:%d tasks:%d tasksx=%d tasksy=%d \n",myrank,ntasks,ntasksx,ntasksy);
    #ifdef WRITE 
        //Print matrix

        int rows=mystarty;
        int columns=mystartx;
        int a=-1;
        for (int i = 0; i < ntasks; i++)
        {
            if (i!=0){
                
                if(i%ntasksy==0)
                {
                    rows=mystarty;
                }
                
            }
            
            
            printf("N(%d,%d)=%d div=%d\n",rows+1,columns+1,kLOC[i],(i+1)%(ntasksy));
            //printf("%d",ntasksy%i);
            if (i!=0){
                if((i+1)%(ntasksy)==0){
                
                    columns++;
                }
            }
            rows++;
            
        } 
    #endif

    free(kLOC);   
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