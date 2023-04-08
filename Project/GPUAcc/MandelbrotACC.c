#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <openacc.h>
#include "Image.h"

#define ac_mat(X,i,j,n) X[(i)*(n)+(j)]
//Posar timing del calcul
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
    double coordX[n]; //Can be arrays, not changing potser no, el heap? es petit
    double coordY[n];
    double* k=(double*)calloc(n*n,sizeof(double)); //Pot ser un static array

    coordX[0]=WindowX[0];
    coordY[0]=WindowY[0];
    //Allocate matrix coordinates: 
    #pragma acc parallel loop
    for (int i = 1; i < n; i++)
    {
        coordX[i]=coordX[0]+i*deltax;
        coordY[i]=coordY[0]+i*deltay;
    }
    
    //Optimized code
    time1=omp_get_wtime();
    //Compute the mandelbrot set
    #pragma acc parallel loop collapse(2)
    for (int px = 0; px < n; ++px)
    {   
        for (int py = 0; py < n; ++py)
        {
            double x0=coordX[px]; //Posar dintre altre for per OPEN*
            double y0=coordY[py];
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
            
            ac_mat(k,py,px,n) = iter-iteration;
            
        }
    }
    time2=omp_get_wtime();
    dub_time = time2 - time1 ;
    printf("Time for n=%d and iter=%d gpu ACC -----> %lf \n", n,iter,dub_time);

    

    #ifdef DEBUG
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
        writeImage("test",k,n,n,ch);
    #endif
    free(k);

}