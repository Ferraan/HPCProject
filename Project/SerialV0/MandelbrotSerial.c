#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "Image.h"
#define ac_mat(X,i,j,n) X[(i)*(n)+(j)]
//Posar timing del calcul
int main(int argc, char *argv[]){
    clock_t time1, time2;
    double dub_time;
    const int n=atoi(argv[1]);  //Up to n=40000 using 16Gb RAM and a little bit of swap
    const int iter=atoi(argv[2]);
    const int ch=1;
    const double WindowX[2]={-2,0.47};
    const double WindowY[2]={-1.12, 1.12};
    const double deltax=(WindowX[1]-WindowX[0])/(double)(n-1);
    const double deltay=(WindowY[1]-WindowY[0])/(double)(n-1);
    double* coordX=(double*)calloc(n,sizeof(double)); //Can be arrays, not changing potser no, el heap? es petit
    double* coordY=(double*)calloc(n,sizeof(double));
    double* k=(double*)calloc(n*n,sizeof(double)); //Pot ser un static array
    double* k2=(double*)calloc(n*n,sizeof(double)); 

    coordX[0]=WindowX[0];
    coordY[0]=WindowY[0];
    //Allocate matrix coordinates: 
    for (int i = 1; i < n; i++)
    {
        coordX[i]=coordX[0]+i*deltax;
        coordY[i]=coordY[0]+i*deltay;
    }
    //Naive code
    time1=clock();
    //Compute the mandelbrot set
    for (int px = 0; px < n; ++px)
    {   
        double cx=coordX[px]; //Posar dintre altre for per OPEN*
        for (int py = 0; py < n; ++py)
        {
            
            double cy=coordY[py];
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
            ac_mat(k,px,py,n) = zz;
            
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
        double x0=coordX[px]; //Posar dintre altre for per OPEN*
        for (int py = 0; py < n; ++py)
        {
            
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
            
            ac_mat(k2,px,py,n) = iter-iteration;
            
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
    writeImage("test",k,n,n,ch);
    writeImage("test2",k2,n,n,ch);

    free(k);
    free(k2);
}