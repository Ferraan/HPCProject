#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "Image.h"
#define ac_mat(X,i,j,n) X[(i)*(n)+(j)]
//Posar timing del calcul
int main(){
    const int n=5000;
    const int iter=1000;
    const int ch=1;
    const double WindowX[2]={-2,1};
    const double WindowY[2]={-1, 1};
    const double deltax=(WindowX[1]-WindowX[0])/(double)(n-1);
    const double deltay=(WindowY[1]-WindowY[0])/(double)(n-1);
    double coordX[n]; //Can be arrays, not changing potser no, el heap? es petit
    double coordY[n];
    double* k=(double*)calloc(n*n,sizeof(double)); //Pot ser un static array


    coordX[0]=WindowX[0];
    coordY[0]=WindowY[0];
    //Allocate matrix coordinates: 
    for (int i = 1; i < n; i++)
    {
        coordX[i]=coordX[0]+i*deltax;
        coordY[i]=coordY[0]+i*deltay;
    }
    

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
            for (int i = 1; i < iter+1; ++i)
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
    writeImage("test",k,n,n,ch);
}