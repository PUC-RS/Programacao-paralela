#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

# define MAXITER 2000

struct complex{
  double real;
  double img;
};

void verifyPointInSet(struct complex complexStruct);

int numoutside = 0;

int main(){
  int i, j;
  double area, error;
  double start, finish;
  double eps = 1.0e-7;
  struct complex z, c;

  start = omp_get_wtime();
#pragma omp parallel for firstprivate(eps) private(c,j)
  for (i=0; i<NPOINTS; i++) {
    for (j = 0; j<NPOINTS; j++) {
      c.real = -2.0+2.5*(double)(i)/(double)(NPOINTS)+eps;
      c.img = 1.125*(double)(j)/(double)(NPOINTS)+eps;
      verifyPointInSet(c);
    }
  }

  finish = omp_get_wtime();

  area=2.0*2.5*1.125*(double)(NPOINTS*NPOINTS-numoutside)/(double)(NPOINTS*NPOINTS);
  error=area/(double)NPOINTS;

  printf("Area of Mandlebrot set = %12.8f +/- %12.8f\n",area,error);
  printf("Time = %12.8f seconds\n",finish-start);
}

void verifyPointInSet(struct complex c){
  double ztemp;
  struct complex z;

  z=c;
      for (int iter=0; iter<MAXITER; iter++){
        ztemp=(z.real*z.real)-(z.img*z.img)+c.real;
        z.img=z.real*z.img*2+c.img;
        z.real=ztemp;
        if ((z.real*z.real+z.img*z.img)>4.0e0) {
        #pragma omp reduction(+:numoutside)
          numoutside++;
          break;
        }
      }
}