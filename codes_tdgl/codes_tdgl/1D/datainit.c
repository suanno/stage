//#include <iostream.h>
//#include <fstream.h>
#include <string.h>
#include <stdio.h>
#include<math.h>
#include<stdlib.h>
#include<omp.h>
#include <time.h>

#define pi  4*atan(1.0)

//pour générer des nombres aléatoires
double randU(double randmin, double randmax)
{
double randU1=0.;
randU1 = randmin*(1-rand()/(double)RAND_MAX)+randmax*rand()/(double)RAND_MAX;
return randU1;
}

int main(){

int i;

int N=10000;
double deca;
double hmoy=0.0;

FILE *fileinit;

srand(time(NULL));
fileinit = fopen("fileL1000init.dat", "w");
#pragma omp parallel for
for (i=0; i<N; i++){
deca=randU(-0.1, 0.1)+hmoy;
fprintf(fileinit, "%.20f\n", deca);
}

fclose(fileinit);

return 0;

}
