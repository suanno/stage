//#include <iostream.h>
//#include <fstream.h>
#include <string.h>
#include <stdio.h>
#include<math.h>
#include<stdlib.h>
#include<omp.h>

#define pi  4*atan(1.0)

//pour générer des nombres aléatoires
float randU(float randmin, float randmax)
{
float randU1=0.;
randU1 = randmin*(1-rand()/(float)RAND_MAX)+randmax*rand()/(float)RAND_MAX;
return randU1;
}

int main(){

int i;
int j;

int N=2500;
double deca;
double hmoy=0.0;

FILE *fileinit;

fileinit = fopen("file0init.dat", "w");
#pragma omp parallel for
for (i=0; i<N; i++){
for (j=0; j<N; j++){
deca=randU(-0.1, 0.1)+hmoy;
fprintf(fileinit, "%.20f\n", deca);
}
}

fclose(fileinit);

}
