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

int N=1000;
double deca;
double hmoy=0.0;
int seed;

FILE *fileinit;

seed = time(NULL);
srand(seed);
fileinit = fopen("fileL10000init.dat", "w");
fprintf(fileinit, "%d\n", seed);
#pragma omp parallel for
for (i=0; i<N; i++){
deca=randU(-0.1, 0.1)+hmoy;
fprintf(fileinit, "%.20f\n", deca);
}

fclose(fileinit);

/*Recreate fileCout of values of C(t) [Progressive
executions of the dynamics will APPEND info]*/
FILE *fileCout;
fileCout = fopen("fileCout.dat", "w");
//fprintf(fileCout, "");
fclose(fileCout);


return 0;

}
