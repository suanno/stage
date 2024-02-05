#include <string.h>
#include <stdio.h>
#include<math.h>

#define pi  4*atan(1.0)

int main(){

int N=10000;
double L=100; // L=100,200,300,400
double dx=0.1; // N/L
double Hm2=1;

/*The scalar field is Phi(x,t).
- x dependance is in the array index
- for time dependance we save current time (t) value "u" and precedent time (t-1) value "up"
u = Phi(t)
up = Phi(t-1)
udt = 
*/
double u[10000];       
double udt[10000];
double up[10002];

double lapu[10000];

double dF[10000];
double uc[10000];
double Lu[10000];

double z;
double decax;
double decau;
double decaoutC=0;
int i;

double d2coef=1.0/pow(dx,2);

int loop;
double dt=0.001; // note: stable upto dt=0.001
double ttime;
double tmin=0;
double tmax=0.01;
double area;
int nloop=(tmax-tmin)/dt;

double C[10000];
for (loop=0;loop<nloop;loop++){
ttime = tmin + (loop+1)*dt;
/*Creating an array of equal width alternated segments of +1 and -1. The width is the half-period of sin(2pi*x)*/
if (sin(2*pi*ttime/1)>=0) C[loop]=1;
else C[loop]=-1;
}

FILE *filerandominit;
FILE *filetdglinit;
FILE *fileCout;

filerandominit = fopen("fileL1000init.dat", "r");
for (i=0; i<N; i++){
fscanf(filerandominit, "%lf\n", &z);
u[i]=z;
}

fclose(filerandominit);

//------------------------------------------

for(loop=0;loop<nloop;loop++) {

ttime = (loop+1)*dt+tmin;

for(i=0; i<N; i++) {
up[i+1] = u[i];
dF[i]=u[i]*u[i]*u[i]-C[loop]*Hm2*u[i];
}
up[0] = up[N];
up[N+1] = up[1];

/*Compute second order second derivative*/
for(i=1; i<N+1; i++) {
lapu[i-1]=d2coef*(up[i+1]+up[i-1]-2*up[i]);
}

for(i=0;i<N;i++) {
Lu[i] = lapu[i] - dF[i];
udt[i] = u[i] + dt*Lu[i];
}

for(i=0;i<N;i++) {
u[i] = udt[i];
}
}

fileCout = fopen("fileCout.dat", "w");
for (loop=0; loop<1000; loop++){
ttime = tmin + (loop+1)*dt;
decaoutC = C[loop];
fprintf(fileCout, "%.5f %.20f\n", ttime, decaoutC);
}
fclose(fileCout);

filetdglinit = fopen("tdgl_init.dat", "w");
for(i=0;i<N;i++) {
decax=i*dx;
decau=udt[i];
fprintf(filetdglinit, "%.2f %.20f\n", decax, decau);
}
fclose(filetdglinit);

return 0;
}

