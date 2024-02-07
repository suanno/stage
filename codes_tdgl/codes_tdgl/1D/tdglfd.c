#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#include <unistd.h>        // chdir
#include <sys/stat.h>      // mkdir


#define pi  4*atan(1.0)

int main(int argc, char* argv[]){

int N=1000;            /*Lattice sites*/
double dx = 0.1;        /*Lattice spacing*/
double L=N*dx;
double Hm2=1;           /* ? */

double dt=0.001; // note: stable upto dt=0.001
double tmin=0;
double tmax=1;  // Suggest "1" because in the fft algorithm tmin = 1 is hardcoded 

/* C(t)=Ampl*sign(sin(2pi/Thalf*t) 
   so it switchd from +Ampl to -Ampl every Thalf time    
*/
double Ampl = 1;
double Thalf = tmax;  // with this choice C will be constantly +1


/* Get inputs from the terminal */
char *ptr;
//printf("argv1 = %lf", atof(argv[1]));
if (argc > 1)
	tmax=strtod(argv[1], &ptr);
if (argc > 2){
  Ampl = strtod(argv[2], &ptr);
}
if (argc > 3){
  Thalf = strtod(argv[3], &ptr);
}
  


/*The scalar field is Phi(x,t) and we call it u(x,t)
- x dependance is in the array index
- for time dependance we save current time (t) value "u" and precedent time (t-1) value "up"
u = Phi(t)
up = Phi(t-1)
udt = dPhi/dt
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
double ttime;

double area;
int nloop=(tmax-tmin)/dt;

double C[nloop];

/*Define the values of C(t) at different times
  We set C(t) switching between +Ampl and -Ampl with the
  period of sin(2pi*t)
*/


for (loop=0;loop<nloop;loop++){
    ttime = tmin + (loop+1)*dt;
    if (sin(pi*ttime/Thalf)>=0) C[loop]=Ampl;
    else C[loop]=-Ampl;
}
/*
for (loop = 0; loop < nloop; loop++){
  C[loop] = 1;
}
*/

/*Read the initial state*/
FILE *filerandominit;
FILE *filetdglinit;
FILE *fileCout;
int seed;

filerandominit = fopen("fileL10000init.dat", "r");
fscanf(filerandominit, "%d\n", &seed);
for (i=0; i<N; i++){
fscanf(filerandominit, "%lf\n", &z);
u[i]=z;
}

fclose(filerandominit);


/*--- Numerical solving A.L. eq ---*/

/*Time loop*/

for(loop=0;loop<nloop;loop++) {

    ttime = (loop+1)*dt+tmin;

    for(i=0; i<N; i++) {
    up[i+1] = u[i];
    dF[i]=u[i]*u[i]*u[i]-C[loop]*Hm2*u[i];
    }

    /*PBC boundaries*/
    up[0] = up[N];
    up[N+1] = up[1];

    /*Compute second order second derivative*/
    for(i=1; i<N+1; i++) {
    lapu[i-1]=d2coef*(up[i+1]+up[i-1]-2*up[i]);
    }

    /*Compute EXPLICIT Euler du/dt and then du (and so u(t+1))*/
    for(i=0;i<N;i++) {
    Lu[i] = lapu[i] - dF[i];
    udt[i] = u[i] + dt*Lu[i];
    }
    for(i=0;i<N;i++) {
    u[i] = udt[i];
    }

    /*Save only state at specific times separated by dt_plot*/

}

/*Save values of C(t) in different times*/
fileCout = fopen("fileCout.dat", "w");
for (loop=0; loop<nloop; loop++){
ttime = tmin + (loop+1)*dt;
decaoutC = C[loop];
fprintf(fileCout, "%.5f %.20f\n", ttime, decaoutC);
}
fclose(fileCout);

/*Save the final state*/
//filetdglinit = fopen("tdgl_init.dat", "w");
filetdglinit = fopen("tdgl_result.dat", "w");
/*Save parameters N, tmax, dx, dt, seed*/
tmax = tmax + tmin; // So we save the TOTAL time of the dynamics
fprintf(filetdglinit, "%d %.2lf %.10lf %.10lf %d %lf %lf\n", N, tmax, dx, dt, seed, Ampl, Thalf);
for(i=0;i<N;i++) {
decax=i*dx;
decau=udt[i];
fprintf(filetdglinit, "%.2f %.20f\n", decax, decau);
}
fclose(filetdglinit);

return 0;
}

