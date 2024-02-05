// Solve the TDGL equation
// Forward integration using implicit Euler scheme
#include <stdio.h>
#include <math.h>
#include <fftw3.h>

#define pi 4*atan(1.0)

int main(){

const int N=10000;
const double dt=0.1;
const double dx=0.1;

double x[10000];
double u[10000];
double ufr[10000];
double ufi[10000];
double udt[10000];
double udtfr[10000];
double udtfi[10000];

double NL[10000];
double NLfr[10000];
double NLfi[10000];

int i;
double decainx=0;
double decainu=0;
double decaoutx=0;
double decaoutu=0;
double decaoutC=0;
double decaoutq2mean=0;
double decainC=0;
double decatime=0;

double ffr[10000];
double qfr[10000];
double d2coef[10000];
double integ_coef[10000];
double q2meannum=0.0;
double q2meandenum=0.0;

FILE *fileinit;
FILE *stateeqn_result;
FILE *fileCout;
//FILE *fileCin;
FILE *fileq2mean;

double ttime=0;
double tmin=1;
double tmax=10;
int nloop=(tmax-tmin)/dt;
int loop;
double q2mean[nloop];
double C[nloop];

//fileCin = fopen("fileC.dat", "r");
//for (loop=0; loop<nloop; loop++){
//fscanf(fileCin, "%lf %lf \n", &decatime, &decainC);
//C[loop]=decainC;
//}
//fclose(fileCin);

for (loop=0;loop<nloop;loop++){
ttime = tmin + (loop+1)*dt;
//C(t_{k+1})
//C[loop]=sin(2*pi*ttime/1);
if (sin(2*pi*ttime/1)>=0) C[loop]=1;
else C[loop]=-1;
}

ttime=0;
loop=0;

char filename[16];

/* TUNG: create plan for FFT*/
fftw_complex *in, *out;
fftw_plan pf,pb;
   
in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
pf = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
pb = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

for (i=1; i<(N/2); i++){
ffr[i]=i;
ffr[N-i]=-i;
}
ffr[0]=0;
ffr[N/2]=N/2;
for (i=0; i<N; i++){
qfr[i]=ffr[i]*2*pi/N;
}

fileinit = fopen("tdgl_init.dat", "r");
for (i=0; i<N; i++){
fscanf(fileinit, "%lf %lf \n", &decainx, &decainu);
x[i]=decainx;
u[i]=decainu;
}
fclose(fileinit);

//------------------------------------

// coefficients for spectral derivatives
for (i=0; i<N; i++){
d2coef[i]=-(qfr[i]/dx)*(qfr[i]/dx);
}

for (loop=0; loop < nloop; loop++){

ttime = tmin + (loop+1)*dt;

for (i=0; i<N; i++){
integ_coef[i]=1-dt*C[loop]-dt*d2coef[i];
}

// transform h
for(i=0; i<N; i++) {
in[i][0]=u[i];
in[i][1]=0.0;
}
fftw_execute(pf); // repeat as needed
for(i=0; i<N; i++) {
ufr[i]=out[i][0];
ufi[i]=out[i][1];
}

q2meannum=0.0;
q2meandenum=0.0;
for(i=0; i<N; i++) {
q2meannum = q2meannum + (qfr[i]/dx)*(qfr[i]/dx)*(ufr[i]*ufr[i]+ufi[i]*ufi[i]);
q2meandenum = q2meandenum + (ufr[i]*ufr[i]+ufi[i]*ufi[i]);
}
q2mean[loop] = q2meannum/q2meandenum;

// ********** begin of main algorithm **********

for (i=0; i<N; i++){
NL[i]=u[i]*u[i]*u[i];
}

for(i=0; i<N; i++) {
in[i][0]=NL[i];
in[i][1]=0.0;
}
fftw_execute(pf); // repeat as needed
for(i=0; i<N; i++) {
NLfr[i]=out[i][0];
NLfi[i]=out[i][1];
}

for (i=0; i<N; i++){
// implicit Euler scheme
udtfr[i]=(ufr[i]-dt*NLfr[i])/integ_coef[i];
udtfi[i]=(ufi[i]-dt*NLfi[i])/integ_coef[i];
}

// ********** end of main algorithm **********

for(i=0; i<N; i++) {
in[i][0]=udtfr[i];
in[i][1]=udtfi[i];
}
fftw_execute(pb); // repeat as needed
for(i=0; i<N; i++) {
udt[i]=out[i][0]/N;
}

for(i=0; i<N; i++) {
u[i]=udt[i];
}

}

stateeqn_result = fopen("tdgl_result.dat", "w");
for (i=0; i<N; i++){
decaoutx=x[i];
decaoutu=u[i];
fprintf(stateeqn_result, "%.1f %.20f\n", decaoutx, decaoutu);
}
fclose(stateeqn_result);

fileq2mean = fopen("q2mean.dat", "w");
for (loop=0; loop<nloop; loop++){
ttime = tmin + (loop+1)*dt;
decaoutq2mean=q2mean[loop];
fprintf(fileq2mean, "%.2f %.20f\n", ttime, decaoutq2mean);
}
fclose(fileq2mean);

fileCout = fopen("fileCout.dat", "w");
for (loop=0; loop<nloop; loop++){
ttime = tmin + (loop+1)*dt;
decaoutC = C[loop];
fprintf(fileCout, "%.5f %.20f\n", ttime, decaoutC);
}
fclose(fileCout);

fftw_destroy_plan(pf);
fftw_destroy_plan(pb);
fftw_free(in);
fftw_free(out);
	return 0;
}

