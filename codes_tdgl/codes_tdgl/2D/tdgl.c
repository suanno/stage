#include <string.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
//#include </opt/intel/composer_xe_2013_sp1.3.174/mkl/include/fftw/fftw3.h>
#include <fftw3.h>

#define pi  4*atan(1.0)
int main(){

int N=2500;
static double h[2500][2500];
static double hfr[2500][2500];
static double hfi[2500][2500];
static double hdt[2500][2500];
static double hdtfr[2500][2500];
static double hdtfi[2500][2500];
static double coeffr[2500][2500];
static double hc[2500][2500]; // h^3
static double hcfr[2500][2500];
static double hcfi[2500][2500];
static double q2[2500][2500]; // qx^2 + qy^2

static double lapfrcoef[2500][2500];
static double integ_coef[2500][2500];

//static double areaar[2500][2500];
//static double Aarr[2500][2500];
//static double Barr[2500][2500];
//static double Carr[2500][2500];

//double A;
//double B;
//double C;

double ffr[N];
double qfr[N];
double L=1000; // data for x is 0 to L
double dx=L/N;
double z;
double Hm2=1;
int i;
int j;
int k;

double deca;
double decaoutq2mean=0.0;

double q2meannum=0.0;
double q2meandenom=0.0;

int loop;
double dt=0.1; // note: stable upto dt=0.6
int Nst=8;
double savetimes[8]={10,100,1000,2000,4000,6000,8000,10000};
double time;
double tmin=1;
double tmax=10001;
int loops=(tmax-tmin)/dt;
double q2mean[loops];
static double q2meannumar[2500][2500];
static double q2meandenomar[2500][2500];

int q2loop=0;
double qtime=0.0;

double C[loops];
for (loop=0;loop<loops;loop++){
time = tmin + (loop+1)*dt;
//C(t_{k+1})
//C[loop]=sin(2*pi*time/1);
if (sin(2*pi*time/1) >= 0.0) C[loop]=1;
else C[loop]=-1;
}

FILE *fileinit;
FILE *filenonconserve2dh;
FILE *fileq2mean;

char filename[16];

// initialize threads for fft
fftw_init_threads();

fftw_complex *in, *out;
fftw_plan pf, pb;
in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N);
out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N);
// create plans with threads
fftw_plan_with_nthreads(omp_get_max_threads());
// plan for forward transform
pf = fftw_plan_dft_2d(N,N,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
// pf = fftw_plan_dft_r2c_2d(N,N,in,out,FFTW_ESTIMATE);
// plan for backward transform
pb = fftw_plan_dft_2d(N,N,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
// pb = fftw_plan_dft_c2r_2d(N,N,in,out,FFTW_ESTIMATE);
//TUNG: end of adding variables for fft

for (i=1; i<(N/2); i++){
ffr[i]=i;
ffr[N-i]=-i;
}
ffr[0]=0;
ffr[N/2]=N/2;
for (i=0; i<N; i++){
qfr[i]=ffr[i]*2*pi/N;
}


fileinit = fopen("file1.dat", "r");
for (i=0; i<N; i++){
for (j=0; j<N; j++){
fscanf(fileinit, "%lf\n", &z);
h[i][j]=z;
}
}


fclose(fileinit);

// coefficients for implicit scheme
for (i=0; i<N; i++){
for (j=0; j<N; j++){
q2[i][j]=(qfr[i]*qfr[i]+qfr[j]*qfr[j]);
lapfrcoef[i][j]=-(qfr[i]*qfr[i]+qfr[j]*qfr[j])/(dx*dx);
}
}

//------------------------------------------

for(loop=0;loop<loops;loop++) {

time = (loop+1)*dt+tmin;

for (i=0; i<N; i++){
for (j=0; j<N; j++){
integ_coef[i][j]=1-dt*C[loop]-dt*lapfrcoef[i][j];
}
}

// inside big loops
#pragma omp parallel for //seulement pour les grands systèmes
for(i=0; i<N; i++) {
for(j=0; j<N; j++) {
hc[i][j] = h[i][j]*h[i][j]*h[i][j];
}
}

// h to hfr, hfi
#pragma omp parallel for //seulement pour les grands systèmes
for(i=0;i<N;i++) {
for(j=0;j<N;j++) {
in[i*N+j][0] = h[i][j];
in[i*N+j][1] = 0.0;
}
}
// backward transform
fftw_execute(pf);
// copy output vector
#pragma omp parallel for //seulement pour les grands systèmes
for(i=0;i<N;i++) {
for(j=0;j<N;j++) {
hfr[i][j] = out[i*N+j][0];
hfi[i][j] = out[i*N+j][1];
}
}

q2meannum=0.0;
q2meandenom=0.0;
#pragma omp parallel for //seulement pour les grands systèmes
for(i=0; i<N; i++) {
for(j=0; j<N; j++) {
q2meannumar[i][j] = 0.5*((qfr[i]+qfr[j])*(qfr[i]+qfr[j])/(dx*dx))*(hfr[i][j]*hfr[i][j]+hfi[i][j]*hfi[i][j]);
q2meandenomar[i][j] = (hfr[i][j]*hfr[i][j]+hfi[i][j]*hfi[i][j]);
}
}
for(i=0; i<N; i++) {
for(j=0; j<N; j++) {
q2meannum=q2meannum+q2meannumar[i][j];
q2meandenom=q2meandenom+q2meandenomar[i][j];
}
}
q2mean[loop] = q2meannum/q2meandenom;

// hc to hcfr, hcfi
#pragma omp parallel for //seulement pour les grands systèmes
for(i=0;i<N;i++) {
for(j=0;j<N;j++) {
in[i*N+j][0] = hc[i][j];
in[i*N+j][1] = 0.0;
}
}
// backward transform
fftw_execute(pf);
// copy output vector
#pragma omp parallel for //seulement pour les grands systèmes
for(i=0;i<N;i++) {
for(j=0;j<N;j++) {
hcfr[i][j] = out[i*N+j][0];
hcfi[i][j] = out[i*N+j][1];
}
}

// implicit scheme
#pragma omp parallel for //seulement pour les grands systèmes
for (i=0; i<N; i++){
for (j=0; j<N; j++){
if(lapfrcoef[i][j] != 0.0) {
hdtfr[i][j] = (hfr[i][j]-dt*hcfr[i][j])/integ_coef[i][j];
hdtfi[i][j] = (hfi[i][j]-dt*hcfi[i][j])/integ_coef[i][j];
}
}
}

// initialize input vector in for backward transform
#pragma omp parallel for //seulement pour les grands systèmes
for(i=0;i<N;i++) {
for(j=0;j<N;j++) {
in[i*N+j][0] = hdtfr[i][j];
in[i*N+j][1] = hdtfi[i][j];
}
}
// backward transform
fftw_execute(pb);
// copy output vector
#pragma omp parallel for //seulement pour les grands systèmes
for(i=0;i<N;i++) {
for(j=0;j<N;j++) {
hdt[i][j] = out[i*N+j][0]/(N*N);
}
}

#pragma omp parallel for //seulement pour les grands systèmes
for(i=0;i<N;i++) {
for(j=0;j<N;j++) {
h[i][j] = hdt[i][j];
}
}

for(k=0;k<Nst;k++) {
if(time==savetimes[k]) {
sprintf( filename, "file%.0f.dat", time );  
filenonconserve2dh = fopen(filename, "w");
//filenonconserve2dh = fopen("nonconserve2dh.dat", "w");
for(i=0;i<N;i++) {
for(j=0;j<N;j++) {
deca=hdt[i][j];
fprintf(filenonconserve2dh, "%.20f\n", deca);
}
}
fclose(filenonconserve2dh);

fileq2mean = fopen("q2mean.dat", "w");
for (q2loop=0; q2loop<loops; q2loop++){
qtime = tmin + (q2loop+1)*dt;
decaoutq2mean=q2mean[q2loop];
fprintf(fileq2mean, "%.2f %.20f\n", qtime, decaoutq2mean);
}
fclose(fileq2mean);
}
}

//if(fmod(time,15.) ==0) {
//printf("%.1f\n", time);
//}

//if(fmod(time,3000.) ==0) {
//sprintf( filename, "file%.0f.dat", time );  
//filenonconserve2dh = fopen(filename, "w");
//filenonconserve2dh = fopen("nonconserve2dh.dat", "w");
//for(i=0;i<N;i++) {
//for(j=0;j<N;j++) {
//deca=hdt[i][j];
//fprintf(filenonconserve2dh, "%.20f\n", deca);
//}
//}
//fclose(filenonconserve2dh);
//}

}

//------------------------------------------

fftw_destroy_plan(pf);
fftw_destroy_plan(pb);
fftw_free(in);
fftw_free(out);
// clean up threads
fftw_cleanup_threads();
return 0;
}

