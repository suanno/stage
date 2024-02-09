// Solve the TDGL equation
// Forward integration using implicit Euler scheme
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

#define pi 4*atan(1.0)

int main(int argc, char  *argv [ ]){

FILE *fileinit;
int N=1000;
double dx=0.1;

double dt=0.1;
double tmin=1;	//It is the tmax of "tdglfd"
double Deltat = 10;
double tmax = tmin + Deltat;

/* C(t)=Ampl*sign(sin(2pi/Thalf*t) 
   so it switchd from +Ampl to -Ampl every Thalf time    
*/
double Ampl = 1;
double Thalf = tmax;  // with this choice C will be constantly +1

/* READ PARAMETERS */
//fileinit = fopen("tdgl_init.dat", "r");
fileinit = fopen("tdgl_result.dat", "r");
/*First line is for parameters and seed*/
int seed;
fscanf(fileinit, "%d %lf %lf %lf %d %lf %lf\n", &N, &tmin, &dx, &dt, &seed, &Ampl, &Thalf);
fclose(fileinit);

/* Get inputs from the terminal */
char *ptr;
//printf("argv1 = %lf", atof(argv[1]));
if (argc > 1)
	Deltat = strtod(argv[1], &ptr);
if (argc > 2){
  	Ampl = strtod(argv[2], &ptr);
}
if (argc > 3){
  	Thalf = strtod(argv[3], &ptr);
	//printf("\n%lf\n", Thalf);
}
tmax = tmin + Deltat;
/*
if (argc > 2){
  	Ampl = strtod(argv[2], &ptr);
}
if (argc > 3){
  	Thalf = strtod(argv[3], &ptr);
}
*/

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
double decaoutAve=0;
double decaoutq2mean=0;
double decainC=0;
double decatime=0;

double ffr[10000];
double qfr[10000];
double d2coef[10000];
double integ_coef[10000];
double q2meannum=0.0;
double q2meandenum=0.0;

FILE *stateeqn_result;
FILE *fileCout;				/*C(t) values*/
FILE *fileAveout;			/*Space average of u(t) values*/
//FILE *fileCin;
FILE *fileq2mean;

double ttime=0;
int nloop=(Deltat)/dt;
int loop;
double q2mean[nloop];
double C[nloop];
double uAve[nloop];

//fileCin = fopen("fileC.dat", "r");
//for (loop=0; loop<nloop; loop++){
//fscanf(fileCin, "%lf %lf \n", &decatime, &decainC);
//C[loop]=decainC;
//}
//fclose(fileCin);

/*Define value of C(t) in time.
	It switches from +A to -A each T/2 time */
for (loop=0;loop<nloop;loop++){
ttime = tmin + (loop+1)*dt;
if(Thalf > 0){
	C[loop]=Ampl*sin(pi*ttime/Thalf);
	//C(t_{k+1})
	/*
	if (sin(pi*ttime/Thalf)>=0) 
		C[loop]=Ampl;
	else 
		C[loop]=-Ampl;
	*/
	}
	else
		C[loop] = Ampl;
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

/*Define q-space lattice*/
for (i=1; i<(N/2); i++){
ffr[i]=i;
ffr[N-i]=-i;
}
ffr[0]=0;
ffr[N/2]=N/2;
for (i=0; i<N; i++){
qfr[i]=ffr[i]*2*pi/N;
}

/*LOAD INITIAL STATE*/
//fileinit = fopen("tdgl_init.dat", "r");
fileinit = fopen("tdgl_result.dat", "r");
/*First line is for parameters and seed*/
fscanf(fileinit, "%d %lf %lf %lf %d %lf %lf\n", &N, &tmin, &dx, &dt, &seed, &Ampl, &Thalf);
/*Now read the initial (smooth) state*/
for (i=0; i<N; i++){
fscanf(fileinit, "%lf %lf \n", &decainx, &decainu);
if(i==0)
	printf("%lf\n", decainx);
x[i]=decainx;
u[i]=decainu;
}
fclose(fileinit);

//------------------------------------

// coefficients for spectral derivatives (-q^2) array
for (i=0; i<N; i++){
d2coef[i]=-(qfr[i]/dx)*(qfr[i]/dx);
}

for (loop=0; loop < nloop; loop++){

ttime = tmin + (loop+1)*dt;

/*Compute space average*/
uAve[loop] = 0;
for (i = 0; i<N; i++){
	uAve[loop] = uAve[loop] + u[i];
}
uAve[loop] = uAve[loop]/N;
//printf("Ave(%lf) = %lf", ttime, uAve[loop]);
/*Denominator of implicit euler [As equation is non-linear, the numerator will be unusual, but the denominator no]*/
for (i=0; i<N; i++){
integ_coef[i]=1-dt*C[loop]-dt*d2coef[i]; //Note that d2coeff = -q^2 (already contains a minus sign)
}

/*Compute FFT of u(x)->u(q)*/
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

/*???*/
q2meannum=0.0;
q2meandenum=0.0;
for(i=0; i<N; i++) {
q2meannum = q2meannum + (qfr[i]/dx)*(qfr[i]/dx)*(ufr[i]*ufr[i]+ufi[i]*ufi[i]);
q2meandenum = q2meandenum + (ufr[i]*ufr[i]+ufi[i]*ufi[i]);
}
q2mean[loop] = q2meannum/q2meandenum;

// ********** begin of main algorithm **********

/*Compute FFT of u(x)^3 */
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
/*IMPLICIT EULER*/
//implicit Euler scheme
udtfr[i]=(ufr[i]-dt*NLfr[i])/integ_coef[i];	//The numerator has unusual shape because eq. in non-linear
udtfi[i]=(ufi[i]-dt*NLfi[i])/integ_coef[i];
}

// ********** end of main algorithm **********

/*INVERSE FFT*/
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

/*Save the final state*/
stateeqn_result = fopen("tdgl_result.dat", "w");
/*Save parameters N, tmax, dx, dt, seed, Ampl, Thalf*/
fprintf(stateeqn_result, "%d %.2lf %.10lf %.10lf %d %lf %lf\n", N, tmax, dx, dt, seed, Ampl, Thalf);
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

fileCout = fopen("fileCout.dat", "a");
for (loop=0; loop<nloop; loop++){
ttime = tmin + (loop+1)*dt;
decaoutC = C[loop];
fprintf(fileCout, "%.5f %.20f\n", ttime, decaoutC);
}
fclose(fileCout);

fileAveout = fopen("fileAveout.dat", "a");
for (loop=0; loop<nloop; loop++){
ttime = tmin + (loop+1)*dt;
decaoutAve = uAve[loop];
fprintf(fileAveout, "%.5f %.20f\n", ttime, decaoutAve);
}
fclose(fileAveout);

fftw_destroy_plan(pf);
fftw_destroy_plan(pb);
fftw_free(in);
fftw_free(out);
	return 0;
}

