#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>		//pour optarg
#include <string.h>
#include <math.h>
//#include <fftw.h>
#define PI 3.1415926535897932384626433832795
//tf_logspace.c
//double *logspace(int, int, int);
//logspace(d1,d2,n,vec_log);
void logspace(int, int, int, double *);

//linspace(N,fs,vec_lin);
void linspacef(int N, double fs, double * vec_lin);
//void linspacet(unsigned long,double,double *);
void linspacet(int N, double ts, double * vec_lin);
void linspaced(int N, int step, double * vec_lin);
void mylinspace(double base, double limit, unsigned int N, double *vec);


//tf.c
void init_nco(double f_nco, double *t, double *cosine, double *sine, int length);

void calc_IQ(double *x, double *cosine, double *sine, double *rawI,
             double *rawQ, int taille);

//void fft(double *in, double *out, int N);
