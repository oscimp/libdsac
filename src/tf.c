#include "tf.h"
void init_nco(double f_nco, double *t, double *cosine, double *sine, int taille)
{
        int i;
        double dpf=2*PI*f_nco;
        for (i = 0; i < taille; i++) {
                cosine[i] = cos(dpf * t[i] +1.507);
                sine[i] = sin(dpf* t[i]+1.507);
        }

}

void calc_IQ(double *x, double *cosine, double *sine, double *rawI,
             double *rawQ, int taille)
{
        int i;
        for (i = 0; i < taille; i++) {
                rawI[i] = x[i] * cosine[i];
                rawQ[i] = x[i] * sine[i];
        }
}

//////void fft(double *in, double *out, int N);
//void fft(double *in, fftw_complex *out, int N);
//void fft(double *in, fftw_complex *out, int N){
////void fft(double *in, double *out, int N);
////void fft(double *in, double *out, int N){
//fftw_plan plan_forward;
//
////plan_forward = fftw_plan_dft_r2c_1d (N, in, out, FFTW_ESTIMATE);
//
////fftw_execute (plan_forward);
//}

