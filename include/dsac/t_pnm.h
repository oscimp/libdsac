#ifndef _T_PNM_H
#define _T_PNM_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "gnuplot_i.h"

#ifdef __cplusplus /* if C++, specify external linkage to C functions */
extern "C" {
#endif

#define PI 3.1415926535897932384626433832795
#define M_2PI 2.0f*M_PI


double constrainAngle(double x);
double angleConv(double angle);
double angleDiff(double a, double b);
double unwrap(double previousAngle, double newAngle);

void unWrappInt64(int64_t *phi1, int64_t *phi_out, int nb_elem);


void save_data(double *raw1, double *raw2, int N, char f_name[100]);
void save_data_plot(double *raw1, double *raw2, int N, char f_name[100]);
void save_1_vector_data(double *raw, int N, char f_name[100]);
void save_2_vector_data_plot(double *raw1, double *raw2, int N, char f_name[100]);
void hanning_window(double *td_in, double *td_out, int N);
void crop_fft(double *in, double *out, double *fft_avg_result, int N, int crop, double overlp, int repeat);
void fft_avg(double *td_in, double *fd_out, int num, int avg);
void convert_to_Lf(double *freq_domain_liner, double *freq_domain_db_out, int num, double lna_gain, double k_phi, double freq_bin, double window_coef);
void decimate_1_vector(double *in_data, double *out_data, int N, int fir_tap, int decimation_ratio);
void decimate_2_vectors(double *in_time, double *in_data, double *out_time, double *out_data, int N, int fir_tap, int decimation_ratio);
void fft_save_logplot(double *in, int num, double freq_bin, char f_name[100]);
void crop_window_fft_average(double *in, double *Lf_out, double *freq_in, int N, int crop_point, int repeat, double overlap, double lna_gain, double k_phi, double fs, char file_name[100]);


void init_nco_calc_IQ(double *raw, double f_nco, double fs, int n, double *rawI, double *rawQ, double phase_shift);
void init_nco_calc_IQ_debug(double *raw, double f_nco, double fs, int n, double *rawI, double *rawQ, double phase_shift);
void time_gen(double *t, int n, double fs);
void phase_calc(double *I, double *Q, double *phase, int n);
void phase_calc2(const double *I, const double *Q, double *phase,const int n);


/*!
 * \brief Calcul la phase d'un signal démodulé
 * Attends les I et les Q d'un signal démodulé afin d'y appliquer un atan2. De plus
 * la fonction gère la dérive de la phase.
 *
 * \param I Tableau contenant les I
 * \param Q Tableau contenant les Q
 * \param phase Tableau contenant le résultat
 * \param data_length Taille du tableau de données
 */
void phase_compute(const double I[], const double Q[], double phase[], const unsigned long data_length);


void subtract(double *a, double *b, double *out, int n);
void add(double *a, double *b, double *out, int n);
void amplify(double *in, double *out, double gain, int n);
//void phase_calc_no_expand(double *I, double *Q, double *phase, int n);
//void subtract_compensate(double *a, double *b, double *out, int n);




void amplitude_calc(double *I, double *Q, int nob, double vfsr, double *amplitude, int n);
void am_noise_conv(double *am_noise, int n, double carrier_power);
void judgement(int n, int crop, int repeat, int decim, double ovlp, int tap, double fs);
void powerspectrum(double *rawa, double *rawb, double *Lf_out, double *freq_out, int n, int crop, int repeat, double ovlp, double fs, char file_name[100]);
double powerspectrum2(double *rawa, double *rawb, double *Lf_out, double *freq_out, int n, int crop, int repeat, double ovlp, double fs, char file_name[100]);
void powerspectrum3(double * data_arm_A, double *data_arm_B, double *Cospctl, int win_size);
void normalize_psd_dBc(  double * Cospctl , const int win_size , const double fs , double * psd);
void powerSpectralDensity0 (  double * armx , double * army , const int win_size , const double fs , double * psd);
void cwfx(double *rawa, double *rawb, double *Lf_out, double *freq_out, int n, int crop, int repeat, double ovlp, double fs, char file_name[100]);
void cwfx_infinite(double *rawa, double *rawb, double *Lf_out, double *freq_out, int n, int crop, int repeat, double ovlp, double fs, char file_name[100]);
void convert_to_lf(double *Lf_lin, double *Lf_log, double *freq, int n, int corr_factor, char file_name[100]);
void convert_to_am(double *Lf_lin, double *Lf_log, double *freq, double carrier_power, int n, int corr_factor, char file_name[100]);
void convert_to_lf_dBrad2(double *Lf_lin, double *Lf_log, double *freq, int n, int corr_factor, char file_name[100]);
void convert_to_am_dBrad2(double *Lf_lin, double *Lf_log, double *freq, double carrier_power, int n, int corr_factor, char file_name[100]);
void convert_to_lf_dBV2(double *Lf_lin, double *Lf_log, double *freq, int n, int corr_factor, double fsr, int nob, double gain, char file_name[100]);
void hogehoge(double *out, double *in, int n, double average_factor);
double get_nominal_frequency_from_phase(double *phase, int n, double dt, double f_nco);
double mean(double * vec , unsigned long size);
void remove_mean(double * vec , unsigned long size);
double vector_max(double * v , const unsigned long size);

#ifdef __cplusplus
}
#endif

#endif // _T_PNM_H
