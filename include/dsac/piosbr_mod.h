/*!
 * \file piosbr_mod.h
 * \brief Description courte...

 * Description longue...
 */


#ifndef DSAC_PIOSBR_MOD_H
#define DSAC_PIOSBR_MOD_H

#ifdef __cplusplus /* if C++, specify external linkage to C functions */
extern "C" {
#endif

double * allocate_double(int size);
void dprint(double *vec, int size, char *buf);
void read_vector(double *vec, char *filename);
void read_1_vector(double *vec, char *filename, int size);	//timaike mod. Juillet2014
void read_2_vectors(double * vec1, double *vec2, char * filename, int size);
void newline(void);
void adc_raw_2_volt_stereo(double *x, double *y, int n, double vmax, double vmin, int M,
                           char * filename);
void show(double *x, int n, char *buf);
//void showc(double *x, int n, char *buf);

void gplot(double *f_log, char *format_f, double *p_log, char *format_p, int n,
	   char *logmode, char *title);

#ifdef __cplusplus
}
#endif

#endif // DSAC_PIOSBR_MOD_H
