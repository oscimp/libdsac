#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>		//pour optarg
#include <string.h>
#include <math.h>
#include <errno.h>

double * allocate_double(int size);
void dprint(double *vec, int size, char *buf);
void read_vector(double *vec, char *filename);
void read_2_vectors(double * vec1, double *vec2, char * filename, int size);
void newline(void);
void adc_raw_2_volt_stereo(double *x, double *y, int n, double vmax, double vmin, int M,
                           char * filename);
void show(double *x, int n, char *buf);
//void showc(double *x, int n, char *buf);

void gplot(double *, char *, double *, char *, int, char *, char *);
