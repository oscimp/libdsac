#include "logalizit.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>

void logalizit(int N, int n, double fs, int d2, double *f_lin, double *f_log,
	      double *p_lin, double *p_log)
{
	int i_inf = (int)(N * pow(10, d2 - 1) / fs);
	int i_sup = (int)(N * pow(10, d2) / fs);
	int avg = 0;
	int index = 0;
	double power=0.0;
	double f_lin_i=0.0;
	int i;

	if (n > i_sup - i_inf) {
		fprintf(stderr,
			"operation not permitted (logalize function)\n");
		exit(EXIT_FAILURE);
	}

	for (i = i_inf; i < i_sup; i++) {
		f_lin_i = f_lin[i];
		if (fabs(f_lin_i - f_log[index]) <
		    fabs(f_lin_i - f_log[index + 1])) {
			avg++;
			power += p_lin[i];
		} else {
			p_log[index] = power / avg;
			avg = 1;
			power = p_lin[i];
			index++;
		}
	}
}
