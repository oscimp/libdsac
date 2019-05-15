#include "tf.h"
//double *logspace(int d1, int d2, int n)
void logspace(int d1, int d2, int n, double *vec_log)
{

	//double *vec_log;
	//vec_log = (double *)malloc(sizeof(double) * (n * (d2 - d1) + 1));
	//memset(vec_log, 0, sizeof(double) * (n * (d2 - d1) + 1));

	int i;
	double k = 0.0;

	//for (i = 1; i < n * (d2 - d1); i++, k++) {
	for (i = 0; i < n * (d2 - d1); i++, k++) {

		vec_log[i] = pow(10, d1 + k / (n - 1));
//printf("%d\t%lf\n",i,vec_log[i]);
	}

	return;
	//return (vec_log);

}
