#include "tf.h"
// linspacef (for freq, with N/2 points),  linspacet (for time, i.e. with N points) and linspaced
// for integer spaced vector;
// last one shoul be renamed linspacei
//
void linspacef(int N, double fs, double * vec_lin)
{
/*ici il faut que la taille de vec_lin soit de N/2 -> error*/
	int i;
	for (i=0 ; i< N/2  ; i++) {
		vec_lin[i]=(fs*i)/N; //should be used

	}
	return;
}


//void linspacet(unsigned long N, double ts, double * vec_lin)
void linspacet(int N, double ts, double * vec_lin)
{
	int i;
	for (i=0 ; i< N  ; i++) {
		vec_lin[i]=ts*i; //should be used

	}
	return;
}


void linspaced(int N, int step, double * vec_lin)
{
	int i;
	for (i=0 ; i< N  ; i++) {
		vec_lin[i]=step*i; //should be used

	}
	return;
}


/*matlab/octave clone (without isinf,isequal verif)*/
void mylinspace(double base, double limit, unsigned int N, double *vec)
{
        double delta = (limit - base) / ((double)N -1);
        unsigned int k = 0;
        for (k = 0; k < N; k++) {
                vec[k] = base + delta * k;
        }

}
