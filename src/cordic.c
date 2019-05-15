#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <cordic.h>
//#define PI 3.1415926535897932384626433832795
#define PI  3.14159265358979323846

/*
0 	45.00 	0.7854
1 	26.57 	0.4636
2 	14.04 	0.2450
3 	7.13 	0.1244
4 	3.58 	0.0624
5 	1.79 	0.0312
6 	0.90 	0.0160
7 	0.45 	0.0080
8 	0.22 	0.0040
9 	0.11 	0.0020
*/

void computeAtanLutS(int size, short *atanLut, int scaleFactor)
{
	int i;
	for (i = 0; i < size; i++)
		atanLut[i] = (short)(atan(pow(2, -i)) * (float)scaleFactor);
}

void computeAtanLut(int size, int *atanLut, int scaleFactor)
{
	int i;
	for (i = 0; i < size; i++)
		atanLut[i] = (int)(atan(pow(2, -i)) * (float)scaleFactor);
}

void computeAtanLutL(int size, long *atanLut, long scaleFactor)
{
	int i;
	for (i = 0; i < size; i++)
		atanLut[i] = (long)(atan(pow(2, -i)) * (double)scaleFactor);
}

void computeAtanLutLL(int size, long long *atanLut, long long scaleFactor)
{
	int i;
	for (i = 0; i < size; i++)
		atanLut[i] =
		    (long long)(atan(pow(2, -i)) * (double)scaleFactor);
}

void computeAtanLut64(int size, int64_t * atanLut, int64_t scaleFactor)
{
	int i;
	for (i = 0; i < size; i++)
		atanLut[i] = (int64_t) (atan(pow(2, -i)) * (double)scaleFactor);
}

void computeAtanLutf(int size, float *atanLut)
{
	int i;
	for (i = 0; i < size; i++)
		atanLut[i] = atan(pow(2, -i));
}

long computeCordic(long x, long y, int iter, int *atanLut)
{
	int i;
	int sign = -1;
	long ydiv, xdiv;
	long ytmp, xtmp;
	long theta, thetatmp;

	theta = 0;

	for (i = 0; i < iter; i++) {
		if (x < 0 && y < 0)
			sign = 1;
		else if (x < 0 || y < 0)
			sign = -1;
		else
			sign = 1;
		sign = -sign;
		ydiv = y >> i;
		xdiv = x >> i;
		//printf("%d %ld %ld ", i, y, ydiv);

		xtmp = x - (sign * ydiv);
		ytmp = y + (sign * xdiv);
		thetatmp = theta - (sign * atanLut[i]);
		theta = thetatmp;
		x = xtmp;
		y = ytmp;
		//printf("%ld %ld %ld\n", x, y, theta);

	}
	return theta;
}

long long computeCordicll(long long x, long long y, int iter,
			  long long *atanLut, int debug)
{
	int i;
	long long sign = -1, ydiv, xdiv, theta = 0;

	for (i = 0; i < iter; i++) {
		/*if (y<0)
		   sign = -1;
		   else
		   sign = 1; */
		if (x < 0 && y < 0)
			sign = 1;
		else if (x < 0 || y < 0)
			sign = -1;
		else
			sign = 1;
		sign = -sign;

		ydiv = y >> i;
		xdiv = x >> i;
		if (debug)
			printf("%d %Ld(%Ld) %Ld(%Ld) ", i, x, xdiv, y, ydiv);

		x -= sign * ydiv;
		y += sign * xdiv;
		/*if (debug) {
			if (x_orig != 0 && y_orig != 0) {
				if (x == 0 && y == 0)
					printf("%Ld %Ld %Ld %Ld\n", x, y, x_orig, y_orig);
				if (xdiv == 0 && ydiv == 0)
					printf("%Ld %Ld %Ld %Ld\n", x, y, x_orig, y_orig);
			}
		}*/
		theta -= sign * atanLut[i];
		if (debug)
			printf("%Ld %Ld %Ld %Ld\n", x, y, theta, atanLut[i]);

	}
	return theta;
}

float computeCordicf(float x, float y, int iter, float *atanLut)
{
	int i;
	float sign = -1;
	float ydiv, xdiv;
	float ytmp, xtmp;
	float theta, thetatmp;

	theta = 0;

	for (i = 0; i < iter; i++) {
		if (x < 0 && y < 0)
			sign = 1;
		else if (x < 0 || y < 0)
			sign = -1;
		else
			sign = 1;
		sign = -sign;
		ydiv = y * pow(2, -(i));
		xdiv = x * pow(2, -(i));
		//printf("%f %f %f ", y, ydiv, pow(2,-(i)));

		xtmp = x - (sign * ydiv);
		ytmp = y + (sign * xdiv);
		thetatmp = theta - (sign * atanLut[i]);
		theta = thetatmp;
		x = xtmp;
		y = ytmp;
		//printf("%f %f %f\n", x, y, theta);

	}
	return theta;
}

double correctDouble(double value, double x, double y)
{
	double pi_mult = ((double)PI);
	/* special cases */
	if (x == 0 && y == 0)
		return 0;
	if (x == 0) {
		if (y > 0)
			return pi_mult / 2;
		if (y < 0)
			return -pi_mult / 2;
	}
	if (y == 0) {
		if (x > 0)
			return 0;
		if (x < 0)
			return pi_mult;
	}

	/* classic cases */
	if (x < 0 && y < 0)
		value -= PI;
	else if (x < 0 && y >= 0)
		value += PI;
	return value;
}

float correctF(float value, float x, float y)
{
	float pi_mult = ((float)PI);
	/* special cases */
	if (x == 0 && y == 0)
		return 0;
	if (x == 0) {
		if (y > 0)
			return pi_mult / 2;
		if (y < 0)
			return -pi_mult / 2;
	}
	if (y == 0) {
		if (x > 0)
			return 0;
		if (x < 0)
			return pi_mult;
	}

	/* classic cases */
	if (x < 0 && y < 0)
		value -= PI;
	else if (x < 0 && y >= 0)
		value += PI;
	return value;
}

long correctL(long value, float x, float y, int scaleFactor)
{
	long long pi_mult = (long)((double)PI * (double)scaleFactor);
	/* special cases */
	if (x == 0 && y == 0)
		return 0;
	if (x == 0) {
		if (y > 0)
			return pi_mult / 2;
		if (y < 0)
			return -pi_mult / 2;
	}
	if (y == 0) {
		if (x > 0)
			return 0;
		if (x < 0)
			return pi_mult;
	}

	/* classic cases */
	if (x < 0 && y < 0)
		value -= ((double)scaleFactor) * PI;
	else if (x < 0 && y > 0)
		value += ((double)scaleFactor) * PI;
	return value;
}

long long correctLL(long long value, long long x, long long y,
		    long long scaleFactor)
{
	long long ret = value;
	long long pi_mult =
	    (long long)((long double)PI * (long double)scaleFactor);
	/* special cases */
	if (x == 0 && y == 0)
		return 0;
	if (x == 0) {
		if (y > 0)
			return pi_mult / 2;
		if (y < 0)
			return -pi_mult / 2;
	}
	if (y == 0) {
		if (x > 0)
			return 0;
		if (x < 0)
			return pi_mult;
	}

	/* classic cases */
	if (x < 0 && y < 0) {
		ret -= pi_mult;
	} else if (x < 0 && y > 0)
		ret += pi_mult;
	return ret;
}

void saveAtanLutLL(char *filename, int iter, long long *atanLutLL)
{
	int i;
	FILE *fd = fopen(filename, "w+");
	for (i = 0; i < iter; i++)
		fprintf(fd, "%Ld\n", atanLutLL[i]);
	fclose(fd);
}
