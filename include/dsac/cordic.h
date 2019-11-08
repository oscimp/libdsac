/*!
 * \file cordic.h
 * \brief Description courte...

 * Description longue...
 */

#ifndef DSAC_CORDIC_H
#define DSAC_CORDIC_H

#ifdef __cplusplus /* if C++, specify external linkage to C functions */
extern "C" {
#endif

long computeCordic(long x, long y, int iter, int *atanLut);
void computeAtanLut(int size, int *atanLut, int scaleFactor);
void computeAtanLutS(int size, short *atanLut, int scaleFactor);
long long computeCordicll(long long x, long long y, int iter, long long *atanLut, int debug);
void computeAtanLutL(int size, long *atanLut, long scaleFactor);
void computeAtanLutLL(int size, long long *atanLut, long long scaleFactor);
void computeAtanLut64(int size, int64_t *atanLut, int64_t scaleFactor);
void saveAtanLutLL(char *filename, int iter, long long *atanLutLL);
float computeCordicf(float x, float y, int iter, float *atanLut);
void computeAtanLutf(int size, float *atanLut);
double correctDouble(double value, double x, double y);
float correctF(float value, float x, float y);
long correctL(long value, float x, float y, int scaleFactor);
long long correctLL(long long value, long long x, long long y, long long scaleFactor);

#ifdef __cplusplus
}
#endif

#endif // DSAC_CORDIC_H
