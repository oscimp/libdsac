#include "fir.h"

#ifdef DSA_OPENMP
#include <omp.h>
#endif // DSA_OPENMP

void fir_complex_32bits(const int32_t vali[], const int32_t valq[], const size_t data_size,
  const int32_t coeff[], const size_t coeff_size,
  const unsigned int decim, int64_t accumi[], int64_t accumq[])
{
  size_t i, ii, offset;
  int64_t tmp_vali, tmp_valq;

  /* convolution loop */
  for (i=0, offset=0; i< data_size-coeff_size; i+=decim, offset++) {
    /* for each conv */
    accumi[offset] = accumq[offset] = 0;
    for (ii=0; ii<coeff_size; ii++) {
      tmp_vali = (int64_t)coeff[ii] * (int64_t)vali[i+ii];
      accumi[offset] += (int64_t)tmp_vali;
      tmp_valq = (int64_t)coeff[ii] * (int64_t)valq[i+ii];
      accumq[offset] += (int64_t)tmp_valq;
      /*if (i==0 && ii < 100){
        printf("vali[%d]=%d\tcoeff[%d]=%d\ttmp_vali=%d\taccumi[0]=%Ld\n",
                ii,vali[ii],ii,coeff[ii],tmp_vali,accumi[i]);
      }*/
    }
    //if (i == 0) {
    //  printf("%Ld %Ld\n", accumi[offset], accumq[offset]);
    //}
  }
}

void fir_complex_double(const double vali[], const double valq[], const size_t data_size,
  const double coeff[], const size_t coeff_size,
  unsigned long decim,
  double accumi[], double accumq[]) {
  size_t i, ii, offset;
  double tmp_vali, tmp_valq;
  double c = 0.0;
  double y = 0.0;
  double t = 0.0;

  // Improve summation thanks Kahan algorithm
  // https://en.wikipedia.org/wiki/Kahan_summation_algorithm
  /* convolution loop */
  for (i=0, offset=0; i< data_size-coeff_size; i+=decim, offset++) {
    /* for each conv */
    accumi[offset] = 0.0;
    accumq[offset] = 0.0;

    for (ii=0; ii<coeff_size; ii++) {
      // Accum I value
      tmp_vali = coeff[ii] * vali[i + ii];
      y = tmp_vali - c;
      t = accumi[offset] + y;
      c = (t - accumi[offset]) - y;
      accumi[offset] = t;

      // Accum Q value
      tmp_valq = coeff[ii] * valq[i + ii];
      y = tmp_valq - c;
      t = accumq[offset] + y;
      c = (t - accumq[offset]) - y;
      accumq[offset] = t;
    }
  }
}

/*a changer decim en  uint*/
void fir_double(const double x[], const size_t data_size, const double coeff[], const size_t coeff_size,
    const unsigned long decim, double accumx[])
{
#ifdef DSA_OPENMP
  #pragma omp parallel
  {
    size_t i, ii, offset;
    double tmp_valx;
    double c = 0.0;

    // Split the task
    const size_t limit = data_size - coeff_size;
    const size_t nthreads = omp_get_num_threads();
    const size_t ithread = omp_get_thread_num();
    const size_t start_i = ithread * limit / nthreads;
    const size_t finish_i = (ithread + 1) * limit /nthreads;
    const size_t start_offset = ithread * limit / decim / nthreads;

    /* convolution loop */
    for (i=start_i, offset=start_offset; i < finish_i; i+=decim, ++offset) {
      /* for each conv */
      accumx[offset] = 0.0;
      for (ii=0; ii<coeff_size; ii++) {
        tmp_valx = coeff[ii] * x[i+ii];
        double y = tmp_valx - c;
        double t = accumx[offset] + y;
        c = (t - accumx[offset]) - y;
        accumx[offset] = t;
      }
    }
  }
#else // !DSA_OPENMP
  size_t i, ii, offset;
  double tmp_valx;
  double c = 0.0;

  // Improve summation thanks Kahan algorithm
  // https://en.wikipedia.org/wiki/Kahan_summation_algorithm
  /* convolution loop */
  for (i=0, offset=0; i< data_size-coeff_size; i+=decim, offset++) {
    /* for each conv */
    accumx[offset] = 0.0;
    for (ii=0; ii<coeff_size; ii++) {
      tmp_valx = coeff[ii] * x[i+ii];
      double y = tmp_valx - c;
      double t = accumx[offset] + y;
      c = (t - accumx[offset]) - y;
      accumx[offset] = t;
    }
  }
#endif // DSA_OPENMP
}

void fir_float(const float x[], const size_t data_size, const float coeff[], const size_t coeff_size,
    const unsigned long decim, float accumx[])
{
#ifdef DSA_OPENMP
  #pragma omp parallel
  {
    size_t i, ii, offset;
    float tmp_valx;
    float c = 0.0;

    // Split the task
    const size_t limit = data_size - coeff_size;
    const size_t nthreads = omp_get_num_threads();
    const size_t ithread = omp_get_thread_num();
    const size_t start_i = ithread * limit / nthreads;
    const size_t finish_i = (ithread + 1) * limit /nthreads;
    const size_t start_offset = ithread * limit / decim / nthreads;

    /* convolution loop */
    for (i=start_i, offset=start_offset; i < finish_i; i+=decim, ++offset) {
      /* for each conv */
      accumx[offset] = 0.0;
      for (ii=0; ii<coeff_size; ii++) {
        tmp_valx = coeff[ii] * x[i+ii];
        float y = tmp_valx - c;
        float t = accumx[offset] + y;
        c = (t - accumx[offset]) - y;
        accumx[offset] = t;
      }
    }
  }
#else // !DSA_OPENMP
  size_t i, ii, offset;
  float tmp_valx;
  float c = 0.0;

  // Improve summation thanks Kahan algorithm
  // https://en.wikipedia.org/wiki/Kahan_summation_algorithm
  /* convolution loop */
  for (i=0, offset=0; i< data_size-coeff_size; i+=decim, offset++) {
    /* for each conv */
    accumx[offset] = 0.0;
    for (ii=0; ii<coeff_size; ii++) {
      tmp_valx = coeff[ii] * x[i+ii];
      float y = tmp_valx - c;
      float t = accumx[offset] + y;
      c = (t - accumx[offset]) - y;
      accumx[offset] = t;
    }
  }
#endif // DSA_OPENMP
}

void fir_int64(const int64_t x[], const size_t data_size, const int64_t coeff[], const size_t coeff_size,
    const unsigned long decim, int64_t accumx[])
{
  size_t i, ii, offset;
  int64_t tmp_valx;

  /* convolution loop */
  for (i=0, offset=0; i< data_size-coeff_size; i+=decim, offset++) {
    /* for each conv */
    accumx[offset] = 0;
    for (ii=0; ii<coeff_size; ii++) {
      tmp_valx = coeff[ii] * x[i+ii];
      accumx[offset] += tmp_valx;
    }
  }
}

void fird(double *coeff, double *x,  int coeff_size, long data_size,
                int decim, double *accumx)
{
  int i, ii, offset;
  double tmp_valx;

  /* convolution loop */
  for (i=0, offset=0; i< data_size-coeff_size; i+=decim, offset++) {
    /* for each conv */
    accumx[offset] = 0;
    for (ii=0; ii<coeff_size; ii++) {
      tmp_valx = coeff[ii] * x[i+ii];
      accumx[offset] += tmp_valx;
      /*if (i==0 && ii < 100){
        printf("vali[%d]=%d\tcoeff[%d]=%d\ttmp_vali=%d\taccumi[0]=%Ld\n",
                            ii,vali[ii],ii,coeff[ii],tmp_vali,accumi[i]);
      }*/
    }
    //if (i == 0) {
    //      printf("%Ld %Ld\n", accumi[offset], accumq[offset]);
    //}
  }
}
