/*!
 * \file logalizit.h
 * \brief Description courte...

 * Description longue...
 */


#ifndef DSAC_LOGALIZIT_H
#define DSAC_LOGALIZIT_H

#ifdef __cplusplus /* if C++, specify external linkage to C functions */
extern "C" {
#endif

void logalizit(int N, int n, double fs, int d2, double *f_lin, double *f_log,
  double *p_lin, double *p_log);

#ifdef __cplusplus
}
#endif

#endif // DSAC_LOGALIZIT_H
