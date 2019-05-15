#ifndef H_REG_LIN_H
#define H_REG_LIN_H

#ifdef __cplusplus /* if C++, specify external linkage to C functions */
extern "C" {
#endif

#include <stddef.h>

void compute_xBarre_xMxB_xSquare(const double x[], size_t data_size,
                 double *xBarre, double xMxB[], double *xSquare);

void reg_lin(double *a, double *b, const double y[], const double x_moins_x_barre[], size_t tab_size,
         double x_barre, size_t data_div, double sum_div);

#ifdef __cplusplus
}
#endif

#endif // H_REG_LIN_H
