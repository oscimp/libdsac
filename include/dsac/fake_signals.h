#ifndef _FAKE_SIGNALS_H
#define _FAKE_SIGNALS_H

#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "conversions.h"

#ifdef __cplusplus /* if C++, specify external linkage to C functions */
extern "C" {
#endif

/*!
* \brief Various fake signals creation
*/

void create_fake_adc_signal_with_noise(double * xtt , double phi_offset , double fs  , double fc , double power , uint16_t nob, char * int_type, double * fake_signal  , unsigned long size);

#ifdef __cplusplus
}
#endif

#endif // _FAKE_SIGNALS_H

