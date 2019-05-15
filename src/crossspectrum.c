#include "crossspectrum.h"

#include <stdlib.h>

#include <fftw3.h>

void cross_spectrum(const double * const chanel_a[], const double * const chanel_b[], double cross_spectrum[], const size_t nb_elements)
{
	unsigned int y = 0;
	for (y = 0; y < nb_elements; ++y)
	{
		cross_spectrum[y] = ((chanel_a[0][y] * chanel_b[0][y]) + (chanel_a[1][y] * chanel_b[1][y]));
	}
}
