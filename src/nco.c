#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <cordic.h>

short *createNCOLutShort(long nb, int DATA_SIZE)
{
	double step = 2*M_PI/nb;
	double i;
	int ii;
	double cos_v; 
	double gain = pow(2, DATA_SIZE-1)-1;

	short *tab = (short *)malloc(sizeof(short) * nb);
	if (!tab) {
		printf("%s: erreur de malloc\n", __func__);
		return NULL;
	}

	for (ii=0, i=0.0f; ii < nb; ii++,i+=step) {
		cos_v = cos(i); 
		tab[ii] = (short)(cos_v * gain); 
	}
	return tab;
}

uint32_t NCOgenFreq(long freq_ref, long freq_out,
		int phase_accum_size, uint32_t start_value,
		short *nco_tab, int addr_size,
		short *out_cos, short *out_sin, int nb_output_elem)
{
	unsigned int decal = phase_accum_size- addr_size;
	uint32_t mask = (int)(pow(2, decal)-1.0f);
	long max_size = (long) pow(2, phase_accum_size)-1;
	double step = (double)freq_ref/pow(2,phase_accum_size);
	double incrf = ((double)freq_out)/step;
	unsigned long incr = (unsigned long) incrf;
	int i;
	uint32_t cos_addr, sin_addr;
	uint32_t phase = start_value;

	uint32_t sin_offset = (uint32_t)pow(2,(phase_accum_size-2));

	for (i=0; i < nb_output_elem; i++) {
		cos_addr = (~mask & phase) >> decal;
		sin_addr = (~mask & (phase-sin_offset)) >> decal;
		out_cos[i] = nco_tab[cos_addr];
		out_sin[i] = nco_tab[sin_addr];
		phase = (phase + incr) & max_size;
	}
	return phase;
}

