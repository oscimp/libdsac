/*!
 * \file nco.h
 * \brief Description courte...

 * Description longue...
 */


#ifndef DSAC_NCO_H
#define DSAC_NCO_H

#ifdef __cplusplus /* if C++, specify external linkage to C functions */
extern "C" {
#endif

/*!
* \brief Produces an cos array. Warning array must be freed manually
* \param nb: array and period length
* \param DATA_SIZE: scale factor for each value => 2^(DATA_SIZE-1)-1
*/
short *createNCOLutShort(long nb, int DATA_SIZE);

/*!
* \brief Produces an sin and cos array. Simulates NCO behavior
* \param freq_ref: clock frequency used as reference
* \param freq_out: wished output frequency
* \param phase_accum_size: number of bits used for phase accumulator
* \param start_value: start value for phase accumulator. Used for multiple
* 	function calls
* \param nco_tab: array containing set of value for cos waveform
* \param addr_size: number of bits used for LUT address
* \param out_cos: array containing cos out values
* \param out_sin: array containing sin out values
* \param nb_output_elem: arrays length
*/
uint32_t NCOgenFreq(long freq_ref, long freq_out,
        int phase_accum_size, uint32_t start_value,
        short *nco_tab, int addr_size,
        short *out_cos, short *out_sin, int nb_output_elem);

#ifdef __cplusplus
}
#endif

#endif // DSAC_NCO_H
