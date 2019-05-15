#include "fake_signals.h"

void create_fake_adc_signal_with_noise(double * xtt , double phi_offset , double fs  , double fc , double power , uint16_t nob, char * int_type, double * fake_signal  , unsigned long size)
{

                /*signal 0dBm */
                        //s_fake = 0.31623 * cos ( 2 * pi * fc * tt + phi_offset + 2 * pi * fc * xtt');
                        //s_fake_adc = int16(s_fake*2^15);
        unsigned long k;
        double temp;//=0.0;
	double temp0; //=0.0;
        double dpn = 2.0 * M_PI  * fc;
	double vp = dbm2vp(power);
	double pow_nob = pow(2,nob-1);


        for (k=0 ; k<size ; k++){

                temp0 = dpn * (k+1) / fs + phi_offset + dpn * xtt[k];
                temp = vp * cos ( temp0 );
                /*ok for fpga, should be a_16 = round(a_u * 2^15) cf oppenheim 4th p449*/

                fake_signal[k]  = (int16_t)(round(temp  * pow_nob));

        }

}
