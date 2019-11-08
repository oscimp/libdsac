#include "t_pnm.h"

#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>

#include "gnuplot_i.h"
#include "piosbr_mod.h"
#include "file_handler.h"

double constrainAngle(double x)
{
	x = fmod(x + M_PI, M_2PI);
	if (x < 0)
		x += M_2PI;
	return x - M_PI;
}

		    // convert to [-360,360]
double angleConv(double angle)
{
	return fmod(constrainAngle(angle), M_2PI);
}

double angleDiff(double a, double b)
{
	double dif = fmod(b - a + M_PI, M_2PI);
	if (dif < 0)
		dif += M_2PI;
	return dif - M_PI;
}

double unwrap(double previousAngle, double newAngle)
{
	return previousAngle - angleDiff(newAngle, angleConv(previousAngle));
}


void unWrappInt64(int64_t *phi1, int64_t *phi_out, int nb_elem)
{
	int64_t phase_diff;
	int64_t phase_accum=0; //ds fpga sur 2 samples, ici fct^t ~ ?
	int64_t prev1 = phi1[0];
	int64_t res;
	int64_t pi_gain = (int64_t)(2*M_PI*pow(2,24)); //ce gain c'est koi-> arg

	int i;

	for (i=1; i<nb_elem; i++) {
		phase_diff = phi1[i] - prev1;
		res = phi1[i] + phase_accum;
		if ((prev1 < 0) && (phi1[i] >= 0)) { /* neg to pos */
			if (phase_diff > pi_gain/2) {
				res = phi1[i] + phase_accum - (pi_gain);
				phase_accum = phase_accum - (pi_gain);
			}
		} else if ((prev1 >= 0) && (phi1[i] < 0)) { /* pos to neg */
			if (phase_diff < -pi_gain/2) {
				res = phi1[i] + phase_accum + (pi_gain);
				phase_accum = phase_accum + (pi_gain);
			}
		}
		prev1 = phi1[i];
		phi_out[i] = res;
	}
}







//save dat + plot
void save_2_vector_data_plot(double *raw1, double *raw2, int N, char f_name[100])
{
	FILE *fp_out;
	fp_out = fopen(f_name, "w");
	int j = N - 1;
	while(j >= 0){
		fprintf(fp_out, "%lf\t%lf\n", raw1[j], raw2[j]);
		j--;
	}
	fclose(fp_out);
	//gplot(raw1,"%e",raw2,"%e",N+1,"logx",f_name);
}

void save_2_vector_data_plot_e(double *raw1, double *raw2, int N, char f_name[100])
{
	FILE *fp_out;
	fp_out = fopen(f_name, "w");
	int j = N - 1;
	while(j >= 0){
		//fprintf(fp_out, "%.10e\t%e.10\n", raw1[j], raw2[j]/1e10/2/PI);
		fprintf(fp_out, "%e\t%e\n", raw1[j], raw2[j]);
		j--;
	}
	fclose(fp_out);
	//gplot(raw1,"%e",raw2,"%e",N,"nop",f_name);
}

void save_2_vector_data_plot_e_inv(double *raw1, double *raw2, int N, char f_name[100])
{
	FILE *fp_out;
	fp_out = fopen(f_name, "w");
	int j = 0;//N - 1;
	while(j <= N-1){
		fprintf(fp_out, "%.10e\t%.10e\n", raw1[j], raw2[j]/1e10/2/PI);
		//fprintf(fp_out, "%e\t%e\n", raw1[j], raw2[j]);
		j++;
	}
	fclose(fp_out);
	//gplot(raw1,"%e",raw2,"%e",N,"nop",f_name);
}

void save_1_vector_data(double *raw, int N, char f_name[100])
{
	FILE *fp_out;
	fp_out = fopen(f_name, "w");
	int j = N - 1;
	while(j >= 0){
		//fprintf(fp_out, "%d\t%e\n", j, raw[j]);
		fprintf(fp_out, "%.10e\n", raw[j]);
		j--;
	}
	fclose(fp_out);
	//gplot(raw1,"%lf",raw2,"%lf",N,"nop",f_name);
}

//window function
void hanning_window(double *window_in_data, double *window_out_data, int N)
{
	int i;
	double *w = (double *)malloc(sizeof(double) * N);
	double S1;
	double S2;

	for(i = 0; i < N; i++) w[i] = 0;

	//--create hanning window--
	for (i = 0; i < N ; i++){
	w[i]=0.5*(1.0-cos((2.0*(double)i*PI)/((double)N-1.0)));
	//printf("%e\n",w[i]);
	}

	//--compute normalization sums--
	for (i = 0; i < N ; i++){
	S1+=w[i];
	S2+=pow(w[i],2);
	}
	//printf("%e\t%e\n", S1/N, S2/N);

	//--window multiply--
	for (i = 0; i < N; i++){
	window_out_data[i] = w[i] * window_in_data[i];
	//window_out_data[i] = window_in_data[i];	//rectangular
	}
	free(w);
}


void hft95_window(double *window_in_data, double *window_out_data, int N)
{
	int i;
	double *w = (double *)malloc(sizeof(double) * N);
	double S1;
	double S2;

	for(i = 0; i < N; i++) w[i] = 0;

	//--create hanning window--
	for (i = 0; i < N ; i++){
		w[i] = 1 - 1.9383379*cos(((2.0*(double)i*PI) / (double)N)) + 1.3045202*cos(2*((2.0*(double)i*PI) / (double)N)) - 0.4028270*cos(3*((2.0*(double)i*PI) / (double)N)) + 0.0350665*cos(4*((2.0*(double)i*PI) / (double)N));
	//printf("%e\n",w[i]);
	}

	//--compute normalization sums--
	for (i = 0; i < N ; i++){
	S1+=w[i];
	S2+=pow(w[i],2);
	}
	//printf("%e\t%e\n", S1/N, S2/N);

	//--window multiply--
	for (i = 0; i < N; i++){
	window_out_data[i] = w[i] * window_in_data[i];
	//window_out_data[i] = window_in_data[i];	//rectangular
	}
	free(w);
}



//cropping & FFT
void crop_fft(double *in, double *out, double *fft_avg_result, int N, int crop, double overlp, int repeat)
{
	int crop_point_end = crop;
	int a = 0;
	int aa;
	int repeat_counter = 0;

	while(repeat_counter < repeat){
		aa = crop - 1 ;
//		printf("crop point start = %8d\tcrop point end = %8d\n", a, crop_point_end - 1);
		while(a < crop_point_end){
			out[aa] = in[N-1-a];
			a++;
			aa--;
		}

//		save_2_vector_data_plot(xaxis_debug, out, crop, "temp.cropped_temp.dat"); //before hann.win. plot

		hanning_window(out, out, crop);
//		save_2_vector_data_plot(xaxis_debug, out, crop, "temp.cropped_temp_temp.dat"); //after hann.win. plot
		fft_avg(out, fft_avg_result, crop, repeat);

		crop_point_end = crop * (1 - overlp) + a;

		a = crop_point_end - crop;

		repeat_counter++;
		//printf("crop point start = %8d\tcrop point end = %8d\n", a + 1, crop_point_end);
	}
}

//FFT module
void fft_avg(double *time_domain_in, double *freq_domain_out, int num, int avg)
{
	int nc =(num / 2) + 1;
	double *in;
	fftw_complex *out;
	fftw_plan plan_forward;

	in = fftw_malloc(sizeof(double) * num);
	out = fftw_malloc(sizeof(fftw_complex) *nc);

	plan_forward = fftw_plan_dft_r2c_1d(num, time_domain_in, out, FFTW_ESTIMATE);
	fftw_execute(plan_forward);

	int y;
	for (y = 0; y < nc; y++){
		freq_domain_out[y] += (2*2*(pow(out[y][0],2)+pow(out[y][1],2)) /pow(num, 2)) / avg;
	}

	fftw_destroy_plan(plan_forward);
	fftw_free(in);
	fftw_free(out);
}

void convert_to_Lf(double *freq_domain_liner, double *freq_domain_db_out, int num, double lna_gain, double k_phi, double freq_bin, double window_coef)
{
	int i;
	for(i = 0; i < num; i++){
		freq_domain_db_out[i] = -10*log10(2) + 10*log10(freq_domain_liner[i]) - lna_gain -20*log10(k_phi) - 10*log10(freq_bin) - 10*log10(window_coef);
	}
}


void decimate_1_vector(double *in_data, double *out_data, int N, int fir_tap, int decimation_ratio)
{
	int di = 0;
	int dii = 0;

	while(di < ((N-fir_tap) / decimation_ratio) + 1){
		out_data[di] = in_data[dii];
		di++;
		dii = dii +decimation_ratio;
	}
}



void decimate_2_vectors(double *in_time, double *in_data, double *out_time, double *out_data, int N, int fir_tap, int decimation_ratio)
{
	int di = 0;
	int dii = 0;

	while(di < ((N-fir_tap) / decimation_ratio) + 1){
		out_time[di] = in_time[dii];
		out_data[di] = in_data[dii];
		di++;
		dii = dii +decimation_ratio;
	}
}


/*
void fft_save_logplot(double *in, int num, double freq_bin, char f_name[100])
{
	int i;
	double *freq = (double *)malloc(sizeof(double) * num);
	FILE *fp_out;
	fp_out = fopen(f_name, "w");

	for(i = 0; i < (num / 2) + 1; i++){
		freq[i] = freq_bin * i;
		fprintf(fp_out, "%e\t%e\n", freq[i], in[i]);
	}
	fclose(fp_out);
	gplot(freq, "%lf", in, "%lf", num, "logx", f_name);
}
*/

void crop_window_fft_average(double *in, double *Lf_out, double *freq_out, int N, int crop_point, int repeat, double overlap, double lna_gain, double k_phi, double fs, char file_name[100])
{
	double freq_bin = 1 / (crop_point / fs);
	double *temp_cropped_data = (double *)malloc(sizeof(double) * crop_point);
	double *temp_fft_avg_result = (double *)malloc(sizeof(double) * crop_point);
	//double *temp_cropped_data = allocate_double(crop_point);
	//double *temp_fft_avg_result = allocate_double(crop_point);

	int i;
	for(i = 0; i < crop_point; i++){
		temp_cropped_data[i] = 0;
		temp_fft_avg_result[i] = 0;
	}

	crop_fft(in, temp_cropped_data, temp_fft_avg_result, N, crop_point, overlap, repeat);
	convert_to_Lf(temp_fft_avg_result, Lf_out, crop_point, lna_gain, k_phi, freq_bin, 1.5);

	//int i;
	FILE *fp_out;
	fp_out = fopen(file_name, "w");

	for(i = 0; i < (crop_point / 2) + 1; i++){
		freq_out[i] = freq_bin * i;
		fprintf(fp_out, "%e\t%e\n", freq_out[i], Lf_out[i]);
	}
	fclose(fp_out);

	free(temp_cropped_data);
	free(temp_fft_avg_result);

	//gplot(freq_out, "%lf", Lf_out, "%lf", crop_point / 2 + 1, "logx", file_name);

}





/*
void init_nco_calc_IQ(double *raw, double f_nco, double fs, int n, double *rawI, double *rawQ, double phase_shift)
{
	int i;
	double dt = 1/fs;
	double t = 0;
	double cosine = 0;
	double sine = 0;
	double *tutu=allocate_double(n);
	for(i = 0; i < n; i++){
		t = dt * i;
		cosine = cos(2*PI*f_nco * t + phase_shift);
		sine = sin(2*PI*f_nco * t + phase_shift);

		rawI[i] = raw[i] * cosine;
		rawQ[i] = raw[i] * sine;
		tutu[i] = cosine;
	}
	double *toto=allocate_double(n);time_gen(toto,n,1);
	gplot(toto,"%f",tutu,"%f",128,"nop","cosine");
	//gplot(axis, "%e", rawa, "%e", N/10, "nop", "rawa");
}
*/

/*n should be long*/
void init_nco_calc_IQ(double *raw, double f_nco, double fs, int n, double *rawI, double *rawQ, double phase_shift)
{
	int i;
	double dpnt = 2. * PI * f_nco / fs;

	for(i = 0; i < n; i++){
		rawI[i] = raw[i] * cos ( dpnt * i + phase_shift );
		rawQ[i] = raw[i] * sin ( dpnt * i + phase_shift );
	}

        /*debug : si on veut plotter ou sortir les vecteurs nco*/
	////	double dpnt = 2. * PI * f_nco / fs;
	////	double * cosine = allocate_double(n);
	////	double * sine = allocate_double(n);
	////	double * xaxis  = allocate_double(n);
	////
	////	for(i = 0; i < n; i++){
	////		cosine[i] = cos ( dpnt * i + phase_shift );
	////		sine[i] = sin ( dpnt * i + phase_shift );
	////		xaxis[i]= i;
	////	}

	////	/*dangerous if successive calls !*/
	////	/*use init_nco_calc_IQ_debug for this*/
	////	//choose your crop here
	////	int crop = 128;
	////	gplot(xaxis,"%f",cosine,"%f",crop,"nop","cosine");
	////	gplot(xaxis,"%f",sine,"%f",crop,"nop","cosine");
	////	//storeDouble("cosine",crop,cosine);
	////	//storeDouble("sine",crop,sine);
	////	//storeDouble("xaxis",crop,xaxis);
}

void init_nco_calc_IQ_debug(double *raw, double f_nco, double fs, int n, double *rawI, double *rawQ, double phase_shift)
{
	/* process data as usual */
	init_nco_calc_IQ(raw,f_nco, fs, n, rawI, rawQ, phase_shift);

	int i;

        /*debug : si on veut plotter ou sortir les vecteurs nco*/
	double dpnt = 2. * PI * f_nco / fs;
	double * cosine = allocate_double(n);
	double * sine = allocate_double(n);
	double * xaxis  = allocate_double(n);

	for(i = 0; i < n; i++){
		cosine[i] = cos ( dpnt * (i+1) + phase_shift );
		sine[i] = sin ( dpnt * (i+1) + phase_shift );
		xaxis[i]= i;
	}

	/*dangerous if successive calls !*/
	//choose your crop here
	int crop = n;
	//gplot(xaxis,"%f",cosine,"%f",crop,"nop","cosine");
	//gplot(xaxis,"%f",sine,"%f",crop,"nop","cosine");
	storeDouble("fakecosine",crop,cosine);
	storeDouble("fakesine",crop,sine);
	//storeDouble("xaxis",crop,xaxis);

}

void time_gen(double *t, int n, double fs)
{
        int i;
        //double dt = 1/fs;
        for(i = 0; i < n; i++)
                //t[i] = dt * i;
		t[i] = (double)i / fs;
}

void phase_calc(double *I, double *Q, double *phase, int n)
{
        int i;
        int j = 0;

        for(i = 0; i < n; i++){
                //phase[i] = atan(I[i] / Q[i]);
		//phase[i] = atan2(I[i], Q[i]);
		phase[i] = atan2(Q[i], I[i]);
	}

	// prosess of phase unwrapping (this algorithm is very slow)
        for(i = 1; i < n; i++){
                if(phase[i] - phase[i-1] > PI/1){// if you use atan(): PI/2
                        j = i;
                        while(j < n){
                                phase[j] = phase[j] - 2*PI;// if you use atan(): -1*PI
                                j++;
                        }
                }
                if(phase[i] - phase[i-1] < -PI/1){// if you use atan(): PI/2
                        j = i;
                        while(j < n){
                                phase[j] = phase[j] + 2*PI;// if you use atan: 1*PI
                                j++;
                        }
                }
        }
}


void phase_calc2(const double *I, const double *Q, double *phase,const int n)
{
	int i;
	for(i = 0; i < n; i++){
		phase[i] = atan2(Q[i], I[i]);
	}

	double temp;
	temp = phase[0];

	for(i = 1; i < n; i++){
		temp = unwrap(temp, phase[i]);
		phase[i] = temp;
	}
}

void phase_compute(const double I[], const double Q[], double phase[], const unsigned long data_length) {
    // Alias de fonction !
    phase_calc2(I, Q, phase, data_length);
}



/*n should be long*/
void subtract(double *a, double *b, double *out, int n)
{
        int i;
        for(i = 0; i < n; i++)
                        out[i] = a[i] - b[i];
}

void add(double *a, double *b, double *out, int n)
{
        int i;
        for(i = 0; i < n; i++)
                out[i] = a[i] + b[i];
}

void amplify(double *in, double *out, double gain, int n)
{
        int i;
        for(i = 0; i < n; i++)
                out[i] = gain * in[i];
}

/* These functions are no longer suported
void phase_calc_no_expand(double *I, double *Q, double *phase, int n)
{
        // ATTENTION: restricted "-PI/2 < phase[i] < PI/2"
        int i;
        for(i = 0; i < n; i++)
                phase[i] = atan(I[i]/ Q[i]);
}

void subtract_compensate(double *a, double *b, double *out, int n)
{
        // ATTENTION: if use this function, conbine with FUNCTION "phase_calc_no_expand()"
        int i;
        out[0] = a[0] - b[0];
        for(i = 0; i < n; i++){
                if(((a[i] - b[i]) - out[i-1]) > 1)
                        out[i] = a[i] - b[i] - PI;
                else
                        out[i] = a[i] - b[i];

                }
        for(i = 1; i < n; i++){
                if((out[i] - out[i-1]) < -1)
                        out[i] += PI;
        }
}
*/

// this func. calculates amplitude noise
// unit of amplitude noise is dBV
// 0dBV = 13dBm @50ohm
// you can get carrier_power P0[dBm] by using func. of powerspectrum2()
void amplitude_calc(double *I, double *Q, int nob, double vfsr, double *amplitude, int n)
{
	int i;
	double  vq = vfsr / (pow(2,nob) -1);
	for(i = 0; i < n; i++)
		//amplitude[i] = (sqrt(pow(I[i], 2) + pow(Q[i], 2)) * 2) / pow(2,16) * v_max;// unit -> [V]
		/*not checked*/
		amplitude[i] = ( sqrt(pow(I[i], 2) + pow(Q[i], 2)) * 2.0 ) / vq; // unit -> [V]
}

// this func. converts dBV -> dBc
// use after cwfx().
void am_noise_conv(double *am_noise, int n, double carrier_power)
{
	int i;
	double temp;
	for(i =0; i < n; i++){
		temp = am_noise[i];
		am_noise[i] = temp + 10*log10(20) - carrier_power;
		//am_noise[i] = temp + 13.0103 - carrier_power;
	}
}

void judgement(int n, int crop, int repeat, int decim, double ovlp, int tap, double fs)
{
        printf("\n");
        printf("read points: %d\n crop points: %d\n repeat: %d\n decim ratio: %d\n overlap: %2.lf\n",
                        n, crop, repeat, decim, ovlp*100.0);

        int i;// counter
        double dec_n[9];// before decimating of decade n
        double repeat_n[9];// minimun repeatable value, i.e. "no overlap"
        double bin_n[9];// freq. bin
        double fh_n[9];// half of Nq freq.
        double chk_n[9];// chk_n[] * crop: actually need points for each decade
        char *judge_n[9];// O: ok /  N: no

        printf("points of pre-decimated: %d\n\n",(int)((n - tap) / decim));
        printf("decade\tpoints\trep_min\t judge\t\tf_bin[Hz]\tf_max[Hz]\n");

        dec_n[8] = n;
        for(i = 7; i > 1; i--){
                dec_n[i] = (int)((dec_n[i+1] - tap) / decim);
                repeat_n[i] = dec_n[i] / (crop * repeat/repeat);

                chk_n[i] = 1+(1-ovlp) * (repeat - 1);
                if(chk_n[i] <= repeat_n[i])
                        judge_n[i] = "O";
                else
                        judge_n[i] = "N";

                fh_n[i] = fs / pow(decim, 8 - i) / 2;
                bin_n[i] = fs / pow(decim, 8 - i) / crop;
                printf("dec%d %9d\t%4d\t > %.0lf %c  \t%8.1lf\t%.2e\n",i, (int)dec_n[i],
                                (int)(repeat_n[i]), chk_n[i]*pow(10,i-3), *judge_n[i], bin_n[i], fh_n[i]);
        }
        printf("\n");
}


void powerspectrum(double *rawa, double *rawb, double *Lf_out, double *freq_out, int n, int crop, int repeat, double ovlp, double fs, char file_name[100])
{
	double freq_bin = 1 / (crop / fs);
	double *cropped_cha = (double *)malloc(sizeof(double) * crop);
	double *cropped_chb = (double *)malloc(sizeof(double) * crop);
	double *temp_fft_xs_result = (double *)malloc(sizeof(double) * crop);
	double *x_spctl = (double *)malloc(sizeof(double) * crop);

	int i, k, y;
	int nc = (crop / 2) + 1;// FFT result points

	fftw_complex *Gyx = fftw_malloc(sizeof(fftw_complex) * nc);// complex cross-spectrum

	double *time = (double *)malloc(sizeof(double) * crop);// use for debug

	// init
	for(i = 0; i < crop; i++){
		cropped_cha[i] = 0;
		cropped_chb[i] = 0;
		time[i] = i * freq_bin;// use for debug
	}

	// crop -> window
	int crop_point_end = crop;
	int a = 0;
	int aa;
	int repeat_counter = 0;
	double *Cospctl = (double *)malloc(sizeof(double) * nc);
	double *Quadspctl = (double *)malloc(sizeof(double) * nc);

	// init
	for(k = 1; k < nc; k++){
		x_spctl[k] = 0;// cross-spectrum
		Cospctl[k] = 0;// co-spectrum
		Quadspctl[k] = 0;// quad-spectrum
	}

	while(repeat_counter < repeat){
		aa = crop - 1;
		while(a < crop_point_end){
			cropped_cha[aa] = rawa[n-1-a];
			cropped_chb[aa] = rawb[n-1-a];
			//printf("%lf %lf %d %d\n", cropped_cha[aa], cropped_chb[aa], a, aa);
			a++;
			aa--;
		}
		hft95_window(cropped_cha, cropped_cha, crop);
		hft95_window(cropped_chb, cropped_chb, crop);
		//gplot(time,"%f", cropped_cha, "%f", crop, "nop", "temp.cropped_cha.dat");
		//gplot(time,"%f", cropped_chb, "%f", crop, "nop", "temp.cropped_cha.dat");

		//start FFT
		fftw_complex *fftout_cha;
		fftw_complex *fftout_chb;
		fftw_plan plan_forward;

		fftout_cha = fftw_malloc(sizeof(fftw_complex) * nc);
		fftout_chb = fftw_malloc(sizeof(fftw_complex) * nc);

		plan_forward = fftw_plan_dft_r2c_1d(crop, cropped_cha, fftout_cha, FFTW_ESTIMATE);
		fftw_execute(plan_forward);
		fftw_destroy_plan(plan_forward);

		plan_forward = fftw_plan_dft_r2c_1d(crop, cropped_chb, fftout_chb, FFTW_ESTIMATE);
		fftw_execute(plan_forward);
		fftw_destroy_plan(plan_forward);


		// calc cross-spectrum
		for(y = 0; y < nc; y++){
			Gyx[y][0] = ((fftout_cha[y][0] * fftout_chb[y][0]) + (fftout_cha[y][1] * fftout_chb[y][1]));
			Gyx[y][1] = ((fftout_chb[y][0] * fftout_cha[y][1]) - (fftout_cha[y][0] * fftout_chb[y][1]));

			Cospctl[y] += Gyx[y][0];
			Quadspctl[y] += Gyx[y][1];
		}
		for(y = 0; y < nc; y++){

			x_spctl[y] =  2 * 2 / pow(crop, 2) * Cospctl[y] / repeat;

		}

		fftw_free(fftout_cha);
		fftw_free(fftout_chb);

		//gplot(time,"%f", x_spctl, "%f", nc, "logx", "temp.x-spectl.dat");

		crop_point_end = crop * (1-ovlp) + a;
		a = crop_point_end - crop;

		repeat_counter++;
	}

	// conversion power spectrum
	for(y = 0; y < nc; y++){
		//x_spctl[y] = (10*log10(x_spctl[y]) +20*log10(2.5 / pow(2,16))+ 10*log10(1000 / 50 / 2));
		Lf_out[y] = (10*log10(x_spctl[y]) +20*log10(2.5 / pow(2,16))+ 10*log10(1000 / 50 / 2));
		freq_out[y] = y * freq_bin;
	}
	save_2_vector_data_plot(freq_out, Lf_out, nc, file_name);
	//gplot(time,"%f", x_spctl, "%f", nc, "logx", file_name);

	free(cropped_cha);
	free(cropped_chb);
	free(temp_fft_xs_result);
	free(x_spctl);
	free(Cospctl);
	free(Quadspctl);
	fftw_free(Gyx);
}

/*
crop = 16384;
for i = 1:crop
        hft95(i) = 1 - 1.9383379*cos(((2.0*i*pi) / crop)) + 1.3045202*cos(2*((2.0*i*pi) / crop)) - 0.4028270*cos(3*((2.0*i*pi) / crop)) + 0.0350665*cos(4*((2.0*i*pi) / crop));
endfor

%ps_s_hann = s(1:crop) .* hanning(crop)';
ps_s_hann = s(1:crop) .* (hft95(crop)/max(hft95(crop)))';
ps_S = fft(ps_s_hann);
ps_PS = ps_S .* conj (ps_S);
ps_PS_norm1 = 2 * 2 * ps_PS(1:crop) / (crop^2);
ps_PS_hardware_normalized = 10*log10(ps_PS_norm1(2:crop/2)) +20*log10(2.5 / 2^16) + 10*log10(1000 / 50 / 2) ;
%ps_PS_hardware_normalized = 10*log10(ps_PS_norm1(2:crop/2)) +20*log10(2.5 / 2^16) + 10*log10(1000 / 50 / 2) -10*log10(1.5);
freq_ps = fs * (2:crop/2) / crop;
carrier_power = max(ps_PS_hardware_normalized)
figure(4);semilogx(freq_ps,ps_PS_hardware_normalized);
max(ps_PS_hardware_normalized)
*/
/*cette fonction est la fonction correcte (l'ad9652), elle doit être réécrite autrement*/
/*en fait ne correspond pas trop au code ci-dessus*/
double powerspectrum2(double *rawa, double *rawb, double *Lf_out, double *freq_out, int n, int crop, int repeat, double ovlp, double fs, char file_name[100])
{
//nb ici repeat = n/N avec n le np points, N la fft_size (=le crop)
	double freq_bin = 1 / (crop / fs);
	printf("freq_bin=%lf\n",freq_bin);
	double *cropped_cha = (double *)malloc(sizeof(double) * crop);
	double *cropped_chb = (double *)malloc(sizeof(double) * crop);
	double *temp_fft_xs_result = (double *)malloc(sizeof(double) * crop);
	double *x_spctl = (double *)malloc(sizeof(double) * crop);
	double peak_power=0.0;
	double peak_freq=0.0;

	int i, k, y;
	int nc = (crop / 2) + 1;// FFT result points

	fftw_complex *Gyx = fftw_malloc(sizeof(fftw_complex) * nc);// complex cross-spectrum

	double *time = (double *)malloc(sizeof(double) * crop);// use for debug

	// init
	for(i = 0; i < crop; i++){
		cropped_cha[i] = 0;
		cropped_chb[i] = 0;
		time[i] = i * freq_bin;// use for debug
	}

	// crop -> window
	int crop_point_end = crop;
	int a = 0;
	int aa;
	int repeat_counter = 0;
	double *Cospctl = (double *)malloc(sizeof(double) * nc);
	double *Quadspctl = (double *)malloc(sizeof(double) * nc);

	// init
	for(k = 1; k < nc; k++){
		x_spctl[k] = 0;// cross-spectrum
		Cospctl[k] = 0;// co-spectrum
		Quadspctl[k] = 0;// quad-spectrum
	}

	while(repeat_counter < repeat){
		aa = crop - 1;
//printf("crop= %d repeat = %d repeat_counter= %d aa= %d  a= %d crop_point_end= %d\n",crop,repeat,repeat_counter,aa,a,crop_point_end);
		while(a < crop_point_end){
			cropped_cha[aa] = rawa[n-1-a];
			cropped_chb[aa] = rawb[n-1-a];
			//printf("%lf %lf %d %d\n", cropped_cha[aa], cropped_chb[aa], a, aa);
			a++;
			aa--;
		}
		hft95_window(cropped_cha, cropped_cha, crop);
		hft95_window(cropped_chb, cropped_chb, crop);
		//gplot(time,"%f", cropped_cha, "%f", crop, "nop", "temp.cropped_cha.dat");
		//gplot(time,"%f", cropped_chb, "%f", crop, "nop", "temp.cropped_cha.dat");

		//start FFT
		fftw_complex *fftout_cha;
		fftw_complex *fftout_chb;
		fftw_plan plan_forward;

		fftout_cha = fftw_malloc(sizeof(fftw_complex) * nc);
		fftout_chb = fftw_malloc(sizeof(fftw_complex) * nc);

		plan_forward = fftw_plan_dft_r2c_1d(crop, cropped_cha, fftout_cha, FFTW_ESTIMATE);
		fftw_execute(plan_forward);
		fftw_destroy_plan(plan_forward);

		plan_forward = fftw_plan_dft_r2c_1d(crop, cropped_chb, fftout_chb, FFTW_ESTIMATE);
		fftw_execute(plan_forward);
		fftw_destroy_plan(plan_forward);


		// calc cross-spectrum
		for(y = 0; y < nc; y++){
			Gyx[y][0] = ((fftout_cha[y][0] * fftout_chb[y][0]) + (fftout_cha[y][1] * fftout_chb[y][1]));
			Gyx[y][1] = ((fftout_chb[y][0] * fftout_cha[y][1]) - (fftout_cha[y][0] * fftout_chb[y][1]));

			Cospctl[y] += Gyx[y][0];
			Quadspctl[y] += Gyx[y][1];
		}
		for(y = 0; y < nc; y++){

			x_spctl[y] =  2 * 2 / pow(crop, 2) * Cospctl[y] / repeat;

		}

		fftw_free(fftout_cha);
		fftw_free(fftout_chb);

		//gplot(time,"%f", x_spctl, "%f", nc, "logx", "temp.x-spectl.dat");

		crop_point_end = crop * (1-ovlp) + a;
		a = crop_point_end - crop;

		repeat_counter++;
	}

	// conversion power spectrum
	for(y = 0; y < nc; y++){
		//x_spctl[y] = (10*log10(x_spctl[y]) +20*log10(2.5 / pow(2,16))+ 10*log10(1000 / 50 / 2));
		//Lf_out[y] = (10*log10(x_spctl[y]) +20*log10(1.41 / pow(2,14))+ 10*log10(1000 / 50 / 2));
		//Lf_out[y] = (10*log10(x_spctl[y]) +20*log10(2.5 / pow(2,16))+ 10*log10(1000 / 50 / 2));
		//we need to pass vfsr and nob, here setup for ad9652
		Lf_out[y] = (10*log10(x_spctl[y]) +20*log10(2.5 / pow(2,16))+ 10*log10(1000 / 50 / 2));
		freq_out[y] = y * freq_bin;
	}
//max/mean
	double temp_power = Lf_out[10]; //cut DC, ret max(Lf_out)
	for(y = 11; y < nc; y++){
		if(temp_power < Lf_out[y]){
			temp_power = Lf_out[y];
			peak_freq = freq_out[y];
		}
	}
	peak_power = temp_power;
	printf("peak power: %lfdBm @%.2eHz\n", peak_power, peak_freq);

	/*here the original function is to save the power spectrum in a file*/
	/*this should be avoided here*/
	save_2_vector_data_plot(freq_out, Lf_out, nc, file_name);
	/**/
	//these may be useful
	//gplot(time,"%f", Lf_out, "%f", nc, "logx", file_name);
	//gplot(time,"%f", x_spctl, "%f", nc, "logx", file_name);

	free(cropped_cha);
	free(cropped_chb);
	free(temp_fft_xs_result);
	free(x_spctl);
	free(Cospctl);
	free(Quadspctl);
	fftw_free(Gyx);

	return(peak_power);
}


/*powerspectrum3!*/
/*renommer  cette fonction, ce n'est pas le spectre en puissance*/
void powerspectrum3(double * data_arm_A, double *data_arm_B, double *Cospctl, int win_size )
{
/*caution : only cos spectrum is given here*/
/* void powerspectrum3cplx(double * data_arm_A, double *data_arm_B, double *Cospctl, double *Quadspctl , int win_size )
would do the trick  */
        int nc = win_size/2+1;
        int k=0;

        fftw_complex *fftout_cha;
        fftw_complex *fftout_chb;

        fftw_plan plan_forward;

        fftout_cha = fftw_malloc(sizeof(fftw_complex) * nc);
        fftout_chb = fftw_malloc(sizeof(fftw_complex) * nc);

        plan_forward = fftw_plan_dft_r2c_1d(win_size, data_arm_A, fftout_cha, FFTW_ESTIMATE);
        fftw_execute(plan_forward);
        fftw_destroy_plan(plan_forward);

        plan_forward = fftw_plan_dft_r2c_1d(win_size, data_arm_B, fftout_chb, FFTW_ESTIMATE);
        fftw_execute(plan_forward);
        fftw_destroy_plan(plan_forward);

        for (k=0 ; k<nc;k++){
                Cospctl[k] = ((fftout_cha[k][0] * fftout_chb[k][0]) + (fftout_cha[k][1] * fftout_chb[k][1]));
        }
}


/*normalize*/
void normalize_psd_dBc(  double * Cospctl , const int win_size , const double fs , double * psd){

	int k=0;
	double hann_norm = 1.7609125905568124;
	double freq_bin_norm = 10.0 * log10(fs/(double)win_size);
	double PS=0.0;

        for (k=0 ; k<win_size/2 ; k++){
		PS = 2.0 * 2.0 * 2.0 * fabs(Cospctl[k] / (  pow(win_size,2)));
		psd[k] = 10.0*log10(PS) - hann_norm - freq_bin_norm;
        }

}
/*peak_power*/

/*caution*/
void powerSpectralDensity0 (  double * armx , double * army , const int win_size , const double fs , double * psd)
{
	int k=0;
	//not for  the moment, may  induce NaN
	/* hann_norm = 10.0 * log10(1.5) */
	double hann_norm = 1.7609125905568124;
	double freq_bin_norm = 10.0 * log10(fs/(double)win_size);
	double * P = allocate_double(win_size);
	double * PS = allocate_double(win_size/2);

	double * armx_han = allocate_double(win_size);
	double * army_han = allocate_double(win_size);

	hanning_window(armx,armx_han,win_size);
	hanning_window(army,army_han,win_size);

	powerspectrum3(armx_han,army_han,P,win_size);
	//powerspectrum3(armx,army,P,win_size);
	//storeDouble("P",3,P);

        for (k=0 ; k<win_size/2 ; k++){
                //PS[k]=2.0*fabs(P[k])/(2.0*(double)win_size*fs);
                //PS[k]=2.0*fabs(P[k])/((double)win_size*fs);
                //psd[k]=10.0*log10(PS[k]) + hann_norm; //!!c'est -
		PS[k] = 2.0 * 2.0 * 2.0 * fabs(P[k] / (  pow(win_size,2)));
		psd[k] = 10.0*log10(PS[k]) - hann_norm - freq_bin_norm;
        }
	//storeDouble("psd",win_size/2,psd);

	free(P);
	free(PS);

}

// Crop->Window->FFT->X_spectrum
// rawa and rawb is phase data in time domain
// ovlp is enabled when repeat >= 2 (if repeat == 1, ovlp is ignored)
// If you use this function as measuring of ADC noise floor, you must subtract 3dB after calculation.
void cwfx(double *rawa, double *rawb, double *Lf_out, double *freq_out, int n, int crop, int repeat, double ovlp, double fs, char file_name[100])
{
	double freq_bin = 1 / (crop / fs);
	double *cropped_cha = (double *)malloc(sizeof(double) * crop);
	double *cropped_chb = (double *)malloc(sizeof(double) * crop);
	double *temp_fft_xs_result = (double *)malloc(sizeof(double) * crop);
	double *x_spctl = (double *)malloc(sizeof(double) * crop);

	int i, k, y;
	int nc = (crop / 2) + 1;// FFT result points

	fftw_complex *Gyx = fftw_malloc(sizeof(fftw_complex) * nc);// complex cross-spectrum

	//double *time = (double *)malloc(sizeof(double) * crop);// use for debug

	// init
	for(i = 0; i < crop; i++){
		cropped_cha[i] = 0;
		cropped_chb[i] = 0;
	//	time[i] = i * freq_bin;// use for debug
	}

	// crop -> window
	int crop_point_end = crop;
	int a = 0;
	int aa;
	int repeat_counter = 0;
	double *Cospctl = (double *)malloc(sizeof(double) * nc);
	double *Quadspctl = (double *)malloc(sizeof(double) * nc);

	// init
	for(k = 1; k < nc; k++){
		x_spctl[k] = 0;// cross-spectrum
		Cospctl[k] = 0;// co-spectrum
		Quadspctl[k] = 0;// quad-spectrum
	}

	while(repeat_counter < repeat){
		aa = crop - 1;
		while(a < crop_point_end){
			cropped_cha[aa] = rawa[n-1-a];
			cropped_chb[aa] = rawb[n-1-a];
			//printf("%lf %lf %d %d\n", cropped_cha[aa], cropped_chb[aa], a, aa);
			a++;
			aa--;
		}
		hanning_window(cropped_cha, cropped_cha, crop);
		hanning_window(cropped_chb, cropped_chb, crop);
		//gplot(time,"%f", cropped_cha, "%f", crop, "nop", "temp.cropped_cha.dat");
		//gplot(time,"%f", cropped_chb, "%f", crop, "nop", "temp.cropped_cha.dat");

		//start FFT
		fftw_complex *fftout_cha;
		fftw_complex *fftout_chb;
		fftw_plan plan_forward;

		fftout_cha = fftw_malloc(sizeof(fftw_complex) * nc);
		fftout_chb = fftw_malloc(sizeof(fftw_complex) * nc);

		plan_forward = fftw_plan_dft_r2c_1d(crop, cropped_cha, fftout_cha, FFTW_ESTIMATE);
		fftw_execute(plan_forward);
		fftw_destroy_plan(plan_forward);

		plan_forward = fftw_plan_dft_r2c_1d(crop, cropped_chb, fftout_chb, FFTW_ESTIMATE);
		fftw_execute(plan_forward);
		fftw_destroy_plan(plan_forward);


		// calc cross-spectrum
		for(y = 0; y < nc; y++){
			Gyx[y][0] = ((fftout_cha[y][0] * fftout_chb[y][0]) + (fftout_cha[y][1] * fftout_chb[y][1]));
			Gyx[y][1] = ((fftout_chb[y][0] * fftout_cha[y][1]) - (fftout_cha[y][0] * fftout_chb[y][1]));

			Cospctl[y] += Gyx[y][0];
			Quadspctl[y] += Gyx[y][1];
		}
		for(y = 0; y < nc; y++){
			// (1/0.5): coef. of hanning window
			//x_spctl[y] = (1/0.5) * 2 * 2 / pow(crop, 2) * sqrt(pow(Cospctl[y], 2) + pow(Quadspctl[y], 2)) / repeat;// |Sxy|
			//x_spctl[y] = (1/0.5) * 2 * 2 / pow(crop, 2) * sqrt(pow(Cospctl[y], 2) + pow(Quadspctl[y], 2)) * cos(atan(Quadspctl[y] / Cospctl[y])) / repeat;// |Re{Sxy}|
			//x_spctl[y] = (1/0.5) * 2 * 2 / pow(crop, 2) * sqrt(pow(Cospctl[y], 2)) / repeat;// |Re{Sxy}|
			x_spctl[y] = (1/0.5) * 2 * 2 / pow(crop, 2) * Cospctl[y] / repeat;// Re{Sxy}
		}

		fftw_free(fftout_cha);
		fftw_free(fftout_chb);

		//gplot(time,"%f", x_spctl, "%f", nc, "logx", "temp.x-spectl.dat");

		// overlap process //
		crop_point_end = crop * (1-ovlp) + a;
		a = crop_point_end - crop;

		repeat_counter++;
	}

	// conversion L(f)
	for(y = 0; y < nc; y++){
		Lf_out[y] = (-10*log10(2) + 10*log10(x_spctl[y]) - 10*log10(freq_bin) - 10*log10(1.5));
		freq_out[y] = y * freq_bin;
	}
	save_2_vector_data_plot(freq_out, Lf_out, nc, file_name);
	//gplot(time,"%f", x_spctl, "%f", nc, "logx", file_name);

	free(cropped_cha);
	free(cropped_chb);
	free(temp_fft_xs_result);
	free(x_spctl);
	free(Cospctl);
	free(Quadspctl);
	fftw_free(Gyx);
}
















void cwfx_infinite(double *rawa, double *rawb, double *Lf_out, double *freq_out, int n, int crop, int repeat, double ovlp, double fs, char file_name[100])
{
	double freq_bin = 1 / (crop / fs);
	double *cropped_cha = (double *)malloc(sizeof(double) * crop);
	double *cropped_chb = (double *)malloc(sizeof(double) * crop);
	double *temp_fft_xs_result = (double *)malloc(sizeof(double) * crop);
	double *x_spctl = (double *)malloc(sizeof(double) * crop);

	int i, k, y;
	int nc = (crop / 2) + 1;// FFT result points

	fftw_complex *Gyx = fftw_malloc(sizeof(fftw_complex) * nc);// complex cross-spectrum

	double *time = (double *)malloc(sizeof(double) * crop);// use for debug

	// init
	for(i = 0; i < crop; i++){
		cropped_cha[i] = 0;
		cropped_chb[i] = 0;
		time[i] = i * freq_bin;// use for debug
	}

	// crop -> window
	int crop_point_end = crop;
	int a = 0;
	int aa;
	int repeat_counter = 0;
	double *Cospctl = (double *)malloc(sizeof(double) * nc);
	//double *Quadspctl = (double *)malloc(sizeof(double) * nc);

	// init
	for(k = 1; k < nc; k++){
		x_spctl[k] = 0;// cross-spectrum
		Cospctl[k] = 0;// co-spectrum
		//Quadspctl[k] = 0;// quad-spectrum
	}

	while(repeat_counter < repeat){
		aa = crop - 1;
		while(a < crop_point_end){
			cropped_cha[aa] = rawa[n-1-a];
			cropped_chb[aa] = rawb[n-1-a];
			//printf("%e %e %d %d\n", cropped_cha[aa], cropped_chb[aa], a, aa);
			a++;
			aa--;
		}
		printf("\n");
		//gplot(time,"%f", cropped_cha, "%f", crop, "nop", "temp.cropped_cha.dat");
		hanning_window(cropped_cha, cropped_cha, crop);
		hanning_window(cropped_chb, cropped_chb, crop);
		//bug gplot violation memoire
		//gplot(time,"%f", cropped_cha, "%f", crop, "nop", "temp.cropped_cha.dat");
		//gplot(time,"%f", cropped_chb, "%f", crop, "nop", "temp.cropped_cha.dat");

		//start FFT
		fftw_complex *fftout_cha;
		fftw_complex *fftout_chb;
		fftw_plan plan_forward;

		fftout_cha = fftw_malloc(sizeof(fftw_complex) * nc);
		fftout_chb = fftw_malloc(sizeof(fftw_complex) * nc);

		plan_forward = fftw_plan_dft_r2c_1d(crop, cropped_cha, fftout_cha, FFTW_ESTIMATE);
		fftw_execute(plan_forward);
		fftw_destroy_plan(plan_forward);

		plan_forward = fftw_plan_dft_r2c_1d(crop, cropped_chb, fftout_chb, FFTW_ESTIMATE);
		fftw_execute(plan_forward);
		fftw_destroy_plan(plan_forward);


		// calc cross-spectrum
		for(y = 0; y < nc; y++){
			Gyx[y][0] = ((fftout_cha[y][0] * fftout_chb[y][0]) + (fftout_cha[y][1] * fftout_chb[y][1]));
//			Gyx[y][1] = ((fftout_chb[y][0] * fftout_cha[y][1]) - (fftout_cha[y][0] * fftout_chb[y][1]));

			Cospctl[y] += Gyx[y][0];
//			Quadspctl[y] += Gyx[y][1];
		}
		for(y = 0; y < nc; y++){
			// (1/0.5): coef. of hanning window
			//x_spctl[y] = (1/0.5) * 2 * 2 / pow(crop, 2) * sqrt(pow(Cospctl[y], 2) + pow(Quadspctl[y], 2)) / repeat;// |Sxy|
			//x_spctl[y] = (1/0.5) * 2 * 2 / pow(crop, 2) * sqrt(pow(Cospctl[y], 2) + pow(Quadspctl[y], 2)) * cos(atan(Quadspctl[y] / Cospctl[y])) / repeat;// |Re{Sxy}|
			//x_spctl[y] = (1/0.5) * 2 * 2 / pow(crop, 2) * sqrt(pow(Cospctl[y], 2)) / repeat;// |Re{Sxy}|
			x_spctl[y] = (1/0.5) * 2 * 2 / pow(crop, 2) * Cospctl[y] / repeat;// Re{Sxy}
		}

		fftw_free(fftout_cha);
		fftw_free(fftout_chb);

		//gplot(time,"%f", x_spctl, "%f", nc, "logx", "temp.x-spectl.dat");

		// overlap process //
		crop_point_end = crop * (1-ovlp) + a;
		a = crop_point_end - crop;

		repeat_counter++;
	}


	// conversion L(f)
	for(y = 0; y < nc; y++){
		//Lf_out[y] = (-10*log10(2) + 10*log10(x_spctl[y]) - 10*log10(freq_bin) - 10*log10(1.5));
		Lf_out[y] = x_spctl[y];
		freq_out[y] = y * freq_bin;
	}
	save_2_vector_data_plot_e(freq_out, Lf_out, nc, file_name);
	//gplot(time,"%f", x_spctl, "%f", nc, "logx", file_name);

	free(cropped_cha);
	free(cropped_chb);
	free(temp_fft_xs_result);
	free(x_spctl);
	free(Cospctl);
	//free(Quadspctl);
	fftw_free(Gyx);
}

void convert_to_lf(double *Lf_lin, double *Lf_log, double *freq, int n, int corr_factor, char file_name[100])// calc. L(f) dBc / Hz
{
	int y;
	double temp;
	double temp2;
	for(y = 0; y < n; y++){
		//temp = (-10*log10(2) + 10*log10(Lf_lin[y]) - 10*log10(freq[1] - freq[0]) - 10*log10(1.5));
		temp = Lf_lin[y];
		//temp2 = (-10*log10(2) + 10*log10(temp/corr_factor) - 10*log10(freq[1] - freq[0]) - 10*log10(1.5));
		temp2 = (10*log10(temp/corr_factor) - 10*log10(freq[1] - freq[0]) - 10*log10(1.5));
		Lf_log[y] = temp2;
		//printf("%e\n", Lf_log[y]);
	}
	//save_2_vector_data_plot(freq, Lf_log, n, file_name);
}

void convert_to_lf_dBrad2(double *Lf_lin, double *Lf_log, double *freq, int n, int corr_factor, char file_name[100])// calc. Sphi dBrad^2 / Hz
{
	int y;
	double temp;
	double temp2;
	for(y = 0; y < n; y++){
		//temp = (-10*log10(2) + 10*log10(Lf_lin[y]) - 10*log10(freq[1] - freq[0]) - 10*log10(1.5));
		temp = Lf_lin[y];
		temp2 = (10*log10(temp/corr_factor) - 10*log10(freq[1] - freq[0]) - 10*log10(1.5));
		Lf_log[y] = temp2;
		//printf("%e\n", Lf_log[y]);
	}
	save_2_vector_data_plot(freq, Lf_log, n, file_name);
}

void convert_to_lf_dBV2(double *Lf_lin, double *Lf_log, double *freq, int n, int corr_factor, double fsr, int nob, double gain, char file_name[100])// calc. dBV^2 / Hz
{
	int y;
	double temp;
	double temp2;
	for(y = 0; y < n; y++){
		//temp = (-10*log10(2) + 10*log10(Lf_lin[y]) - 10*log10(freq[1] - freq[0]) - 10*log10(1.5));
		temp = Lf_lin[y];
		temp2 = 20*log10(fsr/pow(2, nob)) + (10*log10(temp/corr_factor) - 10*log10(freq[1] -
freq[0]) - 10*log10(1.5)) - 20*log10(gain); //-3 = single channel must be added into the gain
		Lf_log[y] = temp2;
		//printf("%e\n", Lf_log[y]);
	}
	//save_2_vector_data_plot(freq, Lf_log, n, file_name);
}

void convert_to_am(double *Lf_lin, double *Lf_log, double *freq, double carrier_power, int n, int corr_factor, char file_name[100])
{
	int y;
	double temp;
	double temp2;
	for(y = 0; y < n; y++){
		//temp = (-10*log10(2) + 10*log10(Lf_lin[y]) - 10*log10(freq[1] - freq[0]) - 10*log10(1.5));
		temp = Lf_lin[y];
		temp2 = (-10*log10(2) + 10*log10(temp/corr_factor) - 10*log10(freq[1] - freq[0]) - 10*log10(1.5)) + 10*log10(20) - carrier_power;
		Lf_log[y] = temp2;
		//printf("%e\n", Lf_log[y]);
	}
	save_2_vector_data_plot(freq, Lf_log, n, file_name);
}

void convert_to_am_dBrad2(double *Lf_lin, double *Lf_log, double *freq, double carrier_power, int n, int corr_factor, char file_name[100])
{
	int y;
	double temp;
	double temp2;
	for(y = 0; y < n; y++){
		//temp = (-10*log10(2) + 10*log10(Lf_lin[y]) - 10*log10(freq[1] - freq[0]) - 10*log10(1.5));
		temp = Lf_lin[y];
		temp2 = (10*log10(temp/corr_factor) - 10*log10(freq[1] - freq[0]) - 10*log10(1.5)) + 10*log10(20) - carrier_power;
		Lf_log[y] = temp2;
		//printf("%e\n", Lf_log[y]);
	}
	save_2_vector_data_plot(freq, Lf_log, n, file_name);
}

void hogehoge(double *out, double *in, int n, double average_factor)
{
	int y;
	//double temp;
	for(y = 0; y < n; y++){
		out[y] = out[y] + in[y] / average_factor;
	}
}

double get_nominal_frequency_from_phase(double *phase, int n, double dt, double f_nco)
{
        // phase phi(t) [rad] -> instantaneous frequency f(t) [Hz] //
        int k;
        double f0=0.0;
        double tau0 = dt;
        double *df = (double *)malloc(sizeof(double) * n);
        for(k = 0; k < n-1; ++k){
                df[k] = (phase[n-1-k] - phase[n-2-k]) / (2*PI*tau0);
                f0 += df[k] / (n - 1);
        }
        free(df);
        f0 = f_nco - f0;
        printf("frequency %.12e[Hz]\n", f0);
        return(f_nco - f0);
}

double mean(double * vec , unsigned long size){

	unsigned long k;
	double accum = 0.0;

	for ( k=0 ; k < size ; k++){
		accum += vec[k];
	}

	return accum  / (double) k ;

}

void remove_mean(double * vec , unsigned long size){

	double mean_vec = mean(vec,size);

	unsigned long k;

	for ( k = 0 ; k < size ; k++ ) {
		vec[k] -= mean_vec;
	}

}


double vector_max(double * v , const unsigned long size){

	unsigned long k;
	double max=.0;

	for ( k = 0 ; k < size ; k++ ) {
		if ( v[k] > max) max =  v[k];
	}

return max;
}
