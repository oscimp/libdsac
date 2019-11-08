#include "file_handler.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>

int storeComplex64Bits(char *filename, int dataSize, long long *dataI, long long *dataQ)
{
	int i;
	FILE *res_fd;

	res_fd = fopen(filename,"w");

	for (i = 0; i < dataSize; i++) {
		fprintf(res_fd, "%Ld %Ld\n", dataI[i], dataQ[i]);
	}

	fflush(res_fd);
	fclose(res_fd);

	return EXIT_SUCCESS;
}

int storeSimple64Bits(char *filename, int dataSize, long long *data)
{
	int i;
	FILE *res_fd;

	res_fd = fopen(filename,"w");

	for (i = 0; i < dataSize; i++) {
		fprintf(res_fd, "%Ld\n", data[i]);
	}

	fflush(res_fd);
	fclose(res_fd);

	return EXIT_SUCCESS;
}

int readContent(short *tab, char *filename, long nb, int offset)
{
	int i, tmp;
	FILE *fd = fopen(filename, "r");

	if (fd == NULL) {
		printf("erreur d'ouverture de %s\n", filename);
		return -1;
	}

	if (offset != 0){
		for (i = 0; i < offset; i++)
			fscanf(fd, "%d\n", &tmp);
	}
	for (i = 0; i < nb; i++) {
		if (fscanf(fd, "%d\n", &tmp) == EOF)
			return -1;
		tab[i] = (short)tmp;
	}
	fclose(fd);
	return 0;
}

int storeDouble(char  * filename , int dataSize , double * data){
        FILE *fd;
        int i=0;

        fd = fopen(filename,"w");

        if (fd == NULL) {
                printf("error while opening %s\n", filename);
                return -1;
        }

        for (i = 0; i < dataSize; i++) {
                fprintf(fd, "%24.16lf\n", data[i]);
        }

        fflush(fd);
        fclose(fd);

	printf("stored: %s\n",filename);

        return EXIT_SUCCESS;
}

int storeInt64(char  * filename , int64_t * data , unsigned long dataSize){
        FILE *fd;
        unsigned long i=0;

        fd = fopen(filename,"w");

        if (fd == NULL) {
                printf("error while opening %s\n", filename);
                return -1;
        }

        for (i = 0; i < dataSize; i++) {
                fprintf(fd, "%lli\n", (long long int)data[i]);
        }

        fflush(fd);
        fclose(fd);

	printf("stored: %s\n",filename);

        return EXIT_SUCCESS;
}

int store8Int64(char  * filename , int64_t * a , int64_t * b, int64_t * c, int64_t * d, int64_t * e, int64_t * f, int64_t * g, int64_t * h ,unsigned long dataSize){
        FILE *fd;
        unsigned long i=0;

        fd = fopen(filename,"w");

        if (fd == NULL) {
                printf("error while opening %s\n", filename);
                return -1;
        }

        for (i = 0; i < dataSize; i++) {
                fprintf(fd, "%lli\t%lli\t%lli\t%lli\t%lli\t%lli\t%lli\t%lli\n", (long long int)a[i],(long long int)b[i],(long long int)c[i],(long long int)d[i],(long long int)e[i],(long long int)f[i],(long long int)g[i],(long long int)h[i]);
                //fprintf(fd, "%Ld\n", (int16_t)data[i]);
        }

        fflush(fd);
        fclose(fd);

	printf("stored: %s\n",filename);

        return EXIT_SUCCESS;
}

int store8Int32(char  * filename , int32_t * a , int32_t * b, int32_t * c, int32_t * d, int32_t * e, int32_t * f, int32_t * g, int32_t * h ,unsigned long dataSize){
        FILE *fd;
        unsigned long i=0;

        fd = fopen(filename,"w");

        if (fd == NULL) {
                printf("error while opening %s\n", filename);
                return -1;
        }

        for (i = 0; i < dataSize; i++) {
                fprintf(fd, "%li\t%li\t%li\t%li\t%li\t%li\t%li\t%li\n", (long int)a[i],(long int)b[i],(long int)c[i],(long int)d[i],(long int)e[i],(long int)f[i],(long int)g[i],(long int)h[i]);
                //fprintf(fd, "%Ld\n", (int16_t)data[i]);
        }

        fflush(fd);
        fclose(fd);

	printf("stored: %s\n",filename);

        return EXIT_SUCCESS;
}

int store8Int32fromDouble(char  * filename , double * a , double * b, double * c, double * d, double * e, double * f, double * g, double * h ,unsigned long dataSize){
        FILE *fd;
        unsigned long i=0;

        fd = fopen(filename,"w");

        if (fd == NULL) {
                printf("error while opening %s\n", filename);
                return -1;
        }

        for (i = 0; i < dataSize; i++) {
                fprintf(fd, "%li\t%li\t%li\t%li\t%li\t%li\t%li\t%li\n", (long int)a[i],(long int)b[i],(long int)c[i],(long int)d[i],(long int)e[i],(long int)f[i],(long int)g[i],(long int)h[i]);
                //fprintf(fd, "%Ld\n", (int16_t)data[i]);
        }

        fflush(fd);
        fclose(fd);

	printf("stored: %s\n",filename);

        return EXIT_SUCCESS;
}

int readDataInt32to2int16d(char * filename, long nb_elements, double *adc1, double *adc2)
{
        long i = 0;
        int32_t data;
        FILE * fd = fopen (filename,"rb");
        if ( fd == NULL){
                printf ("Error: Open %s failed -- %d\n", filename, errno);
                return EXIT_FAILURE;
        }
        for( i =0 ; i< nb_elements ; i++){

                if (fread (&data, sizeof (int32_t), 1, fd) != (size_t) 1){
                        if (!feof (fd)){
                                printf ("Error: fread failed -- %d\n", errno);
                                return -1;
                        }
                        else break;
                }

                adc1[i] = (int16_t)(0xffff&(data));
                adc2[i] = (int16_t)(0xffff&(data>>16));
        }

        fclose(fd);
return EXIT_SUCCESS;
}

int readDataInt64to4int16d(char * filename, long nb_elements, double *adc1, double *adc2 , double *adc3, double *adc4)
{
        long i = 0;
        int64_t data;
	int32_t ch12;
	int32_t ch34;
        FILE * fd = fopen (filename,"rb");
        if ( fd == NULL){
                printf ("Error: Open %s failed -- %d\n", filename, errno);
                return EXIT_FAILURE;
        }
        for( i =0 ; i< nb_elements ; i++){

                if (fread (&data, sizeof (int64_t), 1, fd) != (size_t) 1){
                        if (!feof (fd)){
                                printf ("Error: fread failed -- %d\n", errno);
                                return -1;
                        }
                        else break;
                }

		ch12 = (int32_t)(0xFFFFFFFF&(data));
		ch34 = (int32_t)(0xFFFFFFFF&(data>>32));

                adc1[i] = (int16_t)(0xffff&(ch12));
                adc2[i] = (int16_t)(0xffff&(ch12>>16));
		adc3[i] = (int16_t)(0xffff&(ch34));
		adc4[i] = (int16_t)(0xffff&(ch34>>16));
        }

        fclose(fd);
return EXIT_SUCCESS;
}
