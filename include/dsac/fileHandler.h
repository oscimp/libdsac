#ifndef __FILE_HANDLER__
#define __FILE_HANDLER__ 

int readContent(short *tab, char *filename2, long nb, int offset);
int storeComplex64Bits(char *filename, int dataSize, long long *dataI, long long *dataQ);
int storeSingle64Bits(char *filename, int dataSize, long long *data);
int storeDouble(char  * filename , int dataSize , double * data); //order should be changed
int storeInt64(char  * filename , int64_t * data , unsigned long dataSize);
int store8Int64(char  * filename , int64_t * a , int64_t * b, int64_t * c, int64_t * d, int64_t * e, int64_t * f, int64_t * g, int64_t * h,unsigned long dataSize);
int store8Int32(char  * filename , int32_t * a , int32_t * b, int32_t * c, int32_t * d, int32_t * e, int32_t * f, int32_t * g, int32_t * h ,unsigned long dataSize);
int store8Int32fromDouble(char  * filename , double * a , double * b, double * c, double * d, double * e, double * f, double * g, double * h ,unsigned long dataSize);
int readDataInt32to2int16d(char * filename, long nb_elements, double *adc1, double *adc2);
int readDataInt64to4int16d(char * filename, long nb_elements, double *adc1, double *adc2 , double *adc3, double *adc4);

#endif /* __FILE_HANDLER__ */
