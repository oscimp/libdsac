#include "piosbr_mod.h"
/*exit: marche pas pour test unitaire (return NULL;) mais l'utilisateur doit tester la fonction (if) */
double * allocate_double(int size){
double *vec=(double *)calloc(size, sizeof(double));
        if (vec == NULL) {
		perror("(dsalib:allocate_double)");
                exit(EXIT_FAILURE);
        }

	return(vec);
}


void dprint(double *vec, int size, char *buf)
{
	int i;
	if (strcmp(buf,"none")!=1){
	for (i = 1; i < size + 1; i++)
	printf("%24.16lf\n", vec[i]);
	}
	else
	for (i = 1; i < size + 1; i++)
		printf("%s[%d]=%24.16lf\n", buf, i, vec[i]);
}

void read_vector(double *vec, char *filename)
{
	int i=0;
	FILE *fd;
	fd = fopen(filename, "r");
	while (fscanf(fd, "%lf", &vec[i]) != EOF)
		i++;
	fclose(fd);
	return;
}

void read_1_vector(double * vec, char * filename, int size)//timaike mod. Juillet 2014
{
	int i=size;
	FILE *fd;
	fd=fopen(filename,"r");
	while (i != 0){
		i--;
		fscanf(fd, "%lf", &vec[i]);
		}
	fclose(fd);
	return;
}

void read_2_vectors(double * vec1, double *vec2, char * filename, int size)
{
	int i=size;
	FILE *fd;
	fd=fopen(filename,"r");
	while (i != 0){
		i--;
		fscanf(fd, "%lf %lf", &vec1[i],&vec2[i]);
		}
	fclose(fd);
	return;
}

void newline(void) {printf("\n");}


void adc_raw_2_volt_stereo(double *x, double *y, int n, double vmax,
                           double vmin, int M, char *filename)
{
        FILE *fd;
        int i = 0;
        int rawx, rawy;
        double vlsb = (vmax - vmin) / pow(2, M);
        double offset = (double)(pow(2, M) / 2);
        printf("vlsb=%lf\n", vlsb);
        if ((fd = fopen(filename, "r")) == NULL) {
                fprintf(stderr, "could'nt open the file %s\n", filename);
                exit(EXIT_FAILURE);
        }
        //while (fscanf(fd, "%d %d", rawx,rawy) != EOF){
        while (i != n) {
                fscanf(fd, "%d %d\n", &rawx, &rawy);
                x[i] = ((double)rawx - offset) * vlsb;
                y[i] = ((double)rawy - offset) * vlsb;
                //printf("%d\t%d-->%lf\t%d-->%lf\n",i,rawx,x[i],rawy,y[i]);
                i++;
        }
        fclose(fd);
        return;
}


void show(double *x, int n, char *buf)
{
        int i;
        printf("\n ---%s-----\n", buf);
        for (i = 0; i < n; i++) {
                printf("%d\t%.24lf\n", i, x[i]);
        }

}

//void showc(double *x, int n, char *buf)
//{
//        int i;
//        printf("\n ---%s-----\n", buf);
//        for (i = 0; i < n; i++) {
//                printf("%d\t%.24lf\t%.24lf\n", i, x[i][0],x[i][1]);
//        }
//
//}
void gplot(double *f_log, char *format_f, double *p_log, char *format_p, int n,
	   char *logmode, char *title)
{
	FILE *tmpfile;
	int i;
	char format[20];
	sprintf(format, "%s\t%s\n", format_f, format_p);
//store result in file
	tmpfile = fopen("yourfile.dat", "w");
	if (strcmp(logmode, "logx") == 0) {
		for (i = 0; i < n - 1; i++) {
			//fprintf(tmpfile,"%24.16e\t%24.16e\n",f_log[i],p_log[i]);
			fprintf(tmpfile, format, f_log[i], p_log[i]);
		}

		fclose(tmpfile);
	} else {
		for (i = 0; i < n; i++) {
			//fprintf(tmpfile,"%24.16e\t%24.16e\n",f_log[i],p_log[i]);
			fprintf(tmpfile, format, f_log[i], p_log[i]);
		}

		fclose(tmpfile);

	}

	FILE *pipe = popen("gnuplot -persist", "w");
	/*il faut pouvoir choisir l'option (par ex GNUPLOT_TERMINAL*/
	fprintf(pipe, "set t x11\n");
	//fprintf(pipe, "set terminal 'wxt' size 1024,512\n");
	//fprintf(pipe, "set mouse\n");
	if (strcmp(logmode, "logx") == 0) {
		fprintf(pipe, "set logsc x\n");
		fprintf(pipe, "set mytics 2\n");
		fprintf(pipe, "set mxtics 10\n");
		fprintf(pipe, "set grid mxtics lt 0 lw 2 lc rgb '#00FFFF'\n");
		fprintf(pipe, "set grid mytics lt 0 lw 2 lc rgb '#00FFFF'\n");
	}
	//fprintf(pipe, "plot 'yourfile.dat' using 1:2 lt 1 pt 2 w l title '%s'\n", title);
	//for x11
	fprintf(pipe, "plot 'yourfile.dat' using 1:2 lt 3 pt 2 w l title '%s'\n", title);

	fclose(pipe);
}
