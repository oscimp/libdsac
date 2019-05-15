#include <stdio.h>
#include <stdlib.h>

// #define TAB_SIZE 8
// #define V2 1

void compute_xBarre_xMxB_xSquare(const double x[], size_t data_size,
				 double *xBarre, double xMxB[], double *xSquare)
{
	size_t i;
	*xBarre = 0;
	*xSquare = 0;
	for (i = 0; i < data_size; i++)
		*xBarre += (double)x[i];
	*xBarre /= data_size;

	for (i = 0; i < data_size; i++) {
		xMxB[i] = (double)x[i] - *xBarre;
		*xSquare += xMxB[i] * xMxB[i];
	}
	*xSquare /= data_size;
}

void reg_lin(double *a, double *b, const double y[], const double x_moins_x_barre[], size_t tab_size,
	     double x_barre, size_t data_div, double sum_div)
{
	size_t i;
	double y_barre = 0;
	double ty;
	double sum_fin = 0;
	/* first part direct on acquisition */
	for (i = 0; i < tab_size; i++)
		y_barre += y[i];
	y_barre /= data_div;

	/* second part, data reloading */
	for (i = 0; i < tab_size; i++) {
		ty = y[i] - y_barre;
		sum_fin += x_moins_x_barre[i] * ty;
	}
	sum_fin /= data_div;

	ty = sum_fin / sum_div;
	*a = ty;
	*b = y_barre - (ty * x_barre);
}



#if 0
void reg_lin(float *a, float *b, float *y, float *x_moins_x_barre, int tab_size,
	     float x_barre, int data_div, float sum_div)
{
	int i;
	float y_barre = 0;
	float ty;
	float sum_fin = 0;
	/* first part direct on acquisition */
	for (i = 0; i < tab_size; i++)
		y_barre += y[i];
	y_barre /= data_div;

	/* second part, data reloading */
	for (i = 0; i < tab_size; i++) {
		ty = y[i] - y_barre;
		sum_fin += x_moins_x_barre[i] * ty;
	}
	sum_fin /= data_div;

	ty = sum_fin / sum_div;
	*a = ty;
	*b = y_barre - (ty * x_barre);
}

void compute_xBarre_xMxB_xSquare(int *x, int data_size,
				 float *xBarre, float *xMxB, float *xSquare)
{
	int i;
	*xBarre = 0;
	*xSquare = 0;
	for (i = 0; i < data_size; i++)
		*xBarre += (float)x[i];
	*xBarre /= data_size;

	for (i = 0; i < data_size; i++) {
		xMxB[i] = (float)x[i] - *xBarre;
		*xSquare += xMxB[i] * xMxB[i];
	}
	*xSquare /= data_size;

}

int main(void)
{
	int x[TAB_SIZE];
	float y[TAB_SIZE];
	int i;
	float a, b;

	float x_barre;
	float x_moins_x_barre[TAB_SIZE];

	float xSquare;

	FILE *data = fopen("a.dat", "r");
	if (data == NULL) {
		printf("erreur d'ouverture du fichier\n");
		return EXIT_FAILURE;
	}

	for (i = 0; i < TAB_SIZE; i++)
		fscanf(data, "%d\t%f\n", (int *)&x[i], &y[i]);
	fclose(data);

	compute_xBarre_xMxB_xSquare(x, TAB_SIZE, &x_barre, x_moins_x_barre,
				    &xSquare);
	printf("%f %f\n", x_barre, xSquare);
	reg_lin(&a, &b, y, x_moins_x_barre, TAB_SIZE, x_barre, TAB_SIZE,
		xSquare);

	printf("\n %f %f\n", a, b);
#if 0
	double y_barre;
	double x_m_x_b[TAB_SIZE];
	double y_m_y_b[TAB_SIZE];
	double sum_square_x_m_x_b = 0;
	double sum_y_m_y_b = 0.f;
	double sum_fin = 0.f;
	y_barre /= (float)TAB_SIZE;
	printf("%lf %lf\n", x_barre, y_barre);

	for (i = 0; i < TAB_SIZE; i++) {
		x_m_x_b[i] = x[i] - x_barre;
		y_m_y_b[i] = y[i] - y_barre;
		sum_square_x_m_x_b += (x_m_x_b[i] * x_m_x_b[i]);
		//sum_x_m_x_b += (x_m_x_b[i]);
		//sum_square_y_m_y_b += (y_m_y_b[i] * y_m_y_b[i]);
		sum_y_m_y_b += y_m_y_b[i];
		sum_fin += x_m_x_b[i] * y_m_y_b[i];
		printf("%lf %lf\n", x_m_x_b[i], y_m_y_b[i]);
	}
	sum_y_m_y_b /= (float)TAB_SIZE;
	sum_square_x_m_x_b /= (float)TAB_SIZE;
	sum_fin /= (float)TAB_SIZE;
	printf("\n");
	printf("%.24e %.24e %lf\n", sum_y_m_y_b, sum_square_x_m_x_b, sum_fin);
	printf("\n");

	a = sum_fin / (sum_square_x_m_x_b);
	b = y_barre - (a * x_barre);
	printf("%lf %lf\n", a, b);
	printf("\n");
	printf("%f %f %f %f %f\n", sum_fin, sum_square_x_m_x_b,
	       y_barre, a, x_barre);
#endif
	return EXIT_SUCCESS;
}
#endif