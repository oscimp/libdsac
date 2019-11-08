#include "conversions.h"

#include <stdlib.h>
#include <math.h>

double dbm2vp(double x){
        /*sqrt(50.0*10^(-3+0.1*x))*sqrt(2);*/

        //printf("dbm2vp(%lf)= %f\n",x,10.0 * sqrt( pow(10,(-3+0.1*x)) ));
        return (10.0 * sqrt( pow(10,(-3.0+0.1*x)) ));
}
