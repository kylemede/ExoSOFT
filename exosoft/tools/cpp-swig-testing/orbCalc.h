#if !defined(ORBCALC_H)
#define ORBCALC_H 1
#include<string>
#include <math.h>
#include <stdio.h>
#include <iostream>

double sum2d_double(double *xx, int nx, int ny); 
void modAndReturn2d_double(double *yy,int nx, int ny);
double modAndReturnDouble(double x);
void speedTestSWIG(int ntimes);
#endif
