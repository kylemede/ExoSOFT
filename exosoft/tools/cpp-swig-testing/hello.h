#if !defined(HELLO_H)
#define HELLO_H 1
#include<string>
#include <math.h>
#include <stdio.h>
#include <iostream>

typedef double (*FuncPtr)(double x, double y);
double cadd(double x, double y);
int selection(bool flag);
int fact(int n);
double csin(double i);
void speak(std::string const& str);
double sum_scalar(double *x, int n, int y); 
double sum_1d(int nx,double *x, int ny,double *y);
double sum2d(double *x, int nx, int ny); 
double sumfunc(int nx,double *x, int ny, double *y, double (*op) (double,double));
double sum2d_double(double **x, int nx, int ny); 
#endif
