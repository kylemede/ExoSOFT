/* file: orbit.i */
%module orbit
%{
#include "orbit.h"
#include "string"
%}

%{
#define SWIG_FILE_WITH_INIT
%}
/*include numpy typemaps*/
%include "numpy.i"
/*initialize module*/
%init %{
import_array();
%}
/*typemaps for different non-arrays*/
//%apply double *OUTPUT {double, *result};
/*typemaps for different arrays*/
//pass in 1D array of dynamic size, that can NOT change in CPP func
%apply (double* IN_ARRAY1, int DIM1) {(double *x, int x_n)};
//pass in 2D array of dynamic size, that can NOT change in CPP func
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double *xx, int xx_nx, int xx_ny)};
//pass in 1D array of dynamic size, that CAN change in CPP func
%apply (double* INPLACE_ARRAY1,int DIM1) {(double *p, int p_n)};
//pass in 2D array of dynamic size, that CAN change in CPP func
%apply (double* INPLACE_ARRAY2,int DIM1, int DIM2) {(double *yy, int yy_nx, int yy_ny)};
//pass in a 2D array of dynamic size, that CAN change in CPP func, along with a 1D array that CAN change
%apply (double* INPLACE_ARRAY2,int DIM1, int DIM2, double* INPLACE_ARRAY1, int DIM1) {(double *yy, int yy_nx, int yy_ny, double *p, int p_n)};
//%apply (double* IN_ARRAY2, int DIM1, int DIM2,double* INPLACE_ARRAY2,int DIM1, int DIM2){(double *xx,int xx_nx, int xx_ny,double *zz,int zz_nx, int zz_ny)};
%include "std_string.i"
%include "orbit.h"
