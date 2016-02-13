/* file: orbCalc2.i */
%module orbCalc2
%{
#include "orbCalc2.h"
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
%apply (double* IN_ARRAY1, int DIM1) {(double *x, int n)};
//pass in 2D array of dynamic size, that can NOT change in CPP func
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double *xx, int nx, int ny)};
//pass in 2D array of dynamic size, that CAN change in CPP func
%apply (double* INPLACE_ARRAY2,int DIM1, int DIM2) {(double *yy, int nx, int ny)};
//%apply (double* IN_ARRAY2, int DIM1, int DIM2,double* INPLACE_ARRAY2,int DIM1, int DIM2){(double *xx,int xx_nx, int xx_ny,double *zz,int zz_nx, int zz_ny)};
%include "std_string.i"
%include "orbCalc2.h"
