/* file: postctools.i */
%module postctools
%{
#include "postctools.h"
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
%apply (double* IN_ARRAY1, int DIM1) {(double *par, int par_n)};
//pass in 2D array of dynamic size, that can NOT change in CPP func
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double *zz, int zz_nx, int zz_ny)};
%include "std_string.i"
%include "postctools.h"
