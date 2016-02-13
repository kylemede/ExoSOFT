/* file: hello.i */
%module hello
%{
#include "hello.h"
#include "string"
%}

%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
import_array();
%}
%apply (double* IN_ARRAY1, int DIM1) {(double *x, int n)};
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double *x, int nx, int ny)};
%apply (int DIM1, double* IN_ARRAY1){(int nx, double* x), (int ny, double* y)};
%include "std_string.i"
%include "hello.h"
