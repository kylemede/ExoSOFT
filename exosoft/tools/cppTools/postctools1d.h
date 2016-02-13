//@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
#if !defined(POSTCTOOLS1D_H)
#define POSTCTOOLS1D_H 1
#include<string>
#include <math.h>
#include <stdio.h>
#include <iostream>

class PostCtools1D{
public:
	//internal params
	int i,i_last,j,n;
	double s,ep;
	double ave,var,sum,Sx,Sxx,x;
	double varALL,halfVarALL,numCorrLengths,corrLengthsTotal;
	//input data ary
	double* data;
	int data_nx;
	//funcs
	void loadParamData(double *z, int z_nx);
	void sumCalc(int startPoint,int lastPoint);
	void meanCalc(int startPoint, int lastPoint);
	double varianceCalc(int startPoint, int lastPoint);
	double corrLenCalc();
};
#endif
