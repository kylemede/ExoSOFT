//@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
#if !defined(POSTCTOOLS_H)
#define POSTCTOOLS_H 1
#include<string>
#include <math.h>
#include <stdio.h>
#include <iostream>

class PostCtools{
public:
	//internal params
	int i,i_last,j,n,colNum;
	double s,ep;
	double ave,var,sum,Sx,Sxx,x;
	double varALL,halfVarALL,numCorrLengths,corrLengthsTotal;
	//input data ary
	double* data;
	int data_nx,data_ny;
	//funcs
	void loadParamData(double *zz, int zz_nx, int zz_ny);
	void sumCalc(int startPoint,int lastPoint);
	void meanCalc(int startPoint, int lastPoint);
	double varianceCalc(int startPoint, int lastPoint);
	double corrLenCalc(int parNum);
};
#endif
