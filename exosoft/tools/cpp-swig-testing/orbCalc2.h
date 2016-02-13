#if !defined(ORBCALC2_H)
#define ORBCALC2_H 1
#include<string>
#include <math.h>
#include <stdio.h>
#include <iostream>

class OrbCalcObj{
public:
    //obj variables
    double testDouble;
    double* dataAry;
    int dataAry_x;
    int dataAry_y;
    
    //funcs
    void loadData(double *xx, int nx, int ny);
    void calculator(double  *yy, int nx, int ny);
};
double testFunc(double t);


#endif
