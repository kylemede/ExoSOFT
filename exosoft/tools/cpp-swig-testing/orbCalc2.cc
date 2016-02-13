#include "orbCalc2.h"

double testFunc(double t){
    std::cout<<"\nInside testFunc"<<std::endl;
    t = t*300;
    return t;
}

void OrbCalcObj::loadData(double *xx, int nx, int ny){
    std::cout<<"\nInside loadData function"<<std::endl;
    dataAry = xx;
    dataAry_x = nx;
    dataAry_y = ny;
    std::cout<<"data loaded!"<<std::endl;
};

void OrbCalcObj::calculator(double  *yy, int nx, int ny){
    /*
    The calculator function to perform the primary 
    orbit calculations for the C++ OrbCalcObj object.
    */
    std::cout<<"\nInside OrbCalcObj calculator function"<<std::endl;
    std::cout<<"testDouble = "<<testDouble<<std::endl;
    testDouble = testDouble*10.0;
    std::cout<<"testDouble = "<<testDouble<<std::endl;
    std::cout<<"For STATIC dataAry:"<<std::endl;
    for (int i=0; i<dataAry_x; i++){
		for (int j=0;j<dataAry_y;j++){
	        std::cout<<"[i,j = ["<<i<<","<<j<<"] = "<<dataAry[j+i*dataAry_y]<<std::endl;
		}
	}
	std::cout<<"For INPLACE yy:"<<std::endl;
	int ints;
	double ints2;
	for (int i=0; i<nx; i++){
		for (int j=0;j<ny;j++){
		ints = i+j;
		std::cout<<"ints = "<<ints<<std::endl;
		ints2 = (double)ints;
		std::cout<<"ints2 = "<<ints2<<std::endl;
		yy[j+i*ny]=ints2*1.0+1.0;
	        std::cout<<"[i,j] = ["<<i<<","<<j<<"] = "<<yy[j+i*ny]<<std::endl;
		}
	}
	std::cout<<"testFunc provided = "<<testDouble<<std::endl;
	double r = testFunc(testDouble);
	std::cout<<"testFunc returned = "<<r<<std::endl;
	
}


