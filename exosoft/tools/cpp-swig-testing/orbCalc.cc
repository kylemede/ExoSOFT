#include "orbCalc.h"


double sum2d_double(double *xx, int nx, int ny){
	//locations in numpy 2d arrays are broadcast to a 1D vector in cpp
	//thus: x[i][j] in Python => x[j+i*ny] on the cpp side!!
	double xsum=0;
	for (int i=0; i<nx; i++){
		for (int j=0;j<ny;j++){
		xsum+=xx[j+i*ny];
	        std::cout<<"[i,j = ["<<i<<","<<j<<"] = "<<xx[j+i*ny]<<std::endl;
		}
	}
	return xsum;
}
void modAndReturn2d_double(double *yy,int nx, int ny){
	//locations in numpy 2d arrays are broadcast to a 1D vector in cpp
	//thus: x[i][j] in Python => x[j+i*ny] on the cpp side!!
	std::cout<<"inside modAndReturn2d_double func"<<std::endl;
	for (int i=0; i<nx; i++){
		for (int j=0;j<ny;j++){
		int ints = i+j;
		std::cout<<"ints = "<<ints<<std::endl;
		double ints2 = (double)ints;
		std::cout<<"ints2 = "<<ints2<<std::endl;
		yy[j+i*ny]=ints2*1.0+1.0;
	        std::cout<<"[i,j] = ["<<i<<","<<j<<"] = "<<yy[j+i*ny]<<std::endl;
		}
	}
}

double modAndReturnDouble(double x){
	x+=10.0;
	return x;

}
void speedTestSWIG(int ntimes){
	
	for(int i=0;i<ntimes;i++){
		sin(1.0);	
	}
}	

/*
int main(){
	speak("hello world");
}*/
