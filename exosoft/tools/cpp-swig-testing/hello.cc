#include "hello.h"

double cadd(double x, double y){
	return x+y;
}

void speak(std::string const& str)
{
	std::cout << "Rober says: " << str << std::endl;
}

int selection(bool flag)
{
	if(flag){
		return 1;
	} else {
		return 0;
	}
}

int fact(int n){
	if (n<=1) return 1;
	else return n*fact(n-1);
}

double csin(double i){
	return sin(i);
}

double sum_scalar(double *x, int n, int y){
	double xsum=0;
	for (int i=0; i<n; i++){
		xsum+=x[i]+y;
	}
	return xsum;

}

double sum_1d(int nx, double *x,  int ny, double *y){
	double xsum=0;
	for (int i=0; i<nx; i++){
		xsum+=x[i]+y[i];
	}
	return xsum;

}

double sum2d(double *x, int nx, int ny){
	double xsum=0;
	for (int i=0; i<nx; i++){
		for (int j=0;j<ny;j++){
		xsum+=x[j+i*ny];
		}
	}
	return xsum;
}

double sumfunc(int nx, double *x, int ny, double *y, double (*op) (double,double)){
	double xsum=0;
	printf("nx=%d,ny=%d",nx,ny);
	for (int i=0; i<nx; i++){
		xsum+=op(x[i],y[i]);
	}
	return xsum;

}

/*this would always fail in python, at least till now */
double sum2d_double(double **x, int nx, int ny){
	double xsum=0;
	for (int i=0; i<nx; i++){
		for (int j=0;j<ny;j++){
		xsum+=x[i][j];
		}
	}
	return xsum;
}
/*
int main(){
	speak("hello world");
}*/
