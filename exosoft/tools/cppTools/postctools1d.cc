//@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
#include "postctools1d.h"


void PostCtools1D::loadParamData(double *z, int z_nx){
	//std::cout<<"in loadParamData"<<std::endl;
	data = z;
	data_nx = z_nx;
	//std::cout<<"data_nx = "<<data_nx<<std::endl;
};

void PostCtools1D::sumCalc(int startPoint,int lastPoint)
{
	sum = 0.0;
	//loop through all data points to get total value
	for (j=startPoint;j<(lastPoint+1);j++){
		if ((isnan(sum))||(sum>1e200))
			break;
		else
			sum+=data[j];
	}
};

void PostCtools1D::meanCalc(int startPoint, int lastPoint)
{
	/**
	 * A very simple function that simply adds up the total of all the elements
	 * in a vector<double> and divides it by the number of elements to give
	 * the mean value.  The precision of the calculation is checked to make sure
	 * the result is sufficiently precise.
	 */
	ave = 0;
	sumCalc(startPoint,lastPoint);
	ave = sum/(double(lastPoint-startPoint+1));
};

double PostCtools1D::varianceCalc(int startPoint, int lastPoint)
{
	/**
	 * This will calculate the "bias-corrected sample variance"
	 * and uses an advanced "corrected two-pass algorithm" to reduce roundoff error,
	 * a modified version of function on pg 728 of Numerical Recipes 3rd.
	 */
	n=lastPoint;
	meanCalc(startPoint,lastPoint);
	var=ep=0.0;
	if (n>1){
		//loop through all points to load up sums needed for "corrected two-pass algorithm", eq 14.1.8 pg 724
		for (j=startPoint;j<n;j++)
		{
			s=data[j]-ave;
			ep+=s;
			var+=s*s;
		}
		var=(var-ep*ep/n)/(n-1);
	}
	return var;
};

double PostCtools1D::corrLenCalc(){
	/**
	 * Calculates the average correlation length
	 * of the input parameter's data in a boxcar style.
	 * The correlation length follows that described in Tegmark2004.
	 * ie. The step in the chain where the variance is half the total chain's variance.
	 * We perform this repeatedly across the chain to give a more accurate average.
	 */
	varALL = varianceCalc(0,data_nx);
	halfVarALL = varALL/2.0;
	numCorrLengths=0.0;
	corrLengthsTotal=0.0;
	i_last=0;
	Sx=0.0;
	Sxx=0.0;
	for (i=0;i<data_nx;i++){
		x=data[i];
		Sx+=x;
		Sxx+=x*x;
		var = (Sxx/double(i-i_last+1))-(Sx/double(i-i_last+1))*(Sx/double(i-i_last+1));
		if (var>halfVarALL){
			corrLengthsTotal+=double(i-i_last+1);
			i_last=i+1;
			numCorrLengths+=1.0;
			Sx=0;
			Sxx=0;
		}
	}
	return corrLengthsTotal/numCorrLengths;
};
