//@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
#include "orbit.h"

double testFunc(double t){
    std::cout<<"\nInside testFunc"<<std::endl;
    t = t*300;
    return t;
};

void Orbit::anomalyCalc(double ecc, double T, double Tc,double P, double epoch){
	//------------------
	//Calculate TA and E
	//Remember that in RV, there is a phase shift due to the Tc!=T, that doesn't exist in DI.
	//------------------
	//std::cout<<"\necc = "<<ecc<<", T = "<<T<<", Tc = "<<Tc<<", P = "<<P<<", epoch = "<<epoch<<std::endl;
	M = (2.0*pi*(epoch-2.0*T+Tc))/(P*daysPerYear);
	multiples = (double)((int)(M/(2.0*pi)));
	M -= multiples*(2.0*pi);//shift into [-360,360]
	if (M<0)
		M+=2.0*pi;//shift into [0,360]
	//If e=0 (circular), then theta and E are both equal to M.
	thetaRV=M;
	E=M;
	if (ecc!=0){
		//for RV
		if ((M!=0)and(M!=(2.0*pi))){
			Eprime = M+ecc*sin(M)+((ecc*ecc)/(2.0*M))*sin(2.0*M);
			newtonCount = 0;
			while ( (fabs(E-Eprime)>1.0e-10)&&(newtonCount<50) ){
				E = Eprime;
				Eprime = E-((E-ecc*sin(E)-M)/(1.0-ecc*cos(E)));
				newtonCount +=1;
			}
			//double check it satisfies the original equation
			if (fabs((E-ecc*sin(E))-M)>1.0e-5){
				std::cout<<"PROBLEM!! resulting E from Newton's loop isn't within error limit!!!"<<std::endl;
				if (true){
					std::cout<<"M = "<<M <<std::endl;
					std::cout<<"e = "<<ecc<<std::endl;
					std::cout<<"T = "<<T<<std::endl;
					std::cout<<"Tc = "<<Tc<<std::endl;
					std::cout<<"P = "<<P<<std::endl;
					std::cout<<"Eprime = "<<Eprime <<"\n" <<std::endl;
				}
			}
			//std::cout<<"E RV [deg] = "<<E*(180.0/pi)<<std::endl;
			thetaPrime = acos((cos(E)-ecc)/(1.0-ecc*cos(E)));
			if (E>pi)
				thetaPrime = 2.0*pi-thetaPrime;
			thetaRV = thetaPrime;
			//std::cout<<"theta RV [deg] = "<<thetaRV*(180.0/pi)<<std::endl;
		}
		//for DI
		if (T!=Tc){
			M = (2.0*pi*(epoch-T))/(P*daysPerYear);
			multiples = (double)((int)(M/(2.0*pi)));
			M -= multiples*(2.0*pi);//shift into [-360,360]
			if ((M!=0)and(M!=(2.0*pi))){
				Eprime = M+ecc*sin(M)+((ecc*ecc)/(2.0*M))*sin(2.0*M);
				newtonCount = 0;
				while ( (fabs(E-Eprime)>1.0e-10)&&(newtonCount<50) ){
					E = Eprime;
					Eprime = E-((E-ecc*sin(E)-M)/(1.0-ecc*cos(E)));
					newtonCount +=1;
				}
				//double check it satisfies the original equation
				if (fabs((E-ecc*sin(E))-M)>1.0e-5){
					std::cout<<"PROBLEM!! resulting E from Newton's loop isn't within error limit!!!"<<std::endl;
					if (true){
						std::cout<<"M = "<<M <<std::endl;
						std::cout<<"e = "<<ecc<<std::endl;
						std::cout<<"T = "<<T<<std::endl;
						std::cout<<"Tc = "<<Tc<<std::endl;
						std::cout<<"P = "<<P<<std::endl;
						std::cout<<"Eprime = "<<Eprime <<"\n" <<std::endl;
					}
				}
			}
		}
	}
	EDI = E;
	//std::cout<<"\nin anomaly calc: EDI = "<<EDI<<", thetaRV = "<<thetaRV<<std::endl;//$$$$$$$$$$$$$$$$$$
};

void Orbit::loadStaticVars(double omegaoffsetDI,double omegaoffsetRV,bool lowEcc_in,bool PASA_in){
	omegaOffsetDI = omegaoffsetDI;
	omegaOffsetRV = omegaoffsetRV;
	lowEcc = lowEcc_in;
	PASA = PASA_in;
//	if (lowEcc)
//		std::cout<<"lowEcc is true"<<std::endl;
//	if (lowEcc==false)
//		std::cout<<"lowEcc is false"<<std::endl;
};

void Orbit::loadRealData(double *xx, int xx_nx, int xx_ny){
	if (false)
		std::cout<<"\nInside loadData function"<<std::endl;
    dataRealAry = xx;
    dataRealAry_nx = xx_nx;
    dataRealAry_ny = xx_ny;
    if (false)
    	std::cout<<"data loaded!"<<std::endl;
};

void Orbit::loadConstants(double Grav_in,double pi_in,double KGperMsun_in, double daysPerYear_in,double secPerYear_in,double MperAU_in){
	Grav = Grav_in;
	pi = pi_in;
	KGperMsun = KGperMsun_in;
	daysPerYear = daysPerYear_in;
	secPerYear = secPerYear_in;
	MperAU = MperAU_in;
	if (false)
		std::cout<<"constants loaded!"<<std::endl;
};

void Orbit::convertParsFromRaw(double *p, int p_n){
	/*
	 Convert sqrt(e)sin(omega)&sqrt(e)cos(omega) => e and omega if required.
	 Else, do nothing.
	 */
	if (lowEcc){
		//Only calc if non-circular
		if ((p[4]!=0)and(p[9]!=0)){
			//std::cout<<"converting params to use forms"<<std::endl;
			e = p[4]*p[4]+p[9]*p[9];
			omega = (180.0/pi)*atan2(p[4],p[9]);
			//if (omega<0)
			//	omega+=360.0;
			p[4] = e;
			p[9] = omega;
		}
	}
};

void Orbit::convertParsToRaw(double *p, int p_n){
	/*
	Convert e and omega => sqrt(e)sin(omega)&sqrt(e)cos(omega) if required.
	Else, do nothing.
	 */
	if (lowEcc){
		//std::cout<<"converting params to Raw forms"<<std::endl;
		e = p[4];
		omega = p[9];
		//Only calc if non-circular
		if (e!=0){
			p[4] = sqrt(e)*sin((pi/180.0)*omega);
			p[9] = sqrt(e)*cos((pi/180.0)*omega);
		}
	}
};

void Orbit::calculate(double *yy, int yy_nx, int yy_ny, double *p, int p_n){
    /*
    The calculator function to perform the primary
    orbit calculations for the C++ Orbit object.
    */
	dataModelAry=yy;
	dataModelAry_nx=yy_nx;
	dataModelAry_ny=yy_ny;
	params = p;
	params_n = p_n;
	bool verbose=false;
//	if (verbose)
//		std::cout<<"\nInside Orbit calculator function"<<std::endl;//$$$$$$$$$$$$$$$$$$$$$$$$$
//	if (verbose){//$$$$$$$$$$$$$$$$$$$$$$$$$
//		std::cout<<"real data inside Orbit c++ code:"<<std::endl;//$$$$$$$$$$$$$$$$$$$$$$$$$
//		for (int i=0; i<dataRealAry_nx; i++){//$$$$$$$$$$$$$$$$$$$$$$$$$
//			std::cout<<"[";//$$$$$$$$$$$$$$$$$$$$$$$$$
//			for (int j=0;j<dataRealAry_ny;j++){//$$$$$$$$$$$$$$$$$$$$$$$$$
//				std::cout<<dataRealAry[j+i*dataRealAry_ny]<<", ";//$$$$$$$$$$$$$$$$$$$$$$$$$
//			}//$$$$$$$$$$$$$$$$$$$$$$$$$
//			std::cout<<"]"<<std::endl;//$$$$$$$$$$$$$$$$$$$$$$$$$
//		}//$$$$$$$$$$$$$$$$$$$$$$$$$
//	}

	//--------------------------------------------------------------------------
	//------------------------------- 	START REAL CALCULATE STEPS -------------
	//--------------------------------------------------------------------------
	//calc those that are static for each epoch
	atot=0;
	if (params[0]!=0)
		atot =pow(((params[7]*params[7]*secPerYear*secPerYear*Grav*KGperMsun*(params[0]+params[1]))/(4.0*pi*pi)),(1.0/3.0));
	params[10]=atot/MperAU;
	if ((dataRealAry[6]<1e6)&&(params[12]==0)){
		//NOTE: both forms were tested and showed to be the same to the 1e-11 level.
		//      Thus, the same to the limit of rounding errors.  Using semi-major version as it is simpler
		//Semi-major axis version
		K = ((2.0*pi*(atot/(1.0+(params[0]/params[1])))*sin(params[8]*(pi/180.0)))/(params[7]*secPerYear*pow((1.0-params[4]*params[4]),(1.0/2.0))));
		//Masses version
		//K = pow((2.0*pi*Grav)/(params[7]*secPerYear),(1.0/3.0))*((params[1]*KGperMsun)/pow((params[0]+params[1])*KGperMsun,(2.0/3.0)))*(sin(params[8]*(pi/180.0))/pow((1.0-params[4]*params[4]),(1.0/2.0)));
		params[12]=K;
	}
	//Get the model version of each omega and shift into [0,360]
	omegaDI = params[9]+omegaOffsetDI;
	omegaRV = params[9]+omegaOffsetRV;
	if (omegaDI>360.0)
		omegaDI-=360.0;
	if (omegaDI<0.0)
		omegaDI+=360;
	if (omegaRV>360.0)
		omegaRV-=360.0;
	if (omegaRV<0)
		omegaRV+=360.0;
	//Calculate Tc <-> T if needed
	if (params[5]!=params[6]){
		//if T=Tc already, do nothing.
		if ((params[4]==0)||(omegaRV==90.0)||(omegaRV==270.0)){
			//Circular, so just set equal.
			if (params[6]==0)
				params[6]=params[5];
			else
				params[5]=params[6];
		}
		else{
			ta = pi/2.0 - omegaRV*(pi/180.0);
			halfE = atan2(sqrt(1.0-params[4])*sin(ta/2.0),sqrt(1.0+params[4])*cos(ta/2.0));
			mTTc = 2.0*halfE-params[4]*sin(2.0*halfE);
			deltaT = (mTTc*params[7]*daysPerYear)/(2.0*pi);
			if (params[6]==0)
				params[6] = params[5]+deltaT;
			else
				params[5] = params[6]-deltaT;
			if (false)
				std::cout<<"T = "<<params[5]<<", params[9] = "<<params[9]<<", ta = "<<ta*(180.0/pi)<<", halfE = "<< halfE*(180.0/pi)<<", mTTc = "<<mTTc*(180.0/pi)<<", deltaT = "<<deltaT<<", Tc = "<<params[6]<<std::endl;
		}
	}
	//start loop over each epoch of data
	for (i=0;i<dataModelAry_nx; i++){
		if (verbose)//$$$$$$$$$$$$$$$$$$$$$$$$$
			std::cout<<"\nepoch "<<dataRealAry[0+i*dataRealAry_ny]<<":"<<std::endl;//$$$$$$$$$$$$$$$$$$$$$$$$$
		//Calculate true anomaly for RV and eccentric anomaly for DI
		anomalyCalc(params[4],params[5],params[6],params[7],dataRealAry[0+i*dataRealAry_ny]);
		//std::cout<<"in calc: EDI = "<<EDI<<", thetaRV = "<<thetaRV<<std::endl;//$$$$$$$$$$$$$$$$$$
		//--------------------------
		//Calculate RV
		//--------------------------
		if (dataRealAry[6+i*dataRealAry_ny]<1e6)
			dataModelAry[2+i*dataModelAry_ny]=params[12]*(cos(thetaRV+omegaRV*(pi/180.0))+params[4]*cos(omegaRV*(pi/180.0)))+params[13+int(dataRealAry[7+i*dataRealAry_ny])];
		else
			dataModelAry[2+i*dataModelAry_ny]=0.0;
		if (verbose){//$$$$$$$$$$$$$$$$$$$$$$$$$
			std::cout<<"RV = "<<dataModelAry[2+i*dataModelAry_ny] <<std::endl;//$$$$$$$$$$$$$$$$$$$$$$$$$
		}//$$$$$$$$$$$$$$$$$$$$$$$$$
		//--------------------------
		//Calculate x,y or PA,SA
		//--------------------------
		if ((dataRealAry[2+i*dataRealAry_ny]<1e6)&&(dataRealAry[4+i*dataRealAry_ny]<1e6)){
			// calculate all the Thiele-Innes constants in ["]
			A = ((atot/MperAU)*(params[2]/1000.0))*(cos(params[3]*(pi/180.0))*cos(omegaDI*(pi/180.0))-sin(params[3]*(pi/180.0))*sin(omegaDI*(pi/180.0))*cos(params[8]*(pi/180.0)));
			B = ((atot/MperAU)*(params[2]/1000.0))*(sin(params[3]*(pi/180.0))*cos(omegaDI*(pi/180.0))+cos(params[3]*(pi/180.0))*sin(omegaDI*(pi/180.0))*cos(params[8]*(pi/180.0)));
			F = ((atot/MperAU)*(params[2]/1000.0))*(-cos(params[3]*(pi/180.0))*sin(omegaDI*(pi/180.0))-sin(params[3]*(pi/180.0))*cos(omegaDI*(pi/180.0))*cos(params[8]*(pi/180.0)));
			G = ((atot/MperAU)*(params[2]/1000.0))*(-sin(params[3]*(pi/180.0))*sin(omegaDI*(pi/180.0))+cos(params[3]*(pi/180.0))*cos(omegaDI*(pi/180.0))*cos(params[8]*(pi/180.0)));
			// The coordinates of the unit orbital ellipse in the true plane (Binnendijk)
			X = cos(EDI)-params[4];
			Y = sqrt(1.0-params[4]*params[4])*sin(EDI);
			// Calculate the predicted x&y in ["], or PA[deg], SA["]
			//KEY NOTE: x_TH-I = y_plot = North = Dec = A*X +F*Y
			//          y_TH-I = x_plot = East = RA = B*X +G*Y
			RA = B*X +G*Y;
			Dec = A*X +F*Y;
			//check which and then store RA/Dec or PA/SA accordingly
			if (PASA){
				PA = atan2(RA,Dec);
				if (PA<0.0)
					PA+=2.0*pi;
				SA = sqrt(RA*RA+Dec*Dec);
				dataModelAry[0+i*dataModelAry_ny] = PA*(180.0/pi);
				dataModelAry[1+i*dataModelAry_ny] = SA;
			}
			else{
				dataModelAry[0+i*dataModelAry_ny] = RA;
				dataModelAry[1+i*dataModelAry_ny] = Dec;
			}
		}
		else{
			//std::cout<<"x real = "<< dataRealAry[1+i*dataRealAry_nx]<<", y real = "<< dataRealAry[3+i*dataRealAry_nx]<<std::endl;//$$$$$$$$$$$$$$$$$$$$$$$$$
			dataModelAry[0+i*dataModelAry_ny] = dataModelAry[1+i*dataModelAry_ny] = 0.0;
		}
		if (verbose){//$$$$$$$$$$$$$$$$$$$$$$$$$
			std::cout<<"x (or PA) = "<<dataModelAry[0+i*dataModelAry_ny] <<std::endl;//$$$$$$$$$$$$$$$$$$$$$$$$$
			std::cout<<"y (or SA) = "<<dataModelAry[1+i*dataModelAry_ny] <<std::endl;//$$$$$$$$$$$$$$$$$$$$$$$$$
		}//$$$$$$$$$$$$$$$$$$$$$$$$$
	}
	if (verbose){//$$$$$$$$$$$$$$$$$$$$$$$$$
		std::cout<<"\nModel data (after calculating all model epochs) inside Orbit c++ code:"<<std::endl;//$$$$$$$$$$$$$$$$$$$$$$$$$
		for (i=0; i<dataModelAry_nx; i++){//$$$$$$$$$$$$$$$$$$$$$$$$$
			std::cout<<"[";//$$$$$$$$$$$$$$$$$$$$$$$$$
			for (j=0;j<dataModelAry_ny;j++){//$$$$$$$$$$$$$$$$$$$$$$$$$
				std::cout<<dataModelAry[j+i*dataModelAry_ny]<<", ";//$$$$$$$$$$$$$$$$$$$$$$$$$
			}//$$$$$$$$$$$$$$$$$$$$$$$$$
			std::cout<<"]"<<std::endl;//$$$$$$$$$$$$$$$$$$$$$$$$$
		}//$$$$$$$$$$$$$$$$$$$$$$$$$
	}


}


















