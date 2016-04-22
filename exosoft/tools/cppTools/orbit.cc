//@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
#include "orbit.h"

double testFunc(double t){
    std::cout<<"\nInside testFunc"<<std::endl;
    t = t*300;
    return t;
};

void Orbit::anomalyCalc(double epoch){
	/*
	Calculate TA and E.
	Remember that in RV, there is occasionally a phase shift due to the Tc!=To,
	that doesn't exist in DI.
	*/
	M = (2.0*pi*(epoch-2.0*To+Tc))/(P*daysPerYear);
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
				if (warningsOn){
					std::cout<<"PROBLEM!! resulting E from Newton's loop isn't within error limit!!!"<<std::endl;
					if (true){
						std::cout<<"M = "<<M <<std::endl;
						std::cout<<"e = "<<ecc<<std::endl;
						std::cout<<"To = "<<To<<std::endl;
						std::cout<<"Tc = "<<Tc<<std::endl;
						std::cout<<"P = "<<P<<std::endl;
						std::cout<<"Eprime = "<<Eprime <<"\n" <<std::endl;
					}
				}
			}
			thetaPrime = acos((cos(E)-ecc)/(1.0-ecc*cos(E)));
			if (E>pi)
				thetaPrime = 2.0*pi-thetaPrime;
			thetaRV = thetaPrime;
		}
		//for DI
		if (To!=Tc){
			M = (2.0*pi*(epoch-To))/(P*daysPerYear);
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
					if (warningsOn){
						std::cout<<"PROBLEM!! resulting E from Newton's loop isn't within error limit!!!"<<std::endl;
						if (true){
							std::cout<<"M = "<<M <<std::endl;
							std::cout<<"e = "<<ecc<<std::endl;
							std::cout<<"To = "<<To<<std::endl;
							std::cout<<"Tc = "<<Tc<<std::endl;
							std::cout<<"P = "<<P<<std::endl;
							std::cout<<"Eprime = "<<Eprime <<"\n" <<std::endl;
						}
					}
				}
			}
		}
	}
	EDI = E;
};

void Orbit::loadStaticVars(double omegaoffsetDI,double omegaoffsetRV,bool lowEcc_in,bool PASA_in){
	omegaOffsetDI = omegaoffsetDI;
	omegaOffsetRV = omegaoffsetRV;
	lowEcc = lowEcc_in;
	PASA = PASA_in;
};

void Orbit::loadRealData(double *xx, int xx_nx, int xx_ny){
    dataRealAry = xx;
    dataRealAry_nx = xx_nx;
    dataRealAry_ny = xx_ny;
};

void Orbit::loadConstants(double Grav_in,double pi_in,double KGperMsun_in, double daysPerYear_in,double secPerYear_in,double MperAU_in){
	Grav = Grav_in;
	pi = pi_in;
	KGperMsun = KGperMsun_in;
	daysPerYear = daysPerYear_in;
	secPerYear = secPerYear_in;
	MperAU = MperAU_in;
};

void Orbit::convertParsFromRaw(double *p, int p_n){
	/*
	 Convert sqrt(e)sin(omega)&sqrt(e)cos(omega) => e and omega if required.
	 Else, do nothing.
	 */
	if (lowEcc){
		//Only calc if non-circular
		if ((p[4]!=0)and(p[9]!=0)){
			ecc = p[4]*p[4]+p[9]*p[9];
			omega = (180.0/pi)*atan2(p[4],p[9]);
			p[4] = ecc;
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
		ecc = p[4];
		omega = p[9];
		//Only calc if non-circular
		if (ecc!=0){
			p[4] = sqrt(ecc)*sin((pi/180.0)*omega);
			p[9] = sqrt(ecc)*cos((pi/180.0)*omega);
		}
	}
};

void Orbit::NewtonWarningsOn(bool warningsOn_in){
	/*
	 At the beginning of SA with lowEcc mode, there can be some errors
	 during the newton's method to calculate the anomalies.  This occurs as
	 when SA has yet to find a good starting point and still has bad values
	 for sqrt(e)sin(omega) and sqrt(e)cos(omega).  These warnings
	 are off by default and need to be turned on after this beginning phase.

	 warningOn_in=true will turn them back on.
	 */
	warningsOn = warningsOn_in;
};

void Orbit::parsAryToVariables(){
	/*
	A function to take array of parameters used by simulator.py into more
	cleanly named parameters for use inside orbit.
	 */
	m1 = params[0];
	m2 = params[1];
	parallax = params[2];
	Omega = params[3];
	ecc = params[4];
	To = params[5];
	Tc = params[6];
	P = params[7];
	inc = params[8];
	omega = params[9];
	K = params[12];
};

void Orbit::variablesToParsAry(){
	/*
	A function to push the variables calculated during orbit back into the
	array used in simulator.py.
	 */
	params[5] = To;
	params[6] = Tc;
	params[10] = atotAU;
	params[12] = K;
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
	//First push parameters array into clearly named variables
	parsAryToVariables();
	//--------------------------------------------------------------------------
	//------------------------------- 	START REAL CALCULATE STEPS -------------
	//--------------------------------------------------------------------------
	//calc those that are static for all epochs
	atot=0;
	atotAU=0;
	if (m1!=0)
		atot =pow(((P*P*secPerYear*secPerYear*Grav*KGperMsun*(m1+m2))/(4.0*pi*pi)),(1.0/3.0));
	atotAU = atot/MperAU;
	if ((dataRealAry[6]<1e6)&&(K==0)){
		//NOTE: both forms were tested and showed to be the same to the 1e-11 level.
		//      Thus, the same to the limit of rounding errors.
		//      Using semi-major version as it is simpler
		//Semi-major axis version
		K = ((2.0*pi*(atot/(1.0+(m1/m2)))*sin(inc*(pi/180.0)))/(P*secPerYear*pow((1.0-ecc*ecc),(1.0/2.0))));
		//Masses version
		//K = pow((2.0*pi*Grav)/(P*secPerYear),(1.0/3.0))*((m2*KGperMsun)/pow((m1+m2)*KGperMsun,(2.0/3.0)))*(sin(inc*(pi/180.0))/pow((1.0-ecc*ecc),(1.0/2.0)));
	}
	//Get the model version of each omega and shift into [0,360]
	//This is required as the RV values are for the primary and thus omega+pi
	//compared to the value used for the secondary in the DI data.
	omegaDI = omega+omegaOffsetDI;
	omegaRV = omega+omegaOffsetRV;
	if (false){
		if (omegaDI>360.0)
			omegaDI-=360.0;
		if (omegaDI<0.0)
			omegaDI+=360;
		if (omegaRV>360.0)
			omegaRV-=360.0;
		if (omegaRV<0)
			omegaRV+=360.0;
	}
	//Calculate Tc <-> T if needed.  Note, only relevant to RV data.
	if (To!=Tc){
		//if T=Tc already, do nothing.
		if ((ecc==0)||(omegaRV==90.0)||(omegaRV==270.0)){
			//Circular, so just set equal.
			if (Tc==0)
				Tc=To;
			else
				To=Tc;
		}
		else{
			ta = (pi/2.0)-omegaRV*(pi/180.0);
			halfE = atan2(sqrt(1.0-ecc)*sin(ta/2.0),sqrt(1.0+ecc)*cos(ta/2.0));
			mTTc = 2.0*halfE-ecc*sin(2.0*halfE);
			deltaT = (mTTc*P*daysPerYear)/(2.0*pi);
			if (Tc==0)
				Tc = To+deltaT;
			else
				To = Tc-deltaT;
		}
	}
	//start loop over each epoch of data
	for (i=0;i<dataModelAry_nx; i++){
		//Calculate true anomaly for RV and eccentric anomaly for DI
		anomalyCalc(dataRealAry[0+i*dataRealAry_ny]);
		//--------------------------
		//Calculate RV
		//--------------------------
		if (dataRealAry[6+i*dataRealAry_ny]<1e6)
			dataModelAry[2+i*dataModelAry_ny]=K*(cos(thetaRV+omegaRV*(pi/180.0))+ecc*cos(omegaRV*(pi/180.0)))+params[13+int(dataRealAry[7+i*dataRealAry_ny])];
		else
			dataModelAry[2+i*dataModelAry_ny]=0.0;
		//--------------------------
		//Calculate x,y or PA,SA
		//--------------------------
		if ((dataRealAry[2+i*dataRealAry_ny]<1e6)&&(dataRealAry[4+i*dataRealAry_ny]<1e6)){
			// calculate all the Thiele-Innes constants in ["]
			A = (atotAU*(parallax/1000.0))*(cos(Omega*(pi/180.0))*cos(omegaDI*(pi/180.0))-sin(Omega*(pi/180.0))*sin(omegaDI*(pi/180.0))*cos(inc*(pi/180.0)));
			B = (atotAU*(parallax/1000.0))*(sin(Omega*(pi/180.0))*cos(omegaDI*(pi/180.0))+cos(Omega*(pi/180.0))*sin(omegaDI*(pi/180.0))*cos(inc*(pi/180.0)));
			F = (atotAU*(parallax/1000.0))*(-cos(Omega*(pi/180.0))*sin(omegaDI*(pi/180.0))-sin(Omega*(pi/180.0))*cos(omegaDI*(pi/180.0))*cos(inc*(pi/180.0)));
			G = (atotAU*(parallax/1000.0))*(-sin(Omega*(pi/180.0))*sin(omegaDI*(pi/180.0))+cos(Omega*(pi/180.0))*cos(omegaDI*(pi/180.0))*cos(inc*(pi/180.0)));
			// The coordinates of the unit orbital ellipse in the true plane (Binnendijk)
			X = cos(EDI)-ecc;
			Y = sqrt(1.0-ecc*ecc)*sin(EDI);
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
		else
			dataModelAry[0+i*dataModelAry_ny] = dataModelAry[1+i*dataModelAry_ny] = 0.0;
	}
	//push parameters back into array
	variablesToParsAry();
}
//EOF
