#include "PSI1D_Sputtering.h"


#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <sstream>

using namespace std;


// Constructor
PSI1D_Sputtering::PSI1D_Sputtering(
    double an1_in,
    double am1_in,
    double an2_in,
    double am2_in,
    double es_in,
    double density_solid_in )
{


    an1 = an1_in;
    am1 = am1_in;
    an2 = an2_in;
    am2 = am2_in;
    es = es_in;
    n = density_solid_in;

    init();
}

PSI1D_Sputtering::~PSI1D_Sputtering()
{

}





// Calculates the PSI1D for a given Collisions object
void PSI1D_Sputtering::init()
{


    //ionflag -> flag for light/heavy ion.
    //ionflag = 0 => light ion sputtering.
    //ionflag = 1 => heavy ion sputtering.
    if (am1<=4.0)
    	ionflag = 0;
    else
    	ionflag = 1;

    Mratio = am2 / am1;
    Q = 1.633 * pow( an1,(2.0/3.0) ) * pow( an2,(2.0/3.0) ) * pow( ( pow(an1,(2.0/3.0)) +  pow(an2,(2.0/3.0)) ),(1.0/3.0) )
    	* pow( am1,(5.0/6.0) ) * pow( am2,(1.0/6.0) ) / (am1+ am2) * (0.15 + 0.05 * Mratio) / (1 + 0.05 * pow( Mratio,(1.6)))
    	/ pow(es,(2.0/3.0));

    eth = ( 7.0 * pow( Mratio,(-0.54) ) + 0.15 * pow( Mratio,(1.12) ) ) * es;
    mu = 4.0 * am1 * am2 / pow( (am1+am2),2.0 );
    eth1 = es / ( mu * (1-mu) ) ;
    etf = 30.74 * (am1+am2)/am2  * an1 * an2 * pow( ( pow(an1,(2.0/3.0)) +  pow(an2,(2.0/3.0)) ),(1.0/2.0) ) ;
    aL = 0.4685 * pow( ( pow(an1,(2.0/3.0)) +  pow(an2,(2.0/3.0)) ),(-1.0/2.0) );

}




// from Shuyu Dai' c code
double PSI1D_Sputtering::phy_sput_yield(double ke, double theta)
{
    double	 reducedE, stopcs, yldphy;
    double	 f,q,anu,ang_opt,eta,fs,psi,Sigma,t,angcntrb;

	reducedE = ke / etf;
	stopcs = ( 0.5 * log( 1.0+1.2288 * reducedE ) ) / ( reducedE + 0.1728 * sqrt(reducedE) + 0.008 * pow( reducedE, (0.1504) ) );


	if (ionflag==0)
	{
		f = ( 0.94 - 1.33e-3 * Mratio ) * sqrt(es);
		q = sqrt(es/mu/ke);
		anu = aL * pow( n,(1.0/3.0) ) * sqrt(1.0/2.0/reducedE/q);
		ang_opt = 90.0 - 57.29*anu;		//57.29 = 180 / pi
	}
	else
	//For a really good quantitative study you have to input the value of fs from the above reference. I am
	//circumventing the need of making a big database by using fig-3 from the following reference:
	//Institute of Plasma Physics, Nagoya University report - IPPJ-AM-26, Japan.
	{
		if (Mratio<=2.0) fs = 1.75;
		if (2.0<=Mratio && Mratio<4.0) fs = 1.5;
		if (4.0<=Mratio && Mratio<6.0) fs = 1.25;
		if (6.0<=Mratio && Mratio<8.0) fs = 1.0;
		if (Mratio>=8.0) fs = 0.8;
	//The if .. then ... else below is introduced to avoid the blowing up of f when eta = 0.0d0 (Xavier Bonnin suggested
	//this during inclusion of this routine in SOLPS-B2.5).
		eta = 1.0 - sqrt(eth/ke);
		if (eta>0.0)
			f = min(10.0,fs*(1.0+2.5*(1-eta)/eta));
		else
			f = fs;

		psi = pow( ( aL * n),(3.0/2.0) ) * sqrt(an1 * an2 * pow( ( pow(an1,(2.0/3.0)) +  pow(an2,(2.0/3.0)) ),(-1.0/2.0) )/ke);
		ang_opt = 90.0 - 286.0 * pow(psi, 0.45);
	}

	if(ang_opt<0.0)
		ang_opt = 0.0;

	Sigma = cos(ang_opt*1.7453293e-2);	//1.7453293e-2 = pi/180
	t = 1.0/cos(theta*1.7453293e-2);
	angcntrb = exp( f*( Sigma*(1-t)+log(t) ) );


	yldphy = Q * stopcs * ( 1.0 - pow( (eth/ke), (2.0/3.0) ) ) * pow( ( 1.0 - (eth/ke) ), 2.0) * angcntrb;
	if (yldphy<0.0) yldphy = 0.0;

	return yldphy;
}
