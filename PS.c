#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#ifndef max
#define max(x, y)       (((x) > (y)) ? (x) : (y))
#endif

#ifndef min
#define min(x, y)       (((x) < (y)) ? (x) : (y))
#endif

#define points		101
#define DENSITY_C 0.11286  //*number density---unit-- C/A**3 ---*

//theta -> angle (in degrees) with normal to target of incident particle.
//eo	-> energy of the incident particle (eV)
//Z1	-> nuclear charge of incident atom.
//Z2	-> nuclear charge of target atom.
//M1	-> atomic mass of incident atom (amu).  M2	-> atomic mass of target atom (amu).
//es	-> surface binding energy (heat of sublimation) of target (eV).
//tgdns -> target density (gms/cc)
//yldphy -> The physical sputtering yield.


#define	a_angle 0.0
#define	b_angle 90.0
#define	a_energy 000.0
#define	b_energy 1000.0


int		ionflag;
float	Z1,Z2,M1,M2,es,n ;
float	Q,eth,eth1,mu,etf,aL,Mratio ;
float	phy_sput_yield(float , float );


void main()
{

//void physput(theta,eo,z1,z2,am1,am2,yldphy)
int		i,j;	  
float	theta,eo,yldphy,	*energy, *angle;
char	substrate;
FILE	*diff_angle, *diff_energy, *test,*input,*output;

int nangle,ncolumn,ndata;
float *data,yield_avg;



//read the data=========================================
nangle=90;
ncolumn=8;
ndata=nangle*ncolumn;
i=0;

data=(float *)calloc(nangle*ncolumn,sizeof(float));
input = fopen("pangledist10_b5.txt","r");


//Energy
eo=100.1643753;
//Angle
theta= 45.0;



while(!feof(input))
{
	fscanf(input,"%f",&data[i++]);

}

for(j=0;j<nangle;j++)
{
	for(i=0;i<ncolumn;i++)
	{
		printf("%15.8f",data[j*ncolumn+i]);
	}
	printf("\n");
	printf("\n");
}

//end:read the data=========================================


/***************************************************************/

energy = (float *)calloc(points,sizeof(float));
angle  = (float *)calloc(points,sizeof(float));

diff_angle = fopen("diff_angle.dat","w");
diff_energy = fopen("diff_energy.dat","w");
test = fopen("test.dat","w");

/***************************************************************/

substrate = 'c';

//incident particle



Z1=1.0;		//D
M1=2.0;



/***************************************************************/


if(substrate== 'C' || substrate== 'c')
{
	n = DENSITY_C ;
	es = 7.428;
	Z2 = 6.0;
	M2 = 12.0;
}


/***************************************************************/
 
if (theta<0.0 || theta>=90.0) 
	printf("Input error; Correct so that 0 <= theta < 90");

Mratio = M2/M1;

Q = 1.633 * pow( Z1,(2.0/3.0) ) * pow( Z2,(2.0/3.0) ) * pow( ( pow(Z1,(2.0/3.0)) +  pow(Z2,(2.0/3.0)) ),(1.0/3.0) )
	* pow( M1,(5.0/6.0) ) * pow( M2,(1.0/6.0) ) / (M1+ M2) * (0.15 + 0.05 * Mratio) / (1 + 0.05 * pow( Mratio,(1.6))) 
	/ pow(es,(2.0/3.0));

eth = ( 7.0 * pow( Mratio,(-0.54) ) + 0.15 * pow( Mratio,(1.12) ) ) * es;
mu = 4.0 * M1 * M2 / pow( (M1+M2),2.0 );
eth1 = es / ( mu * (1-mu) ) ;
etf = 30.74 * (M1+M2)/M2  * Z1 * Z2 * pow( ( pow(Z1,(2.0/3.0)) +  pow(Z2,(2.0/3.0)) ),(1.0/2.0) ) ;


fprintf(test,"Mratio=%e\netf=%e\neth=%e\neth1=%e\nQ=%e\n",Mratio,etf,eth,eth1,Q);



//ionflag -> flag for light/heavy ion.
//ionflag = 0 => light ion sputtering.
//ionflag = 1 => heavy ion sputtering.
if (M1<=4.0) 
	ionflag = 0;
else
	ionflag = 1;

aL = 0.4685 * pow( ( pow(Z1,(2.0/3.0)) +  pow(Z2,(2.0/3.0)) ),(-1.0/2.0) );

/***************************************************************/

yldphy = phy_sput_yield( eo, theta );

fprintf(test,"Physical sputtering yield=%e\n\n\n",yldphy);


//calculate the average phy sput yield==================================================
//======================================================================================
output = fopen("output.txt","w");

yield_avg=0.0;
for(j=0;j<nangle;j++)
{
	yield_avg+=phy_sput_yield(eo,data[j*ncolumn]-0.5)*data[j*ncolumn+1];
	//fprintf(output,"%f\n",yield_avg);
}

fprintf(output,"%f",yield_avg);




for(i=0; i<points; i++)
{
	energy[i] = a_energy + (b_energy - a_energy)/(points-1) * i ;
	
	yldphy = phy_sput_yield( energy[i], theta );

	fprintf(diff_energy,"%e %e\n",energy[i],yldphy);

/***************************************************************/

	angle[i] = a_angle + (b_angle - a_angle)/(points-1) * i ;
	
	yldphy = phy_sput_yield( eo, angle[i] );

	fprintf(diff_angle,"%e %e\n",angle[i],yldphy);

/***************************************************************/

}

}





float phy_sput_yield(float eo, float theta )
{

float	 reducedE, stopcs, yldphy;
float	 f,q,anu,ang_opt,eta,fs,psi,Sigma,t,angcntrb;

	reducedE = eo / etf;
	stopcs = ( 0.5 * log( 1.0+1.2288 * reducedE ) ) / ( reducedE + 0.1728 * sqrt(reducedE) + 0.008 * pow( reducedE, (0.1504) ) );


	if (ionflag==0) 
	{
		f = ( 0.94 - 1.33e-3 * Mratio ) * sqrt(es);	
		q = sqrt(es/mu/eo);
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
		eta = 1.0 - sqrt(eth/eo);
		if (eta>0.0) 
			f = min(10.0,fs*(1.0+2.5*(1-eta)/eta));
		else
			f = fs;
	
		psi = pow( ( aL * n),(3.0/2.0) ) * sqrt(Z1 * Z2 * pow( ( pow(Z1,(2.0/3.0)) +  pow(Z2,(2.0/3.0)) ),(-1.0/2.0) )/eo);
		ang_opt = 90.0 - 286.0 * pow(psi, 0.45);	
	}

	if(ang_opt<0.0)
		ang_opt = 0.0;

	Sigma = cos(ang_opt*1.7453293e-2);	//1.7453293e-2 = pi/180
	t = 1.0/cos(theta*1.7453293e-2);
	angcntrb = exp( f*( Sigma*(1-t)+log(t) ) );


	yldphy = Q * stopcs * ( 1.0 - pow( (eth/eo), (2.0/3.0) ) ) * pow( ( 1.0 - (eth/eo) ), 2.0) * angcntrb;
	if (yldphy<0.0) yldphy = 0.0;
	
	return(yldphy);

}