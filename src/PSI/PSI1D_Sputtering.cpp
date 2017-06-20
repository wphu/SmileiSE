#include "PSI1D_Sputtering.h"
#include "SmileiMPI.h"
#include "Field2D.h"
#include "H5.h"


#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <sstream>

using namespace std;


// Constructor
PSI1D_Sputtering::PSI1D_Sputtering(
    PicParams& params,
    SmileiMPI* smpi,
    vector<Species*>& vecSpecies,
    unsigned int psi_species1,
    unsigned int psi_species2,
    string psiPosition,
    double emitTemperature
):
PSI1D(params, smpi)
{
    species1 = psi_species1;
    species2 = psi_species2;
    psiPos = psiPosition;
    emitTemp = emitTemperature;

    const_e = params.const_e;

    init(vecSpecies);
}

PSI1D_Sputtering::~PSI1D_Sputtering()
{

}



// Calculates the PSI1D for a given Collisions object
void PSI1D_Sputtering::performPSI(PicParams& params, SmileiMPI* smpi, vector<Species*>& vecSpecies, int itime, ElectroMagn* fields)
{
    // the angle of particle velocity with the surface normal
    double theta;
    // kinetic energy_ion
    double ke;
    double v_square, v_magnitude;
    // sputtering probability
    double pSput;
    int iDim;
    Species   *s1, *s2;
    Particles *p1, *p2;


    s1 = vecSpecies[species1];
    s2 = vecSpecies[species2];
    p1 = &(s1->psi_particles);
    p2 = &(s2->psi_particles);


    iDim = 0;
    nPartEmit = 0;
    int nPart = p1->size();
    for(unsigned int iPart = 0; iPart < nPart; iPart++)
    {
        if( p1->position(iDim,iPart) < smpi->getDomainLocalMin(iDim) || p1->position(iDim,iPart) > smpi->getDomainLocalMax(iDim) ) {
            v_square = pow(p1->momentum(0,iPart),2) + pow(p1->momentum(1,iPart),2) + pow(p1->momentum(2,iPart),2);
            theta = abs( p1->momentum(0,iPart) ) / sqrt( v_square );
            theta *= ( 180.0 / params.const_pi );
            ke = 0.5 * s1->species_param.mass * v_square;
            //ke *= params.norm_temperature;
            pSput = phy_sput_yield( theta, ke/const_e );
            double ran_p = (double)rand() / RAND_MAX;
            if( pSput > ran_p ) {
                nPartEmit++;
            }
        }
    };

    if( smpi->isWestern() || smpi->isEastern() )
    {
        emit(params, vecSpecies);
        s2->insert_particles_to_bins(new_particles, count_of_particles_to_insert_s2);
        new_particles.clear();
    }

}




void PSI1D_Sputtering::emit(PicParams& params, vector<Species*>& vecSpecies)
{
    Species   *s1;
    // Here species2 is sputtered
    s1 = vecSpecies[species2];

    new_particles.initialize(nPartEmit, params);
    if(psiPos == "left"){
        count_of_particles_to_insert_s2.front() = nPartEmit;
        for(int iPart=0; iPart<nPartEmit; iPart++)
        {
            new_particles.position(0,iPart)=(((double)rand() / RAND_MAX))*params.cell_length[0]*posOffset;
            new_particles.position_old(0,iPart) = new_particles.position(0,iPart);

            double ran;
            do {
                ran = (double)rand() / RAND_MAX;
            }
            while (ran == 0.0);
            // initialize using the Maxwell distribution function in x-direction
            double psm = sqrt(2.0 * const_e * emitTemp / s1->species_param.mass) * sqrt(-log(ran));
            double theta = M_PI*(double)rand() / RAND_MAX;
            double phi   = 2.0 * M_PI*(double)rand() / RAND_MAX;
            new_particles.momentum(0,iPart) = abs( psm*sin(theta)*cos(phi) );
            new_particles.momentum(1,iPart) = 0.0;
            new_particles.momentum(2,iPart) = 0.0;

            new_particles.al_imp(0,iPart) = 0.0;
            new_particles.al_imp(1,iPart) = 0.0;
            new_particles.al_imp(2,iPart) = 0.0;
            new_particles.au_imp(0,iPart) = 0.0;
            new_particles.au_imp(1,iPart) = 0.0;
            new_particles.au_imp(2,iPart) = 0.0;

            new_particles.weight(iPart) = s1->species_param.weight;
            new_particles.charge(iPart) = s1->species_param.charge;
        }
    }
    else if(psiPos == "right"){
        count_of_particles_to_insert_s2.back() = nPartEmit;
        for(int iPart=0; iPart<nPartEmit; iPart++)
        {
           new_particles.position(0,iPart)=params.cell_length[0]*params.n_space_global[0] - (((double)rand() / RAND_MAX))*params.cell_length[0]*posOffset;
           new_particles.position_old(0,iPart) = new_particles.position(0,iPart);

           double ran;
           do {
               ran = (double)rand() / RAND_MAX;
           }
           while (ran == 0.0);
           // initialize using the Maxwell distribution function in x-direction
           double psm = sqrt(2.0 * const_e * emitTemp / s1->species_param.mass) * sqrt(-log(ran));
           double theta = M_PI*(double)rand() / RAND_MAX;
           double phi   = 2.0 * M_PI*(double)rand() / RAND_MAX;
           new_particles.momentum(0,iPart) = -abs( psm*sin(theta)*cos(phi) );
           new_particles.momentum(1,iPart) = 0.0;
           new_particles.momentum(2,iPart) = 0.0;

           new_particles.al_imp(0,iPart) = 0.0;
           new_particles.al_imp(1,iPart) = 0.0;
           new_particles.al_imp(2,iPart) = 0.0;
           new_particles.au_imp(0,iPart) = 0.0;
           new_particles.au_imp(1,iPart) = 0.0;
           new_particles.au_imp(2,iPart) = 0.0;

           new_particles.weight(iPart) = s1->species_param.weight;
           new_particles.charge(iPart) = s1->species_param.charge;
       }
    }
    else {
        ERROR("no such psiPos: " << psiPos);
    }
}



// Calculates the PSI1D for a given Collisions object
void PSI1D_Sputtering::init(vector<Species*>& vecSpecies)
{
    Species   *s1, *s2;
    Particles *p1, *p2;


    s1 = vecSpecies[species1];
    s2 = vecSpecies[species2];
    p1 = &(s1->psi_particles);
    p2 = &(s2->psi_particles);

    an1 = s1->species_param.atomic_number;
    am1 = s1->species_param.atomic_mass;
    an2 = s2->species_param.atomic_number;
    am2 = s2->species_param.atomic_mass;
    es = s2->species_param.surface_binding_energy;
    n = s2->species_param.density_solid;

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
