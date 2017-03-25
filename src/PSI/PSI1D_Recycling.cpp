#include "PSI1D_Recycling.h"
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
PSI1D_Recycling::PSI1D_Recycling(
    PicParams& params,
    SmileiMPI* smpi,
    unsigned int psi_species1,
    unsigned int psi_species2,
    string psiPosition,
    double emitTemperature,
    double recycling_factor_temp
):
PSI1D(params, smpi)
{
    species1            = psi_species1;
    species2            = psi_species2;
    psiPos              = psiPosition;
    emitTemp            = emitTemperature;
    recycling_factor    = recycling_factor_temp;

}

PSI1D_Recycling::~PSI1D_Recycling()
{

}



// Calculates the PSI1D for a given Collisions object
void PSI1D_Recycling::performPSI(PicParams& params, SmileiMPI* smpi, vector<Species*>& vecSpecies, int itime, ElectroMagn* fields)
{
    // the angle of particle velocity with the surface normal
    double theta;
    // kinetic energy_ion
    double ke;
    double v_square, v_magnitude;
    // sputtering probability
    double pSput;
    int iDim;
    double nPartEmit_double;
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
        if( psiPos == "left" && p1->position(iDim,iPart) < 0.0 )
        {
            nPartEmit++;
        }
        else if( psiPos == "right" && p1->position(iDim,iPart) > params.sim_length[0] )
        {
            nPartEmit++;
        }
    };
    nPartEmit_double = nPartEmit * recycling_factor;
    nPartEmit = nPartEmit_double;
    nPartEmit_rem += ( nPartEmit_double - nPartEmit );
    if(nPartEmit_rem >= 1.0)
    {
        nPartEmit++;
        nPartEmit_rem = nPartEmit_rem - 1.0;
    }

    if( smpi->isWestern() || smpi->isEastern() ) {
        emit(params, vecSpecies, species2);
        s2->insert_particles_to_bins(new_particles, count_of_particles_to_insert_s2);
        new_particles.clear();
    };

}




void PSI1D_Recycling::emit(PicParams& params, vector<Species*>& vecSpecies, unsigned int species_emit)
{
    Species   *s1;
    s1 = vecSpecies[species_emit];


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
            new_particles.momentum(0,iPart) = abs( psm*sin(theta) );
            new_particles.momentum(1,iPart) = psm*cos(theta)*sin(phi);
            new_particles.momentum(2,iPart) = psm*cos(theta)*cos(phi);

            new_particles.weight(iPart) = s1->species_param.weight;
            new_particles.charge(iPart) = s1->species_param.charge;
        }
    }
    else if(psiPos == "right"){
        count_of_particles_to_insert_s2.back() = nPartEmit;
        for(int iPart=0; iPart<nPartEmit; iPart++)
        {
           new_particles.position(0,iPart)=params.sim_length[0] - (((double)rand() / RAND_MAX))*params.cell_length[0]*posOffset;
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
           new_particles.momentum(0,iPart) = -abs( psm*sin(theta) );
           new_particles.momentum(1,iPart) = psm*cos(theta)*sin(phi);
           new_particles.momentum(2,iPart) = psm*cos(theta)*cos(phi);

           new_particles.weight(iPart) = s1->species_param.weight;
           new_particles.charge(iPart) = s1->species_param.charge;
       }
    }
    else {
        ERROR("no such psiPos: " << psiPos);
    }
}
