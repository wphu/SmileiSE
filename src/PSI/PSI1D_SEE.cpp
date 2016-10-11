#include "PSI1D_SEE.h"
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
PSI1D_SEE::PSI1D_SEE(
    PicParams& params,
    SmileiMPI* smpi,
    unsigned int psi_species1,
    unsigned int psi_species2,
    string psiPosition,
    double emitTemperature,
    double SEEYield
):
PSI1D(params, smpi),
SEEYield(SEEYield)
{
    species1 = psi_species1;
    species2 = psi_species2;
    psiPos = psiPosition;
    emitTemp = emitTemperature;
}

PSI1D_SEE::~PSI1D_SEE()
{

}



// Calculates the PSI1D
void PSI1D_SEE::performPSI(PicParams& params, SmileiMPI* smpi, vector<Species*>& vecSpecies, int itime, ElectroMagn* fields)
{
    // the angle of particle velocity with the surface normal
    double theta;
    // kinetic energy_ion
    double ke;
    double v_square, v_magnitude;
    // sputtering probability
    double pSput;
    double nPartEmit_temp;
    int iDim;
    Species   *s1, *s2;
    Particles *p1, *p2;


    s1 = vecSpecies[species1];
    s2 = vecSpecies[species2];
    p1 = &(s1->particles);
    p2 = &(s2->particles);


    iDim = 0;
    nPartEmit = 0;
    nPartEmit_temp = 0.0;
    int nPart = s1->indexes_of_particles_to_exchange_per_thd[0].size();
    for(unsigned int iPart = 0; iPart < nPart; iPart++)
    {
        if( p1->position(iDim,iPart) < smpi->getDomainLocalMin(iDim) || p1->position(iDim,iPart) > smpi->getDomainLocalMax(iDim) ) {
            nPartEmit_temp += SEEYield;
        }
    };
    nPartEmit = nPartEmit_temp;
    //SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);

    // PSIs usually create new particles, insert new particles to the end of particles, no matter the boundary is left or right
    // not affect the indexes_of_particles_to_exchange before exchanging particles using MPI
    if( smpi->isWestern() || smpi->isEastern() ) {
        emit(params, vecSpecies, species2);
        unsigned int iPart = s1->getNbrOfParticles();
        new_particles.cp_particles(nPartEmit, *p1, iPart);
        new_particles.clear();
        unsigned int ibin = s1->bmin.size();
        s1->bmax[ibin] += nPartEmit;
    };
}


void PSI1D_SEE::emit(PicParams& params, vector<Species*>& vecSpecies, unsigned int species_emit){
    Species   *s1;
    s1 = vecSpecies[species_emit];


    new_particles.initialize(nPartEmit, params);
    if(psiPos == "left"){
         for(int iPart=0; iPart<nPartEmit; iPart++)
         {
            new_particles.position(0,iPart)=(((double)rand() / RAND_MAX))*params.cell_length[0]*posOffset;
            new_particles.position_old(0,iPart) = new_particles.position(0,iPart);

            // initialize using the Maxwell distribution function in x-direction
            double psm = sqrt(2.0 * emitTemp / s1->species_param.mass) * sqrt(-log((double)rand() / RAND_MAX));
            double theta = M_PI*(double)rand() / RAND_MAX;
            double phi   = 2.0 * M_PI*(double)rand() / RAND_MAX;
            new_particles.momentum(0,iPart) = abs( psm*sin(theta)*cos(phi) );
            new_particles.momentum(1,iPart) = 0.0;
            new_particles.momentum(2,iPart) = 0.0;

            new_particles.weight(iPart) = weight_const;
            new_particles.charge(iPart) = s1->species_param.charge_profile.profile;
        }
    }
    else if(psiPos == "right"){
        for(int iPart=0; iPart<nPartEmit; iPart++)
        {
           new_particles.position(0,iPart)=params.cell_length[0]*params.n_space_global[0] - (((double)rand() / RAND_MAX))*params.cell_length[0]*posOffset;
           new_particles.position_old(0,iPart) = new_particles.position(0,iPart);

           // initialize using the Maxwell distribution function in x-direction
           double psm = sqrt(2.0 * emitTemp / s1->species_param.mass) * sqrt(-log((double)rand() / RAND_MAX));
           double theta = M_PI*(double)rand() / RAND_MAX;
           double phi   = 2.0 * M_PI*(double)rand() / RAND_MAX;
           new_particles.momentum(0,iPart) = -abs( psm*sin(theta)*cos(phi) );
           new_particles.momentum(1,iPart) = 0.0;
           new_particles.momentum(2,iPart) = 0.0;

           new_particles.weight(iPart) = weight_const;
           new_particles.charge(iPart) = s1->species_param.charge_profile.profile;
       }
    }
    else {
        ERROR("no such psiPos: " << psiPos);
    }
}
