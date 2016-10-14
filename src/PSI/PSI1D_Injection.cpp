#include "PSI1D_Injection.h"
#include "SmileiMPI_Cart1D.h"
#include "Field1D.h"
#include "H5.h"


#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <sstream>

using namespace std;


// Constructor
PSI1D_Injection::PSI1D_Injection(
    PicParams& params,
    SmileiMPI* smpi,
    string psi_emitKind,
    unsigned int psi_species1,
    string psiPosition,
    unsigned int nPartEmit,
    double emitTemperature,
    double emitJ,
    double weight_const,
    double emitOffset,
    double a_FN,
    double b_FN,
    double work_function,
    string psi_relSpecies
):
PSI1D (params, smpi),
nPartEmit (nPartEmit),
emitJ (emitJ),
weight_const (weight_const),
emitOffset (emitOffset),
a_FN (a_FN),
b_FN (b_FN),
work_function (work_function)

{
    emitKind = psi_emitKind;
    species1 = psi_species1;
    relSpecies = psi_relSpecies;
    psiPos =psiPosition,
    emitTemp = emitTemperature,
    dt_ov_dx = params.timestep / params.cell_length[0];

    ySqrt_factor = pow(params.const_e, 3.0) / (4.0 * params.const_pi * params.const_ephi0 * work_function *work_function);
    a_factor = a_FN * params.const_e * params.const_e / work_function;
    b_factor = -b_FN * pow(work_function, 1.5) / params.const_e;

    if(weight_const == 0.0) {
        weight_const = nominalDensity * pow(params.cell_length[0], 3) / nomPtclsPerCell;
    }
}

PSI1D_Injection::~PSI1D_Injection()
{

}



// Calculates the PSI for a given PSI object
void PSI1D_Injection::performPSI(PicParams& params, SmileiMPI* smpi, vector<Species*>& vecSpecies, int itime, ElectroMagn* fields)
{
    Species   *s1;
    Particles *p1;

    s1 = vecSpecies[species1];
    p1 = &(s1->particles);

    SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);

    if(emitKind == "fieldEmit"){
        // field emission is calculated using Fowler-Nordheim formulae
        // from "modelling vacuum arcs: from plasma initiation to surface interactions"
        Field1D* Ex1D = static_cast<Field1D*>(fields->Ex_);
        emitField = (*Ex1D)(0);
        //emitField *= params.norm_efield;
        emitJ = (a_factor * emitField * emitField / t_y2(emitField)) /
                exp(b_factor * v_y(emitField) / emitField);
        //emitJ /= params.norm_j;
        nPartEmit = emitJ * dt_ov_dx * weight_const;
    }
    else if(emitKind == "relEmit"){
        PSI1D_Injection *relPsi_injection = static_cast<PSI1D_Injection*>(relPsi);
        nPartEmit = relPsi_injection->nPartEmit * relEmit_factor;
    }

    // PSIs usually create new particles, insert new particles to the end of particles, no matter the boundary is left or right
    // not affect the indexes_of_particles_to_exchange before exchanging particles using MPI
    if( smpi1D->isWestern() || smpi1D->isEastern() ) {
        emit(params, vecSpecies);
        unsigned int iPart = s1->getNbrOfParticles();
        emit_particles.cp_particles(nPartEmit, *p1, iPart);
        emit_particles.clear();
        unsigned int ibin = s1->bmin.size();
        s1->bmax[ibin] += nPartEmit;
    }
}


void PSI1D_Injection::emit(PicParams& params, vector<Species*>& vecSpecies){
    Species   *s1;
    s1 = vecSpecies[species1];


    emit_particles.initialize(nPartEmit, params);
    if(psiPos == "left"){
         for(int iPart=0; iPart<nPartEmit; iPart++)
         {
            emit_particles.position(0,iPart)=(((double)rand() / RAND_MAX))*params.cell_length[0]*emitOffset;
            emit_particles.position_old(0,iPart) = emit_particles.position(0,iPart);

            // initialize using the Maxwell distribution function in x-direction
            double psm = sqrt(2.0 * emitTemp / s1->species_param.mass) * sqrt(-log((double)rand() / RAND_MAX));
            double theta = M_PI*(double)rand() / RAND_MAX;
            double phi   = 2.0 * M_PI*(double)rand() / RAND_MAX;
            emit_particles.momentum(0,iPart) = abs( psm*sin(theta)*cos(phi) );
            emit_particles.momentum(1,iPart) = 0.0;
            emit_particles.momentum(2,iPart) = 0.0;

            emit_particles.weight(iPart) = weight_const;
            emit_particles.charge(iPart) = s1->species_param.charge_profile.profile;
        }
    }
    else if(psiPos == "right"){
        for(int iPart=0; iPart<nPartEmit; iPart++)
        {
           emit_particles.position(0,iPart)=params.cell_length[0]*params.n_space_global[0] - (((double)rand() / RAND_MAX))*params.cell_length[0]*emitOffset;
           emit_particles.position_old(0,iPart) = emit_particles.position(0,iPart);

           // initialize using the Maxwell distribution function in x-direction
           double psm = sqrt(2.0 * emitTemp / s1->species_param.mass) * sqrt(-log((double)rand() / RAND_MAX));
           double theta = M_PI*(double)rand() / RAND_MAX;
           double phi   = 2.0 * M_PI*(double)rand() / RAND_MAX;
           emit_particles.momentum(0,iPart) = -abs( psm*sin(theta)*cos(phi) );
           emit_particles.momentum(1,iPart) = 0.0;
           emit_particles.momentum(2,iPart) = 0.0;

           emit_particles.weight(iPart) = weight_const;
           emit_particles.charge(iPart) = s1->species_param.charge_profile.profile;
       }
    }
    else {
        ERROR("no such emitPos: " << psiPos);
    }

}
