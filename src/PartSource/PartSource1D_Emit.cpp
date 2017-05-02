#include "PartSource1D_Emit.h"
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
PartSource1D_Emit::PartSource1D_Emit(
    PicParams&      params,
    SmileiMPI*      smpi,
    string          emit_emitKind,
    unsigned int    emit_species1,
    string          emitPosition,
    int             emitNumber,
    double          emitTemperature,
    double          emitJ,
    double          emitFlux,
    double          emitOffset,
    double          a_FN,
    double          b_FN,
    double          work_function,
    string          emit_relSpecies
):
PartSource1D    (params, smpi),
emitNumber      (emitNumber),
emitJ           (emitJ),
emitFlux        (emitFlux),
emitOffset      (emitOffset),
a_FN            (a_FN),
b_FN            (b_FN),
work_function   (work_function)

{
    emitKind    = emit_emitKind;
    species1    = emit_species1;
    relSpecies  = emit_relSpecies;
    emitPos     = emitPosition,
    emitTemp    = emitTemperature,
    dt_ov_dx    = params.timestep / params.cell_length[0];
    dt          = params.timestep;
    YZArea      = 1.0;

    ySqrt_factor    = pow(params.const_e, 3.0) / (4.0 * params.const_pi * params.const_ephi0 * work_function *work_function);
    a_factor        = a_FN * params.const_e * params.const_e / work_function;
    b_factor        = -b_FN * pow(work_function, 1.5) / params.const_e;

    emitStep    = 1.0 + emitNumber * params.species_param[species1].weight * params.cell_length[0] / ( emitFlux * params.timestep );
    emitRem     = emitFlux * emitStep * params.timestep / ( params.species_param[species1].weight * params.cell_length[0] ) - emitNumber;
    emitRemTot  = 0.0;
    MESSAGE("emitStep = "<<emitStep);

    count_of_particles_to_insert_s1.resize(params.n_space[0]);
    for(int i = 0; i < count_of_particles_to_insert_s1.size(); i++)
    {
        count_of_particles_to_insert_s1[i] = 0;
    }
}

PartSource1D_Emit::~PartSource1D_Emit()
{

}



// Calculates the PSI for a given PSI object
void PartSource1D_Emit::emitLoad(PicParams& params, SmileiMPI* smpi, vector<Species*>& vecSpecies, int itime, ElectroMagn* fields)
{
    Species   *s1;
    Particles *p1;

    if(itime%emitStep == 0)
    {
        double emitNumber_temp;
        emitRemTot += emitRem;
        emitNumber_temp = emitNumber;
        if(emitRemTot > 1.0)
        {
            emitNumber_temp = emitNumber + 1;
            emitRemTot -= 1.0;
        }


        s1 = vecSpecies[species1];
        p1 = &(s1->particles);

        SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);

        for(int ibin=0; ibin<count_of_particles_to_insert_s1.size(); ibin++)
        {
            count_of_particles_to_insert_s1[ibin] = 0;
        }

        if(emitKind == "fieldEmit"){
            // field emission is calculated using Fowler-Nordheim formulae
            // from "modelling vacuum arcs: from plasma initiation to surface interactions"
            Field1D* Ex1D = static_cast<Field1D*>(fields->Ex_);
            emitField = (*Ex1D)(0);
            //emitField *= params.norm_efield;
            emitJ = (a_factor * emitField * emitField / t_y2(emitField)) /
                    exp(b_factor * v_y(emitField) / emitField);
            //emitJ /= params.norm_j;
            nPartEmit = emitJ * dt_ov_dx / (weight_const*s1->species_param.charge);
        }
        else if(emitKind == "relEmit"){
            PartSource1D_Emit *relPartSource_Emit = static_cast<PartSource1D_Emit*>(relPartSource);
            nPartEmit = relPartSource_Emit->nPartEmit * relEmit_factor;
        }
        else if(emitKind == "regular"){
            nPartEmit = emitNumber_temp;
            //MESSAGE("Injected number: "<<nPartEmit<<"  "<<emitJ<<"  "<<dt_ov_dx<<" "<<weight_const<<"  "<<s1->species_param.charge);
        }

        // PSIs usually create new particles, insert new particles to the end of particles, no matter the boundary is left or right
        // not affect the indexes_of_particles_to_exchange before exchanging particles using MPI
        if( emitPos=="left" && smpi1D->isWestern() || emitPos=="right" && smpi1D->isEastern() ) {
            //MESSAGE("Befor particle number: "<<s1->getNbrOfParticles());
            emit(params, vecSpecies);
            s1->insert_particles_to_bins(new_particles, count_of_particles_to_insert_s1);
            new_particles.clear();
            //MESSAGE("Now particle number: "<<s1->getNbrOfParticles());
            //MESSAGE("Now particle capacity: "<<s1->getParticlesCapacity());

        }
    }

}

void PartSource1D_Emit::emit(PicParams& params, vector<Species*>& vecSpecies){
    Species   *s1;
    Particles *p1;
    s1 = vecSpecies[species1];

    new_particles.initialize(nPartEmit, params);
    if(emitPos == "left"){
        count_of_particles_to_insert_s1.front() = nPartEmit;
        for(int iPart=0; iPart<nPartEmit; iPart++)
        {
            new_particles.position(0,iPart)=(((double)rand() / RAND_MAX))*params.cell_length[0]*emitOffset;
            new_particles.position_old(0,iPart) = new_particles.position(0,iPart);

            double ran;
            do {
                ran = (double)rand() / RAND_MAX;
            }
            while (ran == 0.0);
            // Velocity magnitude: from Maxwell velocity distribution
            // The angle between velocity of emitted particle and the surface normal: cosine
            //      cos(alpha) = sqrt(random number 0-1)
            // The azimuthal angle is uniformly distributed on the interval [0 2pi]
            double psm = sqrt(2.0 * params.const_e * emitTemp / s1->species_param.mass) * sqrt(-log(ran));
            double cosAlpha = sqrt((double)rand() / RAND_MAX);
            double sinAlpha = sqrt(1.0 - cosAlpha * cosAlpha);
            double phi   = 2.0 * M_PI*(double)rand() / RAND_MAX;

            new_particles.momentum(0,iPart) = abs( psm * cosAlpha );
            new_particles.momentum(1,iPart) = psm * sinAlpha * cos(phi);
            new_particles.momentum(2,iPart) = psm * sinAlpha * sin(phi);

            new_particles.weight(iPart) = s1->species_param.weight;
            new_particles.charge(iPart) = s1->species_param.charge;
        }
    }
    else if(emitPos == "right"){
        count_of_particles_to_insert_s1.back() = nPartEmit;
        for(int iPart=0; iPart<nPartEmit; iPart++)
        {
           new_particles.position(0,iPart) = params.sim_length[0] - (((double)rand() / RAND_MAX))*params.cell_length[0]*emitOffset;
           new_particles.position_old(0,iPart) = new_particles.position(0,iPart);

           double ran;
           do {
               ran = (double)rand() / RAND_MAX;
           }
           while (ran == 0.0);
           // Velocity magnitude: from Maxwell velocity distribution
           // The angle between velocity of emitted particle and the surface normal: cosine
           //      cos(alpha) = sqrt(random number 0-1)
           // The azimuthal angle is uniformly distributed on the interval [0 2pi]
           double psm = sqrt(2.0 * params.const_e * emitTemp / s1->species_param.mass) * sqrt(-log(ran));
           double cosAlpha = sqrt((double)rand() / RAND_MAX);
           double sinAlpha = sqrt(1.0 - cosAlpha * cosAlpha);
           double phi   = 2.0 * M_PI*(double)rand() / RAND_MAX;

           new_particles.momentum(0,iPart) = -abs( psm * cosAlpha );
           new_particles.momentum(1,iPart) = psm * sinAlpha * cos(phi);
           new_particles.momentum(2,iPart) = psm * sinAlpha * sin(phi);

           new_particles.weight(iPart) = s1->species_param.weight;
           new_particles.charge(iPart) = s1->species_param.charge;
       }
    }
    else {
        ERROR("no such emitPos: " << emitPos);
    }

}
