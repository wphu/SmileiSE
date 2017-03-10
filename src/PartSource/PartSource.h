/*
PartSource class
*/

#ifndef PARTSOURCE_H
#define PARTSOURCE_H

#include <vector>
#include <string>

#include "PicParams.h"
#include "SmileiMPI.h"
#include "Species.h"

using namespace std;
class PartSource
{

public:
    //! Constructor for PSI between two species
    PartSource(PicParams& params, SmileiMPI* smpi)
    {
            count_of_particles_to_insert_s1.resize(params.n_space[0]);
            count_of_particles_to_insert.resize(params.n_space[0]);
            numPart_in_each_bin.resize(params.n_space[0]);

            new_particles.initialize(0, params);
            mean_velocity.resize(3, 0.0);
    };
    virtual ~PartSource(){};

    unsigned int species1;
    Particles new_particles;
    vector<int> count_of_particles_to_insert_s1;

    vector<int> numPart_in_each_bin;
    vector<int> count_of_particles_to_insert;
    vector<int> indexes_of_particles_to_erase;

    // Particle number in one cell, different from particle number in one bin for 2D case
    int numPart_in_each_cell;

    //! Method called in the main smilei loop to apply PartSource at each timestep
    virtual void emitLoad(PicParams&, SmileiMPI* smpi, std::vector<Species*>&,int, ElectroMagn* ){};

    vector<double> mean_velocity;

    string PartSource_type;

    // emit kind, regular or fieldEmit for injection PSI
    string emitKind;
    // load Kind: "dn" or "nT"
    string loadKind;

    // relevant PSI, emitting number of two species may be relevant
    // such as nPartEmit(A) = relCoff * nPartEmit(B)
    PartSource *relPartSource;
    string relSpecies;
    virtual void setRelPartSource(PartSource* relevantPartSource)
    {
        relPartSource = relevantPartSource;
    }

    int n_PartSource;

    double weight_const;

    int everyTime;

private:

    double const_e;

};

#endif
