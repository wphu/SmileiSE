/*
PSI class
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
    PartSource(PicParams&, SmileiMPI*){};
    virtual ~PartSource(){};


    //! Identification number of the PSI object
    int n_PartSource;

    // emit kind, regular or fieldEmit for injection PSI
    string emitKind;

    // PSI position : only left and right for 1D case
    string emitPos;

    // relevant PSI, emitting number of two species may be relevant
    // such as nPartEmit(A) = relCoff * nPartEmit(B)
    PartSource *relPartSource;
    string relSpecies;

    // position offset of injected or sputtered particles
    double posOffset;
    // the energy/temperature of the new particles
    double emitTemp;
    double weight_const;

    void setRelPartSource(PartSource* relevantPartSource)
    {
        relPartSource = relevantPartSource;
    }

    //! Group of the species numbers that are associated for PSI.
    //> actually, each species gourp only contains one species for PSI
    //> for PSI_Injection, only species1 is used;
    //> for sputtering and secondary electron emission, species1 is the incident particle.
    unsigned int species1, species2;

    //! Method called in the main smilei loop to apply PSI at each timestep
    virtual void emitLoad(PicParams&, SmileiMPI* smpi, std::vector<Species*>&,int, ElectroMagn* ){};

    Particles new_particles;
    vector<int> count_of_particles_to_insert_s1;
    double const_e;


private:

};

#endif
