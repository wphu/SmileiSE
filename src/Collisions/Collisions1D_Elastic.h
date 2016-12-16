/*
Collisions1D class
*/

#ifndef COLLISIONS1D_ELASTIC_H
#define COLLISIONS1D_ELASTIC_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "Collisions1D.h"
#include "H5.h"

class Collisions1D_Elastic : public Collisions1D
{

public:
    //! Constructor for Collisions1D between two species
    Collisions1D_Elastic(PicParams&,std::vector<Species*>&,SmileiMPI*,unsigned int,std::vector<unsigned int>,std::vector<unsigned int>,double,bool,int);
    ~Collisions1D_Elastic();

    double cross_section(double ke);

    //! Method called in the main smilei loop to apply collisions at each timestep
    void collide(PicParams&, SmileiMPI* smpi, std::vector<Species*>&,int);

private:
    //>the ionization threshold energy
    double energy_ion;

};


#endif
