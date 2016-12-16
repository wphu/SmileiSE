/*
Collisions2D_Elastic class
*/

#ifndef COLLISIONS2D_ELASTIC_H
#define COLLISIONS2D_ELASTIC_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "Collisions2D.h"
#include "H5.h"

class Collisions2D_Elastic : public Collisions2D
{

public:
    //! Constructor for Collisions2D between two species
    Collisions2D_Elastic(PicParams&,std::vector<Species*>&,SmileiMPI*,unsigned int,std::vector<unsigned int>,std::vector<unsigned int>,double,bool,int);
    ~Collisions2D_Elastic();

    double cross_section(double ke);

    //! Method called in the main smilei loop to apply collisions at each timestep
    void collide(PicParams&, SmileiMPI* smpi, std::vector<Species*>&,int);

private:
    //>the ionization threshold energy
    double energy_ion;

};


#endif
