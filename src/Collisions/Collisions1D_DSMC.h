/*
Collisions1D class
*/

#ifndef COLLISIONS1D_DSMC_H
#define COLLISIONS1D_DSMC_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "Collisions1D.h"
#include "H5.h"

class Collisions1D_DSMC : public Collisions1D
{

public:
    //! Constructor for Collisions1D between two species
    Collisions1D_DSMC(PicParams&,std::vector<Species*>&,SmileiMPI*,unsigned int,std::vector<unsigned int>,std::vector<unsigned int>,double,bool,int);
    ~Collisions1D_DSMC();

    double cross_section(double ke);

    //! Method called in the main smilei loop to apply collisions at each timestep
    void collide(PicParams&,std::vector<Species*>&,int);


private:
    inline double scatter_particles(Particles* particle1, int iPart1, Particles* particle2, int iPart2);

};


#endif
