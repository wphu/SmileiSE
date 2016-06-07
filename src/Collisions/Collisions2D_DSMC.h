/*
Collisions2D_DSMC class
*/

#ifndef COLLISIONS2D_DSMC_H
#define COLLISIONS2D_DSMC_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "Collisions2D.h"
#include "H5.h"

class Collisions2D_DSMC : public Collisions2D
{

public:
    //! Constructor for Collisions2D between two species
    Collisions2D_DSMC(PicParams&,std::vector<Species*>&,SmileiMPI*,unsigned int,std::vector<unsigned int>,std::vector<unsigned int>,double,bool,int);
    ~Collisions2D_DSMC();

    void scatter_particles(Particles* particle1, int iPart1, Particles* particle2, int iPart2);
    double evalSigma(double cr);
    inline double relative_velocity(Particles* particle1, int iPart1, Particles* particle2, int iPart2)
    {
        double rv;
        rv = sqrt( pow((particle1->momentum(0, iPart1) - particle2->momentum(0, iPart2)), 2)
                +  pow((particle1->momentum(1, iPart1) - particle2->momentum(1, iPart2)), 2)
                +  pow((particle1->momentum(2, iPart1) - particle2->momentum(2, iPart2)), 2)
                );
        return rv;
    };


    //! Method called in the main smilei loop to apply collisions at each timestep
    void collide(PicParams&,std::vector<Species*>&,int);

private:


};


#endif
