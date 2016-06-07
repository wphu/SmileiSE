/*

Collisions class - Frederic Perez - 03/2015

This is based on the work described here
http://dx.doi.org/10.1063/1.4742167

Binary collisions, between macro-particles, are treated according to a scheme
first described by Nanbu (http://dx.doi.org/10.1103/PhysRevE.55.4642).

To include collisions in the simulations, add a block in the input file,
similar to the following:

# COLLISIONS
# species1    = list of strings, the names of the first species that collide
# species2    = list of strings, the names of the second species that collide
#               (can be the same as species1)
# coulomb_log = float, Coulomb logarithm. If negative or zero, then automatically computed.
Collisions(
	species1 = ["ion1"],
	species2 = ["electron1"],
	coulomb_log = 2.0
)

Several collision types can be defined. For each type, add a group "Collisions()".

*/

#ifndef COLLISIONS_H
#define COLLISIONS_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "H5.h"

class Collisions
{

public:
    //! Constructor for Collisions between two species
    Collisions(){};
    virtual ~Collisions(){};

    //! Identification number of the Collisions object
    int n_collisions;

    //! Group of the species numbers that are associated for Collisions.
    //> species_group1 (2,3) may include more than one species only for coulomb collision and DSMC,
    //> species_group1 (2,3) only include one species for other collision.
    std::vector<unsigned int> species_group1, species_group2, species_group3;

    //! True if collisions inside a group of species, False if collisions between different groups of species
    bool intra_collisions;

    //! Number of timesteps between each dump of collisions debugging
    int debug_every;

    //! Method called in the main smilei loop to apply collisions at each timestep
    virtual void collide(PicParams&,std::vector<Species*>&,int){};

    int totbins;
    int start;






private:




};


#endif
