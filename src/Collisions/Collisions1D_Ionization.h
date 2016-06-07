/*
Collisions1D_Ionization class
*/

#ifndef COLLISIONS_IONIZATION_H
#define COLLISIONS_IONIZATION_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "Collisions1D.h"
#include "H5.h"

using namespace std;

class Collisions1D_Ionization : public Collisions1D
{

public:
    //! Constructor for Collisions1D between two species
    Collisions1D_Ionization(PicParams&,std::vector<Species*>&,SmileiMPI*,unsigned int,std::vector<unsigned int>,std::vector<unsigned int>,double,bool,int);
    ~Collisions1D_Ionization();

    double cross_section(double ke);
    void calculate_scatter_velocity(double ke, double v_magnitude, double mass1, double mass2,
                                    vector<int>& momentum_unit, vector<int>& momentum_temp);


    //! Method called in the main smilei loop to apply collisions at each timestep
    void collide(PicParams&,std::vector<Species*>&,int);

private:
    //>the ionization threshold energy
    double energy_ion;

};


#endif
