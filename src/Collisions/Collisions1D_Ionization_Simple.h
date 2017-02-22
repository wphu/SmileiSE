/*
Collisions1D_Ionization_Simple class
*/

#ifndef COLLISIONS1D_IONIZATION_SIMPLE_H
#define COLLISIONS1D_IONIZATION_SIMPLE_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "Collisions1D.h"
#include "H5.h"

using namespace std;

class Collisions1D_Ionization_Simple : public Collisions1D
{

public:
    //! Constructor for Collisions1D between two species
    Collisions1D_Ionization_Simple(PicParams&,std::vector<Species*>&,SmileiMPI*,unsigned int,std::vector<unsigned int>,std::vector<unsigned int>,std::vector<unsigned int>,double,bool,int);
    ~Collisions1D_Ionization_Simple();

    double cross_section(double ke);
    void calculate_scatter_velocity(double ke, double v_magnitude, double mass1, double mass2,
                                    vector<double>& momentum_unit, vector<double>& momentum_temp);


    //! Method called in the main smilei loop to apply collisions at each timestep
    void collide(PicParams&, SmileiMPI* smpi, ElectroMagn* fields, std::vector<Species*>&,int);

private:
    //>the ionization threshold energy
    double energy_ion;

};


#endif
