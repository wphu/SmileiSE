/*
Collisions1D_Ionization class
*/

#ifndef COLLISIONS1D_RECOMBINATION_TB_H
#define COLLISIONS1D_RECOMBINATION_TB_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "Collisions1D.h"
#include "H5.h"

using namespace std;

class Particles;

class Collisions1D_Recombination_TB : public Collisions1D
{

public:
    //! Constructor for Collisions1D between two species
    Collisions1D_Recombination_TB(  PicParams&,std::vector<Species*>&,SmileiMPI*,unsigned int,
                                    std::vector<unsigned int>,std::vector<unsigned int>,std::vector<unsigned int>,
                                    string );
    ~Collisions1D_Recombination_TB();


    // ke1 and ke2 are the energys of the two electrons
    double cross_section(double ke1, double ke2);
    void calculate_scatter_velocity(double v_magnitude, double mass1, double mass2,
                                    vector<double>& momentum_unit, vector<double>& momentum_temp);


    //! Method called in the main smilei loop to apply collisions at each timestep
    void collide(PicParams&, SmileiMPI* smpi, ElectroMagn* fields, std::vector<Species*>&,int);

    // get the maximum value of crossSection*velocity
    double maxCV(Particles* particles, double eMass);

private:
    //>the ionization threshold energy
    double energy_ionization_threshold;

};


#endif
