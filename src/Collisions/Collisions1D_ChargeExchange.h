/*
Collisions1D class
*/

#ifndef COLLISIONS1D_CHARGEEXCHANGE_H
#define COLLISIONS1D_CHARGEEXCHANGE_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "Collisions1D.h"
#include "H5.h"

class Collisions1D_ChargeExchange : public Collisions1D
{

public:
    //! Constructor for Collisions1D between two species
    Collisions1D_ChargeExchange(PicParams&,std::vector<Species*>&,SmileiMPI*,unsigned int,
                                std::vector<unsigned int>,std::vector<unsigned int>,
                                string);
    ~Collisions1D_ChargeExchange();

    double cross_section(double ke);

    //! Method called in the main smilei loop to apply collisions at each timestep
    void collide(PicParams&, SmileiMPI* smpi, ElectroMagn* fields, std::vector<Species*>&,int);

    // get the maximum value of crossSection*velocity
    double maxCV(Species* s1, Species* s2);

private:
    //>the ionization threshold energy
    double energy_ion;

};


#endif
