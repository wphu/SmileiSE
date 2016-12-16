/*
Collisions2D_ChargeExchange class
*/

#ifndef COLLISIONS2D_CHARGEEXCHANGE_H
#define COLLISIONS2D_CHARGEEXCHANGE_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "Collisions2D.h"
#include "H5.h"


class Collisions2D_ChargeExchange : public Collisions2D
{

public:
    //! Constructor
    Collisions2D_ChargeExchange(PicParams&,std::vector<Species*>&,SmileiMPI*,unsigned int,
                                std::vector<unsigned int>,std::vector<unsigned int>,
                                string);
    ~Collisions2D_ChargeExchange();

    double cross_section(double ke);
    //! Method called in the main smilei loop to apply collisions at each timestep
    void collide(PicParams&, SmileiMPI* smpi, std::vector<Species*>&,int);

private:
    //>the ionization threshold energy
    double energy_ion;

};


#endif
