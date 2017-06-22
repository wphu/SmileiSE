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


    // interplate cross section with energy (eV)
    double interpCrossSection(double energy){
        int low = 0;
        int high = crossSection[0].size() - 1;
        int mid = 0;

        if(energy <= crossSection[0][low] )
        {
            return crossSection[1][low];
        }
        else if(energy >= crossSection[0][high])
        {
            low = high -1;
            double dEnergy_inv = 1.0 / (crossSection[0][low] - crossSection[0][high]);
            return crossSection[1][high] * (crossSection[0][low] - energy) * dEnergy_inv +
                   crossSection[1][low] * (energy - crossSection[0][high]) * dEnergy_inv;
        }

        while(low <= high){
            mid = (low + high) / 2;
            if(crossSection[0][mid] < energy){
                low = mid + 1;
            }
            else if(crossSection[0][mid] > energy){
                high = mid - 1;
            }
            else {
                return crossSection[1][mid];
            }
        }
        // now low is 1 larger than high
        double dEnergy_inv = 1.0 / (crossSection[0][low] - crossSection[0][high]);
        return crossSection[1][high] * (crossSection[0][low] - energy) * dEnergy_inv +
               crossSection[1][low] * (energy - crossSection[0][high]) * dEnergy_inv;
    };



private:
    //>the ionization threshold energy
    double energy_ion;

};


#endif
