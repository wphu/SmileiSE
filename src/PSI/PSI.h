/*
PSI class
*/

#ifndef PSI_H
#define PSI_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "H5.h"

class PSI
{

public:
    //! Constructor for PSI between two species
    PSI(){};
    virtual ~PSI(){};

    //! Group of the species numbers that are associated for PSI.
    //> for PSI_Injection, only species1 is used;
    //> for sputtering and secondary electron emission, species1 is the incident particle.
    unsigned int species1, species2;


    //! Method called in the main smilei loop to apply PSI at each timestep
    virtual void performPSI(PicParams&,std::vector<Species*>&,int){};


private:

};

#endif
