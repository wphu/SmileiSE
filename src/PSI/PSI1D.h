/*
PSI1D class
*/

#ifndef PSI1D_H
#define PSI1D_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "H5.h"
#include "PSI.h"

class PSI1D : public PSI
{

public:
    //! Constructor for PSI between two species
    PSI1D(){};
    virtual ~PSI1D(){};


    //! Method called in the main smilei loop to apply collisions at each timestep
    virtual void performPSI(PicParams&,std::vector<Species*>&,int){};

private:

};


#endif
