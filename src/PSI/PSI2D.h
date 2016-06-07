/*
PSI2D class
*/

#ifndef PSI2D_H
#define PSI2D_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "H5.h"
#include "PSI.h"

class PSI2D : public PSI
{

public:
    //! Constructor for PSI between two species
    PSI2D(){};
    virtual ~PSI2D(){};




    //! Method called in the main smilei loop to apply PSI at each timestep
    virtual void performPSI(PicParams&,std::vector<Species*>&,int){};





private:



};


#endif
