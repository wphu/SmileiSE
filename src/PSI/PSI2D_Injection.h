
#ifndef PSI2D_INJECTION_H
#define PSI2D_INJECTION_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "PSI2D.h"
#include "H5.h"

class PSI2D_Injection : public PSI2D
{

public:
    //! Constructor for Collisions between two species
    PSI2D_Injection();
    ~PSI2D_Injection();



    //! Method called in the main smilei loop to apply PSI at each timestep
    void performPSI(PicParams&,std::vector<Species*>&,int);

private:


};


#endif
