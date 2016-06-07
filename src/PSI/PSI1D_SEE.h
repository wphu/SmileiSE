
#ifndef PSI1D_SEE_H
#define PSI1D_SEE_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "PSI1D.h"
#include "H5.h"

class PSI1D_SEE : public PSI1D
{

public:
    //! Constructor for Collisions between two species
    PSI1D_SEE();
    ~PSI1D_SEE();



    //! Method called in the main smilei loop to apply collisions at each timestep
    void performPSI(PicParams&,std::vector<Species*>&,int);

private:


};


#endif
