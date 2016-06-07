
#ifndef PSI2D_SPUTTERING_H
#define PSI2D_SPUTTERING_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "PSI2D.h"
#include "H5.h"

class PSI2D_Sputtering : public PSI2D
{

public:
    //! Constructor
    PSI2D_Sputtering();
    ~PSI2D_Sputtering();



    //! Method called in the main smilei loop to apply PSI at each timestep
    void performPSI(PicParams&,std::vector<Species*>&,int);

private:


};


#endif
