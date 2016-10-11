/*
PSI2D class
*/

#ifndef PSI2D_H
#define PSI2D_H

#include <vector>
#include "PSI.h"

class PSI2D : public PSI
{

public:
    //! Constructor for PSI between two species
    PSI2D(PicParams& params, SmileiMPI* smpi) : PSI(params, smpi){};
    virtual ~PSI2D(){};




    //! Method called in the main smilei loop to apply PSI at each timestep
    virtual void performPSI(PicParams&,std::vector<Species*>&,int){};





private:



};


#endif
