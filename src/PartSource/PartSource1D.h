/*
PSI1D class
*/

#ifndef PARTSOURCE1D_H
#define PARTSOURCE1D_H

#include <vector>
#include "PartSource.h"

class PartSource1D : public PartSource
{

public:
    //! Constructor for PSI between two species
    PartSource1D(PicParams& params, SmileiMPI* smpi) : PartSource(params, smpi) {};
    virtual ~PartSource1D(){};

    // sputtered particle number
    int nPartEmit;

    //! Method called in the main smilei loop to apply collisions at each timestep
    virtual void emitLoad(PicParams&,std::vector<Species*>&,int, ElectroMagn*){};

private:

};


#endif
