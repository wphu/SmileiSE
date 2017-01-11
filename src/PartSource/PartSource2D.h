/*
PartSource2D class
*/

#ifndef PARTSOURCE2D_H
#define PARTSOURCE2D_H

#include <vector>
#include "PartSource.h"

class PartSource2D : public PartSource
{

public:
    //! Constructor for PSI between two species
    PartSource2D(PicParams& params, SmileiMPI* smpi) : PartSource(params, smpi) {};
    virtual ~PartSource2D(){};

    // sputtered particle number
    int nPartEmit;

    //! Method called in the main smilei loop to apply collisions at each timestep
    virtual void emitLoad(PicParams&,std::vector<Species*>&,int, ElectroMagn*){};

private:

};


#endif
