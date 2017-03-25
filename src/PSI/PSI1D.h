/*
PSI1D class
*/

#ifndef PSI1D_H
#define PSI1D_H

#include <vector>
#include "PSI.h"

class PSI1D : public PSI
{

public:
    //! Constructor for PSI between two species
    PSI1D(PicParams& params, SmileiMPI* smpi) : PSI(params, smpi)
    {
        nPartEmit_rem = 0.0;
    };
    virtual ~PSI1D(){};

    // sputtered particle number
    int nPartEmit;
    double nPartEmit_rem;

    //! Method called in the main smilei loop to apply collisions at each timestep
    virtual void performPSI(PicParams&,std::vector<Species*>&,int, ElectroMagn*){};

private:

};


#endif
