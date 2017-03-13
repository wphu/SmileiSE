
#ifndef PSI1D_RECYCLING_H
#define PSI1D_RECYCLING_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "PSI1D.h"
#include "H5.h"

class PSI1D_Recycling : public PSI1D
{

public:
    //! Constructor for Collisions between two species
    PSI1D_Recycling(
        PicParams& params,
        SmileiMPI* smpi,
        unsigned int psi_species1,
        unsigned int psi_species2,
        string psiPosition,
        double emitTemperature,
        double recycling_factor_temp
    );


    ~PSI1D_Recycling();

    //! Method called in the main smilei loop to apply PSI at each timestep
    void performPSI(PicParams&, SmileiMPI* smpi, std::vector<Species*>&,int, ElectroMagn*);

    // emit particles
    void emit(PicParams&, vector<Species*>&, unsigned int);

    // ===========Parameters for Recycling===========================
    // Recycling factor [0 ~ 1]
    double recycling_factor;


private:


};


#endif
