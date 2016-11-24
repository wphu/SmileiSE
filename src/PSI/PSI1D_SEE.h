
#ifndef PSI1D_SEE_H
#define PSI1D_SEE_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "PSI1D.h"
#include "H5.h"

using namespace std;

class PSI1D_SEE : public PSI1D
{

public:
    //! Constructor for Collisions between two species
    PSI1D_SEE(
        PicParams& params,
        SmileiMPI* smpi,
        unsigned int psi_species1,
        unsigned int psi_species2,
        string psiPosition,
        double emitTemperature,
        double SEEYield
    );
    ~PSI1D_SEE();

    // secondary electron emission yield
    double SEEYield;

    Particles new_particles;

    //! Method called in the main smilei loop to apply collisions at each timestep
    void performPSI(PicParams&, SmileiMPI*,std::vector<Species*>&,int, ElectroMagn*);

    void emit(PicParams& params, vector<Species*>& vecSpecies, unsigned int species_emit);
private:


};


#endif
