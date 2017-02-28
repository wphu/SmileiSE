
#ifndef PSI1D_BACKSCATTERING_H
#define PSI1D_BACKSCATTERING_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "PSI1D.h"
#include "H5.h"

class PSI1D_Backscattering : public PSI1D
{

public:
    //! Constructor for Collisions between two species
    PSI1D_Backscattering(
        PicParams& params,
        SmileiMPI* smpi,
        unsigned int psi_species1,
        unsigned int psi_species2,
        string psiPosition,
        double emitTemperature
    );


    ~PSI1D_Backscattering();

    //! Method called in the main smilei loop to apply PSI at each timestep
    void performPSI(PicParams&, SmileiMPI* smpi, std::vector<Species*>&,int, ElectroMagn*);

    // emit particles
    void emit(PicParams&, vector<Species*>&, unsigned int);

    // parameters for sputtering
    double an1;     // atomic number of incident atomic
    double an2;     // atomic number of target atomic
    double am1;     // atomic mass of incident atomic (amu)
    double am2;     // atomic mass of target atomic (amu)
    double es;      //surface binding energy (heat of sublimation) of target (eV).
                    // es = 8.7 for Carbon
    double ionflag; //ionflag -> flag for light/heavy ion.
                    //ionflag = 0 => light ion sputtering.
                    //ionflag = 1 => heavy ion sputtering.

    double n, Q, eth,eth1, mu, etf, aL, Mratio;

    double rnion( theta,energy,nz1,m1,ne,nz2,nw );
    double reion(theta,energy,nz1,m1,ne,nz2,nw);

private:


};


#endif
