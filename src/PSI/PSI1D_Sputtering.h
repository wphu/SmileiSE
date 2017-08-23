/* ==============================================================
Subroutine to evaluate sputtering yeilds of any
monoatomic target material due to Physical Sputtering
Ref: Subroutines for some plasma surface interaction processes:
     hpysical sputtering, chemical erosion, radiation enhanced
     sublimation, backscattering and thermal evaporation.

!!! 戴舒宇和桑超峰的程序，靶板密度的单位不一样，桑超峰的是g/cm^3，戴书宇的是atoms/Ai^3 （Ai = 10^-10 m）
!!! 此处是改自戴舒宇的程序，所以密度单位要转换成　atoms/Ai^3
================================================================*/



#ifndef PSI1D_SPUTTERING_H
#define PSI1D_SPUTTERING_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "PSI1D.h"
#include "H5.h"

class PSI1D_Sputtering : public PSI1D
{

public:
    //! Constructor for Collisions between two species
    PSI1D_Sputtering(
        PicParams& params,
        SmileiMPI* smpi,
        vector<Species*>& vecSpecies,
        unsigned int psi_species1,
        unsigned int psi_species2,
        string psiPosition,
        double emitTemperature
    );


    ~PSI1D_Sputtering();

    //! Method called in the main smilei loop to apply PSI at each timestep
    void performPSI(PicParams&, SmileiMPI* smpi, std::vector<Species*>&,int, ElectroMagn*);

    // emit particles
    void emit(PicParams&, vector<Species*>&);

    // parameters for sputtering
    double an1;     // atomic number of incident atomic
    double an2;     // atomic number of target atomic
    double am1;     // atomic mass of incident atomic (amu)
    double am2;     // atomic mass of target atomic (amu)
    double es;      //surface binding energy (heat of sublimation) of target (eV).
                    // W: 11.75(old value 8.68) C: 7.41
    double ionflag; //ionflag -> flag for light/heavy ion.
                    //ionflag = 0 => light ion sputtering.
                    //ionflag = 1 => heavy ion sputtering.
    double n;       // number density, unit: atoms/Ai^3 （Ai = 10^-10 m）
                    // C: 0.11286 (2.248 g/cm^3), W: (19.35 g/cm^3)

    double Q, eth,eth1, mu, etf, aL, Mratio;

    void init(std::vector<Species*>&);
    double phy_sput_yield(double ke, double theta);

private:


};


#endif
