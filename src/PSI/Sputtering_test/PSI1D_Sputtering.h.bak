/* ==============================================================
Subroutine to evaluate sputtering yeilds of any
monoatomic target material due to Physical Sputtering
Ref: Subroutines for some plasma surface interaction processes:
     hpysical sputtering, chemical erosion, radiation enhanced
     sublimation, backscattering and thermal evaporation.

!!! 戴舒宇和桑超峰的程序输出结果不一样，此处是改自戴舒宇的程序，但桑超峰的程序是和文献给的程序一样
================================================================*/



#ifndef PSI1D_SPUTTERING_H
#define PSI1D_SPUTTERING_H

#include <vector>



class PSI1D_Sputtering
{

public:
    //! Constructor for Collisions between two species
    PSI1D_Sputtering(
    double an1_in,
    double am1_in,
    double an2_in,
    double am2_in,
    double es_in,
    double density_solid_in );


    ~PSI1D_Sputtering();


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
    double n;       // number density---unit-- C/A**3 ---*

    double Q, eth,eth1, mu, etf, aL, Mratio;

    void init();
    double phy_sput_yield(double ke, double theta);

private:


};


#endif
