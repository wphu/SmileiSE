/* ==============================================================
Subroutine to evaluate empirical formulas for number and energy
Backscattering coefficients of light ions incident on Elemental
and Compound targets
Ref: Subroutines for some plasma surface interaction processes:
     hpysical sputtering, chemical erosion, radiation enhanced
     sublimation, backscattering and thermal evaporation.
================================================================*/
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
    void emit(PicParams&, vector<Species*>&);

    // parameters for sputtering
    // the name is from original fortran77 code
    int nz1;         // atomic number of incident atomic
    int m1;          // atomic mass of incident atomic (amu)
    int ne;             // number of constituent elements in the target.
    vector<int> nz2;	// array for atomic numbers of the constituents.
    vector<int> nw;     // array for relative numbers of the constituents.
    //	example :   for Tritium ions of 1 KeV energy incident on the TiO2
    //		        ( Titanium Dioxide ) target, Energy=1000 , nz1=1 ,
    //		        m1=3 , ne=2 , nz2(1)=22 , nw(1)=1 , nz2(2)=8 and
    //		        nw(2)=2 .



    // rnion and reion are number and energy backscattering coefficients, respectively
    // theta is angle of incidence in degree
    // ke is the incident kinetic energy in eV
    void scatter(double &rnion, double &reion, double theta, double ke);

    double rangen(double eps, double energy);
    double expint(double x);




    // ===========data table alt:======================================
    double  a1t[95] =
    {
        // z1 = 1 ( for hydrogen , deuterium and tritium )
        1.0072766,2.0135536,3.1055011,
        // z1 = 2 ( for helium-3 and helium-4 )
        3.0149325,4.0015059,
        // z1 = 3-10
        6.941,9.01218,10.811,12.011,14.0067,15.9994,18.998403,20.179,
        // z1 = 11 - 20
        22.98977,24.305,26.98154,28.0855,30.97376,32.066,35.453,39.948,
        39.0893,40.078,
        // z1 = 21 - 30
        44.95591,47.88,50.9415,51.9961,54.9380,55.847,58.9332,58.69,
        63.546,65.39,
        // z1 = 31 - 40
        69.723,72.59,74.9216,78.96,79.904,83.80,85.4678,87.62,88.9059,
        91.224,
        // z1 = 41 - 50
        92.9064,95.94,99.0,101.07,102.9055,106.42,107.8682,112.41,114.82,
        118.710,
        // z1 = 51 - 60
        121.75,127.60,126.9045,131.29,132.9054,137.33,138.9055,140.12,
        140.9077,144.24,
        // z1 = 61 - 70
        145.0,150.36,151.96,157.25,158.9254,162.50,164.9304,167.26,
        168.9342,173.04,
        // z1 = 71 - 80
        174.967,178.49,180.9479,183.85,186.207,190.2,192.22,195.08,
        196.9655,200.59,
        //z1 = 81 - 90
        204.383,207.2,208.9804,210.0,210.0,222.0,223.0,226.0,227.0,232.0,
        // z12 = 91 - 92
        231.0,238.0289
    };

    // ===========data table a2t =========================================
    double a2t[92] =
    {
        //z2 = 1 - 10
        1.00794,4.002602,6.941,9.01218,10.811,12.011,14.0067,15.9994,
        18.998403,20.179,
        //z2 = 11 - 20
        22.98977,24.305,26.98154,28.0855,30.97376,32.066,35.453,39.948,
        39.0893,40.078,
        // z2 = 21 - 30
        44.95591,47.88,50.9415,51.9961,54.9380,55.847,58.9332,58.69,
        63.546,65.39,
        // z2 = 31 - 40
        69.723,72.59,74.9216,78.96,79.904,83.80,85.4678,87.62,88.9059,
        91.224,
        // z2 = 41 - 50
        92.9064,95.94,99.0,101.07,102.9055,106.42,107.8682,112.41,114.82,
        118.710,
        // z2 = 51 - 60
        121.75,127.60,126.9045,131.29,132.9054,137.33,138.9055,140.12,
        140.9077,144.24,
        // z2 = 61 - 70
        145.0,150.36,151.96,157.25,158.9254,162.50,164.9304,167.26,
        168.9342,173.04,
        // z2 = 71 - 80
        174.967,178.49,180.9479,183.85,186.207,190.2,192.22,195.08,
        196.9655,200.59,
        //z2 = 81 - 90
        204.383,207.2,208.9804,210.0,210.0,222.0,223.0,226.0,227.0,232.0,
        //z2 = 91 - 92
        231.0,238.0289
    };

    // ===============data table d :======================================
    double  d[92] =
    {
        //z2 = 1 - 10
        0.9341e0,0.6693e0,0.6654e0,0.9712e0,1.007e0,1.024e0,1.111e0,
        0.9699e0,0.7357e0,0.6842e0,
        // z2 = 11 - 20
        0.8769e0,1.290e0,1.395e0,1.378e0,1.063e0,1.123e0,1.632e0,1.839e0,
        1.642e0,1.749e0,
        // z2 = 21 - 30
        1.638e0,1.523e0,1.396e0,1.236e0,1.072e0,1.083e0,0.9624e0,1.085e0,
        1.125e0,1.277e0,
        // z2 = 31 - 40
        1.525e0,1.675e0,1.601e0,1.762e0,1.679e0,1.914e0,1.696e0,1.884e0,
        1.9e0,1.993e0,
        // z2 = 41 - 50
        2.039e0,1.894e0,2.001e0,1.795e0,1.738e0,1.534e0,1.644e0,1.698e0,
        1.816e0,1.866e0,
        // z2 = 51 - 60
        2.181e0,2.027e0,2.240e0,2.384e0,2.108e0,2.283e0,2.321e0,2.159e0,
        2.1e0,2.042e0,
        // z2 = 61 - 70
        1.986e0,1.932e0,1.879e0,1.931e0,1.779e0,1.578e0,1.492e0,1.448e0,
        1.406e0,1.365e0,
        // z2 = 71 - 80
        1.394e0,1.431e0,1.348e0,1.3e0,1.477e0,1.439e0,1.403e0,1.269e0,
        1.376e0,1.22e0,
        // z2 = 81 - 90
        1.336e0,1.504e0,1.683e0,1.739e0,1.751e0,1.744e0,1.959e0,2.115e0,
        2.167e0,2.170e0,
        // z2 = 91 - 92
        2.084e0,2.050e0
    };

    inline double fneps(double e, double z1, double a1, double z2, double a2 )
    {
        double f = sqrt( pow(z1, 2.0/3.0) + pow(z2, 2.0/3.0) );
        return 0.032534 * e / ( z1 * z2 * (1.0 + a1/a2) * f );
    }

    inline double rne(double th, double r0, double as1, double as2)
    {
        return r0 + (1.0 - r0) / ( 1.0 + as1 / pow( tan(th), (2.0*as2) ) );
    }



private:


};


#endif
