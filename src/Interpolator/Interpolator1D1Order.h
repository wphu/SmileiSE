#ifndef INTERPOLATOR1D1ORDER_H
#define INTERPOLATOR1D1ORDER_H


#include "Interpolator1D.h"

#include "Field1D.h"
//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 1d3v simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator1D1Order : public Interpolator1D
{

public:
    Interpolator1D1Order(PicParams&, SmileiMPI*);
    ~Interpolator1D1Order(){};

    void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc);
    void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc);

    inline double compute( double* coeff, Field1D* f, int idx) {
    	double interp_res =  coeff[0] * (*f)(idx)   + coeff[1] * (*f)(idx+1);
    	return interp_res;
    };

private:
    // Last prim index computed
    int ip_;
    // Last dual index computed
    int id_;
    // Interpolation coefficient on Prim grid
    double coeffp_[2];
    // Interpolation coefficient on Dual grid
    double coeffd_[2];


};//END class

#endif
