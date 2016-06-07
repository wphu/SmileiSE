#ifndef INTERPOLATOR2D1ORDER_H
#define INTERPOLATOR2D1ORDER_H


#include "Interpolator2D.h"
#include "Field2D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 1d3v simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator2D1Order : public Interpolator2D
{

public:
    Interpolator2D1Order(PicParams&, SmileiMPI*);
    ~Interpolator2D1Order(){};

    void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc);
    void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc);
    inline double compute( double* coeffx, double* coeffy, Field2D* f, int idx, int idy) {
	double interp_res(0.);
	for (int iloc=0 ; iloc<2 ; iloc++) {
	    for (int jloc=0 ; jloc<2 ; jloc++) {
		interp_res += coeffx[iloc] * coeffy[jloc] * (*f)(idx+iloc,idy+jloc);
	    }
	}
	return interp_res;
    };

private:
    // Last prim index computed
    int ip_, jp_;
    // Interpolation coefficient on Prim grid
    double coeffxp_[2], coeffyp_[2];
    // Interpolation coefficient on Dual grid
    double coeffxd_[2], coeffyd_[2];


};//END class

#endif
