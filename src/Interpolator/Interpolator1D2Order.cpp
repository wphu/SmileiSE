#include "Interpolator1D2Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field1D.h"
#include "Particles.h"

using namespace std;

Interpolator1D2Order::Interpolator1D2Order(PicParams &params, SmileiMPI *smpi) : Interpolator1D(params, smpi) {
    dx_inv_ = 1.0/params.cell_length[0];
}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd Order Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator1D2Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc)
{

    // Variable declaration
    double xjn, xjmxi, xjmxi2;

    // Static cast of the electromagnetic fields
    Field1D* Ex1D     = static_cast<Field1D*>(EMfields->Ex_);
    Field1D* Ey1D     = static_cast<Field1D*>(EMfields->Ey_);
    Field1D* Ez1D     = static_cast<Field1D*>(EMfields->Ez_);
    Field1D* Bx1D_m   = static_cast<Field1D*>(EMfields->Bx_m);
    Field1D* By1D_m   = static_cast<Field1D*>(EMfields->By_m);
    Field1D* Bz1D_m   = static_cast<Field1D*>(EMfields->Bz_m);


    // Particle position (in units of the spatial-step)
    xjn    = particles.position(0, ipart)*dx_inv_;


    // --------------------------------------------------------
    // Interpolate the fields from the Primal grid : Ey, Ez, Bx
    // --------------------------------------------------------
    ip_      = round(xjn);      // index of the central point
    xjmxi  = xjn -(double)ip_;  // normalized distance to the central node
    xjmxi2 = pow(xjmxi,2);      // square of the normalized distance to the central node

    // 2nd order interpolation on 3 nodes
    coeffp_[0] = 0.5 * (xjmxi2-xjmxi+0.25);
    coeffp_[1] = (0.75-xjmxi2);
    coeffp_[2] = 0.5 * (xjmxi2+xjmxi+0.25);

    ip_ -= index_domain_begin;

    (*ELoc).y = compute(coeffp_, Ey1D,   ip_);  
    (*ELoc).z = compute(coeffp_, Ez1D,   ip_);  
    (*BLoc).x = compute(coeffp_, Bx1D_m, ip_);  

    // --------------------------------------------------------
    // Interpolate the fields from the Dual grid : Ex, By, Bz
    // --------------------------------------------------------
    id_      = round(xjn+0.5);        // index of the central point
    xjmxi  = xjn - (double)id_ +0.5;  // normalized distance to the central node
    xjmxi2 = pow(xjmxi,2);            // square of the normalized distance to the central node

    // 2nd order interpolation on 3 nodes
    coeffd_[0] = 0.5 * (xjmxi2-xjmxi+0.25);
    coeffd_[1] = (0.75-xjmxi2);
    coeffd_[2] = 0.5 * (xjmxi2+xjmxi+0.25);

    id_ -= index_domain_begin;

    (*ELoc).x = compute(coeffd_, Ex1D,   id_);  
    (*BLoc).y = compute(coeffd_, By1D_m, id_);  
    (*BLoc).z = compute(coeffd_, Bz1D_m, id_);  

}//END Interpolator1D2Order

void Interpolator1D2Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc)
{
    // Interpolate E, B
    // Compute coefficient for ipart position
    (*this)(EMfields, particles, ipart, ELoc, BLoc);

    // Static cast of the electromagnetic fields
    Field1D* Jx1D     = static_cast<Field1D*>(EMfields->Jx_);
    Field1D* Jy1D     = static_cast<Field1D*>(EMfields->Jy_);
    Field1D* Jz1D     = static_cast<Field1D*>(EMfields->Jz_);
    Field1D* Rho1D    = static_cast<Field1D*>(EMfields->rho_);
    
    // --------------------------------------------------------
    // Interpolate the fields from the Primal grid : Jy, Jz, Rho
    // --------------------------------------------------------
    (*JLoc).y = compute(coeffp_, Jy1D,  ip_);  
    (*JLoc).z = compute(coeffp_, Jz1D,  ip_);  
    (*RhoLoc) = compute(coeffp_, Rho1D, ip_);    

    // --------------------------------------------------------
    // Interpolate the fields from the Dual grid : Jx
    // --------------------------------------------------------
    (*JLoc).x = compute(coeffd_, Jx1D,  id_);  
    
}
