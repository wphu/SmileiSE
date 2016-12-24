#include "Interpolator1D1Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field1D.h"
#include "Particles.h"

using namespace std;

Interpolator1D1Order::Interpolator1D1Order(PicParams &params, SmileiMPI *smpi) : Interpolator1D(params, smpi) {
    dx_inv_ = 1.0/params.cell_length[0];
}

// ---------------------------------------------------------------------------------------------------------------------
// 1nd Order Interpolation of the fields at a the particle position (2 nodes are used)
// only has Electric field, no Magnetic field
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator1D1Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc)
{

    // Variable declaration
    double xjn, xjmxi, xjmxi2;

    // Static cast of the electromagnetic fields
    Field1D* Ex1D     = static_cast<Field1D*>(EMfields->Ex_);
    Field1D* Ey1D     = static_cast<Field1D*>(EMfields->Ey_);
    Field1D* Ez1D     = static_cast<Field1D*>(EMfields->Ez_);
    //Field1D* Bx1D_m   = static_cast<Field1D*>(EMfields->Bx_m);
    //Field1D* By1D_m   = static_cast<Field1D*>(EMfields->By_m);
    //Field1D* Bz1D_m   = static_cast<Field1D*>(EMfields->Bz_m);


    // Particle position (in units of the spatial-step)
    xjn    = particles.position(0, ipart)*dx_inv_;


    // --------------------------------------------------------
    // Interpolate the fields from the Primal grid : Ey, Ez, Bx
    // --------------------------------------------------------
    ip_      = floor(xjn);      // index of the central point
    xjmxi  = xjn -(double)ip_;  // normalized distance to the central node

    // 1nd order interpolation on 2 nodes
    coeffp_[0] = 1.0 - xjmxi;
    coeffp_[1] = xjmxi;

    ip_ -= index_domain_begin;

    (*ELoc).x = compute(coeffp_, Ex1D,   ip_);
    (*ELoc).y = compute(coeffp_, Ey1D,   ip_);
    (*ELoc).z = compute(coeffp_, Ez1D,   ip_);
    //(*BLoc).x = compute(coeffp_, Bx1D_m, ip_);
    //(*BLoc).y = compute(coeffp_, By1D_m, ip_);
    //(*BLoc).z = compute(coeffp_, Bz1D_m, ip_);



}//END Interpolator1D1Order




// ---------------------------------------------------------------------------------------------------------------------
// 1nd Order Interpolation of the fields at a the particle position (2 nodes are used)
// There are electric fields and magnetic fields
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator1D1Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc)
{

    // Variable declaration
    double xjn, xjmxi, xjmxi2;

    // Static cast of the electromagnetic fields
    Field1D* Ex1D       = static_cast<Field1D*>(EMfields->Ex_);
    Field1D* Ey1D       = static_cast<Field1D*>(EMfields->Ey_);
    Field1D* Ez1D       = static_cast<Field1D*>(EMfields->Ez_);
    Field1D* Bx1D       = static_cast<Field1D*>(EMfields->Bx_);
    Field1D* By1D       = static_cast<Field1D*>(EMfields->By_);
    Field1D* Bz1D       = static_cast<Field1D*>(EMfields->Bz_);


    // Particle position (in units of the spatial-step)
    xjn    = particles.position(0, ipart)*dx_inv_;


    // --------------------------------------------------------
    // Interpolate the fields from the Primal grid : Ey, Ez, Bx
    // --------------------------------------------------------
    ip_      = floor(xjn);      // index of the central point
    xjmxi  = xjn -(double)ip_;  // normalized distance to the central node

    // 1nd order interpolation on 2 nodes
    coeffp_[0] = 1.0 - xjmxi;
    coeffp_[1] = xjmxi;

    ip_ -= index_domain_begin;

    (*ELoc).x = compute(coeffp_, Ex1D, ip_);
    (*ELoc).y = compute(coeffp_, Ey1D, ip_);
    (*ELoc).z = compute(coeffp_, Ez1D, ip_);
    (*BLoc).x = compute(coeffp_, Bx1D, ip_);
    (*BLoc).y = compute(coeffp_, By1D, ip_);
    (*BLoc).z = compute(coeffp_, Bz1D, ip_);



}//END Interpolator1D1Order







void Interpolator1D1Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc)
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
