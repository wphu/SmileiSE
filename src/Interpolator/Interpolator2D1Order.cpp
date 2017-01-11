#include "Interpolator2D1Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field2D.h"
#include "Particles.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Creator for Interpolator2D1Order
// ---------------------------------------------------------------------------------------------------------------------
Interpolator2D1Order::Interpolator2D1Order(PicParams &params, SmileiMPI *smpi) : Interpolator2D(params, smpi)
{

    dx_inv_ = 1.0/params.cell_length[0];
    dy_inv_ = 1.0/params.cell_length[1];

}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd Order Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator2D1Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc)
{
    // Static cast of the electromagnetic fields
    Field2D* Ex2D = static_cast<Field2D*>(EMfields->Ex_);
    Field2D* Ey2D = static_cast<Field2D*>(EMfields->Ey_);
    Field2D* Ez2D = static_cast<Field2D*>(EMfields->Ez_);
    Field2D* Bx2D = static_cast<Field2D*>(EMfields->Bx_);
    Field2D* By2D = static_cast<Field2D*>(EMfields->By_);
    Field2D* Bz2D = static_cast<Field2D*>(EMfields->Bz_);


    // Normalized particle position
    double xpn = particles.position(0, ipart)*dx_inv_;
    double ypn = particles.position(1, ipart)*dy_inv_;


    // Indexes of the central nodes
    ip_ = floor(xpn);
    jp_ = floor(ypn);


    // Declaration and calculation of the coefficient for interpolation
    double delta;

    delta   = xpn - (double)ip_;
    coeffxp_[0] = 1.0 - delta;
    coeffxp_[1] = delta;


    delta   = ypn - (double)jp_;
    coeffyp_[0] = 1.0 - delta;
    coeffyp_[1] = delta;

    //!\todo CHECK if this is correct for both primal & dual grids !!!
    // First index for summation
    ip_ = ip_ - i_domain_begin;
    jp_ = jp_ - j_domain_begin;


    // -------------------------
    // Interpolation of Ex^(d,p)
    // -------------------------
    (*ELoc).x =  compute( coeffxp_, coeffyp_, Ex2D, ip_, jp_);

    // -------------------------
    // Interpolation of Ey^(p,d)
    // -------------------------
    (*ELoc).y = compute( coeffxp_, coeffyp_, Ey2D, ip_, jp_);

    // -------------------------
    // Interpolation of Ez^(p,p)
    // -------------------------
    (*ELoc).z = compute( coeffxp_, coeffyp_, Ez2D, ip_, jp_);

    // -------------------------
    // Interpolation of Bx^(p,d)
    // -------------------------
    (*BLoc).x = compute( coeffxp_, coeffyp_, Bx2D, ip_, jp_);

    // -------------------------
    // Interpolation of By^(d,p)
    // -------------------------
    (*BLoc).y = compute( coeffxp_, coeffyp_, By2D, ip_, jp_);

    // -------------------------
    // Interpolation of Bz^(d,d)
    // -------------------------
    (*BLoc).z = compute( coeffxp_, coeffyp_, Bz2D, ip_, jp_);

} // END Interpolator2D1Order

void Interpolator2D1Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc)
{

}
