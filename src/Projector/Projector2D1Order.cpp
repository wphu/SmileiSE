#include "Projector2D1Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field2D.h"
#include "Particles.h"
#include "Tools.h"
#include "SmileiMPI_Cart2D.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Projector2D1Order
// ---------------------------------------------------------------------------------------------------------------------
Projector2D1Order::Projector2D1Order (PicParams& params, SmileiMPI* smpi) : Projector2D(params, smpi)
{
    SmileiMPI_Cart2D* smpi2D = static_cast<SmileiMPI_Cart2D*>(smpi);

    dx_inv_   = 1.0/params.cell_length[0];
    dx_ov_dt  = params.cell_length[0] / params.timestep;
    dy_inv_   = 1.0/params.cell_length[1];
    dy_ov_dt  = params.cell_length[1] / params.timestep;

    one_third = 1.0/3.0;

    i_domain_begin = smpi2D->getCellStartingGlobalIndex(0);
    j_domain_begin = smpi2D->getCellStartingGlobalIndex(1);

    DEBUG("cell_length "<< params.cell_length[0]);

}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Projector2D1Order
// ---------------------------------------------------------------------------------------------------------------------
Projector2D1Order::~Projector2D1Order()
{
}


//! Below, in this order :
//!   Project global current densities (EMfields->Jx_/Jy_/Jz_), not used
//!   Projection by species
//!   Project global current charge
//!   Project local current densities (sort)
//!   Project global current densities (ionize)




// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities (EMfields->Jx_/Jy_/Jz_), not used
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D1Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, double gf)
{

} // END Project global current densities, not used


// ---------------------------------------------------------------------------------------------------------------------
//!   Projection by species
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D1Order::operator() (Field* Jx, Field* Jy, Field* Jz, Field* rho, Particles &particles, int ipart, double gf)
{

}//END Projection by species



// ---------------------------------------------------------------------------------------------------------------------
//! Project global current charge
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D1Order::operator() (Field* rho, Particles &particles, int ipart, double weight)
{

    //Static cast of the total charge density
    Field2D* rho2D  = static_cast<Field2D*>(rho);

    //Declaration of local variables
    double delta, delta2;
    double rho_p = weight;   // charge density of the macro-particle
    double Sx[2], Sy[2];             // projection coefficient arrays

    //Locate particle on the primal grid & calculate the projection coefficients
    double       xpn = particles.position(0, ipart) * dx_inv_;  // normalized distance to the first node
    int ic  = floor(xpn);                   // index of the central node
    delta  = xpn - (double)ic;                       // normalized distance to the nearest grid point
    Sx[0]  = 1.0-delta;
    Sx[1]  = delta;

    double       ypn = particles.position(1, ipart) * dy_inv_;  // normalized distance to the first node
    int jc   = floor(ypn);                  // index of the central node
    delta  = ypn - (double)jc;                       // normalized distance to the nearest grid point
    Sy[0]  = 1.0-delta;
    Sy[1]  = delta;

    //cout << "Pos = " << particles.position(0, ipart) << " - i global = " << i << " - i local = " << i-index_domain_begin <<endl;

    int i = ic-i_domain_begin; // index of first point for projection in x
    int j = jc-j_domain_begin; // index of first point for projection in y

    // 2nd order projection for the total charge density
    for (unsigned int iloc=0 ; iloc<2 ; iloc++) {
        for (unsigned int jloc=0 ; jloc<2 ; jloc++) {
            (*rho2D)(i+iloc,j+jloc) += Sx[iloc]*Sy[jloc]*rho_p;
        }
    }

} // END Project global current charge



// ---------------------------------------------------------------------------------------------------------------------
//! Project global current charge
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D1Order::operator() (Field* rho, Particles &particles, int ipart)
{

    //Static cast of the total charge density
    Field2D* rho2D  = static_cast<Field2D*>(rho);

    //Declaration of local variables
    double delta, delta2;
    double rho_p = particles.weight(ipart);   // charge density of the macro-particle
    double Sx[2], Sy[2];             // projection coefficient arrays

    //Locate particle on the primal grid & calculate the projection coefficients
    double       xpn = particles.position(0, ipart) * dx_inv_;  // normalized distance to the first node
    int ic  = floor(xpn);                   // index of the central node
    delta  = xpn - (double)ic;                       // normalized distance to the nearest grid point
    Sx[0]  = 1.0-delta;
    Sx[1]  = delta;

    double       ypn = particles.position(1, ipart) * dy_inv_;  // normalized distance to the first node
    int jc   = floor(ypn);                  // index of the central node
    delta  = ypn - (double)jc;                       // normalized distance to the nearest grid point
    Sy[0]  = 1.0-delta;
    Sy[1]  = delta;

    //cout << "Pos = " << particles.position(0, ipart) << " - i global = " << i << " - i local = " << i-index_domain_begin <<endl;

    int i = ic-i_domain_begin; // index of first point for projection in x
    int j = jc-j_domain_begin; // index of first point for projection in y

    // 2nd order projection for the total charge density
    for (unsigned int iloc=0 ; iloc<2 ; iloc++) {
        for (unsigned int jloc=0 ; jloc<2 ; jloc++) {
            (*rho2D)(i+iloc,j+jloc) += Sx[iloc]*Sy[jloc]*rho_p;
        }
    }

} // END Project global current charge


// ---------------------------------------------------------------------------------------------------------------------
//! Project local current densities (sort)
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D1Order::operator() (double* Jx, double* Jy, double* Jz, double* rho, Particles &particles, int ipart, double gf, unsigned int bin, unsigned int b_dim1)
{

} // END Project local current densities (sort)


// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities (ionize)
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D1Order::operator() (Field* Jx, Field* Jy, Field* Jz, Particles &particles, int ipart, LocalFields Jion)
{
    ERROR("Projection of ionization current not yet defined for 2D 2nd order");

} // END Project global current densities (ionize)
