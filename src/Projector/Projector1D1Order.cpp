#include "Projector1D1Order.h"

#include <cmath>
#include <iostream>
#include <iomanip>

#include "ElectroMagn.h"
#include "Field1D.h"
#include "Particles.h"
#include "Tools.h"
#include "SmileiMPI_Cart1D.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Projector1D1Order
// ---------------------------------------------------------------------------------------------------------------------
Projector1D1Order::Projector1D1Order (PicParams& params, SmileiMPI* smpi) : Projector1D(params, smpi)
{
    SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);

    dx_inv_  = 1.0/params.cell_length[0];
    dx_ov_dt = params.cell_length[0] / params.timestep;

    index_domain_begin = smpi1D->getCellStartingGlobalIndex(0);

    DEBUG("cell_length "<< params.cell_length[0]);

}

Projector1D1Order::~Projector1D1Order()
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
void Projector1D1Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, double gf)
{

} // END Project global current densities, not used


// ---------------------------------------------------------------------------------------------------------------------
//!   Projection by species
// ---------------------------------------------------------------------------------------------------------------------
void Projector1D1Order::operator() (Field* Jx, Field* Jy, Field* Jz, Field* rho, Particles &particles, int ipart, double gf)
{

} // END Projection by species


// ---------------------------------------------------------------------------------------------------------------------
//! Project global current charge
// ---------------------------------------------------------------------------------------------------------------------
void Projector1D1Order::operator() (Field* rho, Particles &particles, int ipart)
{
    Field1D* rho1D  = static_cast<Field1D*>(rho);


    //Declaration of local variables
    int i;
    double xjn,xjmxi;
    double rho_j = particles.charge(ipart)*particles.weight(ipart);  // charge density of the macro-particle

    //Locate particle on the grid
    xjn    = particles.position(0, ipart) * dx_inv_;  // normalized distance to the first node
    i      = floor(xjn);                   // index of the central node
    xjmxi  = xjn - (double)i;              // normalized distance to the nearest grid point

    //cout << "Pos = " << particles.position(0, ipart) << " - i global = " << i << " - i local = " << i-index_domain_begin <<endl;

    i -= index_domain_begin;

    // 1nd order projection for the total density
    //#pragma omp atomic
    (*rho1D)( i  )  += (1.0 - xjmxi)    * rho_j ;
    //#pragma omp atomic
    (*rho1D)( i+1)  += xjmxi             * rho_j;

} // END Project global current charge

// ---------------------------------------------------------------------------------------------------------------------
//! Project local current densities (sort)
// ---------------------------------------------------------------------------------------------------------------------
void Projector1D1Order::operator() (double* Jx, double* Jy, double* Jz, double* rho, Particles &particles, int ipart, double gf, unsigned int bin, unsigned int b_dim0)
{


} // END Project local current densities (sort)

// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities (ionize)
// ---------------------------------------------------------------------------------------------------------------------
void Projector1D1Order::operator() (Field* Jx, Field* Jy, Field* Jz, Particles &particles, int ipart, LocalFields Jion)
{


} // END Project global current densities (ionize)
