#include "Collisions1D_Ionization.h"
#include "SmileiMPI.h"
#include "Field2D.h"
#include "H5.h"


#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <sstream>

using namespace std;


// Constructor
Collisions1D_Ionization::Collisions1D_Ionization(PicParams& param, vector<Species*>& vecSpecies, SmileiMPI* smpi,
                       unsigned int n_collisions,
                       vector<unsigned int> species_group1,
                       vector<unsigned int> species_group2,
                       double coulomb_log,
                       bool intra_collisions,
                       int debug_every)
{

    n_collisions    = (n_collisions    );
    species_group1  = (species_group1  );
    species_group2  = (species_group2  );
    coulomb_log     = (coulomb_log     );
    intra_collisions= (intra_collisions);
    debug_every     = (debug_every     );
    start           = (0               );



    // Calculate total number of bins
    int nbins = vecSpecies[0]->bmin.size();
    totbins = nbins;
    //MPI_Allreduce( smpi->isMaster()?MPI_IN_PLACE:&totbins, &totbins, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
    MPI_Reduce( smpi->isMaster()?MPI_IN_PLACE:&totbins, &totbins, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);

    // if debug requested, prepare hdf5 file



}

Collisions1D_Ionization::~Collisions1D_Ionization()
{

}


// Calculates the collisions for a given Collisions1D object
void Collisions1D_Ionization::collide(PicParams& params, vector<Species*>& vecSpecies, int itime)
{


}



//>the method is eqution (11) from the ref: a Monte Carlo collision model for the particle in cell method: applications to
//>argon and oxygen discharges.
//>and the code is transformed from C.F. Sang's fortran code
void Collisions1D_Ionization::calculate_scatter_velocity(double ke, double v_magnitude, double mass1, double mass2,
vector<int>& momentum_unit, vector<int>& momentum_temp)
{
    double up1, up2, up3;
    double r11, r12, r13, r21, r22, r23, r31, r32, r33;
    double mag;

    double ra = (double)rand() / RAND_MAX;
    double costheta = 1.0 - 2.0 * ra;
    double sintheta = sqrt(1.0 - abs(costheta * costheta) );

    ra = (double)rand() / RAND_MAX;
    double pi = 3.1415926;
    double phi = 2.0 * pi * ra;
    double cosphi = cos(phi);
    double sinphi = sin(phi);

    double ve=v_magnitude*sqrt(1.0-2.0*mass1*(1.0-costheta)/mass2);

    r13 = momentum_unit[0];
    r23 = momentum_unit[1];
    r33 = momentum_unit[2];
    if(r33 == 1.0 ){
        up1= 0.;
        up2= 1.;
        up3= 0.;
    }
    else{
        up1= 0.;
        up2= 0.;
        up3= 1.;
    }

    r12 = r23 * up3 - r33 * up2;
    r22 = r33 * up1 - r13 * up3;
    r32 = r13 * up2 - r23 * up1;
    mag = sqrt(r12 * r12 + r22 * r22 + r32 * r32);
    r12 = r12 / mag;
    r22 = r22 / mag;
    r32 = r32 / mag;
    r11 = r22 * r33 - r32 * r23;
    r21 = r32 * r13 - r12 * r33;
    r31 = r12 * r23 - r22 * r13;
    momentum_temp[0] = ve * (r11 * sintheta * cosphi + r12 * sintheta * sinphi + r13 * costheta);
    momentum_temp[1] = ve * (r21 * sintheta * cosphi + r22 * sintheta * sinphi + r23 * costheta);
    momentum_temp[2] = ve * (r31 * sintheta * cosphi + r32 * sintheta * sinphi + r33 * costheta);


}



double Collisions1D_Ionization::cross_section(double ke)
{

}
