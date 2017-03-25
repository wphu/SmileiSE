#include "Collisions1D_Excitation.h"
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
Collisions1D_Excitation::Collisions1D_Excitation(PicParams& params, vector<Species*>& vecSpecies, SmileiMPI* smpi,
    unsigned int n_col,
    vector<unsigned int> sg1,
    vector<unsigned int> sg2,
    string CS_fileName)
    : Collisions1D(params)
{

    n_collisions    = n_col;
    species_group1  = sg1;
    species_group2  = sg2;
    crossSection_fileName = CS_fileName;

    // Calculate total number of bins
    nbins = vecSpecies[0]->bmin.size();
    totbins = nbins;

    //MPI_Allreduce( smpi->isMaster()?MPI_IN_PLACE:&totbins, &totbins, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
    //MPI_Reduce( smpi->isMaster()?MPI_IN_PLACE:&totbins, &totbins, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);

    readCrossSection();
    energy_excitation_threshold = crossSection[0][0];

}

Collisions1D_Excitation::~Collisions1D_Excitation()
{

}


// Calculates the collisions for a given Collisions1D object
void Collisions1D_Excitation::collide(PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies, int itime)
{
    vector<unsigned int> *sg1, *sg2;

    vector<int> index1, index2;
    vector<int> n1, n2;
    vector<double> density1, density2;
    double n1_max, n2_max;
    vector<double> momentum_unit(3, 0.0), momentum_temp(3, 0.0);
    int idNew;
    int totNCollision = 0;
    vector<int> bmin1, bmax1, bmin2, bmax2;
    unsigned int npairs; // number of pairs of macro-particles
    double npairs_double;
    unsigned int i1, i2;
    Species   *s1, *s2;
    Particles *p1, *p2;
    double m1, m2, W1, W2;

    double  sigma_cr, sigma_cr_max, ke1, ke_primary, ke_secondary,
            ran, P_collision, ran_P;
    double  v_square, v_magnitude, v_magnitude_primary, v_magnitude_secondary;


    sg1 = &species_group1;
    sg2 = &species_group2;

    // electons                         atoms or primary ions
    s1 = vecSpecies[(*sg1)[0]];      s2 = vecSpecies[(*sg2)[0]];
    p1 = &(s1->particles);           p2 = &(s2->particles);
    m1 = s1->species_param.mass;     m2 = s2->species_param.mass;
    W1 = p1->weight(0);              W2 = p2->weight(0);
    bmin1 = s1->bmin;                bmin2 = s2->bmin;
    bmax1 = s1->bmax;                bmax2 = s2->bmax;

    n1.resize(nbins);
    density1.resize(nbins);
    n1_max = 0.0;
    for(unsigned int ibin=0 ; ibin<nbins ; ibin++)
    {
        n1[ibin] = bmax1[ibin] - bmin1[ibin];
    }

    n2.resize(nbins);
    density2.resize(nbins);
    n2_max = 0.0;
    for (unsigned int ibin=0 ; ibin<nbins ; ibin++)
    {
        n2[ibin] = bmax2[ibin] - bmin2[ibin];
        density2[ibin] = n2[ibin] * W2;
        if( density2[ibin] > n2_max ) { n2_max = density2[ibin]; };
    }

    totNCollision = 0;
    sigma_cr_max = maxCV(p1, m1);
    for (unsigned int ibin=0 ; ibin<nbins ; ibin++) {
        //MESSAGE("nbins000"<<"  "<<ibin<<"  "<<bmin2[ibin]<<" "<<bmax2[ibin]);
        //>calculate the particle number of species1 in each cell, and the indexs of particles in the cell
        index1.resize( n1[ibin] );

        for(int iPart = 0; iPart < n1[ibin]; iPart++)
        {
            index1[iPart] = bmin1[ibin] + iPart;
        }
        random_shuffle(index1.begin(), index1.end());

        //>calculate the particle number of species2 in each cell, and the indexs of particles in the cell
        index2.resize( n2[ibin] );

        for(int iPart = 0; iPart < n2[ibin]; iPart++)
        {
            index2[iPart] = bmin2[ibin] + iPart;
        }
        random_shuffle(index2.begin(), index2.end());

        //smpi->barrier();
        //MESSAGE("nbinsaaaa"<<"  "<<ibin<<"  "<<n1[ibin]<<" "<<n2[ibin]);
        // Now start the real loop
        // See equations in http://dx.doi.org/10.1063/1.4742167
        // ----------------------------------------------------
        npairs_double = n1[ibin] * (1 - exp(-density2[ibin] * sigma_cr_max * timestep) );
        npairs = npairs_double;
        npairsRem[ibin] += ( npairs_double - npairs );
        if(npairsRem[ibin] > 1)
        {
            npairsRem[ibin] = npairsRem[ibin] - 1.0;
            npairs++;
        }
        //if(npairs > index1.size() || npairs > index2.size()) {ERROR("npairs > index in collisions");}

        //smpi->barrier();
        //MESSAGE("nbins111"<<"  "<<ibin<<"  "<<n1[ibin]<<" "<<n2[ibin]);
        for(int i = 0; i < npairs; i++)
        {
            //MESSAGE("nparis111"<<"  "<<i);
            i1 = index1[i];
            i2 = index2[i];

            v_square = pow(p1->momentum(0,i1),2) + pow(p1->momentum(1,i1),2) + pow(p1->momentum(2,i1),2);
            v_magnitude = sqrt(v_square);
            //>kinetic energy of species1 (electrons)
            ke1 = 0.5 * m1 * v_square;
            ke_primary = ke1 - energy_excitation_threshold * const_e;
            v_magnitude_primary = sqrt( 2.0 * ke_primary / m1 );

            sigma_cr = v_magnitude * interpCrossSection( ke1 / const_e );
            P_collision = sigma_cr / sigma_cr_max;

            // Generate a random number between 0 and 1
            double ran_p = (double)rand() / RAND_MAX;
            if(ran_P < P_collision){
                // Scatter the electrons
                momentum_unit[0] = p1->momentum(0,i1) / v_magnitude;
                momentum_unit[1] = p1->momentum(1,i1) / v_magnitude;
                momentum_unit[2] = p1->momentum(2,i1) / v_magnitude;
                calculate_scatter_velocity(v_magnitude_primary, m1, m2, momentum_unit, momentum_temp);
                p1->momentum(0,i1) = momentum_temp[0];
                p1->momentum(1,i1) = momentum_temp[1];
                p1->momentum(2,i1) = momentum_temp[2];

                totNCollision++;
            }
            //MESSAGE("nparis222"<<"  "<<i);
        }
        //smpi->barrier();
        //MESSAGE("nbins222"<<"  "<<ibin);

    } // end loop on bins
    //smpi->barrier();
    //MESSAGE("aaaa"<<" "<<s1->bmax.back()<<" "<<p1->size());
    // swap lost particles to the end for ionized neutrals

}

double Collisions1D_Excitation::maxCV(Particles* particles, double eMass){
    int nPart = particles->size();
    double v_square;
    double v_magnitude;
    double ke;
    double maxCrossSectionV = 0.0;
    double crossSectionV;

    maxCrossSectionV = 0.0;
    for(unsigned int iPart = 0; iPart < nPart; iPart++)
    {
        v_square = particles->momentum(0,iPart) * particles->momentum(0,iPart) +
                          particles->momentum(1,iPart) * particles->momentum(1,iPart) +
                          particles->momentum(2,iPart) * particles->momentum(2,iPart);
        v_magnitude = sqrt(v_square);
        // ke is energy (eV)
        ke = 0.5 * eMass * v_square / const_e;
        crossSectionV = v_magnitude * interpCrossSection( ke );
        if(crossSectionV > maxCrossSectionV) {maxCrossSectionV = crossSectionV;};
    }
    return maxCrossSectionV;
}

//>the method is eqution (11) from the ref: a Monte Carlo collision model for the particle in cell method: applications to
//>argon and oxygen discharges.
//>and the code is transformed from C.F. Sang's fortran code
void Collisions1D_Excitation::calculate_scatter_velocity( double v_magnitude, double mass1, double mass2,
    vector<double>& momentum_unit, vector<double>& momentum_temp)
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



double Collisions1D_Excitation::cross_section(double ke)
{

}
