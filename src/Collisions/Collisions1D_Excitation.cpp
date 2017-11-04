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
            ran, P_collision;
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

        if(  smpi->getDomainLocalMin(0) + (ibin+1) * params.cell_length[0] < params.region_collision_zoom[0]
          || smpi->getDomainLocalMin(0) + ibin * params.cell_length[0] > params.region_collision_zoom[1] )
        {
            collision_zoom_factor = 1.0;
        }
        else
        {
            collision_zoom_factor = params.collision_zoom_factor;
        }

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

        // Now start the real loop
        // See equations in http://dx.doi.org/10.1063/1.4742167
        // ----------------------------------------------------
        npairs_double = n1[ibin] * (1 - exp(-density2[ibin] * sigma_cr_max * timesteps_collision * timestep * collision_zoom_factor) );
        npairs = npairs_double;
        npairsRem[ibin] += ( npairs_double - npairs );
        if(npairsRem[ibin] > 1)
        {
            npairsRem[ibin] = npairsRem[ibin] - 1.0;
            npairs++;
        }

        if(npairs > n1[ibin])
        {
            cout<<"npairs is larger than the particle number in a cell!!!"<<endl;
            cout<<"npairs, n1 are: "<<npairs<<" "<<n1[ibin]<<endl;
            npairs = n1[ibin];
        }
        if(npairs > n2[ibin])
        {
            cout<<"npairs is larger than the particle number in a cell!!!"<<endl;
            cout<<"npairs, n2 are: "<<npairs<<" "<<n2[ibin]<<endl;
            npairs = n2[ibin];
        }

        for(int i = 0; i < npairs; i++)
        {
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
            if(ran_p < P_collision){
                // Scatter the electrons
                momentum_unit[0] = p1->momentum(0,i1) / v_magnitude;
                momentum_unit[1] = p1->momentum(1,i1) / v_magnitude;
                momentum_unit[2] = p1->momentum(2,i1) / v_magnitude;
                calculate_scatter_velocity(ke_primary/const_e, v_magnitude_primary, m1, m2, momentum_unit, momentum_temp);

                DEBUGEXEC(
                    for(int idirection = 0; idirection < 3; idirection++)
                    {
                        if( isnan(momentum_temp[idirection]) || isinf(momentum_temp[idirection]) )
                        {
                            cout<<"Excitation Error: Species: "<<s1->species_param.species_type<<" old momentum "<<p1->momentum(idirection,i1)<<endl;
                            cout<<"Excitation Error: Species: "<<s1->species_param.species_type<<" new momentum "<<momentum_temp[idirection]<<endl;
                            cout<<"v_square: "<<v_square<<" ke_primary: "<<ke_primary<<" energy_excitation_threshold: "<<energy_excitation_threshold<<endl;
                        }
                    }
                );

                p1->momentum(0,i1) = momentum_temp[0];
                p1->momentum(1,i1) = momentum_temp[1];
                p1->momentum(2,i1) = momentum_temp[2];

                totNCollision++;
            }
        }


    } // end loop on bins


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





double Collisions1D_Excitation::cross_section(double ke)
{

}
