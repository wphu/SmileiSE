#include "Collisions2D_ChargeExchange.h"
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
Collisions2D_ChargeExchange::Collisions2D_ChargeExchange(PicParams& params, vector<Species*>& vecSpecies, SmileiMPI* smpi,
                       unsigned int n_col,
                       vector<unsigned int> sg1,
                       vector<unsigned int> sg2,
                       string CS_fileName)
: Collisions2D(params)
{

    n_collisions    = n_col;
    species_group1  = sg1;
    species_group2  = sg2;
    crossSection_fileName = CS_fileName;




    // Calculate total number of bins
    nbins = vecSpecies[0]->bmin.size();
    totbins = nbins;
    //MPI_Allreduce( smpi->isMaster()?MPI_IN_PLACE:&totbins, &totbins, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
    MPI_Reduce( smpi->isMaster()?MPI_IN_PLACE:&totbins, &totbins, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);




}

Collisions2D_ChargeExchange::~Collisions2D_ChargeExchange()
{

}


// Calculates the collisions for a given Collisions2D object
void Collisions2D_ChargeExchange::collide(PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies, int itime)
{

    unsigned int nbins = vecSpecies[0]->bmin.size(); // number of bins
    vector<unsigned int> *sg1, *sg2, *sg3, *sgtmp;

    vector<vector<int> > index1, index2;
    vector<int> n1, n2;
    vector<double> momentum_unit(3, 0.0), momentum_temp(3, 0.0);
    int idNew;
    int totNCollision = 0;
    vector<int> bmin1, bmax1, bmin2, bmax2, bmin3, bmax3;

    unsigned int nspec1, nspec2; // numbers of species in each group
    unsigned int npart1, npart2; // numbers of macro-particles in each group
    unsigned int npairs; // number of pairs of macro-particles
    vector<unsigned int> np1, np2; // numbers of macro-particles in each species, in each group
    double n12, n123, n223; // densities of particles
    unsigned int i1, i2, ispec1, ispec2;
    Species   *s1, *s2;
    Particles *p1, *p2;
    double m1, m2, m12, W1, W2, qqm, qqm2, gamma1, gamma2, gamma12, gamma12_inv,
           COM_vx, COM_vy, COM_vz, COM_vsquare, COM_gamma,
           term1, term2, term3, term4, term5, term6, coeff1, coeff2, coeff3, coeff4, twoPi,
           vcv1, vcv2, px_COM, py_COM, pz_COM, p2_COM, p_COM, gamma1_COM, gamma2_COM,
           logL, bmin, s, vrel, smax,
           cosX, sinX, phi, sinXcosPhi, sinXsinPhi, p_perp, inv_p_perp,
           newpx_COM, newpy_COM, newpz_COM, U, vcp;

   double  sigma_cr, sigma_cr_max, v_square, v_magnitude, ke1, ke_primary, ke_secondary,
           ran, P_collision, Pi;

    Field2D *smean, *logLmean, *ncol;//, *temperature
    ostringstream name;
    hid_t did;

    sg1 = &species_group1;
    sg2 = &species_group2;

    s1 = vecSpecies[(*sg1)[0]];      s2 = vecSpecies[(*sg2)[0]];
    p1 = &(s1->particles);           p2 = &(s2->particles);
    m1 = s1->species_param.mass;     m2 = s2->species_param.mass;
    W1 = p1->weight(i1);             W2 = p2->weight(i2);
    bmin1 = s1->bmin;                bmin2 = s2->bmin;
    bmax1 = s1->bmax;                bmax2 = s2->bmax;
    // Loop on bins
    for (unsigned int ibin=0 ; ibin<nbins ; ibin++) {



        //>calculate the particle number of species1 in each cell, and the indexs of particles in the cell
        index1.resize(params.n_space[1]);
        n1.resize(params.n_space[1]);
        for(int i = 0; i < index1.size(); i++)
        {
            index1[i].resize(0);
        }

        for(int iPart = bmin1[ibin]; iPart < bmax1[ibin]; iPart++)
        {
            double ypn = p1->position(1, iPart)*dy_inv_;
            int jp_ = floor(ypn);
            jp_ = jp_ - j_domain_begin;
            index1[jp_].push_back(iPart);
        }
        for(int i = 0; i < params.n_space[1]; i++)
        {
            n1[i] = index1.size() / params.cell_volume;
            random_shuffle(index1[i].begin(), index1[i].end());
        }


        //>calculate the particle number of species2 in each cell, and the indexs of particles in the cell
        index2.resize(params.n_space[1]);
        n2.resize(params.n_space[1]);
        for(int i = 0; i < index2.size(); i++)
        {
            index2[i].resize(0);
        }

        for(int iPart = bmin2[ibin]; iPart < bmax2[ibin]; iPart++)
        {
            double ypn = p2->position(1, iPart)*dy_inv_;
            int jp_ = floor(ypn);
            jp_ = jp_ - j_domain_begin;
            index2[jp_].push_back(iPart);
        }
        for(int i = 0; i < params.n_space[1]; i++)
        {
            n2[i] = index2.size() / params.cell_volume;
            random_shuffle(index2[i].begin(), index2[i].end());
        }



        // Now start the real loop
        // See equations in http://dx.doi.org/10.1063/1.4742167
        // ----------------------------------------------------
        for (int iCell = 0; iCell < params.n_space[1]; iCell++)
        {
            npairs = n1[iCell] * (1 - exp(-n2[iCell] * sigma_cr_max) );
            for(int i = 0; i < npairs; i++)
            {
                i1 = i2 = i;

                v_square = pow(p1->momentum(0,i1),2) + pow(p1->momentum(1,i1),2) + pow(p1->momentum(2,i1),2);
                v_magnitude = sqrt(v_square);
                //>kinetic energy of species1 (electrons)
                ke1 = 0.5 * m1 * v_square;
                //>energy_ion  is the ChargeExchange threshold energy
                ke_primary = ke1 - energy_ion;

                //> the energy of the secondary electron
                ke_secondary = 10.0 * tan(ran * atan(ke_primary / 20.0));
                //> the energy of the primary electron
                ke_primary -= ke_secondary;

                sigma_cr = v_magnitude * cross_section(ke1);
                Pi = sigma_cr / sigma_cr_max;
                // Generate a random number between 0 and 1
                double ran_p = (double)rand() / RAND_MAX;

                if(ran_p < Pi){
                    W1              = p1->weight(i1);
                    p1->weight(i1)  = p2->weight(i2);
                    p2->weight(i2)  = W1;
                    totNCollision++;
                }
            }
        }

    } // end loop on bins


}



double Collisions2D_ChargeExchange::cross_section(double ke)
{

}
