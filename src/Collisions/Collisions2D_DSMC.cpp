#include "Collisions2D_DSMC.h"
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
Collisions2D_DSMC::Collisions2D_DSMC(PicParams& params, vector<Species*>& vecSpecies, SmileiMPI* smpi,
                       unsigned int n_collisions,
                       vector<unsigned int> species_group1,
                       vector<unsigned int> species_group2,
                       double coulomb_log,
                       bool intra_collisions,
                       int debug_every)
: Collisions2D(params)
{

    n_collisions    = (n_collisions    );
    species_group1  = (species_group1  );
    species_group2  = (species_group2  );
    coulomb_log     = (coulomb_log     );
    intra_collisions= (intra_collisions);
    debug_every     = (debug_every     );
    start           = (0               );



    // Calculate total number of bins
    nbins = vecSpecies[0]->bmin.size();
    totbins = nbins;
    //MPI_Allreduce( smpi->isMaster()?MPI_IN_PLACE:&totbins, &totbins, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
    MPI_Reduce( smpi->isMaster()?MPI_IN_PLACE:&totbins, &totbins, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);



}

Collisions2D_DSMC::~Collisions2D_DSMC()
{

}





// Calculates the collisions for a given Collisions2D object
void Collisions2D_DSMC::collide(PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies, int itime)
{

    unsigned int nbins = vecSpecies[0]->bmin.size(); // number of bins
    vector<unsigned int> *sg1, *sg2, *sgtmp;
    vector<int> bmin1, bmax1, bmin2, bmax2;

    vector<vector<vector<int> > > index1, index2;
    vector<int> index1_tot, index2_tot;
    vector<int> n1, n2;
    vector<double> momentum_unit(3, 0.0), momentum_temp(3, 0.0);
    int idNew;
    int totNCollision = 0;





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

   double  sigma, sigma_cr, sigma_cr_max, sigma_cr_max_temp, cr, v_square, v_magnitude, ke1, ke_primary, ke_secondary,
           ran, P_collision, ran_P, Pi;

    Field2D *smean, *logLmean, *ncol;//, *temperature
    ostringstream name;
    hid_t did;

    sg1 = &species_group1;
    sg2 = &species_group2;


    // Loop on bins
    for (unsigned int ibin=0 ; ibin<nbins ; ibin++) {


        for (ispec1=0 ; ispec1<nspec1 ; ispec1++) {
            s1 = vecSpecies[(*sg1)[ispec1]];
            bmin1[ispec1] = s1->bmin[ibin];
            bmax1[ispec1] = s1->bmax[ibin];
            np1[ispec1] = bmax1[ispec1] - bmin1[ispec1];
            npart1 += np1[ispec1];
        }



        //>calculate the particle number of species1 in each cell, and the indexs of particles in the cell
        index1.resize(params.n_space[1]);
        np1.resize(nspec1);
        for(int i = 0; i < index1.size(); i++)
        {
            index1[i].resize(nspec1);
            for(int j = 0; j < index1[i].size(); j++)
            {
                index1[i][j].resize(0);
            }
        }


        //> calculate the particle number of each species in each cell, and the indexs of particles in the cell
        for (ispec1=0 ; ispec1<nspec1 ; ispec1++)
        {
            s1 = vecSpecies[(*sg1)[ispec1]];
            p1 = &(s1->particles);
            m1 = s1->species_param.mass;
            W1 = p1->weight(i1);
            bmin1 = s1->bmin;
            bmax1 = s1->bmax;

            for(int iPart = bmin1[ibin]; iPart < bmax1[ibin]; iPart++)
            {
                double ypn = p1->position(1, iPart)*dy_inv_;
                int jp_ = floor(ypn);
                jp_ = jp_ - j_domain_begin;
                index1[jp_][ispec1].push_back(iPart);
            }
        }



        // Now start the real loop
        //> the code is mainly from the website: https://www.particleincell.com/2012/html5-dsmc/
        //> the original code language is javascript
        // ----------------------------------------------------
        for (int iCell = 0; iCell < params.n_space[1]; iCell++)
        {
            sigma_cr_max_temp = 0.0;
            //> calculate the total number of all neutral particles (all species) in the cell
            npart1 = 0;
            for(ispec1=0 ; ispec1<nspec1 ; ispec1++)
            {
                np1[ispec1] = index1[iCell][ispec1].size();
                npart1 += np1[ispec1];
            }
            //> for now, the weights of all neutral species are the same (W1)
            npairs = int (0.5 * npart1 * npart1 * W1 * sigma_cr_max * params.timestep  / params.cell_volume);

            //> for updating sigma_cr_max
            sigma_cr_max_temp = 0;

            for(int i =0; i < npairs; i++)
            {
                //> select the first particle
                i1 = int( ((double)rand() / RAND_MAX) * npart1 );
                //> select the second one, making sure the two particles are unique
                do{
                    i2 = int( ((double)rand() / RAND_MAX) * npart1 );
                }
                while(i1 == i2);


                // find species and index i1 of particle "1"
                for (ispec1=0 ; ispec1<nspec1 ; ispec1++) {
                    if (i1 < np1[ispec1]) break;
                    i1 -= np1[ispec1];
                }
                i1 += bmin1[ispec1];
                // find species and index i2 of particle "2"
                //i2 = index2[i];
                for (ispec2=0 ; ispec2<nspec1 ; ispec2++) {
                    if (i2 < np1[ispec2]) break;
                    i2 -= np1[ispec2];
                }
                i2 += bmin2[ispec2];

                //> only one species group, but one or two particles (the particle is a c++ class)
                s1 = vecSpecies[(*sg1)[ispec1]]; s2 = vecSpecies[(*sg1)[ispec2]];
                p1 = &(s1->particles);           p2 = &(s2->particles);
                m1 = s1->species_param.mass;     m2 = s2->species_param.mass;
                W1 = p1->weight(i1);             W2 = p2->weight(i2);

                //> relative velocity
                cr = relative_velocity(p1, i1, p2, i2);
                sigma = evalSigma(cr);
                sigma_cr = sigma * cr;
                if(sigma_cr > sigma_cr_max_temp){
                    sigma_cr_max_temp = sigma_cr;
                }
                P_collision = sigma_cr / sigma_cr_max;
                double ran_p = (double)rand() / RAND_MAX;

                if(ran_P < P_collision){
                    scatter_particles(p1, i1, p2, i2);
                }

            }


        } //> end loop on cells




    } // end loop on bins

}




void Collisions2D_DSMC::scatter_particles(Particles* particle1, int iPart1, Particles* particle2, int iPart2)
{
    double rv;
    rv = sqrt( pow((particle1->momentum(0, iPart1) - particle2->momentum(0, iPart2)), 2)
            +  pow((particle1->momentum(1, iPart1) - particle2->momentum(1, iPart2)), 2)
            +  pow((particle1->momentum(2, iPart1) - particle2->momentum(2, iPart2)), 2)
            );
}


double Collisions2D_DSMC::evalSigma(double cr)
{


}
