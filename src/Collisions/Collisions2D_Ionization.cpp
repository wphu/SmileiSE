#include "Collisions2D_Ionization.h"
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
Collisions2D_Ionization::Collisions2D_Ionization(PicParams& params, vector<Species*>& vecSpecies, SmileiMPI* smpi,
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

Collisions2D_Ionization::~Collisions2D_Ionization()
{

}


// Calculates the collisions for a given Collisions2D object
void Collisions2D_Ionization::collide(PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies, int itime)
{

    unsigned int nbins = vecSpecies[0]->bmin.size(); // number of bins
    vector<unsigned int> *sg1, *sg2, *sg3, *sgtmp;

    vector<vector<int> > index1, index2;
    vector<int> n1, n2;
    vector<double> density1, density2;
    double n1_max, n2_max;

    vector<double> momentum_unit(3, 0.0), momentum_temp(3, 0.0);
    int idNew;
    int totNCollision = 0;
    vector<int> bmin1, bmax1, bmin2, bmax2, bmin3, bmax3;
    unsigned int npairs; // number of pairs of macro-particles
    unsigned int i1, i2, i3;
    Species   *s1, *s2, *s3;
    Particles *p1, *p2, *p3;
    double m1, m2, m3, m12, W1, W2, W3;

    double  sigma_cr, sigma_cr_max, ke1, ke_primary, ke_secondary,
            ran, P_collision;
    double  v_square, v_magnitude, v_magnitude_primary, v_magnitude_secondary;


    Field2D *smean, *logLmean, *ncol;//, *temperature
    ostringstream name;
    hid_t did;

    sg1 = &species_group1;
    sg2 = &species_group2;
    sg3 = &species_group3;

    s1 = vecSpecies[(*sg1)[0]];      s2 = vecSpecies[(*sg2)[0]];        s3 = vecSpecies[(*sg3)[0]];
    p1 = &(s1->particles);           p2 = &(s2->particles);             p3 = &(s3->particles);
    m1 = s1->species_param.mass;     m2 = s2->species_param.mass;       m3 = s3->species_param.mass;
    W1 = p1->weight(i1);             W2 = p2->weight(i2);               W3 = p3->weight(i3);
    bmin1 = s1->bmin;                bmin2 = s2->bmin;                  bmin3 = s3->bmin;
    bmax1 = s1->bmax;                bmax2 = s2->bmax;                  bmax3 = s3->bmax;


    count_of_particles_to_insert_s1.resize(nbins);
    count_of_particles_to_insert_s3.resize(nbins);
    count_of_particles_to_erase_s2.resize(nbins);
    for(int ibin=0; ibin<nbins; ibin++)
    {
        count_of_particles_to_insert_s1[ibin] = 0;
        count_of_particles_to_insert_s3[ibin] = 0;
        count_of_particles_to_erase_s2[ibin] = 0;
    }

    totNCollision = 0;
    sigma_cr_max = maxCV(p1, m1);
    // Loop on bins
    for (unsigned int ibin=0 ; ibin<nbins ; ibin++) {



        //>calculate the particle number of species1 in each cell, and the indexs of particles in the cell
        index1.resize(params.n_space[1]);
        n1.resize(params.n_space[1]);
        density1.resize(params.n_space[1]);
        n1_max = 0.0;
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
            n1[i] = index1[i].size();
            density1[i] = n1[i] * W1;
            random_shuffle(index1[i].begin(), index1[i].end());
        }


        //>calculate the particle number of species2 in each cell, and the indexs of particles in the cell
        index2.resize(params.n_space[1]);
        n2.resize(params.n_space[1]);
        density2.resize(params.n_space[1]);
        n2_max = 0.0;
        for(int i = 0; i < index2.size(); i++)
        {
            index2[i].resize(0);
        }

        for(int iPart = bmin2[ibin]; iPart < bmax2[ibin]; iPart++)
        {
            double ypn = p2->position(1, iPart)*dy_inv_;
            double jp_ = floor(ypn);
            jp_ = jp_ - j_domain_begin;
            index2[jp_].push_back(iPart);
        }
        for(int i = 0; i < params.n_space[1]; i++)
        {
            n2[i] = index2[i].size();
            density2[i] = n2[i] * W2;
            random_shuffle(index2[i].begin(), index2[i].end());
        }



        // Now start the real loop in cells in y-direction of each bin
        // See equations in http://dx.doi.org/10.1063/1.4742167
        // ----------------------------------------------------
        for (int iCell = 0; iCell < params.n_space[1]; iCell++)
        {
            npairs = n1[iCell] * ( 1 - exp(-density2[iCell] * sigma_cr_max * timestep) );
            for(int i = 0; i < npairs; i++)
            {
                i1 = index1[iCell][i];
                i2 = index2[iCell][i];

                v_square = pow(p1->momentum(0,i1),2) + pow(p1->momentum(1,i1),2) + pow(p1->momentum(2,i1),2);
                v_magnitude = sqrt(v_square);
                //>kinetic energy of species1 (electrons)
                ke1 = 0.5 * m1 * v_square;
                //>energy_ion  is the ionization threshold energy
                ke_primary = ke1 - energy_ion;

                //> the energy of the secondary electron
                ke_secondary = 10.0 * tan(ran * atan(ke_primary / 20.0));
                //> the energy of the primary electron
                ke_primary -= ke_secondary;

                v_magnitude_primary = sqrt( 2.0 * ke_primary / m1 );
                v_magnitude_secondary = sqrt( 2.0 * ke_secondary / m1 );

                sigma_cr = v_magnitude * cross_section(ke1);
                P_collision = sigma_cr / sigma_cr_max;
                // Generate a random number between 0 and 1
                double ran_p = (double)rand() / RAND_MAX;

                if(ran_p < P_collision){
                    count_of_particles_to_erase_s2[ibin]++;
                    //>calculate the scatter velocity of primary electron
                    momentum_unit[0] = p1->momentum(0,i1) / v_magnitude;
                    momentum_unit[1] = p1->momentum(1,i1) / v_magnitude;
                    momentum_unit[2] = p1->momentum(2,i1) / v_magnitude;
                    //MESSAGE("v_magnitude"<<"  "<<v_magnitude_primary<<"  "<<v_magnitude_secondary);
                    //MESSAGE("momentum1"<<" "<<p1->momentum(0, i1)<<"  "<<p1->momentum(1, i1)<<"  "<<p1->momentum(2, i1));
                    calculate_scatter_velocity(v_magnitude_primary, m1, m2, momentum_unit, momentum_temp);
                    p1->momentum(0,i1) = momentum_temp[0];
                    p1->momentum(1,i1) = momentum_temp[1];
                    p1->momentum(2,i1) = momentum_temp[2];

                    //>calculate the scatter velocity of secondary electron
                    calculate_scatter_velocity(v_magnitude_secondary, m1, m2, momentum_unit, momentum_temp);
                    //>create new particle in the end of p1, we should sort_part when all bins are done!!!
                    //MESSAGE("momentum2"<<" "<<p1->momentum(0, i1)<<"  "<<p1->momentum(1, i1)<<"  "<<p1->momentum(2, i1));
                    p1->cp_particle(i1, new_particles1);
                    //MESSAGE("momentum1"<<" "<<new_particles1.momentum(0, idNew)<<"  "<<new_particles1.momentum(1, idNew));

                    idNew = new_particles1.size() - 1;
                    new_particles1.momentum(0, idNew) = momentum_temp[0];
                    new_particles1.momentum(1, idNew) = momentum_temp[1];
                    new_particles1.momentum(2, idNew) = momentum_temp[2];
                    count_of_particles_to_insert_s1[ibin]++;
                    //MESSAGE("momentum3"<<" "<<p1->momentum(0, idNew)<<"  "<<p1->momentum(1, idNew)<<"  "<<p1->momentum(2, idNew));

                    p2->cp_particle(i2, new_particles3);
                    //>create the ionized ion (species3)
                    idNew = new_particles3.size() - 1;
                    new_particles3.charge(idNew) = s3->species_param.charge;
                    count_of_particles_to_insert_s3[ibin]++;

                    indexes_of_particles_to_erase_s2.push_back(i2);
                    totNCollision++;
                }
            }
        }

    } // end loop on bins

    s2->erase_particles_from_bins(indexes_of_particles_to_erase_s2);
    indexes_of_particles_to_erase_s2.clear();

    s1->insert_particles_to_bins(new_particles1, count_of_particles_to_insert_s1);
    s3->insert_particles_to_bins(new_particles3, count_of_particles_to_insert_s3);
    new_particles1.clear();
    new_particles3.clear();

}


double Collisions2D_Ionization::maxCV(Particles* particles, double eMass){
    int nPart = particles->size();
    double v_square;
    double v_magnitude;
    double ke;
    double maxCrossSectionV = 0.0;
    double crossSectionV;

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
void Collisions2D_Ionization::calculate_scatter_velocity(double v_magnitude, double mass1, double mass2,
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



double Collisions2D_Ionization::cross_section(double ke)
{

}
