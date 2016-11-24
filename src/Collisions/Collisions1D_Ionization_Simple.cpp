#include "Collisions1D_Ionization_Simple.h"
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
Collisions1D_Ionization_Simple::Collisions1D_Ionization_Simple(PicParams& param, vector<Species*>& vecSpecies, SmileiMPI* smpi,
                       unsigned int n_collisions,
                       vector<unsigned int> species_group1,
                       vector<unsigned int> species_group2,
                       vector<unsigned int> species_group3,
                       double coulomb_log,
                       bool intra_collisions,
                       int debug_every)
{

    n_collisions    = (n_collisions    );
    species_group1  = (species_group1  );
    species_group2  = (species_group2  );
    species_group3  = (species_group3  );
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

Collisions1D_Ionization_Simple::~Collisions1D_Ionization_Simple()
{

}


// Calculates the collisions for a given Collisions1D object
void Collisions1D_Ionization_Simple::collide(PicParams& params, vector<Species*>& vecSpecies, int itime)
{

    unsigned int nbins = vecSpecies[0]->bmin.size(); // number of bins
    vector<unsigned int> *sg1, *sg2, *sg3;

    vector<int> index1, index2;
    vector<int> n1, n2;
    vector<double> momentum_unit(3, 0.0), momentum_temp(3, 0.0);
    int idNew;
    int totNCollision = 0;
    vector<int> bmin1, bmax1, bmin2, bmax2, bmin3, bmax3;
    unsigned int npairs; // number of pairs of macro-particles
    vector<unsigned int> np1, np2; // numbers of macro-particles in each species, in each group
    double n12, n123, n223; // densities of particles
    unsigned int i1, i2, i3, ispec1, ispec2, ispec3;
    Species   *s1, *s2, *s3;
    Particles *p1, *p2, *p3;
    double m1, m2, m3, m12, W1, W2, W3;

    double  sigma_cr, sigma_cr_max, v_square, v_magnitude, ke1, ke_primary, ke_secondary,
            ran, P_collision, ran_P;

    double *temp=new double[3];
    double *vel=new double[3];
    std::vector<double> max_jutt_cumul;

    vel[0] = 0.0;
    vel[1] = 0.0;
    vel[2] = 0.0;
    temp[0] = 0.0;
    temp[1] = 0.0;
    temp[2] = 0.0;
    max_jutt_cumul.resize(1,0.0);

    Field2D *smean, *logLmean, *ncol;//, *temperature
    ostringstream name;
    hid_t did;

    sg1 = &species_group1;
    sg2 = &species_group2;

    s1 = vecSpecies[(*sg1)[0]];      s2 = vecSpecies[(*sg2)[0]];        s3 = vecSpecies[(*sg3)[0]];
    p1 = &(s1->particles);           p2 = &(s2->particles);             p3 = &(s3->particles);
    m1 = s1->species_param.mass;     m2 = s2->species_param.mass;       m3 = s3->species_param.mass;
    W1 = p1->weight(i1);             W2 = p2->weight(i2);               W3 = p3->weight(i3);
    bmin1 = s1->bmin;                bmin2 = s2->bmin;                  bmin3 = s3->bmin;
    bmax1 = s1->bmax;                bmax2 = s2->bmax;                  bmax3 = s3->bmax;
    // Loop on bins
    for (unsigned int ibin=0 ; ibin<nbins ; ibin++) {

        // Now start the real loop
        // See equations in http://dx.doi.org/10.1063/1.4742167
        // ----------------------------------------------------
        for (int iPart=bmin1[ibin] ; iPart<bmax1[ibin] ; iPart++)
        {
            v_square = pow(p1->momentum(0,iPart),2) + pow(p1->momentum(1,iPart),2) + pow(p1->momentum(2,iPart),2);
            v_magnitude = sqrt(v_square);

            //>kinetic energy of species1 (electrons)
            ke1 = 0.5 * m1 * v_square;
            //>energy_ion  is the ionization threshold energy
            ke_primary = ke1 - energy_ion;

            //> the energy of the secondary electron
            ke_secondary = 10.0 * tan(ran * atan(ke_primary / 20.0));
            //> the energy of the primary electron
            ke_primary -= ke_secondary;

            sigma_cr = v_magnitude * cross_section(ke1);
            P_collision = sigma_cr / sigma_cr_max;
            // Generate a random number between 0 and 1
            double ran_p = (double)rand() / RAND_MAX;

            if(ran_P < P_collision){
                momentum_unit[0] = p1->momentum(0,iPart) / v_magnitude;
                momentum_unit[1] = p1->momentum(1,iPart) / v_magnitude;
                momentum_unit[2] = p1->momentum(2,iPart) / v_magnitude;

                //>calculate the scatter velocity of primary electron
                calculate_scatter_velocity(ke_primary, v_magnitude, m1, m2, momentum_unit, momentum_temp);
                p1->momentum(0,iPart) = momentum_temp[0];
                p1->momentum(1,iPart) = momentum_temp[1];
                p1->momentum(2,iPart) = momentum_temp[2];


                //>calculate the scatter velocity of secondary electron
                calculate_scatter_velocity(ke_secondary, v_magnitude, m1, m2, momentum_unit, momentum_temp);
                //>create new particle in the end of p1, we should sort_part when all bins are done!!!
                p1->create_particle();
                idNew = p1->size() - 1;
                p1->momentum(0, idNew) = momentum_temp[0];
                p1->momentum(1, idNew) = momentum_temp[1];
                p1->momentum(2, idNew) = momentum_temp[2];

                //>create the ionized ion (species3)
                p3->create_particle();
                idNew = p3->size() - 1;

                s3->initMomentum(1,idNew, temp, vel,
                             s2->species_param.initMomentum_type, max_jutt_cumul, params);

                p3->position(0, idNew) = p1->position(0, iPart);
                p3->position(1, idNew) = p1->position(1, iPart);
                p3->position(2, idNew) = p1->position(2, iPart);

                p3->position_old(0, idNew) = p3->position(0, idNew);
                p3->position_old(1, idNew) = p3->position(1, idNew);
                p3->position_old(2, idNew) = p3->position(2, idNew);

                p3->charge(idNew) = p2->charge(0) + 1;
                p3->weight(idNew) = p1->weight(iPart);


                //>swap the particle of p2 to the end, we should erase all the particles
                //>(which are swapped to the end), when all bins are done!!!
                //idNew = p2->size() - 1 - totNCollision;
                //p2->swap_part(i2, idNew);

                totNCollision++;

            }
        } // end loop on particles

    } // end loop on bins

    //>delete ionized neutrals
    idNew = p2->size() - totNCollision;
    p2->erase_particle_trail(idNew);
    //>update the bmax, because created particles are put the end of species1 and species3
    //>then sort_part to move the particles to the bins which they belong to
    s1->bmax.back() += totNCollision;
    s1->sort_part();
    s3->bmax.back() += totNCollision;
    s3->sort_part();

}



//>the method is eqution (11) from the ref: a Monte Carlo collision model for the particle in cell method: applications to
//>argon and oxygen discharges.
//>and the code is transformed from C.F. Sang's fortran code
void Collisions1D_Ionization_Simple::calculate_scatter_velocity(double ke, double v_magnitude, double mass1, double mass2,
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



double Collisions1D_Ionization_Simple::cross_section(double ke)
{

}
