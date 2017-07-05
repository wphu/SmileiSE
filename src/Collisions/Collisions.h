/*

Collisions class - Frederic Perez - 03/2015

This is based on the work described here
http://dx.doi.org/10.1063/1.4742167

Binary collisions, between macro-particles, are treated according to a scheme
first described by Nanbu (http://dx.doi.org/10.1103/PhysRevE.55.4642).

To include collisions in the simulations, add a block in the input file,
similar to the following:

# COLLISIONS
# species1    = list of strings, the names of the first species that collide
# species2    = list of strings, the names of the second species that collide
#               (can be the same as species1)
# coulomb_log = float, Coulomb logarithm. If negative or zero, then automatically computed.
Collisions(
	species1 = ["ion1"],
	species2 = ["electron1"],
	coulomb_log = 2.0
)

Several collision types can be defined. For each type, add a group "Collisions()".

*/

#ifndef COLLISIONS_H
#define COLLISIONS_H

#include <vector>
#include <string>
#include <fstream>


#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "H5.h"

class Particles;

using namespace std;

class Collisions
{

public:
    //! Constructor for Collisions between two species
    Collisions(PicParams &params)
    {
        const_e = params.const_e;
        timestep = params.timestep;
        const_pi = params.const_pi;
        const_ephi0 = params.const_ephi0;

        new_particles1.initialize(0, params);
        new_particles2.initialize(0, params);
        new_particles3.initialize(0, params);
        oversize = params.oversize;

        timesteps_collision = params.timesteps_collision;
    };
    virtual ~Collisions(){};

    //! Identification number of the Collisions object
    int n_collisions;

    //! Group of the species numbers that are associated for Collisions.
    //> species_group1 (2,3) may include more than one species only for coulomb collision and DSMC,
    //> species_group1 (2,3) only include one species for other collision.
    std::vector<unsigned int> species_group1, species_group2, species_group3;

    // record the lost particle indexes because of collision,
    // when all collisions are done, erase all lost particles at the same time
    std::vector<int> indexes_of_particles_to_erase_s1;
    std::vector<int> indexes_of_particles_to_erase_s2;
    std::vector<int> indexes_of_particles_to_erase_s3;

    std::vector<int> count_of_particles_to_erase_s1;
    std::vector<int> count_of_particles_to_erase_s2;
    std::vector<int> count_of_particles_to_erase_s3;


    // record the insert particle indexes because of collision,
    // when all collisions are done, insert all new particles to bins at the same time
    std::vector<int> indexes_of_particles_to_insert_s1;
    std::vector<int> indexes_of_particles_to_insert_s2;
    std::vector<int> indexes_of_particles_to_insert_s3;

    std::vector<int> count_of_particles_to_insert_s1;
    std::vector<int> count_of_particles_to_insert_s2;
    std::vector<int> count_of_particles_to_insert_s3;

    //
    Particles new_particles1;
    Particles new_particles2;
    Particles new_particles3;



    //! True if collisions inside a group of species, False if collisions between different groups of species
    bool intra_collisions;

    //! Number of timesteps between each dump of collisions debugging
    int debug_every;

    //! Method called in the main smilei loop to apply collisions at each timestep
    // relativistic case
    virtual void collide_relativistic(PicParams&, SmileiMPI* smpi, std::vector<Species*>&,int){};

    // non-relativistic case
    virtual void collide(PicParams&, SmileiMPI* smpi, ElectroMagn* fields, std::vector<Species*>&,int){};

    virtual void readCrossSection(){
        ifstream inFile;
        double energy, cross_section;
        crossSection.resize(2);
        inFile.open(crossSection_fileName.c_str());
        while(inFile >> energy && inFile >> cross_section){
            crossSection[0].push_back(energy);
            crossSection[1].push_back(cross_section);
            //crossSection[1].push_back(cross_section);
        }
        inFile.close();
    };
    // interplate cross section with energy (eV)
    double interpCrossSection(double energy){
        int low = 0;
        int high = crossSection[0].size() - 1;
        int mid = 0;

        if(energy <= crossSection[0][low] )
        {
            return 0.0;
        }
        else if(energy >= crossSection[0][high])
        {
            low = high -1;
            double dEnergy_inv = 1.0 / (crossSection[0][low] - crossSection[0][high]);
            return crossSection[1][high] * (crossSection[0][low] - energy) * dEnergy_inv +
                   crossSection[1][low] * (energy - crossSection[0][high]) * dEnergy_inv;
        }

        while(low <= high){
            mid = (low + high) / 2;
            if(crossSection[0][mid] < energy){
                low = mid + 1;
            }
            else if(crossSection[0][mid] > energy){
                high = mid - 1;
            }
            else {
                return crossSection[1][mid];
            }
        }
        // now low is 1 larger than high
        double dEnergy_inv = 1.0 / (crossSection[0][low] - crossSection[0][high]);
        return crossSection[1][high] * (crossSection[0][low] - energy) * dEnergy_inv +
               crossSection[1][low] * (energy - crossSection[0][high]) * dEnergy_inv;
    };

    // Simple method: Calculate electron scattered velocity for electron-neutral collisions
    //>the method is eqution (11) from the ref: a Monte Carlo collision model for the particle in cell method: applications to
    //>argon and oxygen discharges.
    //>and the code is transformed from C.F. Sang's fortran code
    void calculate_scatter_velocity_simple(double v_magnitude, double mass1, double mass2,
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
    };


    // Complex method: Calculate electron scattered velocity for electron-neutral collisions
    void calculate_scatter_velocity(double ke, double v_magnitude, double mass1, double mass2,
                                    vector<double>& momentum_unit, vector<double>& momentum_temp)
    {
        double up1, up2, up3;
        double r11, r12, r13, r21, r22, r23, r31, r32, r33;
        double mag;

        double ra = (double)rand() / RAND_MAX;
        double cosX = ( 2.0 + ke - 2.0 * pow(1.0+ke, ra) ) / ke;
        double sinX = sqrt(1.0 - abs(cosX * cosX) );

        ra = (double)rand() / RAND_MAX;
        double pi = 3.1415926;
        double phi = 2.0 * pi * ra;
        double cosphi = cos(phi);
        double sinphi = sin(phi);
        double sinTheta_inv;

        double ve=v_magnitude*sqrt(1.0-2.0*mass1*(1.0-cosX)/mass2);

        r11 = momentum_unit[0];
        r12 = momentum_unit[1];
        r13 = momentum_unit[2];

        if(r11 == 1.0 || r11 == -1.0)
        {
            up1= 0.;
            up2= 1.;
            up3= 0.;
        }
        else
        {
            up1= 1.;
            up2= 0.;
            up3= 0.;
        }

        double cosTheta = r11 * up1 + r12 * up2 + r13 * up3;
        sinTheta_inv = 1.0 / sqrt(1.0 - cosTheta*cosTheta);

        r21 = sinTheta_inv * ( r12 * up3 - r13 * up2 );
        r22 = sinTheta_inv * ( r13 * up1 - r11 * up3 );
        r23 = sinTheta_inv * ( r11 * up2 - r12 * up1 );


        r31 = sinTheta_inv * ( r22 * r13 - r23 * r12 );
        r32 = sinTheta_inv * ( r23 * r11 - r21 * r13 );
        r33 = sinTheta_inv * ( r21 * r12 - r22 * r11 );


        momentum_temp[0] = ve * (r11 * cosX + r21 * sinX * sinphi + r31 * sinX*cosphi);
        momentum_temp[1] = ve * (r12 * cosX + r22 * sinX * sinphi + r32 * sinX*cosphi);
        momentum_temp[2] = ve * (r13 * cosX + r23 * sinX * sinphi + r33 * sinX*cosphi);
    };





    int totbins;
    int nbins;
    int start;
    double norm_temperature;
    double const_e;
    double const_c;
    double const_pi;
    double const_ephi0;
    double const_h;
    double timestep;
    int timesteps_collision;
    std::vector<unsigned int> oversize;
    std::vector< double > npairsRem;

    string crossSection_fileName;
    vector< vector<double> > crossSection;





private:



};


#endif
