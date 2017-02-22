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

        new_particles1.initialize(0, params);
        new_particles2.initialize(0, params);
        new_particles3.initialize(0, params);
        oversize = params.oversize;
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
            crossSection[1].push_back(cross_section*2.0e4);
            //crossSection[1].push_back(cross_section);
        }
        inFile.close();
    };
    // interplate cross section with energy (eV)
    double interpCrossSection(double energy){
        int low = 0;
        int high = crossSection[0].size() - 1;
        int mid = 0;
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
    std::vector<unsigned int> oversize;

    string crossSection_fileName;
    vector< vector<double> > crossSection;





private:



};


#endif
