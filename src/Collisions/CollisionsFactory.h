#ifndef COLLISIONSFACTORY_H
#define COLLISIONSFACTORY_H

#include "Collisions1D_ChargeExchange.h"
#include "Collisions1D_Coulomb.h"
#include "Collisions1D_DSMC.h"
#include "Collisions1D_Elastic.h"
#include "Collisions1D_Ionization.h"

#include "Collisions2D_ChargeExchange.h"
#include "Collisions2D_Coulomb.h"
#include "Collisions2D_DSMC.h"
#include "Collisions2D_Elastic.h"
#include "Collisions2D_Ionization.h"



#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <sstream>
using namespace std;

class CollisionsFactory {

public:

	// Reads the input file and creates the Collisions objects accordingly
	static vector<Collisions*> create(PicParams& params, InputData &ifile, vector<Species*>& vecSpecies, SmileiMPI* smpi)
	{
	    vector<Collisions*> vecCollisions;

		string collisions_type;
	    vector<string> sg1, sg2;
	    vector<unsigned int> sgroup1, sgroup2;
	    double clog;
	    bool intra, debye_length_required = false;
	    int debug_every;
	    ostringstream mystream;


	    // Loop over each binary collisions group and parse info
	    unsigned int numcollisions=ifile.nComponents("Collisions");
	    for (unsigned int n_collisions = 0; n_collisions < numcollisions; n_collisions++) {

			ifile.extract("collisions_type",collisions_type,"Collisions",n_collisions);
			if(params.geometry == "1d3v" && collisions_type == "coulomb"){
				MESSAGE("Parameters for collisions #" << n_collisions << " :");

		        // Read the input file by searching for the keywords "species1" and "species2"
		        // which are the names of the two species that will collide
		        sg1.resize(0);
		        sg2.resize(0);
		        ifile.extract("species1",sg1,"Collisions",n_collisions);
		        ifile.extract("species2",sg2,"Collisions",n_collisions);

		        // Obtain the lists of species numbers from the lists of species names.
		        sgroup1 = params.FindSpecies(sg1);
		        sgroup2 = params.FindSpecies(sg2);

		        // Each group of species sgroup1 and sgroup2 must not be empty
		        if (sgroup1.size()==0) ERROR("No valid `species1` requested in collisions #" << n_collisions);
		        if (sgroup2.size()==0) ERROR("No valid `species2` requested in collisions #" << n_collisions);

		        // sgroup1 and sgroup2 can be equal, but cannot have common species if they are not equal
		        if (sgroup1 != sgroup2) {
		            for (unsigned int i1=0; i1<sgroup1.size(); i1++) {
		                for (unsigned int i2=0; i2<sgroup2.size(); i2++) {
		                    if (sgroup1[i1] == sgroup2[i2])
		                        ERROR("Unauthorized species (#" << sgroup1[i1]
		                              << ") in collisions #" << n_collisions
		                              << " (inter-collisions must not have a species colliding with itself)");
		                }
		            }
		            intra = false;
		        } else {
		            intra = true;
		        }

		        // Coulomb logarithm (if negative or unset, then automatically computed)
		        clog = 0.; // default
		        ifile.extract("coulomb_log",clog,"Collisions",n_collisions);
		        if (clog <= 0.) debye_length_required = true; // auto coulomb log requires debye length

		        // Number of timesteps between each debug output (if 0 or unset, no debug)
		        debug_every = 0; // default
		        ifile.extract("debug_every",debug_every,"Collisions",n_collisions);

		        // Print collisions parameters
		        mystream.str(""); // clear
		        for (unsigned int rs=0 ; rs<sgroup1.size() ; rs++) mystream << " #" << sgroup1[rs];
		        MESSAGE(1,"First  group of species :" << mystream.str());
		        mystream.str(""); // clear
		        for (unsigned int rs=0 ; rs<sgroup2.size() ; rs++) mystream << " #" << sgroup2[rs];
		        MESSAGE(1,"Second group of species :" << mystream.str());
		        MESSAGE(1,"Coulomb logarithm       : " << clog);
		        MESSAGE(1,"Intra collisions        : " << (intra?"True":"False"));
		        mystream.str(""); // clear
		        mystream << "Every " << debug_every << " timesteps";
		        MESSAGE(1,"Debug                   : " << (debug_every<=0?"No debug":mystream.str()));

		        // Add new Collisions objects to vector
		        vecCollisions.push_back( new Collisions1D_Coulomb(params,vecSpecies,smpi,n_collisions,sgroup1,sgroup2,clog,intra,debug_every) );


			}

			else if(collisions_type == "ionization"){

				// Add new Collisions objects to vector
				//> Three species participate in the ionization collision
		        vecCollisions.push_back( new Collisions1D_Ionization(params,vecSpecies,smpi,n_collisions,sgroup1,sgroup2,clog,intra,debug_every) );

			}
			else if(collisions_type == "DSMC"){
				// Add new Collisions objects to vector
				//> Only one species group participate in the ionization collision, all the particles from
				//> different species are the same
		        vecCollisions.push_back( new Collisions1D_DSMC(params,vecSpecies,smpi,n_collisions,sgroup1,sgroup2,clog,intra,debug_every) );

			}
			else if(collisions_type == "ChargeExchange"){
				// Add new Collisions objects to vector
				//> Only one species group participate in the ionization collision, all the particles from
				//> different species are the same
		        vecCollisions.push_back( new Collisions1D_ChargeExchange(params,vecSpecies,smpi,n_collisions,sgroup1,sgroup2,clog,intra,debug_every) );

			}
			else if(collisions_type == "Elastic"){
				// Add new Collisions objects to vector
				//> Only one species group participate in the ionization collision, all the particles from
				//> different species are the same
		        vecCollisions.push_back( new Collisions1D_Elastic(params,vecSpecies,smpi,n_collisions,sgroup1,sgroup2,clog,intra,debug_every) );

			}
			else {
				//ERROR("no collisions_type match: "<<collisions_type);
			}
	    }

	    // Needs wavelength_SI to be defined
	    if (numcollisions > 0)
	        if (params.wavelength_SI <= 0.)
	            ERROR("The parameter `wavelength_SI` needs to be defined and positive in order to compute collisions");

	    // pass the variable "debye_length_required" into the Collision class

		MESSAGE(1,"aaaaaaa ");
	    return vecCollisions;

	};

};

#endif
