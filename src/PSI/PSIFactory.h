#ifndef PSIFACTORY_H
#define PSIFACTORY_H

#include "PSI1D_Injection.h"
#include "PSI1D_SEE.h"
#include "PSI1D_Sputtering.h"

#include "PSI2D_Injection.h"
#include "PSI2D_SEE.h"
#include "PSI2D_Sputtering.h"

#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <sstream>
using namespace std;

class PSIFactory {

public:

	// Reads the input file and creates the PSI objects accordingly
	static vector<PSI*> create(PicParams& params, InputData &ifile, vector<Species*>& vecSpecies, SmileiMPI* smpi)
	{
	    vector<PSI*> vecPSI;

		string PSI_type;
	    vector<string> sg1, sg2;
	    vector<unsigned int> sgroup1, sgroup2;
	    double clog;
	    bool intra, debye_length_required = false;
	    int debug_every;
	    ostringstream mystream;


	    // Loop over each binary PSI group and parse info
	    unsigned int numPSI=ifile.nComponents("PSI");
	    for (unsigned int n_PSI = 0; n_PSI < numPSI; n_PSI++) {

			ifile.extract("PSI_type",PSI_type,"PSI",n_PSI);
			if(PSI_type == "Injection"){
				MESSAGE("Parameters for PSI #" << n_PSI << " :");

		        // Read the input file by searching for the keywords "species1" and "species2"
		        // which are the names of the two species that will collide
		        sg1.resize(0);
		        sg2.resize(0);
		        ifile.extract("species1",sg1,"PSI",n_PSI);
		        ifile.extract("species2",sg2,"PSI",n_PSI);

		        // Obtain the lists of species numbers from the lists of species names.
		        sgroup1 = params.FindSpecies(sg1);
		        sgroup2 = params.FindSpecies(sg2);

		        // Each group of species sgroup1 and sgroup2 must not be empty
		        if (sgroup1.size()==0) ERROR("No valid `species1` requested in PSI #" << n_PSI);
		        if (sgroup2.size()==0) ERROR("No valid `species2` requested in PSI #" << n_PSI);

		        // sgroup1 and sgroup2 can be equal, but cannot have common species if they are not equal
		        if (sgroup1 != sgroup2) {
		            for (unsigned int i1=0; i1<sgroup1.size(); i1++) {
		                for (unsigned int i2=0; i2<sgroup2.size(); i2++) {
		                    if (sgroup1[i1] == sgroup2[i2])
		                        ERROR("Unauthorized species (#" << sgroup1[i1]
		                              << ") in PSI #" << n_PSI
		                              << " (inter-PSI must not have a species colliding with itself)");
		                }
		            }
		            intra = false;
		        } else {
		            intra = true;
		        }

		        // Injection logarithm (if negative or unset, then automatically computed)
		        clog = 0.; // default
		        ifile.extract("coulomb_log",clog,"PSI",n_PSI);
		        if (clog <= 0.) debye_length_required = true; // auto coulomb log requires debye length

		        // Number of timesteps between each debug output (if 0 or unset, no debug)
		        debug_every = 0; // default
		        ifile.extract("debug_every",debug_every,"PSI",n_PSI);

		        // Print PSI parameters
		        mystream.str(""); // clear
		        for (unsigned int rs=0 ; rs<sgroup1.size() ; rs++) mystream << " #" << sgroup1[rs];
		        MESSAGE(1,"First  group of species :" << mystream.str());
		        mystream.str(""); // clear
		        for (unsigned int rs=0 ; rs<sgroup2.size() ; rs++) mystream << " #" << sgroup2[rs];
		        MESSAGE(1,"Second group of species :" << mystream.str());
		        MESSAGE(1,"Injection logarithm       : " << clog);
		        MESSAGE(1,"Intra PSI        : " << (intra?"True":"False"));
		        mystream.str(""); // clear
		        mystream << "Every " << debug_every << " timesteps";
		        MESSAGE(1,"Debug                   : " << (debug_every<=0?"No debug":mystream.str()));

		        // Add new PSI objects to vector
		        vecPSI.push_back( new PSI1D_Injection() );


			}

			else if(PSI_type == "sputtering"){

				// Add new PSI objects to vector
				//> Three species participate in the ionization collision
		        vecPSI.push_back( new PSI1D_Sputtering() );

			}
			else if(PSI_type == "SEE"){
				// Add new PSI objects to vector
				//> Only one species group participate in the ionization collision, all the particles from
				//> different species are the same
		        vecPSI.push_back( new PSI1D_SEE() );

			}

			else {
				ERROR("no PSI_type match: "<<PSI_type);
			}
	    }

	    // Needs wavelength_SI to be defined
	    if (numPSI > 0)
	        if (params.wavelength_SI <= 0.)
	            ERROR("The parameter `wavelength_SI` needs to be defined and positive in order to compute PSI");


	    return vecPSI;
	};

};

#endif
