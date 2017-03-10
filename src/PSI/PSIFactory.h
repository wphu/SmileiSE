#ifndef PSIFACTORY_H
#define PSIFACTORY_H

#include "PSI1D_SEE.h"
#include "PSI1D_Sputtering.h"
#include "PSI1D_Backscattering.h"
#include "PSI1D_Recycling.h"

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
		string emitKind;
	    vector<string> sg1, sg2;
	    vector<unsigned int> sgroup1, sgroup2;
	    string psiPos;
		double emitTemp;
		double weight_const;
		double emitOffset;
		double a_FN;
		double b_FN;
		double work_function;
		double emitJ;
		double emitFlux;
		unsigned int nPartEmit;
		double recycling_factor;
		string relSpecies;


	    bool intra, debye_length_required = false;
	    int debug_every;
	    ostringstream mystream;


	    // Loop over each binary PSI group and parse info
	    unsigned int numPSI=ifile.nComponents("PSI");
	    for (unsigned int n_PSI = 0; n_PSI < numPSI; n_PSI++) {

			ifile.extract("PSI_type",PSI_type,"PSI",n_PSI);
			if(params.geometry == "1d3v" && PSI_type == "sputtering"){
				MESSAGE("Parameters for PSI #" << n_PSI << " :");


		        ifile.extract("species1",sg1,"PSI",n_PSI);
				ifile.extract("species2",sg2,"PSI",n_PSI);
		        // Obtain the lists of species numbers from the lists of species names.
		        sgroup1 = params.FindSpecies(sg1);
				sgroup2 = params.FindSpecies(sg2);
		        // Each group of species sgroup1 and sgroup2 must not be empty
		        if (sgroup1.size()==0) ERROR("No valid `species1` requested in PSI #" << n_PSI);
				if (sgroup2.size()==0) ERROR("No valid `species2` requested in PSI #" << n_PSI);

				// Injection logarithm (if negative or unset, then automatically computed)
		        psiPos = "left"; // default
		        ifile.extract("psiPos",psiPos,"PSI",n_PSI);

		        // Number of timesteps between each debug output (if 0 or unset, no debug)
		        emitTemp = 0.0; // default
		        ifile.extract("emitTemp",emitTemp,"PSI",n_PSI);

		        vecPSI.push_back( new PSI1D_Sputtering(params, smpi,vecSpecies, sgroup1[0],sgroup2[0], psiPos, emitTemp) );

			}
			else if(params.geometry == "1d3v" && PSI_type == "BackScattering"){
				MESSAGE("Parameters for PSI #" << n_PSI << " :");


		        ifile.extract("species1",sg1,"PSI",n_PSI);
				ifile.extract("species2",sg2,"PSI",n_PSI);
		        // Obtain the lists of species numbers from the lists of species names.
		        sgroup1 = params.FindSpecies(sg1);
				sgroup2 = params.FindSpecies(sg2);
		        // Each group of species sgroup1 and sgroup2 must not be empty
		        if (sgroup1.size()==0) ERROR("No valid `species1` requested in PSI #" << n_PSI);
				if (sgroup2.size()==0) ERROR("No valid `species2` requested in PSI #" << n_PSI);

				// Injection logarithm (if negative or unset, then automatically computed)
		        psiPos = "left"; // default
		        ifile.extract("psiPos",psiPos,"PSI",n_PSI);

		        // Number of timesteps between each debug output (if 0 or unset, no debug)
		        emitTemp = 0.0; // default
		        ifile.extract("emitTemp",emitTemp,"PSI",n_PSI);

		        vecPSI.push_back( new PSI1D_Backscattering(params, smpi,sgroup1[0],sgroup2[0], psiPos, emitTemp) );
			}

			else if(params.geometry == "1d3v" && PSI_type == "Recycling"){
				MESSAGE("Parameters for PSI #" << n_PSI << " :");


		        ifile.extract("species1",sg1,"PSI",n_PSI);
				ifile.extract("species2",sg2,"PSI",n_PSI);
		        // Obtain the lists of species numbers from the lists of species names.
		        sgroup1 = params.FindSpecies(sg1);
		        // Each group of species sgroup1 and sgroup2 must not be empty
		        if (sgroup1.size()==0) ERROR("No valid `species1` requested in PSI #" << n_PSI);

				// Injection logarithm (if negative or unset, then automatically computed)
		        psiPos = "left"; // default
		        ifile.extract("psiPos",psiPos,"PSI",n_PSI);

		        // Number of timesteps between each debug output (if 0 or unset, no debug)
		        emitTemp = 0.0; // default
		        ifile.extract("emitTemp",emitTemp,"PSI",n_PSI);

				recycling_factor = 0.0; // default
				ifile.extract("recycling_factor",recycling_factor,"PSI",n_PSI);

		        vecPSI.push_back( new PSI1D_Recycling(params, smpi,sgroup1[0], psiPos, emitTemp, recycling_factor) );
			}

			else if(params.geometry == "1d3v" && PSI_type == "SEE"){
				MESSAGE("Parameters for PSI #" << n_PSI << " :");


		        ifile.extract("species1",sg1,"PSI",n_PSI);
				ifile.extract("species2",sg2,"PSI",n_PSI);
		        // Obtain the lists of species numbers from the lists of species names.
		        sgroup1 = params.FindSpecies(sg1);
				sgroup2 = params.FindSpecies(sg2);
		        // Each group of species sgroup1 and sgroup2 must not be empty
		        if (sgroup1.size()==0) ERROR("No valid `species1` requested in PSI #" << n_PSI);
				if (sgroup2.size()==0) ERROR("No valid `species2` requested in PSI #" << n_PSI);

				// Injection logarithm (if negative or unset, then automatically computed)
		        psiPos = "left"; // default
		        ifile.extract("psiPos",psiPos,"PSI",n_PSI);

		        // Number of timesteps between each debug output (if 0 or unset, no debug)
		        emitTemp = 0.0; // default
		        ifile.extract("emitTemp",emitTemp,"PSI",n_PSI);

				// secondary electron emission yield
		        double SEEYield = 0.0; // default
		        ifile.extract("SEEYield",SEEYield,"PSI",n_PSI);

		        vecPSI.push_back( new PSI1D_SEE(params, smpi,sgroup1[0],sgroup2[0], psiPos, emitTemp, SEEYield) );

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
