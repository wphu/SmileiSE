#ifndef PSIFACTORY_H
#define PSIFACTORY_H

#include "PSI1D_Injection.h"
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
			if(params.geometry == "1d3v" && PSI_type == "Injection"){
				MESSAGE("Parameters for PSI #" << n_PSI << " :");


		        ifile.extract("species1",sg1,"PSI",n_PSI);
		        // Obtain the lists of species numbers from the lists of species names.
		        sgroup1 = params.FindSpecies(sg1);
		        // Each group of species sgroup1 and sgroup2 must not be empty
		        if (sgroup1.size()==0) ERROR("No valid `species1` requested in PSI #" << n_PSI);


				weight_const = vecSpecies[sgroup1[0]]->species_param.weight; // default

				// Injection logarithm (if negative or unset, then automatically computed)
		        emitKind = "regular"; // default
		        ifile.extract("emitKind",emitKind,"PSI",n_PSI);

		        // Injection logarithm (if negative or unset, then automatically computed)
		        psiPos = "left"; // default
		        ifile.extract("psiPos",psiPos,"PSI",n_PSI);

		        // Number of timesteps between each debug output (if 0 or unset, no debug)
		        emitTemp = 0.0; // default
		        ifile.extract("emitTemp",emitTemp,"PSI",n_PSI);

				// emitOffset
		        ifile.extract("emitOffset",emitOffset,"PSI",n_PSI);

				// Number of timesteps between each debug output (if 0 or unset, no debug)
		        ifile.extract("a_FN",a_FN,"PSI",n_PSI);

				// Number of timesteps between each debug output (if 0 or unset, no debug)
		        ifile.extract("b_FN",b_FN,"PSI",n_PSI);

				// Number of timesteps between each debug output (if 0 or unset, no debug)
		        ifile.extract("work_function",work_function,"PSI",n_PSI);

				// Number of timesteps between each debug output (if 0 or unset, no debug)
				emitJ = 0;		// default
		        if( !ifile.extract("emitJ",emitJ,"PSI",n_PSI) );

				emitFlux = 0;		// default
		        if( !ifile.extract("emitFlux",emitFlux,"PSI",n_PSI) );

				// Number of timesteps between each debug output (if 0 or unset, no debug)
				nPartEmit = 0;		// default
		        if( !ifile.extract("nPartEmit",nPartEmit,"PSI",n_PSI) );

				// Number of timesteps between each debug output (if 0 or unset, no debug)
				relSpecies = "";		// default
		        ifile.extract("relSpecies",relSpecies,"PSI",n_PSI);

		        // Print PSI parameters
		        mystream.str(""); // clear
		        for (unsigned int rs=0 ; rs<sgroup1.size() ; rs++) mystream << " #" << sgroup1[rs];
		        MESSAGE(1,"First  group of species :" << mystream.str());

		        // Add new PSI objects to vector
		        vecPSI.push_back( new PSI1D_Injection(params, smpi, emitKind, sgroup1[0], psiPos, nPartEmit,
								  emitTemp, emitJ, emitFlux, weight_const, emitOffset, a_FN, b_FN, work_function, relSpecies) );

			}

			else if(params.geometry == "1d3v" && PSI_type == "sputtering"){
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


		// set relevant PSI for some Injection PSI

		for (unsigned int n_PSI = 0; n_PSI < numPSI; n_PSI++)
		{
			if(vecPSI[n_PSI]->emitKind == "relEmit"){
				for (unsigned int m_PSI = 0; m_PSI < numPSI; n_PSI++)
				{
					unsigned int ispec = vecPSI[m_PSI]->species1;
					if(vecSpecies[ispec]->species_param.species_type == vecPSI[n_PSI]->relSpecies){
						vecPSI[n_PSI]->setRelPsi(vecPSI[m_PSI]);
					}
				}
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
