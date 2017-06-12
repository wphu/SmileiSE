#ifndef PARTSOURCEFACTORY_H
#define PARTSOURCEFACTORY_H

#include "PartSource1D_Emit.h"
#include "PartSource1D_Load.h"
#include "PartSource2D_Load.h"
#include "InputData.h"
//#include "PartSource2D_Emit.h"

#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <sstream>
using namespace std;



class PartSourceFactory {

public:

	// Reads the input file and creates the PSI objects accordingly
	static vector<PartSource*> create(PicParams& params, InputData &ifile, vector<Species*>& vecSpecies, SmileiMPI* smpi)
	{
	    vector<PartSource*> vecPartSource;

		string PartSource_type;
		string emitKind;
	    vector<string> sg1, sg2;
	    vector<unsigned int> sgroup1, sgroup2;
		vector<double> mean_velocity;
	    string emitPos;
		double emitTemp;
		double weight_const;
		double emitOffset;
		double a_FN;
		double b_FN;
		double work_function;
		double emitJ;
		double emitFlux;
		int emitNumber;
		unsigned int nPartEmit;
		string relSpecies;

		string loadKind;
		int loadNumber;
		int everyTime;
		double loadDensity;
	    double loadTemperature;
		double loadDn;
	    double loadPos_start;
	    double loadPos_end;
		double loadPos_Ystart;
	    double loadPos_Yend;

	    bool intra, debye_length_required = false;
	    int debug_every;
	    ostringstream mystream;

		mean_velocity.resize(3, 0.0);


	    // Loop over each binary PartSource group and parse info
	    unsigned int numPartSource=ifile.nComponents("PartSource");
	    for (unsigned int n_PartSource = 0; n_PartSource < numPartSource; n_PartSource++) {

			ifile.extract("PartSource_type",PartSource_type,"PartSource",n_PartSource);
			if(params.geometry == "1d3v" && PartSource_type == "Emit")
			{
		        ifile.extract("species1",sg1,"PartSource",n_PartSource);
		        // Obtain the lists of species numbers from the lists of species names.
		        sgroup1 = params.FindSpecies(sg1);
		        // Each group of species sgroup1 and sgroup2 must not be empty
		        if (sgroup1.size()==0) ERROR("No valid `species1` requested in PSI #" << n_PartSource);


				weight_const = vecSpecies[sgroup1[0]]->species_param.weight; // default

				// Injection logarithm (if negative or unset, then automatically computed)
		        emitKind = "regular"; // default
		        ifile.extract("emitKind",emitKind,"PartSource",n_PartSource);

		        // Injection logarithm (if negative or unset, then automatically computed)
		        emitPos = "left"; // default
		        ifile.extract("emitPos",emitPos,"PartSource",n_PartSource);

		        // Number of timesteps between each debug output (if 0 or unset, no debug)
		        emitTemp = 0.0; // default
		        ifile.extract("emitTemp",emitTemp,"PartSource",n_PartSource);

				// emitOffset
		        ifile.extract("emitOffset",emitOffset,"PartSource",n_PartSource);

				// Number of timesteps between each debug output (if 0 or unset, no debug)
		        ifile.extract("a_FN",a_FN,"PartSource",n_PartSource);

				// Number of timesteps between each debug output (if 0 or unset, no debug)
		        ifile.extract("b_FN",b_FN,"PartSource",n_PartSource);

				// Number of timesteps between each debug output (if 0 or unset, no debug)
		        ifile.extract("work_function",work_function,"PartSource",n_PartSource);

				// Number of timesteps between each debug output (if 0 or unset, no debug)
				emitJ = 0;		// default
		        if( !ifile.extract("emitJ",emitJ,"PartSource",n_PartSource) );

				emitFlux = 0.0;		// default
				if( !ifile.extract("emitFlux",emitFlux,"PartSource",n_PartSource) );

				// Number of timesteps between each debug output (if 0 or unset, no debug)
				emitNumber = 0;		// default
		        if( !ifile.extract("emitNumber",emitNumber,"PartSource",n_PartSource) );

				// Number of timesteps between each debug output (if 0 or unset, no debug)
				relSpecies = "";		// default
		        ifile.extract("relSpecies",relSpecies,"PartSource",n_PartSource);

		        // Print PSI parameters
		        mystream.str(""); // clear
		        for (unsigned int rs=0 ; rs<sgroup1.size() ; rs++) mystream << " #" << sgroup1[rs];
		        MESSAGE(1,"First  group of species :" << mystream.str());

		        // Add new PSI objects to vector
		        vecPartSource.push_back( new PartSource1D_Emit(params, smpi, emitKind, sgroup1[0], emitPos, emitNumber,
								  emitTemp, emitJ, emitFlux, emitOffset, a_FN, b_FN, work_function, relSpecies) );

			}

			else if(params.geometry == "1d3v" && PartSource_type == "Load")
			{
		        ifile.extract("species1",sg1,"PartSource",n_PartSource);
		        // Obtain the lists of species numbers from the lists of species names.
		        sgroup1 = params.FindSpecies(sg1);
		        // Each group of species sgroup1 and sgroup2 must not be empty
		        if (sgroup1.size()==0) ERROR("No valid `species1` requested in PSI #" << n_PartSource);

				ifile.extract("mean_velocity",mean_velocity,"PartSource",n_PartSource);

		        ifile.extract("loadKind",loadKind,"PartSource",n_PartSource);
				ifile.extract("loadNumber",loadNumber,"PartSource",n_PartSource);

		        loadDensity = 0.0; // default
		        ifile.extract("loadDensity",loadDensity,"PartSource",n_PartSource);

				loadTemperature = 0.0; // default
		        ifile.extract("loadTemperature",loadTemperature,"PartSource",n_PartSource);

				loadDn= 0.0; // default
		        ifile.extract("loadDn",loadDn,"PartSource",n_PartSource);

				loadPos_start = 0.0; // default
				ifile.extract("loadPos_start",loadPos_start,"PartSource",n_PartSource);

				loadPos_end = 0.0; // default
				ifile.extract("loadPos_end",loadPos_end,"PartSource",n_PartSource);

		        // Print PSI parameters
		        mystream.str(""); // clear
		        for (unsigned int rs=0 ; rs<sgroup1.size() ; rs++) mystream << " #" << sgroup1[rs];
		        MESSAGE(1,"First  group of species :" << mystream.str());

		        // Add new PSI objects to vector
		        //vecPartSource.push_back( new PartSource1D_Load(params, smpi, sgroup1[0], loadDensity, loadTemperature, loadPos_start, loadPos_end) );
				vecPartSource.push_back( new PartSource1D_Load(params, smpi, sgroup1[0], mean_velocity, loadKind, loadNumber, loadDn,loadDensity,
				loadTemperature, loadPos_start, loadPos_end) );

			}

			else if(params.geometry == "2d3v" && PartSource_type == "Load")
			{
		        ifile.extract("species1",sg1,"PartSource",n_PartSource);
		        // Obtain the lists of species numbers from the lists of species names.
		        sgroup1 = params.FindSpecies(sg1);
		        // Each group of species sgroup1 and sgroup2 must not be empty
		        if (sgroup1.size()==0) ERROR("No valid `species1` requested in PSI #" << n_PartSource);

				ifile.extract("mean_velocity",mean_velocity,"PartSource",n_PartSource);

		        ifile.extract("loadKind",loadKind,"PartSource",n_PartSource);
				ifile.extract("loadNumber",loadNumber,"PartSource",n_PartSource);
				ifile.extract("everyTime",everyTime,"PartSource",n_PartSource);

		        loadDensity = 0.0; // default
		        ifile.extract("loadDensity",loadDensity,"PartSource",n_PartSource);

				loadTemperature = 0.0; // default
		        ifile.extract("loadTemperature",loadTemperature,"PartSource",n_PartSource);

				loadDn= 0.0; // default
		        ifile.extract("loadDn",loadDn,"PartSource",n_PartSource);

				loadPos_start = 0.0; // default
				ifile.extract("loadPos_start",loadPos_start,"PartSource",n_PartSource);

				loadPos_end = 0.0; // default
				ifile.extract("loadPos_end",loadPos_end,"PartSource",n_PartSource);

				loadPos_Ystart = 0.0; // default
				ifile.extract("loadPos_Ystart",loadPos_Ystart,"PartSource",n_PartSource);

				loadPos_Yend = 0.0; // default
				ifile.extract("loadPos_Yend",loadPos_Yend,"PartSource",n_PartSource);

		        // Print PSI parameters
		        mystream.str(""); // clear
		        for (unsigned int rs=0 ; rs<sgroup1.size() ; rs++) mystream << " #" << sgroup1[rs];
		        MESSAGE(1,"First  group of species :" << mystream.str());

		        // Add new PSI objects to vector
		        //vecPartSource.push_back( new PartSource1D_Load(params, smpi, sgroup1[0], loadDensity, loadTemperature, loadPos_start, loadPos_end) );
				vecPartSource.push_back( new PartSource2D_Load(params, smpi, sgroup1[0], mean_velocity, loadKind, loadNumber, everyTime, loadDn,loadDensity,
				loadTemperature, loadPos_start, loadPos_end, loadPos_Ystart, loadPos_Yend) );

			}

			else {
				ERROR("no PartSource_type match: "<<PartSource_type);
			}
	    }


		// set relevant PSI for some Injection PSI

		for (unsigned int n_PartSource = 0; n_PartSource < numPartSource; n_PartSource++)
		{
			if(vecPartSource[n_PartSource]->emitKind == "relEmit"){
				for (unsigned int m_PartSource = 0; m_PartSource < numPartSource; n_PartSource++)
				{
					unsigned int ispec = vecPartSource[m_PartSource]->species1;
					if(vecSpecies[ispec]->species_param.species_type == vecPartSource[n_PartSource]->relSpecies){
						vecPartSource[n_PartSource]->setRelPartSource(vecPartSource[m_PartSource]);
					}
				}
			}
		}
	    return vecPartSource;
	};

};

#endif
