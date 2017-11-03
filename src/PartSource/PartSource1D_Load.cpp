#include "PartSource1D_Load.h"
#include "SmileiMPI_Cart1D.h"
#include "Field1D.h"
#include "H5.h"


#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <sstream>

using namespace std;


// Constructor
PartSource1D_Load::PartSource1D_Load(
    PicParams&      params,
    SmileiMPI*      smpi,
    unsigned int    load_species1,
    vector<unsigned int>     load_species_dependent,
    vector<double>  mean_vel,
    string          load_kind,
    int             load_number,
    double          load_dn,
    double          load_density,
    double          load_temperature,
    vector<int>     load_timeStep_vector,
    vector<double>  load_temperature_vector,
    double          load_temperature_unLimit_factor,
    vector<double>  load_dn_vector,
    double          load_Pos_start,
    double          load_Pos_end
):
PartSource1D (params, smpi)

{
    dt_ov_dx = params.timestep / params.cell_length[0];
    dt = params.timestep;

    species1        = load_species1;
    mean_velocity   = mean_vel;
    loadKind        = load_kind;
    loadNumber      = load_number;
    loadDn          = load_dn;
    loadDensity     = load_density;
    loadTemperature = load_temperature;
    loadPos_start   = load_Pos_start;
    loadPos_end     = load_Pos_end;

    loadTimeStepVector      = load_timeStep_vector;
    loadTemperatureVector   = load_temperature_vector;
    loadTemperature_upLimit_factor = load_temperature_unLimit_factor;
    loadDnVector            = load_dn_vector;
    species_group_dependent = load_species_dependent;

    YZArea = 1.0;
    loadq = loadDn * loadTemperature;

    // load particles by number_density and Temperature
    if(loadKind == "nT") {
        for(int ibin = 0; ibin < numPart_in_each_bin.size(); ibin++)
        {
            numPart_in_each_bin[ibin] = loadDensity / params.species_param[species1].weight;
            //cout<<"numPart_in_each_bin "<<numPart_in_each_bin[ibin]<<endl;
        }
    }
    else if(loadKind == "dn")
    {
        if(loadTimeStepVector.size() != 0)
        {
            loadTemperature = loadTemperatureVector[0];
            loadDn = loadDnVector[0];
        }

        loadStep = 1.0 + loadNumber * params.species_param[species1].weight / (loadDn * params.timestep);
        loadRem = loadDn * loadStep * params.timestep / params.species_param[species1].weight - loadNumber;
        loadRemTot = 0.0;
        MESSAGE("loadStep = "<<loadStep);
    }

    else if(loadKind == "nq")
    {
        loadStep = 1.0 + loadNumber * params.species_param[species1].weight / (loadDn * params.timestep);
        loadNumber_init = loadDn * loadStep * params.timestep / params.species_param[species1].weight;
        loadRemTot = 0.0;
		loadNumber_heat = loadNumber_init;
		loadNumber = loadNumber_heat;
        MESSAGE("loadStep = "<<loadStep);

        int mpiSize = smpi->getSize();
        if(mpiSize % 2 == 0)
        {
            mpiRank_source_middle = mpiSize / 2;
            index_source_middle = 3;
        }
        else
        {
            mpiRank_source_middle = mpiSize / 2;
            index_source_middle = params.n_space[0] / 2;
        }
        // get the middle position of source region
    }



    // the MPI domain is not in the source region
    if(smpi->getDomainLocalMax(0) <= loadPos_start || smpi->getDomainLocalMin(0) >= loadPos_end) {
        loadPos_start = 0.0;
        loadBin_start = 0;
        loadPos_end = 0.0;
        loadBin_end = 0;
    }
    // the left end of source region is in the MPI domain
    else if(smpi->getDomainLocalMin(0) <= loadPos_start && smpi->getDomainLocalMax(0) > loadPos_start
    && smpi->getDomainLocalMax(0) <= loadPos_end) {
        double temp = (loadPos_start - smpi->getDomainLocalMin(0)) / params.cell_length[0];
        loadBin_start = temp + 0.001;
        if(loadBin_start < 0) { loadBin_start = 0; }
        loadPos_start = smpi->getDomainLocalMin(0) + loadBin_start*params.cell_length[0];
        //numPart_in_each_bin[loadBin_start] /=
        //( params.cell_length[0] / (( smpi->getDomainLocalMin(0) + (loadBin_start+1)*params.cell_length[0] ) - loadPos_start) );
        loadPos_end = smpi->getDomainLocalMax(0);
        loadBin_end = params.n_space[0] - 1;
    }
    // the right end of source region is in the MPI domain
    else if(smpi->getDomainLocalMin(0) < loadPos_end && smpi->getDomainLocalMax(0) >= loadPos_end
    && smpi->getDomainLocalMin(0) >= loadPos_start) {
        loadPos_start = smpi->getDomainLocalMin(0);
        loadBin_start = 0;
        loadBin_end = (loadPos_end - smpi->getDomainLocalMin(0)) / params.cell_length[0];
        if(loadBin_end > params.n_space[0] - 1) { loadBin_end = params.n_space[0] - 1; }
        //numPart_in_each_bin[loadBin_end] /=
        //( params.cell_length[0] / (loadPos_end - ( smpi->getDomainLocalMin(0) + loadBin_end*params.cell_length[0] )) );
        loadPos_end = smpi->getDomainLocalMin(0) + (loadBin_end+1)*params.cell_length[0];
    }
    // the whole MPI domain is in the source region
    else if(smpi->getDomainLocalMin(0) >= loadPos_start && smpi->getDomainLocalMax(0) <= loadPos_end) {
        loadPos_start = smpi->getDomainLocalMin(0);
        loadBin_start = 0;
        loadPos_end = smpi->getDomainLocalMax(0);
        loadBin_end = params.n_space[0] - 1;
    }
    // the whole MPI domain source region is in the bin
    else if(smpi->getDomainLocalMin(0) < loadPos_start && smpi->getDomainLocalMax(0) > loadPos_end) {
        //loadBin_start = (loadPos_start - smpi->getDomainLocalMin(0)) / params.cell_length[0];
        //numPart_in_each_bin[loadBin_start] /=
        //( params.cell_length[0] / (( smpi->getDomainLocalMin(0) + (loadBin_start+1)*params.cell_length[0] ) - loadPos_start) );
        //loadBin_end = (loadPos_end - smpi->getDomainLocalMin(0)) / params.cell_length[0];
        //numPart_in_each_bin[loadBin_end] /=
        //( params.cell_length[0] / (loadPos_end - ( smpi->getDomainLocalMin(0) + loadBin_start*params.cell_length[0] )) );
        double temp = (loadPos_start - smpi->getDomainLocalMin(0)) / params.cell_length[0];
        loadBin_start = temp + 0.001;
        if(loadBin_start < 0) { loadBin_start = 0; }
        loadPos_start = smpi->getDomainLocalMin(0) + loadBin_start*params.cell_length[0];
        loadBin_end = (loadPos_end - smpi->getDomainLocalMin(0)) / params.cell_length[0];
        if(loadBin_end > params.n_space[0] - 1) { loadBin_end = params.n_space[0] - 1; }
        loadPos_end = smpi->getDomainLocalMin(0) + (loadBin_end+1)*params.cell_length[0];
    }

    if(loadBin_end != loadBin_start)
    {
        cout<<"MPI "<<smpi->getRank()<<endl;
        cout<<"LoadBin start and end: "<<loadBin_start<<" "<<loadBin_end<<endl;
        cout<<"LoadPos start and end: "<<loadPos_start<<" "<<loadPos_end<<endl;
    }

    // Parameters for "nq"
    timeStep_checkFor_nq = 0;
    temperature_pre = loadTemperature;
    source_density_pre = 0.0;
    loadTemperature_init = loadTemperature;
    loadTemperature_exceed = 0.0;
	loadTemperature_heat = (loadTemperature_init * loadNumber_init - loadNumber * loadTemperature)
							/ (loadNumber_heat * loadStep);

    // Parameters for "dn"
    nextTimeStep = 0;
    nextTimeStep_index = 1;


}

PartSource1D_Load::~PartSource1D_Load()
{

}



// Calculates the PSI for a given PSI object
void PartSource1D_Load::emitLoad(PicParams& params, SmileiMPI* smpi, vector<Species*>& vecSpecies, int itime, ElectroMagn* fields)
{
    Species   *s1;
    Particles *p1;
    int iPart, nPart;
    vector <double> cell_length;            // cell_length for initialize position, maybe not equal to real cell lenghth
    vector<double> max_jutt_cumul;
    double *indexes=new double[params.nDim_particle];
    double *temp=new double[3];
    double *vel=new double[3];

    // Parameters for "nq"
    double source_density;
    double zoom_factor;
	double loadDn_temp;

    if(loadKind == "nq")
    {
        if(itime%loadStep == 0)
        {
            Field1D* rho1D = static_cast<Field1D*>(fields->rho_s[ species_group_dependent[0] ]);
            source_density = (*rho1D)(index_source_middle);
            smpi->bcast_double(&source_density, 1, mpiRank_source_middle);

			if(source_density < loadDensity && source_density < source_density_pre )
			{
				zoom_factor = (1.0 - 0.5*itime/params.n_time) * (source_density_pre - source_density) / loadDensity;
				loadDn *= (1.0 + zoom_factor);
				cout<<"loadDn111  "<<loadDn<<"  "<<zoom_factor<<endl;
			}
			else if(source_density > loadDensity && source_density > source_density_pre)
			{
				zoom_factor = (1.0 - 0.5*itime/params.n_time) * (source_density - source_density_pre) / loadDensity;
				loadDn_temp = loadDn;
				loadDn *= (1.0 - zoom_factor);
				cout<<"loadDn2222  "<<loadDn<<"  "<<zoom_factor<<endl;
			}


			
			if(loadDn <= 0.0)
			{
				cout<<"PartSource1D, loadDn: "<<loadDn<<endl;
				loadDn = loadDn_temp;
			}
            loadTemperature_exceed = loadq / loadDn;
            source_density_pre = source_density;
            if(loadTemperature_exceed < loadTemperature_init * loadTemperature_upLimit_factor)
            {
                loadTemperature = loadTemperature_exceed;
            }

            double loadNumber_temp = loadDn * loadStep * params.timestep / params.species_param[species1].weight;
            loadRemTot += loadNumber_temp;
			loadNumber = loadRemTot;
			loadRemTot -= loadNumber;
			if(loadNumber < loadNumber_init)
			{
				loadTemperature_heat = (loadTemperature_init * loadNumber_init - loadNumber * loadTemperature)
												 / loadNumber_heat;	
				if(loadTemperature_heat < 0.0)
				{
					loadTemperature_heat = 0.0;
					loadTemperature = loadTemperature_init * loadNumber_init / loadNumber;

				}
			}
			else
			{
				loadTemperature_heat = 0.0;
				loadTemperature = loadTemperature_init * loadNumber_init / loadNumber;
				
			}
			cout<<"loadTemperature_heat  "<<loadTemperature_heat<<endl;
			cout<<"loadTemperature_init  "<<loadTemperature_init<<endl;
			cout<<"loadNumber_init  "<<loadNumber_init<<endl;
			cout<<"loadNumber  "<<loadNumber<<endl;
			cout<<"loadTemperature  "<<loadTemperature<<endl;

            if(loadTemperature/temperature_pre > 2.0 || loadTemperature/temperature_pre < 0.5)
            {
                cout<<"Temperature pre and now: "<<temperature_pre<<" "<<loadTemperature<<endl;
                temperature_pre = loadTemperature;
            }
        }
    }


    if(loadTimeStepVector.size() > 1)
    {
        if(nextTimeStep == 0)
        {
            nextTimeStep = loadTimeStepVector[1];
        }
        if(itime > nextTimeStep)
        {
            loadTemperature = loadTemperatureVector[nextTimeStep_index];
            loadDn = loadDnVector[nextTimeStep_index];
            if(nextTimeStep_index < loadTimeStepVector.size() -1 )
            {
                cout<<"PartSource1D_Load: nextTimeStep_index "<<nextTimeStep_index<<endl;
                nextTimeStep_index++;
                nextTimeStep = loadTimeStepVector[nextTimeStep_index];
            }
            else
            {
                nextTimeStep *= 1000;
            }

            loadStep = 1.0 + loadNumber * params.species_param[species1].weight / (loadDn * params.timestep);
            loadRem = loadDn * loadStep * params.timestep / params.species_param[species1].weight - loadNumber;
        }

    }


    if(loadKind == "nT" && loadBin_end != loadBin_start)
    {
        s1 = vecSpecies[species1];
        p1 = &(s1->particles);

        cell_length.resize(params.nDim_particle);
        max_jutt_cumul.resize(0);
        temp[0] = loadTemperature;
        temp[1] = loadTemperature;
        temp[2] = loadTemperature;
        vel[0] = mean_velocity[0];
        vel[1] = mean_velocity[1];
        vel[2] = mean_velocity[2];

        indexes_of_particles_to_erase.clear();
        for(int ibin = 0; ibin < count_of_particles_to_insert.size(); ibin++ )
        {
            count_of_particles_to_insert[ibin] = 0;
        }

        // erase unnecessary paritcles in source region
        for(int ibin=loadBin_start; ibin<=loadBin_end; ibin++)
        {
            if( s1->bmax[ibin] - s1->bmin[ibin] > numPart_in_each_bin[ibin] ){
                for(int iPart = s1->bmin[ibin] + numPart_in_each_bin[ibin]; iPart < s1->bmax[ibin]; iPart++)
                {
                    indexes_of_particles_to_erase.push_back(iPart);
                }
                count_of_particles_to_insert[ibin] = 0;
            }
            else {
                count_of_particles_to_insert[ibin] = numPart_in_each_bin[ibin] - ( s1->bmax[ibin] - s1->bmin[ibin] );
            }
        }
        s1->erase_particles_from_bins(indexes_of_particles_to_erase);

        // insert lost particles into bins
        new_particles.clear();
        for(int ibin = loadBin_start; ibin <= loadBin_end; ibin++ )
        {
            new_particles.create_particles(count_of_particles_to_insert[ibin]);
        }
        s1->insert_particles_to_bins(new_particles, count_of_particles_to_insert);

        // re-initialize paritcles in source region
        for(int ibin=loadBin_start; ibin<=loadBin_end; ibin++)
        {
            iPart = s1->bmin[ibin];
            nPart = numPart_in_each_bin[ibin];
            if(nPart != s1->bmax[ibin] - s1->bmin[ibin]) {
                for(int i = 0; i < numPart_in_each_bin.size(); i++)
                {
                    cout<<"numPart_in_each_bin "<<numPart_in_each_bin[i]<<" "<<i<<endl;
                }
                cout<<"loadBin start and end: "<<loadBin_start<<" "<<loadBin_end<<" "<<params.n_space[0]<<endl;
                ERROR("nPart not equal to bmax - bmin!! "<<nPart<<" "<<s1->bmax[ibin] - s1->bmin[ibin]<<" "<<ibin);
            }
            cell_length[0] = params.cell_length[0];
            indexes[0] = smpi->getDomainLocalMin(0) + ibin*params.cell_length[0];

            s1->initPosition(nPart, iPart, indexes, params.nDim_particle,
                         cell_length, s1->species_param.initPosition_type);

            s1->initMomentum(nPart,iPart, temp, vel,
                         s1->species_param.initMomentum_type, max_jutt_cumul, params);

            s1->initWeight_constant(nPart, species1, iPart, s1->species_param.weight);
            s1->initCharge(nPart, species1, iPart, s1->species_param.charge);
        }
    }

    else if(loadKind == "dn" && itime%loadStep == 0 && loadBin_end != loadBin_start)
    {
        s1 = vecSpecies[species1];
        p1 = &(s1->particles);

        cell_length.resize(params.nDim_particle);
        max_jutt_cumul.resize(0);
        temp[0] = loadTemperature;
        temp[1] = loadTemperature;
        temp[2] = loadTemperature;
        vel[0] = mean_velocity[0];
        vel[1] = mean_velocity[1];
        vel[2] = mean_velocity[2];

        double loadNumber_temp;
        loadRemTot += loadRem;
        loadNumber_temp = loadNumber;
        if(loadRemTot > 1.0)
        {
            loadNumber_temp = loadNumber + 1;
            loadRemTot -= 1.0;
        }

        for(int ibin = 0; ibin < count_of_particles_to_insert.size(); ibin++ )
        {
            count_of_particles_to_insert[ibin] = 0;
            if(ibin >= loadBin_start && ibin <= loadBin_end)
            {

                count_of_particles_to_insert[ibin] = loadNumber_temp;
            }
        }
        //cout<<"number: "<<loadDensity * params.timestep / s1->species_param.weight<<endl;

        new_particles.clear();
        for(int ibin = loadBin_start; ibin <= loadBin_end; ibin++ )
        {
            new_particles.create_particles(count_of_particles_to_insert[ibin]);
        }
        s1->insert_particles_to_bins(new_particles, count_of_particles_to_insert);

        // re-initialize paritcles in source region
        for(int ibin=loadBin_start; ibin<=loadBin_end; ibin++)
        {
            iPart = s1->bmax[ibin] - count_of_particles_to_insert[ibin];
            nPart = count_of_particles_to_insert[ibin];
            cell_length[0] = params.cell_length[0];
            indexes[0] = smpi->getDomainLocalMin(0) + ibin*params.cell_length[0];

            s1->initPosition(nPart, iPart, indexes, params.nDim_particle,
                         cell_length, s1->species_param.initPosition_type);

            s1->initMomentum(nPart,iPart, temp, vel,
                         s1->species_param.initMomentum_type, max_jutt_cumul, params);

            s1->initWeight_constant(nPart, species1, iPart, s1->species_param.weight);
            s1->initCharge(nPart, species1, iPart, s1->species_param.charge);
        }

    }
    else if(loadKind == "nq" && loadBin_end != loadBin_start)
    {
        s1 = vecSpecies[species1];
        p1 = &(s1->particles);

    		if(itime%loadStep == 0)
    		{
    			cell_length.resize(params.nDim_particle);
    			max_jutt_cumul.resize(0);
    			temp[0] = loadTemperature;
    			temp[1] = loadTemperature;
    			temp[2] = loadTemperature;
    			vel[0] = mean_velocity[0];
    			vel[1] = mean_velocity[1];
    			vel[2] = mean_velocity[2];

    			for(int ibin = 0; ibin < count_of_particles_to_insert.size(); ibin++ )
    			{
    				count_of_particles_to_insert[ibin] = 0;
    				if(ibin >= loadBin_start && ibin <= loadBin_end)
    				{

    					count_of_particles_to_insert[ibin] = loadNumber;
    				}
    			}
    			//cout<<"number: "<<loadDensity * params.timestep / s1->species_param.weight<<endl;

    			new_particles.clear();
    			for(int ibin = loadBin_start; ibin <= loadBin_end; ibin++ )
    			{
    				new_particles.create_particles(count_of_particles_to_insert[ibin]);
    			}
    			s1->insert_particles_to_bins(new_particles, count_of_particles_to_insert);

    			// re-initialize paritcles in source region
    			for(int ibin=loadBin_start; ibin<=loadBin_end; ibin++)
    			{
    				iPart = s1->bmax[ibin] - count_of_particles_to_insert[ibin];
    				nPart = count_of_particles_to_insert[ibin];
    				cell_length[0] = params.cell_length[0];
    				indexes[0] = smpi->getDomainLocalMin(0) + ibin*params.cell_length[0];

    				s1->initPosition(nPart, iPart, indexes, params.nDim_particle,
    							 cell_length, s1->species_param.initPosition_type);

    				s1->initMomentum(nPart,iPart, temp, vel,
    							 s1->species_param.initMomentum_type, max_jutt_cumul, params);

    				s1->initWeight_constant(nPart, species1, iPart, s1->species_param.weight);
    				s1->initCharge(nPart, species1, iPart, s1->species_param.charge);
    			}
				
				// heat paritcles in source region
				for(int ibin=loadBin_start; ibin<=loadBin_end; ibin++)
				{
					iPart = s1->bmin[ibin];
					nPart = s1->bmax[ibin] - s1->bmin[ibin];
					if(nPart > 0)
					{
						//cout<<"heat  "<<loadTemperature_heat*loadNumber_heat/nPart<<endl;
						s1->heat(nPart, nPart, iPart, loadTemperature_heat*loadNumber_heat/nPart, params);
					}
					//else if(loadNumber_heat <= nPart && loadNumber_heat > 0)
					//{
					//	s1->heat(nPart, loadNumber_heat, iPart, loadTemperature_heat, params);
					//}
				}
				
				
    		}




    }
    delete [] indexes;
    delete [] temp;
    delete [] vel;
}
