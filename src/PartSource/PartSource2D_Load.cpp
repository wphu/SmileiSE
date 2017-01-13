#include "PartSource2D_Load.h"
#include "SmileiMPI_Cart2D.h"
#include "Field2D.h"
#include "H5.h"


#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <sstream>

using namespace std;


// Constructor
PartSource2D_Load::PartSource2D_Load(
    PicParams& params,
    SmileiMPI* smpi,
    unsigned int load_species1,
    vector<double> mean_vel,
    string load_kind,
    int load_step,
    int every_time,
    double load_dn,
    double load_density,
    double load_temperature,
    double load_Pos_start,
    double load_Pos_end,
    double load_Pos_Ystart,
    double load_Pos_Yend
):
PartSource2D (params, smpi)

{
    dt_ov_dx = params.timestep / params.cell_length[0];
    dt = params.timestep;

    species1        = load_species1;
    mean_velocity   = mean_vel;
    loadKind        = load_kind;
    loadStep        = load_step;
    everyTime       = every_time;
    loadDn          = load_dn;
    loadDensity     = load_density;
    loadTemperature = load_temperature;
    loadPos_start   = load_Pos_start;
    loadPos_end     = load_Pos_end;
    loadPos_Ystart   = load_Pos_Ystart;
    loadPos_Yend     = load_Pos_Yend;

    YZArea = 1.0;

    // define the x-direction range of source region
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

    // define the y-direction range of source region
    // the MPI domain is not in the source region
    if(smpi->getDomainLocalMax(1) <= loadPos_Ystart || smpi->getDomainLocalMin(1) >= loadPos_Yend) {
        loadPos_Ystart = 0.0;
        loadBin_Ystart = 0;
        loadPos_Yend = 0.0;
        loadBin_Yend = 0;
    }
    // the down end of source region is in the MPI domain
    else if(smpi->getDomainLocalMin(1) <= loadPos_Ystart && smpi->getDomainLocalMax(1) > loadPos_Ystart
    && smpi->getDomainLocalMax(1) <= loadPos_Yend) {
        double temp = (loadPos_Ystart - smpi->getDomainLocalMin(1)) / params.cell_length[1];
        loadBin_Ystart = temp + 0.001;
        if(loadBin_Ystart < 0) { loadBin_Ystart = 0; }
        loadPos_Ystart = smpi->getDomainLocalMin(1) + loadBin_Ystart*params.cell_length[1];
        //numPart_in_each_bin[loadBin_start] /=
        //( params.cell_length[0] / (( smpi->getDomainLocalMin(0) + (loadBin_start+1)*params.cell_length[0] ) - loadPos_start) );
        loadPos_Yend = smpi->getDomainLocalMax(1);
        loadBin_Yend = params.n_space[1] - 1;
    }
    // the up end of source region is in the MPI domain
    else if(smpi->getDomainLocalMin(1) < loadPos_Yend && smpi->getDomainLocalMax(1) >= loadPos_Yend
    && smpi->getDomainLocalMin(1) >= loadPos_Ystart) {
        loadPos_Ystart = smpi->getDomainLocalMin(1);
        loadBin_Ystart = 0;
        loadBin_Yend = (loadPos_Yend - smpi->getDomainLocalMin(1)) / params.cell_length[1];
        if(loadBin_Yend > params.n_space[1] - 1) { loadBin_Yend = params.n_space[1] - 1; }
        //numPart_in_each_bin[loadBin_end] /=
        //( params.cell_length[0] / (loadPos_end - ( smpi->getDomainLocalMin(0) + loadBin_end*params.cell_length[0] )) );
        loadPos_Yend = smpi->getDomainLocalMin(1) + (loadBin_Yend+1)*params.cell_length[1];
    }
    // the whole MPI domain is in the source region
    else if(smpi->getDomainLocalMin(1) >= loadPos_Ystart && smpi->getDomainLocalMax(1) <= loadPos_Yend) {
        loadPos_Ystart = smpi->getDomainLocalMin(1);
        loadBin_Ystart = 0;
        loadPos_Yend = smpi->getDomainLocalMax(1);
        loadBin_Yend = params.n_space[1] - 1;
    }
    // the whole MPI domain source region is in the bin
    else if(smpi->getDomainLocalMin(1) < loadPos_Ystart && smpi->getDomainLocalMax(1) > loadPos_Yend) {
        //loadBin_start = (loadPos_start - smpi->getDomainLocalMin(0)) / params.cell_length[0];
        //numPart_in_each_bin[loadBin_start] /=
        //( params.cell_length[0] / (( smpi->getDomainLocalMin(0) + (loadBin_start+1)*params.cell_length[0] ) - loadPos_start) );
        //loadBin_end = (loadPos_end - smpi->getDomainLocalMin(0)) / params.cell_length[0];
        //numPart_in_each_bin[loadBin_end] /=
        //( params.cell_length[0] / (loadPos_end - ( smpi->getDomainLocalMin(0) + loadBin_start*params.cell_length[0] )) );
        double temp = (loadPos_Ystart - smpi->getDomainLocalMin(1)) / params.cell_length[1];
        loadBin_Ystart = temp + 0.001;
        if(loadBin_Ystart < 0) { loadBin_Ystart = 0; }
        loadPos_Ystart = smpi->getDomainLocalMin(1) + loadBin_Ystart*params.cell_length[1];
        loadBin_Yend = (loadPos_Yend - smpi->getDomainLocalMin(1)) / params.cell_length[1];
        if(loadBin_Yend > params.n_space[1] - 1) { loadBin_Yend = params.n_space[1] - 1; }
        loadPos_Yend = smpi->getDomainLocalMin(1) + (loadBin_Yend+1)*params.cell_length[1];

    }

    // load particles by number_density and Temperature
    numPart_in_each_cell = loadDensity / params.species_param[species1].weight;
    if(loadKind == "nT") {
        for(int ibin = 0; ibin < numPart_in_each_bin.size(); ibin++)
        {
            numPart_in_each_bin[ibin] = (loadBin_Yend - loadBin_Ystart + 1) * numPart_in_each_cell;
            //cout<<"numPart_in_each_bin "<<numPart_in_each_bin[ibin]<<endl;
        }
    }


    emitKind = "none";

}

PartSource2D_Load::~PartSource2D_Load()
{

}



// Calculates the PartSource for a given PartSource object
void PartSource2D_Load::emitLoad(PicParams& params, SmileiMPI* smpi, vector<Species*>& vecSpecies, int itime, ElectroMagn* fields)
{
    Species   *s1;
    Particles *p1;
    int iPart, nPart;
    vector <double> cell_length;            // cell_length for initialize position, maybe not equal to real cell lenghth
    vector<double> max_jutt_cumul;
    double *indexes=new double[params.nDim_particle];
    double *temp=new double[3];
    double *vel=new double[3];

    if(everyTime == 0 && itime > 1) { return; }
    if(loadKind == "nT" && loadBin_end != loadBin_start && loadBin_Yend != loadBin_Ystart) {
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
            for(int iPart = s1->bmin[ibin]; iPart < s1->bmax[ibin]; iPart++)
            {
                if( p1->position(1,iPart) > loadPos_Ystart && p1->position(1,iPart) < loadPos_Yend ) {
                    indexes_of_particles_to_erase.push_back(iPart);
                }
            }
            count_of_particles_to_insert[ibin] = numPart_in_each_bin[ibin];
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
            // Because bin is only in x-direction, y-direction has no bin
            // The "Bin" in loadBin_Ystart and loadBin_Yend is just to be consistent with the x-direction
            for(int j = loadBin_Ystart; j <= loadBin_Yend; j++)
            {
                iPart = s1->bmax[ibin] - (loadBin_Yend - j + 1) * numPart_in_each_cell ;
                nPart = numPart_in_each_cell;
                cell_length[0] = params.cell_length[0];
                cell_length[1] = params.cell_length[1];
                indexes[0] = smpi->getDomainLocalMin(0) + ibin*params.cell_length[0];
                indexes[1] = smpi->getDomainLocalMin(1) + j   *params.cell_length[1];

                s1->initPosition(nPart, iPart, indexes, params.nDim_particle,
                             cell_length, s1->species_param.initPosition_type);

                s1->initMomentum(nPart,iPart, temp, vel,
                             s1->species_param.initMomentum_type, max_jutt_cumul, params);

                s1->initWeight_constant(nPart, species1, iPart, s1->species_param.weight);
                s1->initCharge(nPart, species1, iPart, s1->species_param.charge);
            }
        }
    }
    else if(loadKind == "dn" && itime%loadStep == 0 && loadBin_end != loadBin_start) {
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

        for(int ibin = 0; ibin < count_of_particles_to_insert.size(); ibin++ )
        {
            count_of_particles_to_insert[ibin] = 0;
            if(ibin >= loadBin_start && ibin <= loadBin_end) { count_of_particles_to_insert[ibin] = loadDn * loadStep * params.timestep / s1->species_param.weight; }
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

    delete [] indexes;
    delete [] temp;
    delete [] vel;
}
