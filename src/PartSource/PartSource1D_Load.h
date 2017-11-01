
#ifndef PARTSOURCE1D_LOAD_H
#define PARTSOURCE1D_LOAD_H

#include <vector>
#include "PartSource1D.h"


using namespace std;

class PartSource1D_Load : public PartSource1D
{

public:
    //! Constructor for Collisions between two species
    PartSource1D_Load(
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
    double          load_Pos_end );

    ~PartSource1D_Load();



    //! Method called in the main smilei loop to apply collisions at each timestep
    void emitLoad(PicParams&, SmileiMPI* smpi, std::vector<Species*>&,int, ElectroMagn* );


    // =================Parameters for loading particles=================
    vector<int> loadTimeStepVector;
    vector<double> loadTemperatureVector;
    vector<double> loadDnVector;
    double loadDensity;
    double loadTemperature;
    double loadTemperature_init;
    double loadTemperature_exceed;
    // load density per second [m-3 s-1]
    double loadDn;
    // load heat flux density = loadDn * loadTemperature
    double loadq;
    int loadStep;
    // Number of particles loaded in one cell at one loadStep
    int loadNumber;
    double loadNumber_heat;
    // loadRem = loadStep * loadDn *... - loadNumber
    double loadRem;
    double loadRemTot;

    // Position for loading particles in the current MPI region!!!
    double loadPos_start;
    double loadPos_end;
    int loadBin_start;
    int loadBin_end;

    double loadPos_Ystart;
    double loadPos_Yend;
    int loadBin_Ystart;
    int loadBin_Yend;

    // MPI rank of source middle region
    int mpiRank_source_middle;
    // index of field at source middle;
    int index_source_middle;


private:
    double dt_ov_dx;
    double dt;
    double YZArea;

    // nominalDensity and nomPtclsPerCell is used to set the weight_const
    // weight_cosnt = nominalDensity * CellVolume / nomPtclsPerCell
    double nominalDensity;
    double nomPtclsPerCell;
    // emitting tempreature


    // ========== Some variables for emitLoad method ===================

    // Parameters for "nq"
    int timeStep_checkFor_nq;
    double temperature_pre;
    double source_density_pre;
    double loadTemperature_upLimit_factor;

    // Parameters for "dn"
    int nextTimeStep = 0;
    int nextTimeStep_index = 1;





};


#endif
