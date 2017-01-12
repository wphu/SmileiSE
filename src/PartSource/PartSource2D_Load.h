
#ifndef PARTSOURCE2D_LOAD_H
#define PARTSOURCE2D_LOAD_H

#include <vector>
#include "PartSource2D.h"


using namespace std;

class PartSource2D_Load : public PartSource2D
{

public:
    //! Constructor for Collisions between two species
    PartSource2D_Load(
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
    double load_Pos_Yend);

    ~PartSource2D_Load();



    //! Method called in the main smilei loop to apply collisions at each timestep
    void emitLoad(PicParams&, SmileiMPI* smpi, std::vector<Species*>&,int, ElectroMagn* );

    // emit particles
    void emit(PicParams&, vector<Species*>&);


    //particle number emitted every timestep
    unsigned int nPartEmit;




private:
    double dt_ov_dx;
    double dt;
    double YZArea;

    double weight_const;
    // nominalDensity and nomPtclsPerCell is used to set the weight_const
    // weight_cosnt = nominalDensity * CellVolume / nomPtclsPerCell
    double nominalDensity;
    double nomPtclsPerCell;
    // emitting tempreature


};


#endif