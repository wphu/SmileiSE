
#ifndef PARTSOURCE1D_EMIT_H
#define PARTSOURCE1D_EMIT_H

#include <vector>
#include "PartSource1D.h"


using namespace std;

class PartSource1D_Emit : public PartSource1D
{

public:
    //! Constructor for Collisions between two species
    PartSource1D_Emit(
        PicParams& params,
        SmileiMPI* smpi,
        string emit_emitKind,
        unsigned int emit_species1,
        string emitPosition,
        unsigned int nPartEmit,
        double emitTemperature,
        double emitJ,
        double weight_const,
        double emitOffset,
        double a_FN,
        double b_FN,
        double work_function,
        string emit_relSpecies);

    ~PartSource1D_Emit();



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
    // electric field used for field emit, equal to electric field on the boundary
    double emitField;
    // Emitting current density from boundary
    double emitJ;
    // weight of emitting particles
    double weight_const;
    // nominalDensity and nomPtclsPerCell is used to set the weight_const
    // weight_cosnt = nominalDensity * CellVolume / nomPtclsPerCell
    double nominalDensity;
    double nomPtclsPerCell;
    // emitting tempreature
    double emitOffset;

    // some parameters for field emission
    // from "modelling vacuum arcs: from plasma initiation to surface interactions"
    double a_FN;
    double b_FN;
    double work_function;
    double ySqrt_factor;
    double a_factor;
    double b_factor;

    double t_y2(double Eloc) { return 1.0; };
    double v_y(double Eloc)
    {
        return 0.956 - 1.062*ySqrt_factor*Eloc;
    }

    //paramets for relEmit
    double relEmit_factor;



};


#endif
