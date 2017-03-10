
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
        PicParams&      params,
        SmileiMPI*      smpi,
        string          emit_emitKind,
        unsigned int    emit_species1,
        string          emitPosition,
        int             emitNumber,
        double          emitTemperature,
        double          emitJ,
        double          emitFlux,
        double          emitOffset,
        double          a_FN,
        double          b_FN,
        double          work_function,
        string          emit_relSpecies );

    ~PartSource1D_Emit();



    //! Method called in the main smilei loop to apply collisions at each timestep
    void emitLoad(PicParams&, SmileiMPI* smpi, std::vector<Species*>&,int, ElectroMagn* );

    // emit particles
    void emit(PicParams&, vector<Species*>&);


    //particle number emitted every timestep
    unsigned int nPartEmit;

    // PSI position : only left and right for 1D case
    string emitPos;

    // the energy/temperature of the new particles
    double emitTemp;

    // Emitted particle flux from boundary
    double emitFlux;

    // emitting tempreature
    double emitOffset;

    // The particle number emitted from boundary at one time (emitStep), not one timestep(!!!)
    int emitNumber;
    int emitStep;

    // emitRem = emitStep * emitDn *... - emitNumber
    double emitRem;
    double emitRemTot;

private:
    double dt_ov_dx;
    double dt;
    double YZArea;
    // electric field used for field emit, equal to electric field on the boundary
    double emitField;
    // Emitted current density from boundary
    double emitJ;

    // nominalDensity and nomPtclsPerCell is used to set the weight_const
    // weight_cosnt = nominalDensity * CellVolume / nomPtclsPerCell
    double nominalDensity;
    double nomPtclsPerCell;


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
