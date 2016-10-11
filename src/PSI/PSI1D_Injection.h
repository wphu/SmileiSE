
#ifndef PSI1D_INJECTION_H
#define PSI1D_INJECTION_H

#include <vector>
#include "PSI1D.h"


using namespace std;

class PSI1D_Injection : public PSI1D
{

public:
    //! Constructor for Collisions between two species
    PSI1D_Injection(
        PicParams& params,
        SmileiMPI* smpi,
        string emitKind,
        unsigned int psi_species1,
        string psiPosition,
        unsigned int nPartEmit,
        double emitTemperature,
        double emitJ,
        double weight_const,
        double emitOffset,
        double a_FN,
        double b_FN,
        double work_function,
        string relSpecies);

    ~PSI1D_Injection();



    //! Method called in the main smilei loop to apply collisions at each timestep
    void performPSI(PicParams&, SmileiMPI* smpi, std::vector<Species*>&,int, ElectroMagn* );

    // emit particles
    void emit(PicParams&, vector<Species*>&);


    //particle number emitted every timestep
    unsigned int nPartEmit;

    //! Emitted temporary species
    Particles emit_particles;




private:
    double dt_ov_dx;
    // electric field used for field emit, equal to electric field on the boundary
    double emitField;
    // Emitting current density from boundary
    double emitJ;
    // weight of emitting particles
    double weight_const;
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
