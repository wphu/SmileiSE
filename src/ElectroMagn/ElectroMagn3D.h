#ifndef ELECTROMAGN3D_H
#define ELECTROMAGN3D_H

#include "ElectroMagn.h"

class PicParams;

//! class ElectroMagn3D containing all information on the electromagnetic fields & currents for 2d3v simulations
class ElectroMagn3D : public ElectroMagn
{
public:
    //! Constructor for ElectroMagn3D
    ElectroMagn3D(PicParams &params, InputData &input_data, SmileiMPI* smpi);

    //! Destructor for ElectroMagn3D
    ~ElectroMagn3D();


    //! Method used to solve Maxwell-Ampere equation
    void solveMaxwellAmpere();

    //! Method used to save the Magnetic fields (used to center them)
    void saveMagneticFields();

    //! Method used to center the Magnetic fields (used to push the particles)
    void centerMagneticFields();

    //! Method used to reset/increment the averaged fields
    void incrementAvgFields(unsigned int time_step);

    //! Method used to initialize the total charge densities and currents
    void restartRhoJ();
    //! Method used to initialize the total charge densities and currents of species
    void restartRhoJs(int ispec, bool currents);

    //! Method used to compute the total charge density and currents by summing over all species
    void computeTotalRhoJ();

    // gather time-average fields to process 0 to output
    void gatherAvgFields(SmileiMPI *smpi);

    // gather fields to process 0 to calculate EM fields
    void gatherFields(SmileiMPI *smpi);

    //! \todo Create properties the laser time-profile (MG & TV)

    //! Number of nodes on the primal grid in the x-direction
    unsigned int nx_p;

    //! Number of nodes on the dual grid in the x-direction
    unsigned int nx_d;

    //! Number of nodes on the primal grid in the y-direction
    unsigned int ny_p;

    //! Number of nodes on the dual grid in the y-direction
    unsigned int nz_d;

    //! Number of nodes on the primal grid in the z-direction
    unsigned int nz_p;

    //! Number of nodes on the dual grid in the z-direction
    unsigned int ny_d;

    //! Spatial step dx for 3D3V cartesian simulations
    double dx;

    //! Spatial step dy for 3D3V cartesian simulations
    double dy;

    //! Spatial step dz for 3D3V cartesian simulations
    double dz;

    //! Ratio of the time-step by the spatial-step dt/dx for 3D3V cartesian simulations
    double dt_ov_dx;

    //! Ratio of the time-step by the spatial-step dt/dy for 3D3V cartesian simulations
    double dt_ov_dy;

    //! Ratio of the time-step by the spatial-step dt/dz for 3D3V cartesian simulations
    double dt_ov_dz;

    //! Ratio of the spatial-step by the time-step dx/dt for 3D3V cartesian simulations
    double dx_ov_dt;

    //! Ratio of the spatial-step by the time-step dy/dt for 3D3V cartesian simulations
    double dy_ov_dt;

    //! Ratio of the spatial-step by the time-step dz/dt for 3D3V cartesian simulations
    double dz_ov_dt;

    //! compute Poynting on borders
    void computePoynting();

    //> gather time-average fields to process 0 to output

private:

    //! from smpi is west
    const bool isWestern;

    //! from smpi is east
    const bool isEastern;

    //! from smpi is north
    const bool isSouthern;

    //! from smpi is south
    const bool isNorthern;
};

#endif
