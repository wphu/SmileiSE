#ifndef DIAGNOSTIC_H
#define DIAGNOSTIC_H

#include "PicParams.h"
#include "SmileiMPI.h"

#include <iostream>
#include <vector>

using namespace std;

class Diagnostic {

public :

    Diagnostic(PicParams &params);
    virtual ~Diagnostic() {};

    //! Runs the diag for all patches for local diags.
    virtual void run( SmileiMPI* smpi, vector<Species*>& vecSpecies, ElectroMagn* EMfields, int timestep ) {};

    const unsigned int n_species;

protected :

    // pi * 0.5
    double PI_ov_2;
    int dump_step;
    int avg_step;
    double timestep;
    double const_e;

    // Phi is the angle between the magnetic field and the y-direction
    double sinPhi, cosPhi;

    vector<double> sim_length;

    vector<unsigned int> oversize;
};

#endif
