#ifndef DIAGNOSTIC_H
#define DIAGNOSTIC_H

#include "PicParams.h"
#include "SmileiMPI.h"

#include <iostream>
#include <vector>

using namespace std;

class Diagnostic {

public :

    Diagnostic() {};
    virtual ~Diagnostic() {};


    //! Misc init.
    virtual void init(PicParams& params, SmileiMPI* smpi) {};

    //! Runs the diag for a given patch for global diags.
    virtual void run( int timestep ) {};

    //! Runs the diag for all patches for local diags.
    virtual void run( SmileiMPI* smpi, int timestep ) {};


protected :


};

#endif
