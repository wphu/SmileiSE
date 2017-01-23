#ifndef DIAGNOSTIC1D_H
#define DIAGNOSTIC1D_H

#include "Diagnostic.h"
#include "PicParams.h"
#include "SmileiMPI.h"

class Diagnostic1D : public Diagnostic {

public :

    Diagnostic1D() {};
    virtual ~Diagnostic1D() {};


    //! Misc init.
    virtual void init(PicParams& params, SmileiMPI* smpi) {};

    //! Runs the diag for a given patch for global diags.
    virtual void run( int timestep ) {};

    //! Runs the diag for all patches for local diags.
    virtual void run( SmileiMPI* smpi, int timestep ) {};

	vector< vector<double> > particleFlux;	//particleFlux[ispec][iDirection]
	vector< vector<double> > heatFlux;


protected :


};

#endif
