#ifndef SOLVER_H
#define SOLVER_H

#include "PicParams.h"
#include "SmileiMPI.h"

class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Solver
//  --------------------------------------------------------------------------------------------------------------------
class Solver
{

public:
    //! Creator for Solver
    Solver(PicParams &params) {};
    virtual ~Solver() {};

    //! Overloading of () operator
    virtual void operator()( ElectroMagn* fields)=0;
    virtual void operator()( ElectroMagn* fields, SmileiMPI* smpi)=0;
    virtual void solve_SLU(Field* rho, Field* phi){};
    virtual void finishSLU(){};
    virtual void initSLU_test(){};
    virtual void initSLU(){};
    virtual void solve_Exy(Field* phi, Field* Ex, Field* Ey){};

protected:

};//END class

#endif
