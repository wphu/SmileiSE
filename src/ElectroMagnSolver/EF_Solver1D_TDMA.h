#ifndef EF_SOLVER1D_TDMA
#define EF_SOLVER1D_TDMA

#include "Solver1D.h"
#include "SmileiMPI_Cart1D.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class EF_Solver1D_TDMA : public Solver1D
{

public:
    //! Creator for EF_Solver1D_TDMA
    EF_Solver1D_TDMA(PicParams &params, SmileiMPI* smpi);
    virtual ~EF_Solver1D_TDMA();

    //! Overloading of () operator
    virtual void operator()( ElectroMagn* fields){};
    virtual void operator()( ElectroMagn* fields, SmileiMPI* smpi);

    void solve_TDMA(Field* rho, Field* phi);
    void initTDMA();
    void solve_Ex(Field* phi, Field* Ex);
    //>>>SuperLU parameters

    //> boundary conditions for left and right sides: = 1 is Dirichlet, and = 2 is Neumann
    int bc_x_left, bc_x_right;

    //> boundary field value for Dirichlet boundary condition
    double field_left, field_right;
    //> boundary field gradient value (normal derivative) for Neumann boundary condition
    double field_derivative_left, field_derivative_right;

    //>>>geometry parameters
    int nx;
    double dx;
    double dx_inv_;
    double dx_sq;


protected:

    double *a, *b, *c, *f, *e, *d;

};//END class

#endif
