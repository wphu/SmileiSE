#ifndef EF_SOLVER1D_TDMA_IMP_H
#define EF_SOLVER1D_TDMA_IMP_H

#include <vector>

#include "Solver1D.h"
#include "SmileiMPI_Cart1D.h"
class ElectroMagn;

using namespace std;
//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class EF_Solver1D_TDMA_imp : public Solver1D
{

public:
    //! Creator for EF_Solver1D_TDMA_imp
    EF_Solver1D_TDMA_imp(PicParams &params, SmileiMPI* smpi, int nx_sou_left);
    virtual ~EF_Solver1D_TDMA_imp();

    //! Overloading of () operator
    virtual void operator()( ElectroMagn* fields){};
    virtual void operator()( ElectroMagn* fields, SmileiMPI* smpi);

    void initTDMA(PicParams &params);
    void solve_TDMA_imp(ElectroMagn* fields);
    void solve_Ex(Field* phi, Field* Ex);

    // no source region for electric field
    void initTDMA_org();
    void solve_TDMA_org(Field* rho, Field* phi);
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

    // grid point number of source region in the left side
    int nx_source_left;


protected:

    double *a, *b, *c, *f, *e, *d;
    double dt;
    // the [0,0] term of tensor chi(3,3), equation (5) and (3)
    double chi00;
    vector< double > factor_chi;

};//END class

#endif
