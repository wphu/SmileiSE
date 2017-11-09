// ref: 求解周期性三对角方程组的广义Thomas算法_王兴波


#ifndef EF_SOLVER1D_GENERALTHOMAS_H
#define EF_SOLVER1D_GENERALTHOMAS_H

#include "Solver1D.h"
#include "SmileiMPI_Cart1D.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class EF_Solver1D_GeneralThomas : public Solver1D
{

public:
    //! Creator for EF_Solver1D_GeneralThomas
    EF_Solver1D_GeneralThomas(PicParams &params, SmileiMPI* smpi, int nx_sou_left);
    virtual ~EF_Solver1D_GeneralThomas();

    //! Overloading of () operator
    virtual void operator()( ElectroMagn* fields){};
    virtual void operator()( ElectroMagn* fields, SmileiMPI* smpi);

    void initGeneralThomas();
    void solve_GeneralThomas(Field* rho, Field* phi);
    void solve_Ex(Field* phi, Field* Ex);

    // no source region for electric field
    void initGeneralThomas_org();
    void solve_GeneralThomas_org(Field* rho, Field* phi);
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

    vector<double> x, y, a, b, c, p, u, q, t, h, g;

};//END class

#endif
