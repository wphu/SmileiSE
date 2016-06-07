
#include "EF_Solver1D_TDMA.h"

#include "ElectroMagn.h"
#include "Field1D.h"

EF_Solver1D_TDMA::EF_Solver1D_TDMA(PicParams &params, SmileiMPI* smpi)
    : Solver1D(params)
{
    SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);

    dx = params.cell_length[0];
    dx_inv_ = 1.0 / dx;
    dx_sq = dx * dx;
    nx = params.n_space_global[0]+1;

    if(params.bc_em_type_x[0] == "Dirichlet"){
        bc_x_left = 1;
        field_left = params.bc_em_value_x[0];
    }
    else if(params.bc_em_type_x[0] == "Neumann"){
        bc_x_left = 2;
        field_derivative_left = params.bc_em_value_x[0];
    }

    if(params.bc_em_type_x[1] == "Dirichlet"){
        bc_x_right = 1;
        field_right = params.bc_em_value_x[1];
    }
    else if(params.bc_em_type_x[1] == "Neumann"){
        bc_x_right = 2;
        field_derivative_right = params.bc_em_value_x[1];
    }



    if(smpi1D->isMaster()){
        //grid1D = static_cast<Grid1D*>(grid);
        initTDMA();
    }
    //initSLU();
}

EF_Solver1D_TDMA::~EF_Solver1D_TDMA()
{
}

void EF_Solver1D_TDMA::operator()( ElectroMagn* fields, SmileiMPI* smpi)
{
    SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);
    // Static-cast of the fields

    Field1D* Ex1D           = static_cast<Field1D*>(fields->Ex_);
    Field1D* rho1D          = static_cast<Field1D*>(fields->rho_);
    Field1D* phi1D          = static_cast<Field1D*>(fields->phi_);
    Field1D* rho1D_global   = static_cast<Field1D*>(fields->rho_global);
    Field1D* phi1D_global   = static_cast<Field1D*>(fields->phi_global);
    Field1D* Ex1D_global    = static_cast<Field1D*>(fields->Ex_global);

    smpi1D->barrier();
    smpi1D->gatherRho(rho1D_global, rho1D);

    if(smpi1D->isMaster()){
        solve_TDMA(rho1D_global, phi1D_global);
        //phi1D_global->put_to(0.0);
        solve_Ex(phi1D_global, Ex1D_global);
    }

    //Ex1D_global->put_to(0.0);
    //Ey1D_global->put_to(0.0);


    smpi1D->barrier();

    smpi1D->scatterField(Ex1D_global, Ex1D);
    smpi1D->scatterField(phi1D_global, phi1D);

}


void EF_Solver1D_TDMA::initTDMA()
{
    a = new double[nx];
    b = new double[nx];
    c = new double[nx];
    f = new double[nx];
    e = new double[nx];
    d = new double[nx];
    for(int i =1; i < nx-1; i++)
    {
        a[i] = 1.0;
        b[i] = -2.0;
        c[i] = 1.0;
    }

    if(bc_x_left == 1){
        a[0] = 0.0;
        b[0] = 1.0;
        c[0] = 0.0;
    }
    else if(bc_x_left == 2){
        a[0] = 0.0;
        b[0] = 1.0;
        c[0] = -1.0;
    }

    if(bc_x_right == 1){
        a[nx-1] = 0.0;
        b[nx-1] = 1.0;
        c[nx-1] = 0.0;
    }
    else if(bc_x_right == 2){
        a[nx-1] = -1.0;
        b[nx-1] = 1.0;
        c[nx-1] = 0.0;
    }



}


void EF_Solver1D_TDMA::solve_TDMA(Field* rho, Field* phi)
{
    Field1D* rho1D = static_cast<Field1D*>(rho);
    Field1D* phi1D = static_cast<Field1D*>(phi);

    //> The boundary value can be changed with time
    for(int i = 1; i < nx-1; i++)
    {
        f[i] = -dx_sq * (*rho1D)(i);
    }

    if(bc_x_left == 1){
        f[0] = field_left;
    }
    else if(bc_x_left == 2){
        f[0] = field_derivative_left;
    }

    if(bc_x_right == 1){
        f[nx-1] = field_right;
    }
    else if(bc_x_right == 2){
        f[nx-1] = -field_derivative_right;
    }

    e[0] = c[0] / b[0];
    d[0] = f[0] / b[0];
    for(int i =1; i < nx-1; i++)
    {
        e[i] = c[i] / ( b[i] - a[i] * e[i-1] );
        d[i] = ( f[i] -a[i] * d[i-1] ) / ( b[i] - a[i] * e[i-1] );
    }

    (*phi1D)(nx-1) = ( f[nx-1] - a[nx-1] * d[nx-2] ) / ( b[nx-1] - a[nx-1] * e[nx-2] );
    for(int i = nx-2; i >= 0; i--)
    {
        (*phi1D)(i) = d[i] - e[i] * (*phi1D)(i+1);
    }


}


void EF_Solver1D_TDMA::solve_Ex(Field* phi, Field* Ex)
{
    Field1D* phi1D = static_cast<Field1D*>(phi);
    Field1D* Ex1D = static_cast<Field1D*>(Ex);


    for(int i = 1; i < nx-1; i++)
    {
        (*Ex1D)(i) = - ((*phi1D)(i+1) - (*phi1D)(i-1)) *0.5 * dx_inv_;
    }

    (*Ex1D)(0) = -(-3.0 * (*phi1D)(0) + 4.0 * (*phi1D)(1) - (*phi1D)(2)) *0.5 * dx_inv_;
    (*Ex1D)(nx-1) = -((*phi1D)(nx-3) - 4.0 * (*phi1D)(nx-2) + 3.0 * (*phi1D)(nx-1)) *0.5 * dx_inv_;

}
