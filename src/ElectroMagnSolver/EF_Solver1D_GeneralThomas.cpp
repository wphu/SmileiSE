
#include "EF_Solver1D_GeneralThomas.h"

#include "ElectroMagn.h"
#include "Field1D.h"

EF_Solver1D_GeneralThomas::EF_Solver1D_GeneralThomas(PicParams &params, SmileiMPI* smpi, int nx_sou_left)
    : Solver1D(params)
{
    SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);

    dx = params.cell_length[0];
    dx_inv_ = 1.0 / dx;
    dx_sq = dx * dx;
    nx = params.n_space_global[0]+1;
    nx_source_left = nx_sou_left;

    if(params.bc_em_type_x[0] == "Dirichlet"){
        bc_x_left = 1;
        bc_e_value[0][0] = params.bc_em_value_x[0];
    }
    else if(params.bc_em_type_x[0] == "Neumann"){
        bc_x_left = 2;
        bc_e_derivative[0][0] = params.bc_em_value_x[0];
    }
    else
    {
      bc_x_left = 0;
    }

    if(params.bc_em_type_x[1] == "Dirichlet"){
        bc_x_right = 1;
        bc_e_value[0][1] = params.bc_em_value_x[1];
    }
    else if(params.bc_em_type_x[1] == "Neumann"){
        bc_x_right = 2;
        bc_e_derivative[0][1] = params.bc_em_value_x[1];
    }
    else
    {
      bc_x_right = 0;
    }



    if(smpi1D->isMaster()){
        //grid1D = static_cast<Grid1D*>(grid);
        initGeneralThomas();
    }
    //initSLU();
}

EF_Solver1D_GeneralThomas::~EF_Solver1D_GeneralThomas()
{
}

void EF_Solver1D_GeneralThomas::operator()( ElectroMagn* fields, SmileiMPI* smpi)
{
    SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);
    // Static-cast of the fields

    Field1D* Ex1D           = static_cast<Field1D*>(fields->Ex_);
    Field1D* rho1D          = static_cast<Field1D*>(fields->rho_);
    Field1D* phi1D          = static_cast<Field1D*>(fields->phi_);
    Field1D* rho1D_global   = static_cast<Field1D*>(fields->rho_global);
    Field1D* phi1D_global   = static_cast<Field1D*>(fields->phi_global);
    Field1D* Ex1D_global    = static_cast<Field1D*>(fields->Ex_global);

    if(smpi1D->isMaster()){
        solve_GeneralThomas(rho1D_global, phi1D_global);
        //phi1D_global->put_to(0.0);
        solve_Ex(phi1D_global, Ex1D_global);
    }

    //Ex1D_global->put_to(0.0);
    //Ey1D_global->put_to(0.0);


    smpi1D->barrier();

    smpi1D->scatterField(Ex1D_global, Ex1D);
    smpi1D->scatterField(phi1D_global, phi1D);

}


void EF_Solver1D_GeneralThomas::initGeneralThomas()
{
    x.resize(nx, 0.0);
    y.resize(nx, 0.0);
    a.resize(nx, 0.0);
    b.resize(nx, 0.0);
    c.resize(nx, 0.0);
    p.resize(nx, 0.0);
    u.resize(nx, 0.0);
    q.resize(nx, 0.0);
    t.resize(nx, 0.0);
    h.resize(nx, 0.0);
    g.resize(nx, 0.0);

    for(int i = 0; i < nx; i++)
    {
        a[i] = 1.0;
        b[i] = -2.0;
        c[i] = 1.0;
    }

    if(bc_x_left == 1){
        a[0] = 0.0;
        b[0] = 1.0;
        c[0] = 0.0;

        a[nx-1] = 0.0;
        b[nx-1] = 1.0;
        c[nx-1] = 0.0;
    }
    else if(bc_x_left == 2){
        a[0] = 0.0;
        b[0] = 1.0;
        c[0] = -1.0;

        a[nx-1] = -1.0;
        b[nx-1] = 1.0;
        c[nx-1] = 0.0;
    }


}


void EF_Solver1D_GeneralThomas::solve_GeneralThomas(Field* rho, Field* phi)
{
    Field1D* rho1D = static_cast<Field1D*>(rho);
    Field1D* phi1D = static_cast<Field1D*>(phi);

    //> The boundary value can be changed with time
    for(int i = 0; i < nx; i++)
    {
        y[i] = -dx_sq * const_ephi0_inv * (*rho1D)(i);
    }

    if(bc_x_left == 1){
        y[0] = bc_e_value[0][0];
        y[nx-1] = bc_e_value[0][1];
    }
    else if(bc_x_left == 2){
        y[0] = bc_e_derivative[0][0];
        y[nx-1] = -bc_e_derivative[0][1];
    }


    p[0] = b[0];
    u[0] = y[0] / p[0];
    q[0] = - c[0] / p[0];
    t[0] = - a[0] / p[0];

    for(int i = 1; i < nx - 1; i++)
    {
      p[i] = a[i] * q[i-1] + b[i];
      u[i] = (y[i] - a[i] * u[i-1]) / p[i];
      q[i] = - c[i] / p[i];
      t[i] = - a[i]*t[i-1] / p[i];
    }

    h[nx-1] = (y[nx-1] - a[nx-1]*u[nx-2]) / ( a[nx-1] * ( q[nx-2] + t[nx-2] ) + b[nx-1] );
    g[nx-1] = c[nx-1] / ( a[nx-1] * ( g[nx-2] + t[nx-2] ) + b[nx-1] );
    h[nx-2] = u[nx-2] + ( q[nx-2] + t[nx-2] ) * h[nx-1];
    g[nx-2] = ( q[nx-2] + t[nx-2] ) * g[nx-1];
    for(int i = nx - 3; i > 0; i--)
    {
      h[i] = u[i] + q[i] * h[i+1] + t[i] * h[nx-1];
      g[i] = q[i] * g[i+1] + t[i] * g[nx-1];
    }

    x[0] = ( y[0] - c[0] * h[1] - a[0] * h[nx-1] ) / ( b[0] + c[0] * g[1] + a[0] * g[nx-1] );
    for(int i = 1; i <= nx-1; i++)
    {
      x[i] = h[i] + g[i] * x[0];
      DEBUGEXEC(if (!std::isfinite(x[i])) ERROR("Not finite "<< i << " = " << h[i] << g[i] << x[0]));
    }

    for(int i = 0; i < nx; i++)
    {
        (*phi1D)(i) = x[i];
    }


}


void EF_Solver1D_GeneralThomas::solve_Ex(Field* phi, Field* Ex)
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
