
#include "EF_Solver1D_TDMA_imp.h"

#include "ElectroMagn.h"
#include "Field1D.h"

EF_Solver1D_TDMA_imp::EF_Solver1D_TDMA_imp(PicParams &params, SmileiMPI* smpi, int nx_sou_left)
    : Solver1D(params)
{
    SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);

    dx = params.cell_length[0];
    dx_inv_ = 1.0 / dx;
    dx_sq = dx * dx;
    nx = params.n_space_global[0]+1;
    nx_source_left = nx_sou_left;
    dt = params.timestep;


    if(params.bc_em_type_x[0] == "Dirichlet"){
        bc_x_left = 1;
        bc_e_value[0][0] = params.bc_em_value_x[0];
    }
    else if(params.bc_em_type_x[0] == "Neumann"){
        bc_x_left = 2;
        bc_e_derivative[0][0] = params.bc_em_value_x[0];
    }

    if(params.bc_em_type_x[1] == "Dirichlet"){
        bc_x_right = 1;
        bc_e_value[0][1] = params.bc_em_value_x[1];
    }
    else if(params.bc_em_type_x[1] == "Neumann"){
        bc_x_right = 2;
        bc_e_derivative[0][1] = params.bc_em_value_x[1];
    }



    if(smpi1D->isMaster()){
        //grid1D = static_cast<Grid1D*>(grid);
        initTDMA(params);
    }
    //initSLU();
}

EF_Solver1D_TDMA_imp::~EF_Solver1D_TDMA_imp()
{
}

void EF_Solver1D_TDMA_imp::operator()( ElectroMagn* fields, SmileiMPI* smpi)
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
        solve_TDMA_imp(fields);
        //phi1D_global->put_to(0.0);
        solve_Ex(phi1D_global, Ex1D_global);
    }

    //Ex1D_global->put_to(0.0);
    //Ey1D_global->put_to(0.0);


    smpi1D->barrier();

    smpi1D->scatterField(Ex1D_global, Ex1D);
    smpi1D->scatterField(phi1D_global, phi1D);

}


void EF_Solver1D_TDMA_imp::initTDMA(PicParams &params)
{
    double charge_over_mass_;
    double Omega_square;
    double Omega0_square;
    a = new double[nx - nx_source_left];
    b = new double[nx - nx_source_left];
    c = new double[nx - nx_source_left];
    f = new double[nx - nx_source_left];
    e = new double[nx - nx_source_left];
    d = new double[nx - nx_source_left];
    for(int i =1; i < nx-1-nx_source_left; i++)
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
        a[nx-1-nx_source_left] = 0.0;
        b[nx-1-nx_source_left] = 1.0;
        c[nx-1-nx_source_left] = 0.0;
    }
    else if(bc_x_right == 2){
        a[nx-1-nx_source_left] = -1.0;
        b[nx-1-nx_source_left] = 1.0;
        c[nx-1-nx_source_left] = 0.0;
    }

    factor_chi.resize(params.species_param.size());
    for(int iS = 0; iS < factor_chi.size(); iS++)
    {

        charge_over_mass_ = params.species_param[iS].charge / params.species_param[iS].mass;
        factor_chi[iS] = dt * dt * 0.5 * const_ephi0_inv * charge_over_mass_ * params.species_param[iS].charge;
        // Calculate the [0,0] term of the rotation tensor T
        Omega0_square = pow(0.5 * charge_over_mass_ * dt * params.externB[0], 2);
        Omega_square =  pow(0.5 * charge_over_mass_ * dt * params.externB[0], 2)
                      + pow(0.5 * charge_over_mass_ * dt * params.externB[1], 2)
                      + pow(0.5 * charge_over_mass_ * dt * params.externB[2], 2);
        factor_chi[iS] *= ( (1.0 + Omega0_square) / (1.0 + Omega_square) );
    }


}


void EF_Solver1D_TDMA_imp::solve_TDMA_imp(ElectroMagn* fields)
{
    Field1D* rho1D = static_cast<Field1D*>(fields->rho_global);
    Field1D* phi1D = static_cast<Field1D*>(fields->phi_global);
    double ephi_inv;
    double chi00;
    double charge_over_mass_;

    //> The boundary value can be changed with time
    for(int i = 1; i < nx-1-nx_source_left; i++)
    {
        chi00 = 0.0;
        for(int iS = 0; iS < fields->n_species; iS++)
        {
            chi00 += factor_chi[iS] * (*fields->rho_s_global[iS])(i+nx_source_left);
        }
        ephi_inv = const_ephi0_inv / ( 1.0 + chi00 );
        f[i] = -dx_sq * ephi_inv * (*rho1D)(i+nx_source_left);
    }

    if(bc_x_left == 1){
        f[0] = bc_e_value[0][0];
    }
    else if(bc_x_left == 2){
        f[0] = bc_e_derivative[0][0];
    }

    if(bc_x_right == 1){
        f[nx-1-nx_source_left] = bc_e_value[0][1];
    }
    else if(bc_x_right == 2){
        f[nx-1-nx_source_left] = -bc_e_derivative[0][1];
    }

    e[0] = c[0] / b[0];
    d[0] = f[0] / b[0];
    for(int i =1; i < nx-1-nx_source_left; i++)
    {
        e[i] = c[i] / ( b[i] - a[i] * e[i-1] );
        d[i] = ( f[i] -a[i] * d[i-1] ) / ( b[i] - a[i] * e[i-1] );
    }

    (*phi1D)(nx-1) = ( f[nx-1-nx_source_left] - a[nx-1-nx_source_left] * d[nx-2-nx_source_left] ) / ( b[nx-1-nx_source_left] - a[nx-1-nx_source_left] * e[nx-2-nx_source_left] );
    for(int i = nx-2-nx_source_left; i >= 0; i--)
    {
        (*phi1D)(i+nx_source_left) = d[i] - e[i] * (*phi1D)(i+1+nx_source_left);
    }

    for(int i=0; i<nx_source_left; i++)
    {
        (*phi1D)(i) = (*phi1D)(nx_source_left);
    }


}


void EF_Solver1D_TDMA_imp::solve_Ex(Field* phi, Field* Ex)
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


// no source region for electric field
void EF_Solver1D_TDMA_imp::initTDMA_org()
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

// no source region for electric field
void EF_Solver1D_TDMA_imp::solve_TDMA_org(Field* rho, Field* phi)
{
    Field1D* rho1D = static_cast<Field1D*>(rho);
    Field1D* phi1D = static_cast<Field1D*>(phi);
    //> The boundary value can be changed with time
    for(int i = 1; i < nx-1; i++)
    {
        f[i] = -dx_sq * const_ephi0_inv * (*rho1D)(i);
    }
    if(bc_x_left == 1){
        f[0] = bc_e_value[0][0];
    }
    else if(bc_x_left == 2){
        f[0] = bc_e_derivative[0][0];
    }
    if(bc_x_right == 1){
        f[nx-1] = bc_e_value[0][1];
    }
    else if(bc_x_right == 2){
        f[nx-1] = -bc_e_derivative[0][1];
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
