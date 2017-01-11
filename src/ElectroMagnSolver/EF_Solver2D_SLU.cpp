
#include "EF_Solver2D_SLU.h"

#include "ElectroMagn.h"
#include "Field2D.h"
#include <memory>

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;


EF_Solver2D_SLU::EF_Solver2D_SLU(PicParams& params, Grid* grid, SmileiMPI* smpi):
Solver2D(params)
{
    SmileiMPI_Cart2D* smpi2D = static_cast<SmileiMPI_Cart2D*>(smpi);


    dx = params.cell_length[0];
    dy = params.cell_length[1];
    dxy = dx * dy;

    grid2D = static_cast<Grid2D*>(grid);
    //grid2D = new Grid2D(params);
    if(smpi2D->isMaster()){
        //grid2D = static_cast<Grid2D*>(grid);
        initSLU();
    }
    //initSLU();
}


EF_Solver2D_SLU::~EF_Solver2D_SLU()
{
}

void EF_Solver2D_SLU::operator() ( ElectroMagn* fields )
{

}



void EF_Solver2D_SLU::operator() ( ElectroMagn* fields , SmileiMPI* smpi)
{
    SmileiMPI_Cart2D* smpi2D = static_cast<SmileiMPI_Cart2D*>(smpi);
    // Static-cast of the fields
    Field2D* Ex2D = static_cast<Field2D*>(fields->Ex_);
    Field2D* Ey2D = static_cast<Field2D*>(fields->Ey_);


    Field2D* rho2D           = static_cast<Field2D*>(fields->rho_);
    Field2D* rho2D_global    = static_cast<Field2D*>(fields->rho_global);
    Field2D* phi2D_global    = static_cast<Field2D*>(fields->phi_global);
    Field2D* Ex2D_global    = static_cast<Field2D*>(fields->Ex_global);
    Field2D* Ey2D_global    = static_cast<Field2D*>(fields->Ey_global);

    smpi2D->barrier();
    smpi2D->gatherRho(rho2D_global, rho2D);

    if(smpi2D->isMaster()){
        solve_SLU(rho2D_global, phi2D_global);
        //phi2D_global->put_to(0.0);
        solve_Exy(phi2D_global, Ex2D_global, Ey2D_global);
    }

    //Ex2D_global->put_to(0.0);
    //Ey2D_global->put_to(0.0);


    smpi2D->barrier();
    smpi2D->scatterField(Ex2D_global, Ex2D);
    smpi2D->scatterField(Ey2D_global, Ey2D);
}


void EF_Solver2D_SLU::initSLU_test(){

    SuperMatrix A, L, U, B;
    double   *a, *rhs;
    double   s, u, p, e, r, l;
    int      *asub, *xa;
    int      *perm_r; /* row permutations from partial pivoting */
    int      *perm_c; /* column permutation vector */
    int      nrhs, info, i, m, n, nnz, permc_spec;
    superlu_options_t options;
    SuperLUStat_t stat;
    /* Initialize matrix A. */
    m = n = 5;
    nnz = 12;
    if ( !(a = doubleMalloc(nnz)) ) ABORT("Malloc fails for a[].");
    if ( !(asub = intMalloc(nnz)) ) ABORT("Malloc fails for asub[].");
    if ( !(xa = intMalloc(n+1)) ) ABORT("Malloc fails for xa[].");
    s = 19.0; u = 21.0; p = 16.0; e = 5.0; r = 18.0; l = 12.0;
    a[0] = s; a[1] = l; a[2] = l; a[3] = u; a[4] = l; a[5] = l;
    a[6] = u; a[7] = p; a[8] = u; a[9] = e; a[10]= u; a[11]= r;
    asub[0] = 0; asub[1] = 1; asub[2] = 4; asub[3] = 1;
    asub[4] = 2; asub[5] = 4; asub[6] = 0; asub[7] = 2;
    asub[8] = 0; asub[9] = 3; asub[10]= 3; asub[11]= 4;
    xa[0] = 0; xa[1] = 3; xa[2] = 6; xa[3] = 8; xa[4] = 10; xa[5] = 12;

    /* Create matrix A in the format expected by SuperLU. */
    dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);

    /* Create right-hand side matrix B. */
    nrhs = 1;
    if ( !(rhs = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for rhs[].");
    for (i = 0; i < m; ++i) rhs[i] = 1.0;
    dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);

    if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");
    if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");

    /* Set the default input options. */
    set_default_options(&options);
    options.ColPerm = NATURAL;

    /* Initialize the statistics variables. */
    StatInit(&stat);

    /* Solve the linear system. */
    dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

    printf("LU factorization: dgssvx() returns info %d\n", info);


}


void EF_Solver2D_SLU::initSLU(){



    double* val;
    double* b;
    int* row;
    int* col;


    /* Defaults */
    lwork = 0;
    nrhs  = 1;
    equil = YES;
    u     = 1.0;
    trans = NOTRANS;

    set_default_options(&options);
    options.Equil = equil;
    options.DiagPivotThresh = u;
    options.Trans = trans;


    //>>>structure the matrix A
    val = new double[grid2D->ncp*5];
    b = new double[grid2D->ncp];
    row = new int[grid2D->ncp*5];
    col = new int[grid2D->ncp*5];

    for(int i = 0; i < grid2D->ncp*5; i++)
    {
        val[i] = 0.0;
    }


    //unique_ptr<double[]> val(new double[grid2D->ncp*5]);
    //unique_ptr<double[]> b(new double[grid2D->ncp]);
    //unique_ptr<int[]> row(new int[grid2D->ncp*5]);
    //unique_ptr<int[]> col(new int[grid2D->ncp*5]);


    int i,j,k,ii,ll,kk,v,hu,hd,hr,hl,i_ncp,i_nnz,nz_col,i_val;

    nnz=0;
    ii=0;
    v=0;
    nx=grid2D->nx;
    ny=grid2D->ny;


    for(i=0; i<nx; i++)
    {
        for(j=0; j<ny; j++)
        {
            // normal points in the calculation region
            if(grid2D->bndr_global_2D[i][j]==0) {
                hl = grid2D->numcp_global_2D[i][j] - grid2D->numcp_global_2D[i-1][j];
                hr = grid2D->numcp_global_2D[i+1][j] - grid2D->numcp_global_2D[i][j];
                for(k=0; k<grid2D->ncp; k++){
                    b[k]=0.0;
                }
                b[ii]=-4.0; b[ii-hl]=1.0; b[ii-1]=1.0; b[ii+hr]=1.0; b[ii+1]=1.0;
                nnz=nnz+5;
                for(k=0; k<grid2D->ncp; k++){
                    if(b[k] != 0.0){
                        if(v>=grid2D->ncp*5) cout<<"error"<<v<<endl;
                      val[v] = b [k];
                      row[v] = ii;
                      col[v] = k;
                      v++;
                    }
                }
                ii++;
            }

            // Dirchlet boudnary points
            else if(grid2D->bndr_global_2D[i][j]==1) {
                for(k=0; k<grid2D->ncp; k++){
                    b[k]=0.0;
                }
                b[ii]=1.0;
                nnz++;
                for(k=0; k<grid2D->ncp; k++){
                    if(b[k] != 0.0){
                        if(v>=grid2D->ncp*5) cout<<"error"<<v<<endl;
                        val[v] = b [k];
                        row[v] = ii;
                        col[v] = k;
                        v++;
                    }
                }
                ii++;
            }

            // periodic boudnary points
            else if( grid2D->bndr_global_2D[i][j]==8 && j==0) {
                hu = grid2D->numcp_global_2D[i][ny-1] - grid2D->numcp_global_2D[i][j];
                for(k=0; k<grid2D->ncp; k++){
                    b[k]=0.0;
                }
                b[ii]=1.0; b[ii+hu]=-1.0;
                nnz=nnz+2;
                for(k=0; k<grid2D->ncp; k++){
                    if(b[k] != 0.0){
                        if(v>=grid2D->ncp*5) cout<<"error"<<v<<endl;
                        val[v] = b [k];
                        row[v] = ii;
                        col[v] = k;
                        v++;
                    }
                }
                ii++;
            }

            // periodic boudnary points
            else if ( grid2D->bndr_global_2D[i][j] == 8 && j == ny-1 ) {
                hl = grid2D->numcp_global_2D[i][j] - grid2D->numcp_global_2D[i-1][j];
                hr = grid2D->numcp_global_2D[i+1][j] - grid2D->numcp_global_2D[i][j];
                hd = grid2D->numcp_global_2D[i][ny-1] - grid2D->numcp_global_2D[i][1];
                for(k=0; k<grid2D->ncp; k++){
                b[k]=0.0;
                }
                b[ii]=-4.0; b[ii-hl]=1.0; b[ii-1]=1.0; b[ii+hr]=1.0; b[ii-hd]=1.0;
                nnz=nnz+5;
                for(k=0; k<grid2D->ncp; k++){
                if(b[k] != 0.0){
                    if(v>=grid2D->ncp*5) cout<<"error"<<v<<endl;
                  val[v] = b [k];
                  row[v] = ii;
                  col[v] = k;
                  v++;
                }
                }
                ii++;
            }

            // periodic boudnary points
            else if( grid2D->bndr_global_2D[i][j]==8 && i==0) {
                hr = grid2D->numcp_global_2D[nx-1][j] - grid2D->numcp_global_2D[i][j];
                for(k=0; k<grid2D->ncp; k++){
                    b[k]=0.0;
                }
                b[ii]=1.0; b[ii+hr]=-1.0;
                nnz=nnz+2;
                for(k=0; k<grid2D->ncp; k++){
                    if(b[k] != 0.0){
                        if(v>=grid2D->ncp*5) cout<<"error"<<v<<endl;
                        val[v] = b [k];
                        row[v] = ii;
                        col[v] = k;
                        v++;
                    }
                }
                ii++;
            }

            // periodic boudnary points
            else if ( grid2D->bndr_global_2D[i][j] == 8 && i == nx-1 ) {
                hl = grid2D->numcp_global_2D[i][j] - grid2D->numcp_global_2D[i-1][j];
                hr = grid2D->numcp_global_2D[i][j] - grid2D->numcp_global_2D[1][j];
                for(k=0; k<grid2D->ncp; k++){
                b[k]=0.0;
                }
                b[ii]=-4.0; b[ii-hl]=1.0; b[ii-1]=1.0; b[ii-hr]=1.0; b[ii+1]=1.0;
                nnz=nnz+5;
                for(k=0; k<grid2D->ncp; k++){
                if(b[k] != 0.0){
                    if(v>=grid2D->ncp*5) cout<<"error"<<v<<endl;
                  val[v] = b [k];
                  row[v] = ii;
                  col[v] = k;
                  v++;
                }
                }
                ii++;
            }


        }
    }

//for(int i=0; i<grid2D->ncp*5;i++) cout<<i<<" "<<val[i]<<endl;


    //>>>convert the temp "val wor col" to A (compressed column format, i.e. Harwell-Boeing format)
    //a = new double[nnz];
    //asub = new int[nnz];
    //xa = new int[grid2D->ncp+1];
    if ( !(a = doubleMalloc(nnz)) ) ABORT("Malloc fails for a[].");
    if ( !(asub = intMalloc(nnz)) ) ABORT("Malloc fails for asub[].");
    if ( !(xa = intMalloc(grid2D->ncp+1)) ) ABORT("Malloc fails for xa[].");



    i_ncp=0;
    i_nnz=0;
    nz_col=0;
    i_val=0;

    //>>>scan colomn (total ncp column)
    for( i_ncp=0; i_ncp<grid2D->ncp; i_ncp++){
      //>>>is the first not zero element in the current column
      nz_col=0;
      for( i_nnz=0; i_nnz<nnz; i_nnz++){
        if(col[i_nnz] == i_ncp){
          a[i_val] = val[i_nnz];
          asub[i_val] = row[i_nnz];
          if ( nz_col == 0) {
            xa[i_ncp] = i_val;
            nz_col = 1;
          }
          i_val++;
        }
      }
    }//>>>end scan
    cout<<"adrrr "<<val[0]<<" "<<val[grid2D->ncp*5-2]<<" "<<val[grid2D->ncp*5-1]<<endl;
    delete[] val;
    delete[] b;
    delete[] row;
    delete[] col;

    xa[grid2D->ncp]=nnz;
    cout<<"ncp: "<<grid2D->ncp<<" "<<nnz<<endl;

    m = grid2D->ncp;
    n = grid2D->ncp;
    nrhs = 1;
    //rhsb = new double[m * nrhs];
    //rhsx = new double[m * nrhs];
    if ( !(rhsb = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for rhsb[].");
    if ( !(rhsx = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for rhsx[].");

    dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&B, m, nrhs, rhsb, m, SLU_DN, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&X, m, nrhs, rhsx, m, SLU_DN, SLU_D, SLU_GE);
    xact = doubleMalloc(n * nrhs);
    ldx = n;
    dGenXtrue(n, nrhs, xact, ldx);
    dFillRHS(trans, nrhs, xact, ldx, &A, &B);

    if ( !(etree = intMalloc(n)) ) ABORT("Malloc fails for etree[].");
    if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");
    if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(R = (double *) SUPERLU_MALLOC(A.nrow * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for R[].");
    if ( !(C = (double *) SUPERLU_MALLOC(A.ncol * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for C[].");
    if ( !(ferr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for ferr[].");
    if ( !(berr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for berr[].");

    /* Initialize the statistics variables. */
    StatInit(&stat);

    /* ONLY PERFORM THE LU DECOMPOSITION */
    B.ncol = 0;  /* Indicate not to solve the system */
    cout<<"LU fact "<<endl;
    //dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

    dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
           &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
           &mem_usage, &stat, &info);

    printf("LU factorization: dgssvx() returns info %d\n", info);
    StatFree(&stat);


}


void EF_Solver2D_SLU::solve_SLU(Field* rho, Field* phi){

    Field2D* rho2D = static_cast<Field2D*>(rho);
    Field2D* phi2D = static_cast<Field2D*>(phi);

    //>>>convert Field2D rho to SuperLU right hand side matrix
    int ii;
    ii = 0;
    for ( int i=0; i<nx; i++)
    {
      for ( int j=0; j<ny; j++) {
        if ( grid2D->bndr_global_2D[i][j] == 0 ) {
          rhsb[ii] = - dxy * const_ephi0_inv * (*rho2D)(i,j);
          ii++;
        }
        else if ( grid2D->bndr_global_2D[i][j] == 1) {
          rhsb[ii] = grid2D->bndrVal_global_2D[i][j];
          ii++;
        }
        else if ( grid2D->bndr_global_2D[i][j] == 8 && ( j == 0 || i == 0 )) {
          rhsb[ii] = 0.0;
          ii++;
        }
        else if ( grid2D->bndr_global_2D[i][j] == 8 && ( j == ny-1 || i == nx-1 )) {
          rhsb[ii] = - dxy * const_ephi0_inv * (*rho2D)(i,j);
          ii++;
        }
        else {
        }

      }
    }//>>>end convert



    /* ------------------------------------------------------------
       NOW WE SOLVE THE LINEAR SYSTEM USING THE FACTORED FORM OF A.
       ------------------------------------------------------------*/
    options.Fact = FACTORED; /* Indicate the factored form of A is supplied. */
    B.ncol = nrhs;  /* Set the number of right-hand side */

    /* Initialize the statistics variables. */
    StatInit(&stat);
    dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
           &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
           &mem_usage, &stat, &info);

    //printf("Triangular solve: dgssvx() returns info %d\n", info);

    //>>>convert SuperLU solution X to Field2D phi
    ii=0;
    for ( int i=0; i<nx; i++)
      for ( int j=0; j<ny; j++) {
        if ( grid2D->bndr_global_2D[i][j] == 0 || grid2D->bndr_global_2D[i][j] == 1
        || grid2D->bndr_global_2D[i][j] == 2 || grid2D->bndr_global_2D[i][j] == 8) {
          (*phi2D)(i,j) = rhsx[ii];
          ii++;
        }

        if(grid2D->bndr_global_2D[i][j] == 5) {
            (*phi2D)(i,j) = grid2D->bndrVal_global_2D[i][j];
        }

      }//>>>end convert

    StatFree(&stat);


}


void EF_Solver2D_SLU::solve_Exy(Field* phi, Field* Ex, Field* Ey){

    Field2D* phi2D = static_cast<Field2D*>(phi);
    Field2D* Ex2D = static_cast<Field2D*>(Ex);
    Field2D* Ey2D = static_cast<Field2D*>(Ey);


    for(int j = 0; j < ny; j++)
    {
        for(int i = 1; i < nx-1; i++)
        {
            (*Ex2D)(i,j) = - ((*phi2D)(i+1,j) - (*phi2D)(i-1,j)) / (2.0*dx);
        }

        (*Ex2D)(0,j) = -(-3.0 * (*phi2D)(0,j) + 4.0 * (*phi2D)(1,j) - (*phi2D)(2,j)) / (2.0*dx);
        (*Ex2D)(nx-1,j) = -((*phi2D)(nx-3,j) - 4.0 * (*phi2D)(nx-2,j) + 3.0 * (*phi2D)(nx-1,j)) / (2.0*dx);
    }


    for(int i = 0; i < nx; i++)
    {
        for(int j = 1; j < ny-1; j++)
        {
            (*Ey2D)(i,j) = - ((*phi2D)(i,j+1) - (*phi2D)(i,j-1)) / (2.0*dy);
        }

        (*Ey2D)(i,0) = - (-3.0 * (*phi2D)(i,0) + 4.0 * (*phi2D)(i,1) - (*phi2D)(i,2)) / (2.0*dy);
        (*Ey2D)(i,ny-1) = - ((*phi2D)(i,ny-3) - 4.0 * (*phi2D)(i,ny-2) + 3.0 * (*phi2D)(i,ny-1)) / (2.0*dy);
    }


}




void EF_Solver2D_SLU::finishSLU(){

    SUPERLU_FREE (rhsb);
    SUPERLU_FREE (rhsx);
    SUPERLU_FREE (xact);
    SUPERLU_FREE (etree);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    SUPERLU_FREE (R);
    SUPERLU_FREE (C);
    SUPERLU_FREE (ferr);
    SUPERLU_FREE (berr);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperMatrix_Store(&X);
    if ( lwork == 0 ) {
        Destroy_SuperNode_Matrix(&L);
        Destroy_CompCol_Matrix(&U);
    } else if ( lwork > 0 ) {
        SUPERLU_FREE(work);
    }

}
