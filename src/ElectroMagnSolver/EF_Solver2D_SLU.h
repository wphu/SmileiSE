#ifndef EF_SOLVER2D_SLU_H
#define EF_SOLVER2D_SLU_H

#include "Solver2D.h"
#include "slu_ddefs.h"
#include "Grid2D.h"
#include "Field.h"
#include "SmileiMPI_Cart2D.h"

class ElectroMagn;
//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class EF_Solver2D_SLU : public Solver2D
{

public:
    //! Creator for EF_SOLVER2D_SLU
    EF_Solver2D_SLU(PicParams& params, Grid* grid, SmileiMPI* smpi);
    virtual ~EF_Solver2D_SLU();

    //! Overloading of () operator
    virtual void operator()( ElectroMagn* fields);
    virtual void operator()( ElectroMagn* fields, SmileiMPI* smpi);
    void solve_SLU(Field* rho, Field* phi);
    void finishSLU();
    void initSLU_test();
    void initSLU();
    void solve_Exy(Field* phi, Field* Ex, Field* Ey);
    //>>>SuperLU parameters


    //>>>geometry parameters
    int nx,ny;
    double dx, dy, dxy;

    Grid2D* grid2D;

protected:

    char           equed[1];
    yes_no_t       equil;
    trans_t        trans;
    SuperMatrix    A, L, U;
    SuperMatrix    B, X;
    NCformat       *Astore;
    NCformat       *Ustore;
    SCformat       *Lstore;
    double         *a;
    int            *asub, *xa;
    int            *perm_c; /* column permutation vector */
    int            *perm_r; /* row permutations from partial pivoting */
    int            *etree;
    void           *work;
    int            info, lwork, nrhs, ldx;
    int            m, n, nnz;
    double         *rhsb, *rhsx, *xact;
    double         *R, *C;
    double         *ferr, *berr;
    double         u, rpg, rcond;
    mem_usage_t    mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;
    //>>>end





};//END class

#endif
