#include "SmileiMPI.h"

#include <cmath>
#include <cstring>

#include <iostream>
#include <sstream>

#include "PicParams.h"
#include "Tools.h"

#include "ElectroMagn.h"
#include "Field.h"
#include "Species.h"

using namespace std;

SmileiMPI::SmileiMPI( int* argc, char*** argv )
{
    int mpi_provided;

    MPI_Init_thread( argc, argv, MPI_THREAD_FUNNELED, &mpi_provided );
    if (mpi_provided == MPI_THREAD_SINGLE){
        MESSAGE("openMP not supported");
    }

    SMILEI_COMM_WORLD = MPI_COMM_WORLD;
    MPI_Comm_size( SMILEI_COMM_WORLD, &smilei_sz );
    MPI_Comm_rank( SMILEI_COMM_WORLD, &smilei_rk );

    // make the random seed different
    srand( (unsigned)time(NULL) + smilei_rk * smilei_sz * 10.0 );

}

SmileiMPI::SmileiMPI( SmileiMPI *smpi )
{
    SMILEI_COMM_WORLD = smpi->SMILEI_COMM_WORLD;
    MPI_Comm_size( SMILEI_COMM_WORLD, &smilei_sz );
    MPI_Comm_rank( SMILEI_COMM_WORLD, &smilei_rk );

    oversize = smpi->oversize;
    n_space_global = smpi->n_space_global;
    cell_starting_global_index = smpi->cell_starting_global_index;
    min_local = smpi->min_local;
    max_local = smpi->max_local;

}

SmileiMPI::~SmileiMPI()
{
    int status = 0;
    MPI_Finalized( &status );
    if (!status) MPI_Finalize();

}

void SmileiMPI::bcast( string& val )
{
    int charSize=0;
    if (isMaster()) charSize = val.size()+1;
    MPI_Bcast(&charSize, 1, MPI_INT, 0, SMILEI_COMM_WORLD);

    char tmp[charSize];
    strcpy(tmp, val.c_str());
    MPI_Bcast(&tmp, charSize, MPI_CHAR, 0, SMILEI_COMM_WORLD);

    if (!isMaster()) val=tmp;

}

void SmileiMPI::init( PicParams& params )
{
    oversize.resize(params.nDim_field, 0);
    cell_starting_global_index.resize(params.nDim_field, 0);
    min_local.resize(params.nDim_field, 0.);
    max_local.resize(params.nDim_field, 0.);
    n_space_global.resize(params.nDim_field, 0);
    n_space_global_gather.resize(params.nDim_field, 0);
}


void SmileiMPI::sumRho( ElectroMagn* EMfields )
{
    sumField( EMfields->rho_ );

}

void SmileiMPI::sumRhoJ( ElectroMagn* EMfields )
{
    // sum total charge density and currents
    sumField( EMfields->rho_ );
    sumField( EMfields->Jx_ );
    sumField( EMfields->Jy_ );
    sumField( EMfields->Jz_ );

}
void SmileiMPI::sumRhoJs( ElectroMagn* EMfields, int ispec, bool currents )
{
   // sum density and currents for all species
   sumField( EMfields->rho_s[ispec] );
   if(currents){
       sumField( EMfields->Jx_s[ispec] );
       sumField( EMfields->Jy_s[ispec] );
       sumField( EMfields->Jz_s[ispec] );
   }
}


void SmileiMPI::reduceDoubleVector( double* src, double* des, int n )
{
    MPI_Allreduce( src, des, n, MPI_DOUBLE, MPI_SUM, SMILEI_COMM_WORLD );
}


int SmileiMPI::globalNbrParticles(Species* species) {
    int nParticles(0);
    int locNbrParticles(0);
    locNbrParticles = species->getNbrOfParticles();
    MPI_Reduce( &locNbrParticles, &nParticles, 1, MPI_INT, MPI_SUM, 0, SMILEI_COMM_WORLD );
    return nParticles;
}
