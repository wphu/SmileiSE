#include "Collisions1D_ChargeExchange.h"
#include "SmileiMPI.h"
#include "Field2D.h"
#include "H5.h"


#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <sstream>

using namespace std;


// Constructor
Collisions1D_ChargeExchange::Collisions1D_ChargeExchange(PicParams& param, vector<Species*>& vecSpecies, SmileiMPI* smpi,
                       unsigned int n_collisions,
                       vector<unsigned int> species_group1,
                       vector<unsigned int> species_group2,
                       double coulomb_log,
                       bool intra_collisions,
                       int debug_every)
{

    n_collisions    = (n_collisions    );
    species_group1  = (species_group1  );
    species_group2  = (species_group2  );
    coulomb_log     = (coulomb_log     );
    intra_collisions= (intra_collisions);
    debug_every     = (debug_every     );
    start           = (0               );



    // Calculate total number of bins
    int nbins = vecSpecies[0]->bmin.size();
    totbins = nbins;
    //MPI_Allreduce( smpi->isMaster()?MPI_IN_PLACE:&totbins, &totbins, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
    MPI_Reduce( smpi->isMaster()?MPI_IN_PLACE:&totbins, &totbins, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);

    // if debug requested, prepare hdf5 file



}

Collisions1D_ChargeExchange::~Collisions1D_ChargeExchange()
{

}


// Calculates the collisions for a given Collisions1D object
void Collisions1D_ChargeExchange::collide(PicParams& params, vector<Species*>& vecSpecies, int itime)
{

}



double  Collisions1D_ChargeExchange::cross_section(double ke)
{

}
