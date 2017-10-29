#ifndef SMILEIMPI_H
#define SMILEIMPI_H

#include <string>
#include <vector>

#include <mpi.h>

#include "PicParams.h"
#include "Tools.h"
#include "Array4D.h"

class PicParams;
class Particles;
class Species;
class ElectroMagn;
class Field;
class Grid;


//  --------------------------------------------------------------------------------------------------------------------
//! Class SmileiMPI
//  --------------------------------------------------------------------------------------------------------------------
class SmileiMPI {
    friend class SmileiIO;
public:
    //! Create intial MPI environment
    SmileiMPI( int* argc, char*** argv );
    //! Create MPI environment for the data geometry from
    //! \param smpi the initil MPI environment
    SmileiMPI(SmileiMPI *smpi);
    //! Default creator for SmileiMPI
    SmileiMPI() {};
    //! Destructor for SmileiMPI
    virtual ~SmileiMPI();

    //! Initialize geometry members dimension from
    //! \param params Parameters
    //! @see oversize
    //! @see cell_starting_global_index
    //! @see min_local
    //! @see max_local
    //! @see n_space_global
    void init( PicParams& params );


    //! Create MPI communicator
    virtual void createTopology( PicParams& params ) {};
    //! Echanges particles of Species, list of particles comes frome Species::dynamics
    //! See child classes
    virtual void exchangeParticles(Species* species, int ispec, PicParams& params, int tnum, int iDim) {};

    //virtual MPI_Datatype createMPIparticles( Particles* particles, int nbrOfProp ) {MPI_Datatype type ; return type; }
    virtual MPI_Datatype createMPIparticles( Particles* particles, int nbrOfProp ) {return NULL;}

    //! Create MPI_Datatype to exchange/sum fields on ghost data
    //! See child classes
    virtual void createType( PicParams& params ) {};


    //! Sum rho and densities on 2 x oversize[]
    void sumRho( ElectroMagn* EMfields );
    //! Sum rho and all J on the shared domain between processors
    //! 2 x oversize + 1 ( + 1 if direction is dual )
    void sumRhoJ( ElectroMagn* EMfields );
    //! Sum rho_s and all J_s on the shared domain between processors
    void sumRhoJs( ElectroMagn* EMfields, int ispec, bool currents );

    //! Basic method to sum a field, defined in child class
    virtual void sumField      ( Field* field ) {};

    //! Method to identify the rank 0 MPI process
    inline bool isMaster() {
        return (smilei_rk==0);
    }
    //! Method to synchronize MPI process in the current MPI communicator
    inline void barrier() {
        MPI_Barrier( SMILEI_COMM_WORLD );
    }
    //! Return MPI_Comm_rank
    inline int getRank() {
        return smilei_rk;
    }
    //! Return MPI_Comm_size
    inline int getSize() {
        return smilei_sz;
    }
    //! Return global starting (including oversize, ex : rank 0 returns -oversize) index for direction i
    //! \param i direction
    //! @see cell_starting_global_index
    inline int    getCellStartingGlobalIndex(int i) const {
        return cell_starting_global_index[i];
    }
    //! Return real (excluding oversize) min coordinates (ex : rank 0 retourn 0.) for direction i
    //! @see min_local
    inline double getDomainLocalMin(int i) const {
        return min_local[i];
    }
    //! Return real (excluding oversize) max coordinates for direction i
    //! @see max_local
    inline double getDomainLocalMax(int i) const {
        return max_local[i];
    }

    //! Set global starting index for direction i
    //! @see cell_starting_global_index
    inline int&    getCellStartingGlobalIndex(int i)  {
        return cell_starting_global_index[i];
    }
    //! Set real min coordinate for direction i
    //! @see min_local
    inline double& getDomainLocalMin(int i)  {
        return min_local[i];
    }
    //! Set real max coordinate for direction i
    //! @see max_local
    inline double& getDomainLocalMax(int i)  {
        return max_local[i];
    }


    inline void bcast_double(double* buffer, int N, int root_rank)
    {
        MPI_Bcast(buffer, N, MPI_DOUBLE, root_rank, SMILEI_COMM_WORLD);
    }



    //! Temporary storage of particles Id to exchange, merge from per thread storage
    //! A single communication per direction managed by thread master in OpenMP
    //! @see Species::indexes_of_particles_to_exchange_per_thrd
    std::vector<int>                 indexes_of_particles_to_exchange;

    //! Should be pure virtual, see child classes
    virtual bool isEastern(){WARNING("Problem");return false;}
    //! Should be pure virtual, see child classes
    virtual bool isWestern(){WARNING("Problem");return false;}
    //! Should be pure virtual, see child classes
    virtual bool isSouthern(){WARNING("Problem");return false;}
    //! Should be pure virtual, see child classes
    virtual bool isNorthern(){WARNING("Problem");return false;}


    //! MPI process Id in the current communicator
    int smilei_sz;
    //! Number of MPI process in the current communicator
    int smilei_rk;

    int globalNbrParticles(Species* species);

    // Broadcast a string in current communicator
    void bcast( std::string& val );

    //> gatherRho is to gather charge density and species density, need sum guard cell points;
    //> gatherField is to gather Fields like Ex, Bx and so on, not need sum guard cell points
    virtual void scatterGrid( Grid* grid ){};
    virtual void gatherRho( Field* field_global ,Field* field  ){};
    virtual void gatherField( Field* field_global ,Field* field  ){};
    virtual void scatterField( Field* field_global ,Field* field ){};
    virtual void gatherVDF( Array4D* array_global, Array4D* array ){};


    //! Real (exclunding oversize) global number of cells (res_space x sim_length)
    std::vector<unsigned int> n_space_global;
  	//>>> (inclunding oversize) global number of cells (res_space x sim_length) to MPI_GaterV and MPI_scatter
    std::vector<int> n_space_global_gather;

    //> global dimensions (including ghost grid points) of each process stored in ROOT process for gathering
    //> and scattering, only dims_global_gather[0] is used for 1d.
	int dims_global_gather[2];

	//>>>variables for MPI_GaterV and MPI_scatter (grid, rho, phi)
	int* grid_global_gather;
	double* field_global_gather;

	//>local dimensions of each process stored in ROOT process for gathering and scattering
    //>can not use 2d vector, because the memory of date is not contigous for MPI
	int* dims_gather;
	int* dims_gather_temp;

    std::vector <int> recv_cnt, recv_disp;
    std::vector <int> send_cnt, send_disp;
    std::vector <int> recv_cnt_VDF, recv_disp_VDF, recv_cnt_VDF_temp;
        std::vector <int> send_cnt_VDF;

    virtual void reduceDoubleVector( double* src, double* des, int n);

protected:
    //! Global MPI Communicator
    MPI_Comm SMILEI_COMM_WORLD;

    //! Sort particles to exchange per side (2), contains indexes
    //! Reinitialized per direction
    std::vector<int> buff_index_send[2];
    //! buff_index_recv_sz : number of particles to recv per side (2)
    //! Reinitialized per direction
    int buff_index_recv_sz[2];

    //! Size of ghost data (= 2 x oversize + 1 + 1 if dual direction), depend on :
    //!    - projection/interpolation order
    //!    - rate of particles exchange (to implement)
    std::vector<unsigned int> oversize;
    //! cell_starting_global_index : index of 1st cell of local subdomain in the global domain
    //!     - concerns ghost data
    //!     - "- oversize" on rank 0
    std::vector<int> cell_starting_global_index;
    //! "Real" min limit of local domain (ghost data not concerned)
    //!     - "0." on rank 0
    std::vector<double> min_local;
    //! "Real" max limit of local domain (ghost data not concerned)
    std::vector<double> max_local;

    bool isImplicit;
    bool isSameWeight;

};

#endif
