#ifndef SMILEIMPI_CART3D_H
#define SMILEIMPI_CART3D_H

#include <vector>
#include <string>

#include <mpi.h>

#include "SmileiMPI.h"

class Species;
class Grid3D;

//  --------------------------------------------------------------------------------------------------------------------
//! Class SmileiMPI_Cart3D
//  --------------------------------------------------------------------------------------------------------------------
class SmileiMPI_Cart3D : public SmileiMPI {
public:
    friend class SmileiIO_Cart3D;

    //! Create intial MPI environment
    SmileiMPI_Cart3D( int* argc, char*** argv );
    //! Create MPI environment for the data geometry from
    //! \param smpi the initil MPI environment
    SmileiMPI_Cart3D(SmileiMPI *smpi);
    //! Destructor for SmileiMPI
    virtual ~SmileiMPI_Cart3D();

    //! Create MPI communicator
    virtual void createTopology(PicParams& params);
    //! Echanges particles of Species, list of particles comes frome Species::dynamics
    virtual void exchangeParticles(Species* species, int ispec, PicParams& params, int tnum, int iDim);

    //! Create MPI_Datatype to exchange/sum fields on ghost data
    void createType( PicParams& params );

    //! Create MPI_Datatype to exchange all properties of particle in 1 communication
    MPI_Datatype createMPIparticles( Particles* particles, int nbrOfProp );

    //! Basic method to sum a field
    virtual void sumField      ( Field* field );


    //! Return coordinates in the cartesian MPI communicator
    //! \param i direction
    inline int getProcCoord(int i) {
        return coords_[i];
    }
    //! Return number of MPI process in the cartesian MPI communicator
    //! \param i direction
    inline int getNbrOfProcs(int i) {
        return number_of_procs[i];
    }

    //! Identify western MPI process, for boundary condition
    inline bool isWestern() {
        return ((coords_[0]==0)&&(periods_[0]==0));
    }
    //! Identify eastern MPI process, for boundary condition
    inline bool isEastern() {
        return ((coords_[0]==number_of_procs[0]-1)&&(periods_[0]==0));
    }
    //! Identify southern MPI process, for boundary condition
    inline bool isSouthern() {
        return ((coords_[1]==0)&&(periods_[1]==0));
    }
    //! Identify northern MPI process, for boundary condition
    inline bool isNorthern() {
        return ((coords_[1]==number_of_procs[1]-1)&&(periods_[1]==0));
    }

    //! Identify corner MPI ranks (3D, 2 sides)
    int extrem_ranks[2][2];

    void scatterGrid( Grid* grid );
    void gatherRho( Field* field_global ,Field* field  );
    void gatherField( Field* field_global ,Field* field  );
    void scatterField( Field* field_global ,Field* field );


protected:
    //! 3D Cartesian communicator
    MPI_Comm SMILEI_COMM_3D;
    //! Number of dimensions
    int ndims_;
    //! Number of MPI process per direction in the cartesian topology
    int* number_of_procs;

    //! Array of coordinates in the cartesian topology
    int* coords_;
    //! Periodicity of the geometry
    int* periods_;
    //! Reorder MPI rank (not)
    int reorder_;

    //! Number of neighbors per directions (=2)
    int nbNeighbors_;
    //! Id of neighbors, per direction (up to 3), per side (2)
    int neighbor_[3][2];

    //! MPI_Datatype to exchange [ndims_][iDim=0 prim/dial][iDim=1 prim/dial]
    MPI_Datatype ntype_   [3][2][2];
    //! MPI_Datatype to sum [ndims_][iDim=0 prim/dial][iDim=1 prim/dial]
    MPI_Datatype ntypeSum_[2][2][2];
    //! Buffer for buffered communication
    //int bufsize;
    //void *b;


};

#endif
