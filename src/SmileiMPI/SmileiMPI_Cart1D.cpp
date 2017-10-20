#include "SmileiMPI_Cart1D.h"

#include <cmath>
#include <cstring>
#include <algorithm>
#include <string>
#include <fstream>

#include <mpi.h>
#include "Species.h"

#include "ElectroMagn.h"
#include "Field1D.h"

#include "Tools.h"

#include "Field1D.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// SmileiMPI_Cart1D: creator for Smilei MPI environment in 1D cartesian
// ---------------------------------------------------------------------------------------------------------------------
SmileiMPI_Cart1D::SmileiMPI_Cart1D( int* argc, char*** argv )
: SmileiMPI( argc, argv )
{
}


// ---------------------------------------------------------------------------------------------------------------------
// SmileiMPI_Cart1D: creator for Smilei MPI environment in 1D cartesian
// ---------------------------------------------------------------------------------------------------------------------
SmileiMPI_Cart1D::SmileiMPI_Cart1D( SmileiMPI* smpi)
: SmileiMPI( smpi )
{
    ndims_ = 1;
    number_of_procs  = new int[ndims_];
    coords_  = new int[ndims_];
    periods_  = new int[ndims_];
    reorder_ = 0;

    nbNeighbors_ = 2;

    for (int i=0 ; i<ndims_ ; i++) periods_[i] = 0;
    for (int i=0 ; i<ndims_ ; i++) coords_[i] = 0;
    for (int i=0 ; i<ndims_ ; i++) number_of_procs[i] = 1;

    for (int iDim=0 ; iDim<ndims_ ; iDim++) {
        for (int i=0 ; i<nbNeighbors_ ; i++) {
            neighbor_[iDim][i] = MPI_PROC_NULL;
            buff_index_send[i].resize(0);
            buff_index_recv_sz[i] = 0;
        }
    }

}


// ---------------------------------------------------------------------------------------------------------------------
// SmileiMPI_Cart1D: creator for Smilei MPI environment in 1D cartesian
// ---------------------------------------------------------------------------------------------------------------------
SmileiMPI_Cart1D::~SmileiMPI_Cart1D()
{
    delete [] number_of_procs;
    delete [] periods_;
    delete [] coords_;

    if ( SMILEI_COMM_1D != MPI_COMM_NULL) MPI_Comm_free(&SMILEI_COMM_1D);

}
// ---------------------------------------------------------------------------------------------------------------------
// SmileiMPI_Cart1D: create the topology for Smilei MPI environment in 1D cartesian
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI_Cart1D::createTopology(PicParams& params)
{

    for (unsigned int i=0 ; i<params.nDim_field ; i++) {
        params.n_space_global[i] = round(params.sim_length[i]/params.cell_length[i]);
        MESSAGE("Total number of cells in direction " << i << ": " << params.n_space_global[i]);
    }

    number_of_procs[0] = smilei_sz;

    MESSAGE("MPI Domain decomposition : " << smilei_sz);

    // Geometry periodic in x
    if ( (params.bc_em_type_x[0]=="periodic") || (params.bc_em_type_x[1]=="periodic") ) {
        periods_[0] = 1;
        MESSAGE("Periodic geometry in x-direction");
    }

    MPI_Cart_create( SMILEI_COMM_WORLD, ndims_, number_of_procs, periods_, reorder_, &SMILEI_COMM_1D );
    MPI_Cart_coords( SMILEI_COMM_1D, smilei_rk, ndims_, coords_ );


    for (int iDim=0 ; iDim<ndims_ ; iDim++) {
        MPI_Cart_shift( SMILEI_COMM_1D, iDim, 1, &(neighbor_[iDim][0]), &(neighbor_[iDim][1]) );
        //DEBUG(3, smilei_rk, "Neighbors of process in direction " << iDim << " : " << neighbor_[iDim][0] << " ; " << neighbor_[iDim][1] << " Null :" << MPI_PROC_NULL );
    }


    for (unsigned int i=0 ; i<params.nDim_field ; i++) {

        n_space_global[i] = params.n_space_global[i];
        //> I think this if should be removed
        if ( 1 ) {

            params.n_space[i] = params.n_space_global[i] / number_of_procs[i];
            cell_starting_global_index[i] = coords_[i]*(params.n_space_global[i] / number_of_procs[i]);

            if ( number_of_procs[i]*params.n_space[i] != params.n_space_global[i] ) {
                // Correction on the last MPI process of the direction to use the wished number of cells
                if (coords_[i]==number_of_procs[i]-1) {
                    params.n_space[i] = params.n_space_global[i] - params.n_space[i]*(number_of_procs[i]-1);
                }
            }
        }


        oversize[i] = params.oversize[i] = params.interpolation_order + (params.exchange_particles_each-1);
        if ( params.n_space[i] <= 2*oversize[i] ) {
            WARNING ( "Increase space resolution or reduce number of MPI process in direction " << i );
        }

        // min/max_local : describe local domain in which particles cat be moved
        //                 different from domain on which E, B, J are defined
        min_local[i] = (cell_starting_global_index[i]                  )*params.cell_length[i];
        max_local[i] = (cell_starting_global_index[i]+params.n_space[i])*params.cell_length[i];

        cell_starting_global_index[i] -= params.oversize[i];

    }

    //>>>calculate nspace_global_gather to gather and scatter the Field and Grid data
    dims_global_gather[0] = n_space_global[0]+(1+2*params.oversize[0])*number_of_procs[0];

    field_global_gather = new double[dims_global_gather[0]];

    dims_gather         = new int[smilei_sz];
    dims_gather_temp    = new int[smilei_sz];
    for (unsigned int i=0;i< smilei_sz ; i++)
    {
    	if(i==smilei_rk){
    		dims_gather_temp[i]     = params.n_space[0] + 1 + 2*params.oversize[0];
            //cout<<"dims_gather "<<params.n_space[0]<<" "<<params.number_of_procs[0]<<endl;
    	}
    	else {
    		dims_gather_temp[i]     = 0;
    	}

    }

    MPI_Allreduce(dims_gather_temp, dims_gather, smilei_sz, MPI_INT,MPI_SUM, SMILEI_COMM_1D);

    recv_disp.resize(smilei_sz);
    recv_cnt.resize(smilei_sz);
    send_disp.resize(smilei_sz);
    send_cnt.resize(smilei_sz);

    for(int i = 0; i < smilei_sz; i++)
    {
        recv_cnt[i] = dims_gather[i];
        send_cnt[i] = dims_gather[i];
        if(i == 0){
            recv_disp[i] = 0;
            send_disp[i] = 0;
        }
        else{
            recv_disp[i] = recv_disp[i-1] + recv_cnt[i-1];
            send_disp[i] = send_disp[i-1] + send_cnt[i-1];
        }
    }


    //DEBUG(3, smilei_rk, "n_space = " << params.n_space[0] );


    // -------------------------------------------------------
    // Compute & store the ranks of processes dealing with the
    // corner of the simulation box
    // -------------------------------------------------------

    extrem_ranks[0][0] = 0;
    int rank_min =  0;
    if (coords_[0] == 0) {
        rank_min = smilei_rk;
    }
    MPI_Allreduce(&rank_min, &extrem_ranks[0][0], 1, MPI_INT, MPI_SUM, SMILEI_COMM_1D);
    extrem_ranks[0][1] = 0;
    int rank_max = 0;
    if (coords_[0]==number_of_procs[0]-1) {
        rank_max = smilei_rk;
    }
    MPI_Allreduce(&rank_max, &extrem_ranks[0][1], 1, MPI_INT, MPI_SUM, SMILEI_COMM_1D);


}

void SmileiMPI_Cart1D::exchangeParticles(Species* species, int ispec, PicParams& params,int tnum, int iDim)
{

    Particles &cuParticles = species->particles;
    std::vector<int>* cubmin = &species->bmin;
    std::vector<int>* cubmax = &species->bmax;

    MPI_Status Stat;
    int n_particles;
    int tid;
    int tmp = 0;
    int k=0;
    int i,ii, iPart;
    int n_part_recv, n_part_send;

    // ------------------------------------------------------------------------------
    // Build lists of indexes of particle to exchange per neighbor
    // Computed from indexes_of_particles_to_exchange computed during particles' BC
    // ------------------------------------------------------------------------------

    std::vector< std::vector<int> >* indexes_of_particles_to_exchange_per_thd = &species->indexes_of_particles_to_exchange_per_thd;
    //std::vector<int>                 indexes_of_particles_to_exchange;

#pragma omp single
    {
        indexes_of_particles_to_exchange.clear();
    }
#pragma omp barrier

    for (tid=0 ; tid < tnum ; tid++){
        tmp += ((*indexes_of_particles_to_exchange_per_thd)[tid]).size(); //Compute the position where to start copying
    }

    if (tnum == (int)indexes_of_particles_to_exchange_per_thd->size()-1){ //If last thread
        indexes_of_particles_to_exchange.resize( tmp + ((*indexes_of_particles_to_exchange_per_thd)[tnum]).size());
    }
#pragma omp barrier
    //Copy the list per_thread to the global list
    //One thread at a time (works)
#pragma omp master
    {
        for (tid=0 ; tid < (int)indexes_of_particles_to_exchange_per_thd->size() ; tid++) {
            memcpy(&indexes_of_particles_to_exchange[k], &((*indexes_of_particles_to_exchange_per_thd)[tid])[0],((*indexes_of_particles_to_exchange_per_thd)[tid]).size()*sizeof(int));
            k += ((*indexes_of_particles_to_exchange_per_thd)[tid]).size();
        }
        // All threads together (doesn't work)
        /*if (((*indexes_of_particles_to_exchange_per_thd)[tnum]).size() > 0){
         //cout << "tmp = "<<tmp << endl;
         //cout << "tnum = "<< tnum << endl;
         memcpy(&indexes_of_particles_to_exchange[tmp], &((*indexes_of_particles_to_exchange_per_thd)[tnum])[0],((*indexes_of_particles_to_exchange_per_thd)[tnum]).size()*sizeof(int));
         }*/
        //#pragma omp master
        //{
        sort( indexes_of_particles_to_exchange.begin(), indexes_of_particles_to_exchange.end() );

        n_part_send = indexes_of_particles_to_exchange.size();


        for (i=0 ; i<n_part_send ; i++) {
            iPart = indexes_of_particles_to_exchange[i];
            if      ( cuParticles.position(0,iPart) < min_local[0]) {
                buff_index_send[0].push_back( indexes_of_particles_to_exchange[i] );
            }
            else if ( cuParticles.position(0,iPart) >= max_local[0]) {
                buff_index_send[1].push_back( indexes_of_particles_to_exchange[i] );
            }
        } // END for iPart = f(i)

        Particles partVectorSend[1][2];
        partVectorSend[0][0].initialize(0, params, ispec);
        partVectorSend[0][1].initialize(0, params, ispec);
        Particles partVectorRecv[1][2];
        partVectorRecv[0][0].initialize(0, params, ispec);
        partVectorRecv[0][1].initialize(0, params, ispec);

        /********************************************************************************/
        // Exchange particles
        /********************************************************************************/
        // Loop over neighbors in a direction

        // iDim = 0
        // Send to neighbor_[0][iNeighbor] / Recv from neighbor_[0][(iNeighbor+1)%2] :
        // MPI_COMM_SIZE = 2 :  neighbor_[0][0]  |  Current process  |  neighbor_[0][1]
        // Rank = 0 : iNeighbor = 0 : neighbor_[0][0] = NONE : neighbor_[0][(0+1)%2 = 1
        //            iNeighbor = 1 : neighbor_[0][1] = 1    : neighbor_[0][(1+1)%2 = NONE
        // Rank = 1 : iNeighbor = 0 : neighbor_[0][0] = 0    : neighbor_[0][(0+1)%2 = NONE
        //            iNeighbor = 1 : neighbor_[0][1] = NONE : neighbor_[0][(1+1)%2 = 0

        ///********************************************************************************/
        //// Exchange number of particles to exchange to establish or not a communication
        ///********************************************************************************/

        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            n_part_send = buff_index_send[iNeighbor].size();
            if ( (neighbor_[0][0]!=MPI_PROC_NULL) && (neighbor_[0][1]!=MPI_PROC_NULL) ) {
                //Send-receive
                MPI_Sendrecv( &n_part_send, 1, MPI_INT, neighbor_[0][iNeighbor], 0, &buff_index_recv_sz[(iNeighbor+1)%2], 1, MPI_INT, neighbor_[0][(iNeighbor+1)%2], 0, SMILEI_COMM_1D,&Stat);
            } else if (neighbor_[0][iNeighbor]!=MPI_PROC_NULL) {
                //Send
                MPI_Send( &n_part_send, 1, MPI_INT, neighbor_[0][iNeighbor], 0, SMILEI_COMM_1D);
            } else if (neighbor_[0][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
                //Receive
                MPI_Recv( &buff_index_recv_sz[(iNeighbor+1)%2], 1, MPI_INT, neighbor_[0][(iNeighbor+1)%2], 0, SMILEI_COMM_1D, &Stat);
            }
        }

        /********************************************************************************/
        // Define buffers to exchange buff_index_send[iNeighbor].size();
        /********************************************************************************/
        //! \todo Define this as a main parameter for the code so that it needs not be defined all the time

        /********************************************************************************/
        // Proceed to effective Particles' communications
        /********************************************************************************/

        //Number of properties per particles = nDim_Particles + 3 + 6 + 1 + 1
        int nbrOfProp( 12 );
        MPI_Datatype typePartSend, typePartRecv;

        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            n_part_send = buff_index_send[iNeighbor].size();
            n_part_recv = buff_index_recv_sz[(iNeighbor+1)%2];
            if ( (neighbor_[0][0]!=MPI_PROC_NULL) && (neighbor_[0][1]!=MPI_PROC_NULL) && (n_part_send!=0) && (n_part_recv!=0) ) {
                //Send-receive
                double x_max = params.cell_length[0]*( params.n_space_global[0] );
                for (int iPart=0 ; iPart<n_part_send ; iPart++) {
                    // Enabled periodicity in X
                    if ( ( iNeighbor==0 ) &&  (coords_[0] == 0 ) &&( cuParticles.position(0,buff_index_send[iNeighbor][iPart]) < 0. ) ) {
                        cuParticles.position(0,buff_index_send[iNeighbor][iPart])     += x_max;
                    }
                    else if ( ( iNeighbor==1 ) &&  (coords_[0] == number_of_procs[0]-1 ) && ( cuParticles.position(0,buff_index_send[iNeighbor][iPart]) >= x_max ) ) {
                        cuParticles.position(0,buff_index_send[iNeighbor][iPart])     -= x_max;
                    }
                    cuParticles.cp_particle(buff_index_send[iNeighbor][iPart], partVectorSend[0][iNeighbor]);
                }

                typePartSend = createMPIparticles( &(partVectorSend[0][iNeighbor]), nbrOfProp );

                partVectorRecv[0][(iNeighbor+1)%2].initialize( n_part_recv, params, ispec );
                typePartRecv = createMPIparticles( &(partVectorRecv[0][(iNeighbor+1)%2]), nbrOfProp );

                MPI_Sendrecv(&((partVectorSend[0][iNeighbor      ]).position(0,0)),        1, typePartSend, neighbor_[0][iNeighbor      ], 0,
                             &((partVectorRecv[0][(iNeighbor+1)%2]).position(0,0)),        1, typePartRecv, neighbor_[0][(iNeighbor+1)%2], 0, SMILEI_COMM_1D, &Stat);
                MPI_Type_free( &typePartSend );
                MPI_Type_free( &typePartRecv );

            } else if ( (neighbor_[0][iNeighbor]!=MPI_PROC_NULL) && (n_part_send!=0) ) {
                //Send
                partVectorSend[0][iNeighbor].reserve(n_part_send, 1);
                double x_max = params.cell_length[0]*( params.n_space_global[0] );
                for (int iPart=0 ; iPart<n_part_send ; iPart++) {
                    // Enabled periodicity in X
                    if ( ( iNeighbor==0 ) &&  (coords_[0] == 0 ) &&( cuParticles.position(0,buff_index_send[iNeighbor][iPart]) < 0. ) ) {
                        cuParticles.position(0,buff_index_send[iNeighbor][iPart])     += x_max;
                    }
                    else if ( ( iNeighbor==1 ) &&  (coords_[0] == number_of_procs[0]-1 ) && ( cuParticles.position(0,buff_index_send[iNeighbor][iPart]) >= x_max ) ) {
                        cuParticles.position(0,buff_index_send[iNeighbor][iPart])     -= x_max;
                    }
                    cuParticles.cp_particle(buff_index_send[iNeighbor][iPart], partVectorSend[0][iNeighbor]);
                }
                typePartSend = createMPIparticles( &(partVectorSend[0][iNeighbor]), nbrOfProp );
                MPI_Send( &((partVectorSend[0][iNeighbor]).position(0,0)), 1, typePartSend, neighbor_[0][iNeighbor], 0, SMILEI_COMM_1D);
                MPI_Type_free( &typePartSend );

            } else if ( (neighbor_[0][(iNeighbor+1)%2]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
                //Receive
                partVectorRecv[0][(iNeighbor+1)%2].initialize( buff_index_recv_sz[(iNeighbor+1)%2], params, ispec );
                typePartRecv = createMPIparticles( &(partVectorRecv[0][(iNeighbor+1)%2]), nbrOfProp );
                MPI_Recv( &((partVectorRecv[0][(iNeighbor+1)%2]).position(0,0)), 1, typePartRecv,  neighbor_[0][(iNeighbor+1)%2], 0, SMILEI_COMM_1D, &Stat );
                MPI_Type_free( &typePartRecv );

            }
        }

        /********************************************************************************/
        // Delete Particles included in buff_send/buff_recv
        /********************************************************************************/
        // Push lost particles at the end of bins
        //! \todo For loop on bins, can use openMP here.
        for (unsigned int ibin = 0 ; ibin < (*cubmax).size() ; ibin++ ) {
            //        DEBUG(ibin << " bounds " << (*cubmin)[ibin] << " " << (*cubmax)[ibin]);
            ii = indexes_of_particles_to_exchange.size()-1;
            if (ii >= 0) { // Push lost particles to the end of the bin
                iPart = indexes_of_particles_to_exchange[ii];
                while (iPart >= (*cubmax)[ibin] && ii > 0) {
                    ii--;
                    iPart = indexes_of_particles_to_exchange[ii];
                }
                while (iPart == (*cubmax)[ibin]-1 && iPart >= (*cubmin)[ibin] && ii > 0) {
                    (*cubmax)[ibin]--;
                    ii--;
                    iPart = indexes_of_particles_to_exchange[ii];
                }
                while (iPart >= (*cubmin)[ibin] && ii > 0) {
                    cuParticles.overwrite_part1D((*cubmax)[ibin]-1, iPart );
                    (*cubmax)[ibin]--;
                    ii--;
                    iPart = indexes_of_particles_to_exchange[ii];
                }
                if (iPart >= (*cubmin)[ibin] && iPart < (*cubmax)[ibin]) { //On traite la dernière particule (qui peut aussi etre la premiere)
                    cuParticles.overwrite_part1D((*cubmax)[ibin]-1, iPart );
                    (*cubmax)[ibin]--;
                }
            }
        }
        //Shift the bins in memory
        //Warning: this loop must be executed sequentially. Do not use openMP here.
        for (int unsigned ibin = 1 ; ibin < (*cubmax).size() ; ibin++ ) { //First bin don't need to be shifted
            ii = (*cubmin)[ibin]-(*cubmax)[ibin-1]; // Shift the bin in memory by ii slots.
            iPart = min(ii,(*cubmax)[ibin]-(*cubmin)[ibin]); // Number of particles we have to shift = min (Nshift, Nparticle in the bin)
            if(iPart > 0) cuParticles.overwrite_part1D((*cubmax)[ibin]-iPart,(*cubmax)[ibin-1],iPart);
            (*cubmax)[ibin] -= ii;
            (*cubmin)[ibin] = (*cubmax)[ibin-1];
        }


        // Delete useless Particles
        //Theoretically, not even necessary to do anything as long you use bmax as the end of your iterator on particles.
        //Nevertheless, you might want to free memory and have the actual number of particles
        //really equal to the size of the vector. So we do:
        cuParticles.erase_particle_trail((*cubmax).back());
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


        /********************************************************************************/
        // Clean lists of indexes of particle to exchange per neighbor
        /********************************************************************************/
        for (int i=0 ; i<nbNeighbors_ ; i++)
            buff_index_send[i].clear();
        /********************************************************************************/
        // Copy newly arrived particles back to the vector
        /********************************************************************************/
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {

            n_part_recv = buff_index_recv_sz[(iNeighbor+1)%2];
            if ( (neighbor_[0][(iNeighbor+1)%2]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
                if (iNeighbor == 0) { // Copy particles coming from the right at the end of Particles Array
                    n_particles = species->getNbrOfParticles();
                    partVectorRecv[0][(iNeighbor+1)%2].cp_particles(n_part_recv, cuParticles,n_particles);
                    (*cubmax)[(*cubmax).size()-1] += n_part_recv ;
                } else {// Copy particles coming from the left at the beginning of Particles Array
                    //New particles are inserted at the end of bin 0 instead of begining to minimize data movement.
                    partVectorRecv[0][(iNeighbor+1)%2].cp_particles(n_part_recv, cuParticles,(*cubmax)[0]);
                    (*cubmax)[0] += n_part_recv ;
                    for (unsigned int ibin=1 ; ibin < (*cubmax).size() ; ibin++ ) {
                        (*cubmax)[ibin] += n_part_recv ;
                        (*cubmin)[ibin] = (*cubmax)[ibin-1] ;
                    }
                }
            }
        }
    } // END omp master
    //DEBUG( 2, "\tProcess " << smilei_rk << " : " << species->getNbrOfParticles() << " Particles of species " << ispec );
} // END exchangeParticles


MPI_Datatype SmileiMPI_Cart1D::createMPIparticles( Particles* particles, int nbrOfProp )
{
    MPI_Datatype typeParticlesMPI;


    int nbrOfProp2(nbrOfProp);
    if (particles->isTestParticles) nbrOfProp2++;

    MPI_Aint address[nbrOfProp2];
    MPI_Get_address( &(particles->position(0,0)), &(address[0]) );
    MPI_Get_address( &(particles->momentum(0,0)), &(address[1]) );
    MPI_Get_address( &(particles->momentum(1,0)), &(address[2]) );
    MPI_Get_address( &(particles->momentum(2,0)), &(address[3]) );

    MPI_Get_address( &(particles->al_imp(0,0)), &(address[4]) );
    MPI_Get_address( &(particles->al_imp(1,0)), &(address[5]) );
    MPI_Get_address( &(particles->al_imp(2,0)), &(address[6]) );
    MPI_Get_address( &(particles->au_imp(0,0)), &(address[7]) );
    MPI_Get_address( &(particles->au_imp(1,0)), &(address[8]) );
    MPI_Get_address( &(particles->au_imp(2,0)), &(address[9]) );

    MPI_Get_address( &(particles->weight(0)),     &(address[10]) );
    MPI_Get_address( &(particles->charge(0)),     &(address[11]) );
    //MPI_Get_address( &(particles.position_old(0,0)), &address[6] )
    if (particles->isTestParticles)
        MPI_Get_address( &(particles->id(0)),     &(address[nbrOfProp2-1]) );

    int nbr_parts[nbrOfProp2];
    MPI_Aint disp[nbrOfProp2];
    MPI_Datatype partDataType[nbrOfProp2];

    for (int i=0 ; i<nbrOfProp2 ; i++)
        nbr_parts[i] = particles->size();
    disp[0] = 0;
    for (int i=1 ; i<nbrOfProp2 ; i++)
        disp[i] = address[i] - address[0];
    for (int i=0 ; i<nbrOfProp2 ; i++)
        partDataType[i] = MPI_DOUBLE;
    partDataType[nbrOfProp-1] = MPI_DOUBLE;          // if charege is int, this should be MPI_SHORT
    if (particles->isTestParticles)
        partDataType[nbrOfProp2-1] = MPI_UNSIGNED;

    MPI_Type_struct( nbrOfProp2, &(nbr_parts[0]), &(disp[0]), &(partDataType[0]), &typeParticlesMPI);
    MPI_Type_commit( &typeParticlesMPI );

    return typeParticlesMPI;
} // END createMPIparticles


void SmileiMPI_Cart1D::sumField( Field* field )
{
    std::vector<unsigned int> n_elem = field->dims_;
    Field1D* f1D =  static_cast<Field1D*>(field);

    // Use a buffer per direction to exchange data before summing
    Field1D buf[ nbNeighbors_ ];
    // Size buffer is 2 oversize (1 inside & 1 outside of the current subdomain)
    std::vector<unsigned int> oversize2 = oversize;
    oversize2[0] *= 2;

    oversize2[0] += 1 + f1D->isDual_[0];
    for (int i=0; i<nbNeighbors_ ; i++)  buf[i].allocateDims( oversize2 );

    // istart store in the first part starting index of data to send, then the starting index of data to write in
    // Send point of vue : istart =           iNeighbor * ( n_elem[0]- 2*oversize[0] ) + (1-iNeighbor)       * ( 0 );
    // Rank = 0 : iNeighbor = 0 : send - neighbor_[0][0] = NONE
    //            iNeighbor = 1 : send - neighbor_[0][1] = 1 / istart = ( n_elem[0]- 2*oversize[0] )
    // Rank = 1 : iNeighbor = 0 : send - neighbor_[0][0] = 0 / istart = 0
    //            iNeighbor = 1 : send - neighbor_[0][1] = NONE
    // Recv point of vue : istart = ( (iNeighbor+1)%2 ) * ( n_elem[0]- 2*oversize[0] ) + (1-(iNeighbor+1)%2) * ( 0 );
    // Rank = 0 : iNeighbor = 0 : recv - neighbor_[0][1] = 1 / istart = ( n_elem[0]- 2*oversize[0] )
    //            iNeighbor = 1 : recv - neighbor_[0][0] = NONE
    // Rank = 1 : iNeighbor = 0 : recv - neighbor_[0][1] = NONE
    //            iNeighbor = 1 : recv - neighbor_[0][0] = 0 / istart = 0
    int istart;

    MPI_Status sstat[2];
    MPI_Status rstat[2];
    MPI_Request srequest[2];
    MPI_Request rrequest[2];
    /********************************************************************************/
    // Send/Recv in a buffer data to sum
    /********************************************************************************/
    // Loop over neighbors in a direction
    // Send to neighbor_[0][iNeighbor] / Recv from neighbor_[0][(iNeighbor+1)%2] :
    // See in exchangeParticles()
    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {

        if (neighbor_[0][iNeighbor]!=MPI_PROC_NULL) {
            istart = iNeighbor * ( n_elem[0]- oversize2[0] ) + (1-iNeighbor) * ( 0 );
            MPI_Isend( &(f1D->data_[istart]), oversize2[0], MPI_DOUBLE, neighbor_[0][iNeighbor], 0, SMILEI_COMM_1D, &(srequest[iNeighbor]) );
            //cout << "SUM : " << smilei_rk << " send " << oversize2[0] << " data to " << neighbor_[0][iNeighbor] << " starting at " << istart << endl;
        } // END of Send

        if (neighbor_[0][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
            istart = ( (iNeighbor+1)%2 ) * ( n_elem[0]- oversize2[0] ) + (1-(iNeighbor+1)%2) * ( 0 );
            MPI_Irecv( &( (buf[(iNeighbor+1)%2]).data_[0] ), oversize2[0], MPI_DOUBLE, neighbor_[0][(iNeighbor+1)%2], 0, SMILEI_COMM_1D, &(rrequest[(iNeighbor+1)%2]) );
            //cout << "SUM : " << smilei_rk << " recv " << oversize2[0] << " data to " << neighbor_[0][(iNeighbor+1)%2] << " starting at " << istart << endl;
        } // END of Recv

    } // END for iNeighbor


    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
        if (neighbor_[0][iNeighbor]!=MPI_PROC_NULL ) {
            MPI_Wait( &(srequest[iNeighbor]), &(sstat[iNeighbor]) );
        }
        if (neighbor_[0][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
            MPI_Wait( &(rrequest[(iNeighbor+1)%2]), &(rstat[(iNeighbor+1)%2]) );
        }
    }


    // Synchro before summing, to not sum with data ever sum
    barrier();
    /********************************************************************************/
    // Sum data on each process, same operation on both side
    /********************************************************************************/
    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
        istart = ( (iNeighbor+1)%2 ) * ( n_elem[0]- oversize2[0] ) + (1-(iNeighbor+1)%2) * ( 0 );
        // Using Receiver point of vue
        if (neighbor_[0][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
            //cout << "SUM : " << smilei_rk << " sum " << oversize2[0] << " data from " << istart << endl;
            for (unsigned int i=0 ; i<oversize2[0] ; i++)
                f1D->data_[istart+i] += (buf[(iNeighbor+1)%2])(i);
        }
    } // END for iNeighbor


} // END sumField


void SmileiMPI_Cart1D::gatherRho( Field* field_global ,Field* field  )
{

    int procs_rk;
    int iGlobal, jGlobal;
    int iGlobal_gather;
    int nx;

    Field1D* f1D =  static_cast<Field1D*>(field);

    //for(int i = 0; i < f1D->dims_[0]; i++)
    //{
    //    (*f1D)(i) = i;
    //}


    Field1D* f1D_global =  static_cast<Field1D*>(field_global);
    nx = f1D_global->dims_[0];
    f1D_global->put_to(0.0);
    MPI_Gatherv(f1D->data_, send_cnt[smilei_rk], MPI_DOUBLE, field_global_gather, &recv_cnt[0], &recv_disp[0], MPI_DOUBLE, 0, SMILEI_COMM_1D);

    for(int iProcs = 0; iProcs < number_of_procs[0]; iProcs++)
    {
        procs_rk = iProcs;
        for(int i = 0; i < dims_gather[procs_rk]; i++)
        {
            iGlobal = iProcs * (dims_gather[0] - 2*oversize[0] -1) + i -oversize[0];
            if(iProcs == 0 && i < oversize[0] || iProcs == number_of_procs[0] -1 && i > dims_gather[procs_rk] - 1 - oversize[0]){
                iGlobal = abs((int)f1D_global->dims_[0] - abs(iGlobal) - 1);
            }

            iGlobal_gather = send_disp[procs_rk] + i;
            //if(iGlobal >= ii || jGlobal >= jj) cout<<"error "<<iGlobal<<" "<<iProcs<<" "<<dims_gather[0]<<" "<<oversize[0]<<endl;

            //> the differance between gatherRho and gatherField exists only here
            f1D_global->data_[iGlobal] += field_global_gather[iGlobal_gather];
            //if(f1D_global->data_[iGlobal] != 0.0) cout<<"ereeee"; //<<f1D_global->data_[iGlobal]<<endl;
        }
    }


    f1D_global->data_[0] += f1D_global->data_[nx-1];
    f1D_global->data_[nx-1] = f1D_global->data_[0];

} // END gatherRho



void SmileiMPI_Cart1D::gatherField( Field* field_global ,Field* field  )
{

    int procs_rk;
    int iGlobal, jGlobal;
    int iGlobal_gather;
    int nx;

    Field1D* f1D =  static_cast<Field1D*>(field);
    Field1D* f1D_global =  static_cast<Field1D*>(field_global);
    nx = f1D_global->dims_[0];
    f1D_global->put_to(0.0);
    MPI_Gatherv(f1D->data_, send_cnt[smilei_rk], MPI_DOUBLE, field_global_gather, &recv_cnt[0], &recv_disp[0], MPI_DOUBLE, 0, SMILEI_COMM_1D);

    for(int iProcs = 0; iProcs < number_of_procs[0]; iProcs++)
    {
        procs_rk = iProcs;
        for(int i = 0; i < dims_gather[procs_rk]; i++)
        {
            iGlobal = iProcs * (dims_gather[0] - 2*oversize[0] -1) + i -oversize[0];
            if(iProcs == 0 && i < oversize[0] || iProcs == number_of_procs[0] -1 && i > dims_gather[procs_rk] - 1 - oversize[0]){
                iGlobal = abs((int)f1D_global->dims_[0] - abs(iGlobal) - 1);
            }

            iGlobal_gather = send_disp[procs_rk] + i;
            //if(iGlobal >= ii || jGlobal >= jj) cout<<"error "<<iGlobal<<" "<<iProcs<<" "<<dims_gather[0]<<" "<<oversize[0]<<endl;

            //> the differance between gatherRho and gatherField exists only here
            f1D_global->data_[iGlobal] = field_global_gather[iGlobal_gather];
            //if(f1D_global->data_[iGlobal] != 0.0) cout<<"ereeee"; //<<f1D_global->data_[iGlobal]<<endl;
        }
    }


    f1D_global->data_[0] += f1D_global->data_[nx-1];
    f1D_global->data_[nx-1] = f1D_global->data_[0];

} // END gatherField



void SmileiMPI_Cart1D::scatterField( Field* field_global ,Field* field )
{

    int procs_rk;
    int iGlobal, jGlobal;
    int iGlobal_gather;

    Field1D* f1D =  static_cast<Field1D*>(field);
    Field1D* f1D_global =  static_cast<Field1D*>(field_global);

    int ii;
    ii=f1D_global->dims_[0];

    iGlobal = 0;



    for(int iProcs = 0; iProcs < number_of_procs[0]; iProcs++)
    {
        procs_rk = iProcs;
        for(int i = 0; i < dims_gather[procs_rk]; i++)
        {
            iGlobal = iProcs * (dims_gather[0] - 2*oversize[0] -1) + i -oversize[0];
            if(iProcs == 0 && i < oversize[0] || iProcs == number_of_procs[0] -1 && i > dims_gather[procs_rk] - 1 - oversize[0]){
                iGlobal = abs((int)f1D_global->dims_[0] - abs(iGlobal) - 1);
            }

            iGlobal_gather = send_disp[procs_rk] + i;
            if(iGlobal >= ii) cout<<"error "<<iGlobal<<" "<<iProcs<<" "<<dims_gather[0]<<" "<<oversize[0]<<endl;

            field_global_gather[iGlobal_gather] = f1D_global->data_[iGlobal];

            //if(f1D_global->data_[iGlobal] != 0.0) cout<<"ereeee"; //<<f1D_global->data_[iGlobal]<<endl;
        }
    }
    MPI_Scatterv(field_global_gather, &send_cnt[0], &send_disp[0], MPI_DOUBLE, f1D->data_, recv_cnt[smilei_rk], MPI_DOUBLE, 0, SMILEI_COMM_1D);

} // END scatterField



void SmileiMPI_Cart1D::gatherVDF( Array4D* array_global, Array4D* array )
{
    int procs_rk;
    int iGlobal, lGlobal;
    int iGlobal_gather;
    int i_gather;

    recv_disp_VDF.resize(smilei_sz);
    recv_cnt_VDF.resize(smilei_sz);
    recv_cnt_VDF_temp.resize(smilei_sz);
    send_cnt_VDF.resize(smilei_sz);
    double *array_global_gather = new double[array_global->globalDims_];

    for (unsigned int i=0;i< smilei_sz ; i++)
    {
    	if(i==smilei_rk){
    		recv_cnt_VDF_temp[i]     = array->globalDims_;
    	}
    	else {
    		dims_gather_temp[i]     = 0;
    	}
    }

    MPI_Allreduce(&recv_cnt_VDF_temp[0], &recv_cnt_VDF[0], smilei_sz, MPI_INT,MPI_SUM, SMILEI_COMM_1D);

    for(int i = 0; i < smilei_sz; i++)
    {
        send_cnt_VDF[i] = recv_cnt_VDF[i];
        if(i == 0){
            recv_disp_VDF[i] = 0;
        }
        else{
            recv_disp_VDF[i] = recv_disp_VDF[i-1] + recv_cnt_VDF[i-1];
        }
    }

    vector< vector<int> > dims_, dims_temp;
    dims_.resize(4);
    dims_temp.resize(4);

    for(int iDim = 0; iDim < dims_.size(); iDim++)
    {
        dims_temp[iDim].resize(smilei_sz);
        dims_[iDim].resize(smilei_sz);
        for (unsigned int i=0;i< smilei_sz ; i++)
        {
            if(i==smilei_rk){
                dims_temp[iDim][i] = array->dims_[iDim];
            }
            else {
                dims_gather_temp[i]     = 0;
            }
        }

        MPI_Allreduce(&dims_temp[iDim][0], &dims_[iDim][0], smilei_sz, MPI_INT,MPI_SUM, SMILEI_COMM_1D);
    }



    array_global->put_to(0.0);
    MPI_Gatherv(array->data_, send_cnt_VDF[smilei_rk], MPI_DOUBLE, array_global_gather, &recv_cnt_VDF[0], &recv_disp_VDF[0], MPI_DOUBLE, 0, SMILEI_COMM_1D);

    iGlobal = 0;
    i_gather = 0;
    for(int iProcs = 0; iProcs < number_of_procs[0]; iProcs++)
    {
        procs_rk = iProcs;
        for(int i = 0; i < dims_[0][iProcs]; i++)
        {
            for(int l = 0; l < dims_[3][iProcs]; l++)
            {
                lGlobal = l;
                if( iGlobal >= (array_global->dims_[0]) || lGlobal >= (array_global->dims_[3])) {
                    MESSAGE("gatherVDF warning: "<< iGlobal << lGlobal);
                }
                else {
                    (*array_global)(iGlobal,0,0,lGlobal) = array_global_gather[i_gather];
                    if( (*array_global)(iGlobal,0,0,lGlobal) > 0.0 ) {
                        //MESSAGE( "VDF: "<< (*array_global)(iGlobal,0,0,lGlobal) <<"  " << iGlobal <<"  " << lGlobal );
                    }
                }

                i_gather++;
            }
            iGlobal++;
        }
    }

    delete [] array_global_gather;

}
