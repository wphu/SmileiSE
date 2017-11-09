#include "SmileiMPI_Cart3D.h"

#include <cmath>
#include <string>
#include <algorithm>

#include <mpi.h>

#include "Species.h"
#include "Grid3D.h"

#include "ElectroMagn.h"
#include "Field3D.h"

#include "Tools.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// SmileiMPI_Cart3D: creator for Smilei MPI environment in 3D cartesian geometry
// ---------------------------------------------------------------------------------------------------------------------
SmileiMPI_Cart3D::SmileiMPI_Cart3D( int* argc, char*** argv )
: SmileiMPI( argc, argv )
{
}


// ---------------------------------------------------------------------------------------------------------------------
// SmileiMPI_Cart3D: creator for Smilei MPI environment in 3D cartesian geometry
// ---------------------------------------------------------------------------------------------------------------------
SmileiMPI_Cart3D::SmileiMPI_Cart3D( SmileiMPI* smpi)
: SmileiMPI( smpi )
{
    ndims_ = 2;
    number_of_procs = new int[ndims_];
    coords_  = new int[ndims_];
    periods_  = new int[ndims_];
    reorder_ = 0;

    nbNeighbors_ = 2; // number of neighbor processes per direction

    for (int i=0 ; i<ndims_ ; i++) periods_[i] = 0;
    for (int i=0 ; i<ndims_ ; i++) coords_[i] = 0;
    for (int i=0 ; i<ndims_ ; i++) number_of_procs[i] = 1;

    for (int iDim=0 ; iDim<ndims_ ; iDim++)
        for (int iNeighbors=0 ; iNeighbors<nbNeighbors_ ; iNeighbors++)
            neighbor_[iDim][iNeighbors] = MPI_PROC_NULL;

    for (int i=0 ; i<nbNeighbors_ ; i++) {
        buff_index_send[i].resize(0);
        buff_index_recv_sz[i] = 0;
    }

}


// ---------------------------------------------------------------------------------------------------------------------
// SmileiMPI_Cart3D: creator for Smilei MPI environment in 3D cartesian geometry
// ---------------------------------------------------------------------------------------------------------------------
SmileiMPI_Cart3D::~SmileiMPI_Cart3D()
{
    for (int ix_isPrim=0 ; ix_isPrim<1 ; ix_isPrim++) {
        for (int iy_isPrim=0 ; iy_isPrim<1 ; iy_isPrim++) {
            MPI_Type_free( &ntype_[0][ix_isPrim][iy_isPrim]); //line
            MPI_Type_free( &ntype_[1][ix_isPrim][iy_isPrim]); // column
            MPI_Type_free( &ntype_[2][ix_isPrim][iy_isPrim]); // several lines

            MPI_Type_free( &ntypeSum_[0][ix_isPrim][iy_isPrim]); //line
            MPI_Type_free( &ntypeSum_[1][ix_isPrim][iy_isPrim]); // column
        }
    }

    delete[] number_of_procs;
    delete[]periods_;
    delete[] coords_;


    if ( SMILEI_COMM_3D != MPI_COMM_NULL) MPI_Comm_free(&SMILEI_COMM_3D);
}
// ---------------------------------------------------------------------------------------------------------------------
// SmileiMPI_Cart3D: create the topology for Smilei MPI environment in 3D cartesian geometry
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI_Cart3D::createTopology(PicParams& params)
{

    //! Oversize domain to exchange less particles
    oversize = params.oversize;
    for (unsigned int i=0 ; i<params.nDim_field ; i++) {
        params.n_space_global[i] = round(params.sim_length[i]/params.cell_length[i]);
        MESSAGE("Total number of cells in direction " << i << ": " << params.n_space_global[i]);
    }

    if (params.number_of_procs[0]!=0) {
        for (unsigned int i=0 ; i<params.nDim_field ; i++)
            number_of_procs[i] = params.number_of_procs[i];
        if (number_of_procs[0]*number_of_procs[1]!=smilei_sz) {
            DEBUG(3,"Domain decomposition specified in the namelist don't match with the number of MPI process");
            DEBUG(3,"\tit will be computed to be as square as possible");
            for (unsigned int i=0 ; i<params.nDim_field ; i++)
                params.number_of_procs[i] = 0;
        }
    }
    if (params.number_of_procs[0]==0) {
        double tmp(0.);
        tmp  = params.res_space[0]*params.sim_length[0] / ( params.res_space[1]*params.sim_length[1] );

        number_of_procs[0] = min( smilei_sz, max(1, (int)sqrt ( (double)smilei_sz*tmp*tmp) ) );
        number_of_procs[1] = (int)(smilei_sz / number_of_procs[0]);

        while ( number_of_procs[0]*number_of_procs[1] != smilei_sz ) {
            if (number_of_procs[0]>=number_of_procs[1] ) {
                number_of_procs[0]++;
                number_of_procs[1] = (int)(smilei_sz / number_of_procs[0]);
            }
            else {
                number_of_procs[1]++;
                number_of_procs[0] = (int)(smilei_sz / number_of_procs[1]);
            }
        }

    }
    // Force configuration of MPI domain decomposition
    //number_of_procs[0] = 1;
    //number_of_procs[1] = 16;
    MESSAGE("MPI Domain decomposition : " << smilei_sz << " = " << number_of_procs[0] << " x " << number_of_procs[1]);
    MESSAGE("em. bound. condit. in xmin & xmax directions: "<< params.bc_em_type_x[0] << ", "<< params.bc_em_type_x[1]);
    MESSAGE("em. bound. condit. in ymin & ymax directions: "<< params.bc_em_type_y[0] << ", "<< params.bc_em_type_y[1]);

    // Geometry periodic in x
    if ( (params.bc_em_type_x[0]=="periodic") || (params.bc_em_type_x[1]=="periodic") ) {
        periods_[0] = 1;
        MESSAGE(1,"applied topology for periodic BCs in x-direction");
    }
    // Geometry periodic in y
    if ( (params.bc_em_type_y[0]=="periodic") || (params.bc_em_type_y[1]=="periodic") ) {
        periods_[1] = 1;
        MESSAGE(2,"applied topology for periodic BCs in y-direction");
    }
    MPI_Cart_create( SMILEI_COMM_WORLD, ndims_, number_of_procs, periods_, reorder_, &SMILEI_COMM_3D );
    MPI_Cart_coords( SMILEI_COMM_3D, smilei_rk, ndims_, coords_ );

    for (int iDim=0 ; iDim<ndims_ ; iDim++) {
        MPI_Cart_shift( SMILEI_COMM_3D, iDim, 1, &(neighbor_[iDim][0]), &(neighbor_[iDim][1]) );
        //DEBUG(3,smilei_rk,"Neighbors of process in direction " << iDim << " : " << neighbor_[iDim][0] << " - " << neighbor_[iDim][1]  );
    }
    for (unsigned int i=0 ; i<params.nDim_field ; i++) {
        n_space_global[i] = params.n_space_global[i];
        //if (i!=0) {
        if (1) {
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
            WARNING ( "Increase space resolution or reduce number of MPI process in direction " << i << " "<< params.n_space[i]);
        }

        // min/max_local : describe local domain in which particles cat be moved
        //                 different from domain on which E, B, J are defined
        min_local[i] = (cell_starting_global_index[i]                  )*params.cell_length[i];
        max_local[i] = (cell_starting_global_index[i]+params.n_space[i])*params.cell_length[i];
        //PMESSAGE( 0, smilei_rk, "min_local / mac_local on " << smilei_rk << " = " << min_local[i] << " / " << max_local[i] << " selon la direction " << i );

        cell_starting_global_index[i] -= params.oversize[i];
    }

    //>>>calculate nspace_global_gather to gather and scatter the Field and Grid data
    dims_global_gather[0]=n_space_global[0]+(1+2*params.oversize[0])*number_of_procs[0];
    dims_global_gather[1]=n_space_global[1]+(1+2*params.oversize[1])*number_of_procs[1];

    grid_global_gather= new int[dims_global_gather[0]*dims_global_gather[1]];
    field_global_gather= new double[dims_global_gather[0]*dims_global_gather[1]];

    dims_gather         = new int[2*smilei_sz];
    dims_gather_temp    = new int[2*smilei_sz];
    for (unsigned int i=0;i< smilei_sz ; i++)
    {
    	if(i==smilei_rk){
    		dims_gather_temp[i*2]     = params.n_space[0] + 1 + 2*params.oversize[0];
    		dims_gather_temp[i*2 + 1] = params.n_space[1] + 1 + 2*params.oversize[1];
            //cout<<"dims_gather "<<params.n_space[0]<<" "<<params.n_space[1]<<" "<<params.number_of_procs[0]<<params.number_of_procs[1]<<endl;
    	}
    	else {
    		dims_gather_temp[i*2]     = 0;
    		dims_gather_temp[i*2 + 1] = 0;
    	}

    }

    MPI_Allreduce(dims_gather_temp, dims_gather, 2*smilei_sz, MPI_INT,MPI_SUM, SMILEI_COMM_3D);

    recv_disp.resize(smilei_sz);
    recv_cnt.resize(smilei_sz);
    send_disp.resize(smilei_sz);
    send_cnt.resize(smilei_sz);
    for(int i = 0; i < smilei_sz; i++)
    {
        recv_cnt[i] = dims_gather[i*2] * dims_gather[i*2 + 1];
        send_cnt[i] = dims_gather[i*2] * dims_gather[i*2 + 1];
        if(i == 0){
            recv_disp[i] = 0;
            send_disp[i] = 0;
        }
        else{
            recv_disp[i] = recv_disp[i-1] + recv_cnt[i-1];
            send_disp[i] = send_disp[i-1] + send_cnt[i-1];
        }
    }


    //DEBUG(10,"n_space / rank " << smilei_rk << " = " << params.n_space[0] << " " << params.n_space[1] );

    extrem_ranks[0][0] = 0;
    int rank_min =  0;
    if ( (coords_[0] == 0) && (coords_[1] == 0) )
        rank_min = smilei_rk;
    MPI_Allreduce(&rank_min, &extrem_ranks[0][0], 1, MPI_INT, MPI_SUM, SMILEI_COMM_3D);

    extrem_ranks[0][1] = 0;
    int rank_max = 0;
    if ( (coords_[0]==0) && (coords_[1]==number_of_procs[1]-1) )
        rank_max = smilei_rk;
    MPI_Allreduce(&rank_max, &extrem_ranks[0][1], 1, MPI_INT, MPI_SUM, SMILEI_COMM_3D);

    extrem_ranks[1][0] = 0;
    rank_max = 0;
    if ( (coords_[1]==0) && (coords_[0]==number_of_procs[0]-1) )
        rank_max = smilei_rk;
    MPI_Allreduce(&rank_max, &extrem_ranks[1][0], 1, MPI_INT, MPI_SUM, SMILEI_COMM_3D);

    extrem_ranks[1][1] = 0;
    rank_max = 0;
    if ( (coords_[0]==number_of_procs[0]-1) && (coords_[1]==number_of_procs[1]-1) )
        rank_max = smilei_rk;
    MPI_Allreduce(&rank_max, &extrem_ranks[1][1], 1, MPI_INT, MPI_SUM, SMILEI_COMM_3D);

}

void SmileiMPI_Cart3D::exchangeParticles(Species* species, int ispec, PicParams& params, int tnum, int iDim)
{
    Particles &cuParticles = species->particles;		//hu// cuPaticles: current Particles
    std::vector<int>* cubmin = &species->bmin;
    std::vector<int>* cubmax = &species->bmax;

    std::vector< std::vector<int> >* indexes_of_particles_to_exchange_per_thd = &species->indexes_of_particles_to_exchange_per_thd;
    std::vector<int>                 indexes_of_particles_to_exchange;

#pragma omp master
    {
        /********************************************************************************/
        // Build lists of indexes of particle to exchange per neighbor
        // Computed from indexes_of_particles_to_exchange computed during particles' BC
        /********************************************************************************/
        indexes_of_particles_to_exchange.clear();

        int tmp = 0;
        for (int tid=0 ; tid < (int)indexes_of_particles_to_exchange_per_thd->size() ; tid++)
            tmp += ((*indexes_of_particles_to_exchange_per_thd)[tid]).size();
        indexes_of_particles_to_exchange.resize( tmp );

        int k=0;
        for (int tid=0 ; tid < (int)indexes_of_particles_to_exchange_per_thd->size() ; tid++) {
            for (int ipart = 0 ; ipart < (int) ((*indexes_of_particles_to_exchange_per_thd)[tid]).size() ; ipart++ ) {
                indexes_of_particles_to_exchange[k] =  (*indexes_of_particles_to_exchange_per_thd)[tid][ipart] ;
                k++;
            }
            ((*indexes_of_particles_to_exchange_per_thd))[tid].clear();
        }
        sort( indexes_of_particles_to_exchange.begin(), indexes_of_particles_to_exchange.end() );

        int n_part_send = indexes_of_particles_to_exchange.size();
        int n_part_recv;

        int ii,iPart;
        int n_particles,nmove,lmove;
        int shift[(*cubmax).size()+1];//how much we need to shift each bin in order to leave room for the new particles
        double dbin;

        dbin = params.cell_length[0]*params.clrw; //width of a bin.
        for (unsigned int j=0; j<(*cubmax).size()+1 ;j++){
            shift[j]=0;
        }


        Particles diagonalParticles;
        diagonalParticles.initialize(0,params,ispec);

        for (int i=0 ; i<n_part_send ; i++) {
            iPart = indexes_of_particles_to_exchange[i];
            if ( cuParticles.position(iDim,iPart) < min_local[iDim]) {
                buff_index_send[0].push_back( indexes_of_particles_to_exchange[i] );
            }
            else if ( cuParticles.position(iDim,iPart) >= max_local[iDim]) {
                buff_index_send[1].push_back( indexes_of_particles_to_exchange[i] );
            }
            else if ( !(cuParticles.is_part_in_domain(iPart, this) ) ) {
                // at the end of exchangeParticles, diagonalParticles will be reinjected
                // at the end of cuParticles & indexes_of_particles_to_exchange_per_thd[0] for next iDim
                cuParticles.cp_particle(indexes_of_particles_to_exchange[i], diagonalParticles);
            }
            else { // particle will be deleted, if supp_particle particles still in the domain
            }
        } // END for iPart = f(i)

        Particles partVectorSend[2];
        partVectorSend[0].initialize(0,params,ispec);
        partVectorSend[1].initialize(0,params,ispec);
        Particles partVectorRecv[2];
        partVectorRecv[0].initialize(0,params,ispec);
        partVectorRecv[1].initialize(0,params,ispec);

        /********************************************************************************/
        // Exchange particles
        /********************************************************************************/

        MPI_Status sstat    [2];
        MPI_Status rstat    [2];
        MPI_Request srequest[2];
        MPI_Request rrequest[2];

        /********************************************************************************/
        // Exchange number of particles to exchange to establish or not a communication
        /********************************************************************************/
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            if (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
                n_part_send = (buff_index_send[iNeighbor]).size();
                MPI_Isend( &n_part_send, 1, MPI_INT, neighbor_[iDim][iNeighbor], 0, SMILEI_COMM_3D, &(srequest[iNeighbor]) );
            } // END of Send
            else
                n_part_send = 0;
            if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
                buff_index_recv_sz[(iNeighbor+1)%2] = 0;
                //>MPI_Irecv( &(buff_index_recv_sz[(iNeighbor+1)%2]), 1, MPI_INT, neighbor_[iDim][(iNeighbor+1)%2], 0, SMILEI_COMM_3D, &(rrequest[(iNeighbor+1)%2]) );
                MPI_Recv( &(buff_index_recv_sz[(iNeighbor+1)%2]), 1, MPI_INT, neighbor_[iDim][(iNeighbor+1)%2], 0, SMILEI_COMM_3D, &(rstat[(iNeighbor+1)%2]) );
            }
            else
                buff_index_recv_sz[(iNeighbor+1)%2] = 0;

            barrier();
        }
        barrier();

        /********************************************************************************/
        // Wait for end of communications over number of particles
        /********************************************************************************/
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            if (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
                MPI_Wait( &(srequest[iNeighbor]), &(sstat[iNeighbor]) );
            }
            if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
                //>MPI_Wait( &(rrequest[(iNeighbor+1)%2]), &(rstat[(iNeighbor+1)%2]) );
                if (buff_index_recv_sz[(iNeighbor+1)%2]!=0) {
                  partVectorRecv[(iNeighbor+1)%2].initialize( buff_index_recv_sz[(iNeighbor+1)%2], params, ispec);
                }
            }
        }
        barrier();

        /********************************************************************************/
        // Define buffers to exchange buff_index_send[iNeighbor].size();
        /********************************************************************************/


        /********************************************************************************/
        // Proceed to effective Particles' communications
        /********************************************************************************/

        // Number of properties per particles = nDim_Particles + 3 +6 + 1
        int nbrOfProp = 5;
        if(!isSameWeight)
        {
          // particles->weight(0)
          nbrOfProp++;
        }
        if(isImplicit)
        {
          // particles->al_imp, particles->au_imp
          nbrOfProp += 6;
        }

        MPI_Datatype typePartSend, typePartRecv;

        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {

            // n_part_send : number of particles to send to current neighbor
            n_part_send = (buff_index_send[iNeighbor]).size();
            if ( (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) && (n_part_send!=0) ) {
                double x_max = params.cell_length[iDim]*( params.n_space_global[iDim] );
                for (int iPart=0 ; iPart<n_part_send ; iPart++) {
                    if (periods_[iDim]==1) {
                        // Enabled periodicity
                        if ( ( iNeighbor==0 ) &&  (coords_[iDim] == 0 ) &&( cuParticles.position(iDim,buff_index_send[iNeighbor][iPart]) < 0. ) ) {
                            cuParticles.position(iDim,buff_index_send[iNeighbor][iPart])     += x_max;
                        }
                        else if ( ( iNeighbor==1 ) &&  (coords_[iDim] == number_of_procs[iDim]-1 ) && ( cuParticles.position(iDim,buff_index_send[iNeighbor][iPart]) >= x_max ) ) {
                            cuParticles.position(iDim,buff_index_send[iNeighbor][iPart])     -= x_max;
                        }
                    }
                    cuParticles.cp_particle(buff_index_send[iNeighbor][iPart], partVectorSend[iNeighbor]);
                }

                typePartSend = createMPIparticles( &(partVectorSend[iNeighbor]), nbrOfProp );
                MPI_Isend( &((partVectorSend[iNeighbor]).position(0,0)), 1, typePartSend, neighbor_[iDim][iNeighbor], 0, SMILEI_COMM_3D, &(srequest[iNeighbor]) );
                MPI_Type_free( &typePartSend );

            } // END of Send

            n_part_recv = buff_index_recv_sz[(iNeighbor+1)%2];
            if ( (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
                typePartRecv = createMPIparticles( &(partVectorRecv[(iNeighbor+1)%2]), nbrOfProp );
                //>MPI_Irecv( &((partVectorRecv[(iNeighbor+1)%2]).position(0,0)), 1, typePartRecv,  neighbor_[iDim][(iNeighbor+1)%2], 0, SMILEI_COMM_3D, &(rrequest[(iNeighbor+1)%2]) );
                MPI_Recv( &((partVectorRecv[(iNeighbor+1)%2]).position(0,0)), 1, typePartRecv,  neighbor_[iDim][(iNeighbor+1)%2], 0, SMILEI_COMM_3D, &(rstat[(iNeighbor+1)%2]) );
                MPI_Type_free( &typePartRecv );

            } // END of Recv
            barrier();

        } // END for iNeighbor


        /********************************************************************************/
        // Wait for end of communications over Particles
        /********************************************************************************/
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {

            n_part_send = buff_index_send[iNeighbor].size();
            n_part_recv = buff_index_recv_sz[(iNeighbor+1)%2];

            if ( (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) && (n_part_send!=0) ) {
                MPI_Wait( &(srequest[iNeighbor]), &(sstat[iNeighbor]) );
            }

            if ( (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
                //>MPI_Wait( &(rrequest[(iNeighbor+1)%2]), &(rstat[(iNeighbor+1)%2]) );

                // Extract corner particles, not managed in the following process
                // but reinjected at the end in the main list
                for (int iPart=0 ; iPart<n_part_recv; iPart++ )
                    if ( !(partVectorRecv[(iNeighbor+1)%2]).is_part_in_domain(iPart, this) )
                        (partVectorRecv[(iNeighbor+1)%2]).cp_particle(iPart, diagonalParticles);
                for (int iPart=n_part_recv-1 ; iPart>=0; iPart-- ) {
                    if ( !(partVectorRecv[(iNeighbor+1)%2]).is_part_in_domain(iPart, this) ) {
                        (partVectorRecv[(iNeighbor+1)%2]).erase_particle(iPart);
                        buff_index_recv_sz[(iNeighbor+1)%2]--;
                    }
                }

            }

        }
        barrier();
        /********************************************************************************/
        // Clean lists of indexes of particle to exchange per neighbor
        /********************************************************************************/
        for (int i=0 ; i<nbNeighbors_ ; i++)
            buff_index_send[i].clear();


        /********************************************************************************/
        // Delete Particles included in buff_send/buff_recv
        /********************************************************************************/

        // Push lost particles at the end of bins
        //! \todo For loop on bins, can use openMP here.
        for (unsigned int ibin = 0 ; ibin < (*cubmax).size() ; ibin++ ) {
            //cout << "bounds " << (*cubmin)[ibin] << " " << (*cubmax)[ibin] << endl;
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
                    cuParticles.overwrite_part3D((*cubmax)[ibin]-1, iPart );
                    (*cubmax)[ibin]--;
                    ii--;
                    iPart = indexes_of_particles_to_exchange[ii];
                }
                if (iPart >= (*cubmin)[ibin] && iPart < (*cubmax)[ibin]) { //On traite la dernière particule (qui peut aussi etre la premiere)
                    cuParticles.overwrite_part3D((*cubmax)[ibin]-1, iPart );
                    (*cubmax)[ibin]--;
                }
            }
        }
        //Shift the bins in memory
        //Warning: this loop must be executed sequentially. Do not use openMP here.
        for (int unsigned ibin = 1 ; ibin < (*cubmax).size() ; ibin++ ) { //First bin don't need to be shifted
            ii = (*cubmin)[ibin]-(*cubmax)[ibin-1]; // Shift the bin in memory by ii slots.
            iPart = min(ii,(*cubmax)[ibin]-(*cubmin)[ibin]); // Number of particles we have to shift = min (Nshift, Nparticle in the bin)
            if(iPart > 0) cuParticles.overwrite_part3D((*cubmax)[ibin]-iPart,(*cubmax)[ibin-1],iPart);
            (*cubmax)[ibin] -= ii;
            (*cubmin)[ibin] = (*cubmax)[ibin-1];
        }
        // Delete useless Particles
        //Theoretically, not even necessary to do anything as long you use bmax as the end of your iterator on particles.
        //Nevertheless, you might want to free memory and have the actual number of particles
        //really equal to the size of the vector. So we do:
        cuParticles.erase_particle_trail((*cubmax).back());
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        //********************************************************************************/
        // Copy newly arrived particles back to the vector
        // WARNING: very different behaviour depending on which dimension particles are coming from.
        /********************************************************************************/
        //We first evaluate how many particles arrive in each bin.
        if (iDim==1) {
            //1) Count particles coming from south and north
            for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
                n_part_recv = buff_index_recv_sz[iNeighbor];
                for (int j=0; j<n_part_recv ;j++){
                    ii = int((partVectorRecv[iNeighbor].position(0,j)-min_local[0])/dbin);//bin in which the particle goes.
                    shift[ii+1]++; // It makes the next bins shift.
                }
            }
        }
        if (iDim==0) {
            //2) Add particles coming from west and east
            shift[1] += buff_index_recv_sz[0];//Particles coming from south all go to bin 0 and shift all the other bins.
            shift[(*cubmax).size()] += buff_index_recv_sz[1];//Used only to count the total number of particles arrived.
        }

        //Evaluation of the necessary shift of all bins.
        //Must be done sequentially
        for (unsigned int j=1; j<(*cubmax).size()+1;j++){ //bin 0 is not shifted.Last element of shift stores total number of arriving particles.
            shift[j]+=shift[j-1];
        }
        //Make room for new particles
        //cuParticles.create_particles(shift[(*cubmax).size()]);
        cuParticles.initialize( cuParticles.size()+shift[(*cubmax).size()], params, ispec );

        //Shift bins, must be done sequentially
        for (unsigned int j=(*cubmax).size()-1; j>=1; j--){
            n_particles = (*cubmax)[j]-(*cubmin)[j]; //Nbr of particle in this bin
            nmove = min(n_particles,shift[j]); //Nbr of particles to move
            lmove = max(n_particles,shift[j]); //How far particles must be shifted
            if (nmove>0) cuParticles.overwrite_part3D((*cubmin)[j], (*cubmin)[j]+lmove, nmove);
            (*cubmin)[j] += shift[j];
            (*cubmax)[j] += shift[j];
        }

        //Space has been made now to write the arriving particles into the correct bins
        //iDim == 0  is the easy case, when particles arrive either in first or last bin.
        if (iDim==0) {
            for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
                n_part_recv = buff_index_recv_sz[iNeighbor];
                if ( (neighbor_[0][iNeighbor]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
                    ii = iNeighbor*((*cubmax).size()-1);//0 if iNeighbor=0(particles coming from West) and (*cubmax).size()-1 otherwise.
                    partVectorRecv[iNeighbor].overwrite_part3D(0, cuParticles,(*cubmax)[ii],n_part_recv);
                    (*cubmax)[ii] += n_part_recv ;
                }
            }
        }
        //iDim == 1) is the difficult case, when particles can arrive in any bin.
        if (iDim==1) {
            for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
                n_part_recv = buff_index_recv_sz[iNeighbor];
                if ( (neighbor_[1][iNeighbor]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
                    for(int j=0; j<n_part_recv; j++){
                        ii = int((partVectorRecv[iNeighbor].position(0,j)-min_local[0])/dbin);//bin in which the particle goes.
                        partVectorRecv[iNeighbor].overwrite_part3D(j, cuParticles,(*cubmax)[ii]);
                        (*cubmax)[ii] ++ ;
                    }
                }
            }
        }


        // Inject corner particles at the end of the list, update bmax
        //if (iDim==cuParticles.dimension()-1) cout << "Number of diag particles " << diagonalParticles.size() << endl;
        for (int iPart = 0 ; iPart<diagonalParticles.size() ; iPart++) {
            diagonalParticles.cp_particle(iPart, cuParticles);
            (*indexes_of_particles_to_exchange_per_thd)[0].push_back(cuParticles.size()-1);
            (*cubmax)[(*cubmax).size()-1]++;
        }

    }//end of omp master
} // END exchangeParticles


MPI_Datatype SmileiMPI_Cart3D::createMPIparticles( Particles* particles, int nbrOfProp )
{
    MPI_Datatype typeParticlesMPI;

    int nbrOfProp2(nbrOfProp);
    if (particles->isTestParticles) nbrOfProp2++;
    int iProp = 0;

    MPI_Aint address[nbrOfProp2];
    MPI_Get_address( &(particles->position(0,0)), &(address[0]) );
    MPI_Get_address( &(particles->position(1,0)), &(address[1]) );
    //MPI_Get_address( &(particles->position_old(0,0)), &(address[2]) );
    //MPI_Get_address( &(particles->position_old(1,0)), &(address[3]) );
    MPI_Get_address( &(particles->momentum(0,0)), &(address[2]) );
    MPI_Get_address( &(particles->momentum(1,0)), &(address[3]) );
    MPI_Get_address( &(particles->momentum(2,0)), &(address[4]) );
    iProp = 5;
    if(isImplicit)
    {
      MPI_Get_address( &(particles->al_imp(0,0)), &(address[iProp++]) );
      MPI_Get_address( &(particles->al_imp(1,0)), &(address[iProp++]) );
      MPI_Get_address( &(particles->al_imp(2,0)), &(address[iProp++]) );
      MPI_Get_address( &(particles->au_imp(0,0)), &(address[iProp++]) );
      MPI_Get_address( &(particles->au_imp(1,0)), &(address[iProp++]) );
      MPI_Get_address( &(particles->au_imp(2,0)), &(address[iProp++]) );
    }

    if(!isSameWeight)
    {
      MPI_Get_address( &(particles->weight(0)),     &(address[iProp++]) );
    }

    if (particles->isTestParticles)
        MPI_Get_address( &(particles->id(0)),     &(address[iProp++]) );

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
    partDataType[nbrOfProp-1] = MPI_DOUBLE;
    if (particles->isTestParticles)
        partDataType[nbrOfProp2-1] = MPI_UNSIGNED;

    MPI_Type_struct( nbrOfProp2, &(nbr_parts[0]), &(disp[0]), &(partDataType[0]), &typeParticlesMPI);
    MPI_Type_commit( &typeParticlesMPI );

    return typeParticlesMPI;
} // END createMPIparticles

void SmileiMPI_Cart3D::createType( PicParams& params )
{
    int nx0 = params.n_space[0] + 1 + 2*params.oversize[0];
    int ny0 = params.n_space[1] + 1 + 2*params.oversize[1];
    int clrw = params.clrw;

    // MPI_Datatype ntype_[nDim][primDual][primDual]
    int nx, ny;
    int nline, ncol;
    for (int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++) {
        nx = nx0 + ix_isPrim;
        for (int iy_isPrim=0 ; iy_isPrim<2 ; iy_isPrim++) {
            ny = ny0 + iy_isPrim;
            ntype_[0][ix_isPrim][iy_isPrim] = NULL;
            MPI_Type_contiguous(ny, MPI_DOUBLE, &(ntype_[0][ix_isPrim][iy_isPrim]));    //line
            MPI_Type_commit( &(ntype_[0][ix_isPrim][iy_isPrim]) );
            ntype_[1][ix_isPrim][iy_isPrim] = NULL;
            MPI_Type_vector(nx, 1, ny, MPI_DOUBLE, &(ntype_[1][ix_isPrim][iy_isPrim])); // column
            MPI_Type_commit( &(ntype_[1][ix_isPrim][iy_isPrim]) );
            ntype_[2][ix_isPrim][iy_isPrim] = NULL;
            MPI_Type_contiguous(ny*clrw, MPI_DOUBLE, &(ntype_[2][ix_isPrim][iy_isPrim]));   //clrw lines
            MPI_Type_commit( &(ntype_[2][ix_isPrim][iy_isPrim]) );

            ntypeSum_[0][ix_isPrim][iy_isPrim] = NULL;
            nline = 1 + 2*params.oversize[0] + ix_isPrim;
            MPI_Type_contiguous(nline, ntype_[0][ix_isPrim][iy_isPrim], &(ntypeSum_[0][ix_isPrim][iy_isPrim]));    //line
            MPI_Type_commit( &(ntypeSum_[0][ix_isPrim][iy_isPrim]) );
            ntypeSum_[1][ix_isPrim][iy_isPrim] = NULL;
            ncol  = 1 + 2*params.oversize[1] + iy_isPrim;
            MPI_Type_vector(nx, ncol, ny, MPI_DOUBLE, &(ntypeSum_[1][ix_isPrim][iy_isPrim])); // column
            MPI_Type_commit( &(ntypeSum_[1][ix_isPrim][iy_isPrim]) );

        }
    }

} //END createType


void SmileiMPI_Cart3D::sumField( Field* field )
{
    std::vector<unsigned int> n_elem = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field3D* f3D =  static_cast<Field3D*>(field);


    // Use a buffer per direction to exchange data before summing
    Field3D buf[ndims_][ nbNeighbors_ ];
    // Size buffer is 2 oversize (1 inside & 1 outside of the current subdomain)
    std::vector<unsigned int> oversize2 = oversize;
    oversize2[0] *= 2;
    oversize2[0] += 1 + f3D->isDual_[0];
    oversize2[1] *= 2;
    oversize2[1] += 1 + f3D->isDual_[1];

    for (int iDim=0 ; iDim<ndims_ ; iDim++) {
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            std::vector<unsigned int> tmp(ndims_,0);
            tmp[0] =    iDim  * n_elem[0] + (1-iDim) * oversize2[0];
            tmp[1] = (1-iDim) * n_elem[1] +    iDim  * oversize2[1];
            buf[iDim][iNeighbor].allocateDims( tmp );
        }
    }

    int istart, ix, iy;

    /********************************************************************************/
    // Send/Recv in a buffer data to sum
    /********************************************************************************/
    for (int iDim=0 ; iDim<ndims_ ; iDim++) {

        MPI_Datatype ntype = ntypeSum_[iDim][isDual[0]][isDual[1]];
        //      MPI_Status stat[2];
        //      MPI_Request request[2];
        MPI_Status sstat    [ndims_][2];
        MPI_Status rstat    [ndims_][2];
        MPI_Request srequest[ndims_][2];
        MPI_Request rrequest[ndims_][2];

        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {

            if (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
                istart = iNeighbor * ( n_elem[iDim]- oversize2[iDim] ) + (1-iNeighbor) * ( 0 );
                ix = (1-iDim)*istart;
                iy =    iDim *istart;
                MPI_Isend( &(f3D->data_3D[ix][iy]), 1, ntype, neighbor_[iDim][iNeighbor], 0, SMILEI_COMM_3D, &(srequest[iDim][iNeighbor]) );
            } // END of Send

            if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
                int tmp_elem = (buf[iDim][(iNeighbor+1)%2]).dims_[0]*(buf[iDim][(iNeighbor+1)%2]).dims_[1];
                MPI_Irecv( &( (buf[iDim][(iNeighbor+1)%2]).data_3D[0][0] ), tmp_elem, MPI_DOUBLE, neighbor_[iDim][(iNeighbor+1)%2], 0, SMILEI_COMM_3D, &(rrequest[iDim][(iNeighbor+1)%2]) );
            } // END of Recv

        } // END for iNeighbor


        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            if (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
                MPI_Wait( &(srequest[iDim][iNeighbor]), &(sstat[iDim][iNeighbor]) );
            }
            if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
                MPI_Wait( &(rrequest[iDim][(iNeighbor+1)%2]), &(rstat[iDim][(iNeighbor+1)%2]) );
            }
        }


        // Synchro before summing, to not sum with data ever sum
        // Merge loops, Sum direction by direction permits to not communicate with diagonal neighbors
        barrier();
        /********************************************************************************/
        // Sum data on each process, same operation on both side
        /********************************************************************************/

        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            istart = ( (iNeighbor+1)%2 ) * ( n_elem[iDim]- oversize2[iDim] ) + (1-(iNeighbor+1)%2) * ( 0 );
            int ix0 = (1-iDim)*istart;
            int iy0 =    iDim *istart;
            if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
                for (unsigned int ix=0 ; ix< (buf[iDim][(iNeighbor+1)%2]).dims_[0] ; ix++) {
                    for (unsigned int iy=0 ; iy< (buf[iDim][(iNeighbor+1)%2]).dims_[1] ; iy++)
                        f3D->data_3D[ix0+ix][iy0+iy] += (buf[iDim][(iNeighbor+1)%2])(ix,iy);
                }
            } // END if

        } // END for iNeighbor

        barrier();

    } // END for iDim

} // END sumField



void SmileiMPI_Cart3D::scatterGrid( Grid* grid )
{
    int procs_rk;
    int iGlobal, jGlobal;
    int iGlobal_gather;

    Grid3D* grid3D = static_cast<Grid3D*>(grid);
    //TITLE("scatterGrid1111");
    for(int iProcs = 0; iProcs < number_of_procs[0]; iProcs++)
    {
        for(int jProcs = 0; jProcs < number_of_procs[1]; jProcs++)
        {
            procs_rk = iProcs * number_of_procs[1] + jProcs;
            for(int i = 0; i < dims_gather[procs_rk*2]; i++)
            {
                for(int j = 0; j < dims_gather[procs_rk*2 + 1]; j++)
                {
                    iGlobal = iProcs * (dims_gather[0] - 2*oversize[0] -1) + i -oversize[0];
                    jGlobal = jProcs * (dims_gather[1] - 2*oversize[1] -1) + j -oversize[1];
                    if(iProcs == 0 && i < oversize[0] || iProcs == number_of_procs[0] -1 && i > dims_gather[procs_rk*2] - 1 - oversize[0]){
                        iGlobal = abs((int)grid3D->globalDims_[0] - abs(iGlobal) - 1);
                    }
                    if(jProcs == 0 && j < oversize[1] || jProcs == number_of_procs[1] -1 && j > dims_gather[procs_rk*2+1] - 1 - oversize[1]){
                        jGlobal = abs((int)grid3D->globalDims_[1] - abs(jGlobal) - 1);
                    }
                    iGlobal_gather = send_disp[procs_rk] + i * dims_gather[procs_rk*2+1] + j;
                    //if(iGlobal >= ii || jGlobal >= jj) cout<<"error "<<iGlobal<<" "<<iProcs<<" "<<dims_gather[0]<<" "<<oversize[0]<<endl;

                    grid_global_gather[iGlobal_gather] = grid3D->iswall_global_3D[iGlobal][jGlobal];


                }
            }

        }
    }
    MPI_Scatterv(grid_global_gather, &send_cnt[0], &send_disp[0], MPI_INT, &grid3D->iswall_3D[0][0], recv_cnt[smilei_rk], MPI_INT, 0, SMILEI_COMM_3D);
    //TITLE("scatterGrid");

} // END scatterGrid


void SmileiMPI_Cart3D::gatherRho( Field* field_global ,Field* field  )
{

    int procs_rk;
    int iGlobal, jGlobal;
    int iGlobal_gather;
    int nx, ny;

    Field3D* f3D =  static_cast<Field3D*>(field);
    Field3D* f3D_global =  static_cast<Field3D*>(field_global);
    nx = f3D_global->dims_[0];
    ny = f3D_global->dims_[1];
    f3D_global->put_to(0.0);
    MPI_Gatherv(f3D->data_, send_cnt[smilei_rk], MPI_DOUBLE, field_global_gather, &recv_cnt[0], &recv_disp[0], MPI_DOUBLE, 0, SMILEI_COMM_3D);

    for(int iProcs = 0; iProcs < number_of_procs[0]; iProcs++)
    {
        for(int jProcs = 0; jProcs < number_of_procs[1]; jProcs++)
        {
            procs_rk = iProcs * number_of_procs[1] + jProcs;
            for(int i = 0; i < dims_gather[procs_rk*2]; i++)
            {
                for(int j = 0; j < dims_gather[procs_rk*2 + 1]; j++)
                {
                    iGlobal = iProcs * (dims_gather[0] - 2*oversize[0] -1) + i -oversize[0];
                    jGlobal = jProcs * (dims_gather[1] - 2*oversize[1] -1) + j -oversize[1];
                    if(iProcs == 0 && i < oversize[0] || iProcs == number_of_procs[0] -1 && i > dims_gather[procs_rk*2] - 1 - oversize[0]){
                        iGlobal = abs((int)f3D_global->dims_[0] - abs(iGlobal) - 1);
                    }
                    if(jProcs == 0 && j < oversize[1] || jProcs == number_of_procs[1] -1 && j > dims_gather[procs_rk*2+1] - 1 - oversize[1]){
                        jGlobal = abs((int)f3D_global->dims_[1] - abs(jGlobal) - 1);
                    }
                    iGlobal_gather = send_disp[procs_rk] + i * dims_gather[procs_rk*2+1] + j;
                    //if(iGlobal >= ii || jGlobal >= jj) cout<<"error "<<iGlobal<<" "<<iProcs<<" "<<dims_gather[0]<<" "<<oversize[0]<<endl;

                    f3D_global->data_3D[iGlobal][jGlobal] += field_global_gather[iGlobal_gather];

                    //if(f3D_global->data_3D[iGlobal][jGlobal] != 0.0) cout<<"ereeee"; //<<f3D_global->data_3D[iGlobal][jGlobal]<<endl;


                }
            }

        }
    }

    //> Handle the boundary points and corner points,
    //> this is meaningful for periodic boundary condition,
    //> but not affect the results of other boudary conditions.
    for(int i = 0; i < nx; i++)
    {
        for(int j = 0; j < ny; j++)
        {
            if( i == 0){
                f3D_global->data_3D[i][j] += f3D_global->data_3D[nx-1][j];
            }
            else if(j == 0){
                f3D_global->data_3D[i][j] += f3D_global->data_3D[i][ny-1];
            }
        }

    }
    for(int i = 0; i < nx; i++)
    {
        for(int j = 0; j < ny; j++)
        {
            if( i == nx-1){
                f3D_global->data_3D[i][j] = f3D_global->data_3D[0][j];
            }
            else if(j == ny-1){
                f3D_global->data_3D[i][j] = f3D_global->data_3D[i][0];
            }
        }

    }



} // END gatherRho


void SmileiMPI_Cart3D::gatherField( Field* field_global ,Field* field  )
{

    int procs_rk;
    int iGlobal, jGlobal;
    int iGlobal_gather;
    int nx, ny;

    Field3D* f3D =  static_cast<Field3D*>(field);
    Field3D* f3D_global =  static_cast<Field3D*>(field_global);
    nx = f3D_global->dims_[0];
    ny = f3D_global->dims_[1];
    f3D_global->put_to(0.0);
    MPI_Gatherv(f3D->data_, send_cnt[smilei_rk], MPI_DOUBLE, field_global_gather, &recv_cnt[0], &recv_disp[0], MPI_DOUBLE, 0, SMILEI_COMM_3D);

    for(int iProcs = 0; iProcs < number_of_procs[0]; iProcs++)
    {
        for(int jProcs = 0; jProcs < number_of_procs[1]; jProcs++)
        {
            procs_rk = iProcs * number_of_procs[1] + jProcs;
            for(int i = 0; i < dims_gather[procs_rk*2]; i++)
            {
                for(int j = 0; j < dims_gather[procs_rk*2 + 1]; j++)
                {
                    iGlobal = iProcs * (dims_gather[0] - 2*oversize[0] -1) + i -oversize[0];
                    jGlobal = jProcs * (dims_gather[1] - 2*oversize[1] -1) + j -oversize[1];
                    if(iProcs == 0 && i < oversize[0] || iProcs == number_of_procs[0] -1 && i > dims_gather[procs_rk*2] - 1 - oversize[0]){
                        iGlobal = abs((int)f3D_global->dims_[0] - abs(iGlobal) - 1);
                    }
                    if(jProcs == 0 && j < oversize[1] || jProcs == number_of_procs[1] -1 && j > dims_gather[procs_rk*2+1] - 1 - oversize[1]){
                        jGlobal = abs((int)f3D_global->dims_[1] - abs(jGlobal) - 1);
                    }
                    iGlobal_gather = send_disp[procs_rk] + i * dims_gather[procs_rk*2+1] + j;
                    //if(iGlobal >= ii || jGlobal >= jj) cout<<"error "<<iGlobal<<" "<<iProcs<<" "<<dims_gather[0]<<" "<<oversize[0]<<endl;

                    f3D_global->data_3D[iGlobal][jGlobal] = field_global_gather[iGlobal_gather];

                    //if(f3D_global->data_3D[iGlobal][jGlobal] != 0.0) cout<<"ereeee"; //<<f3D_global->data_3D[iGlobal][jGlobal]<<endl;


                }
            }

        }
    }

    //> Handle the boundary points and corner points,
    //> this is meaningful for periodic boundary condition,
    //> but not affect the results of other boudary conditions.
    for(int i = 0; i < nx; i++)
    {
        for(int j = 0; j < ny; j++)
        {
            if( i == 0){
                f3D_global->data_3D[i][j] += f3D_global->data_3D[nx-1][j];
            }
            else if(j == 0){
                f3D_global->data_3D[i][j] += f3D_global->data_3D[i][ny-1];
            }
        }

    }
    for(int i = 0; i < nx; i++)
    {
        for(int j = 0; j < ny; j++)
        {
            if( i == nx-1){
                f3D_global->data_3D[i][j] = f3D_global->data_3D[0][j];
            }
            else if(j == ny-1){
                f3D_global->data_3D[i][j] = f3D_global->data_3D[i][0];
            }
        }

    }



} // END gatherField



void SmileiMPI_Cart3D::scatterField( Field* field_global ,Field* field )
{

    int procs_rk;
    int iGlobal, jGlobal;
    int iGlobal_gather;

    Field3D* f3D =  static_cast<Field3D*>(field);
    Field3D* f3D_global =  static_cast<Field3D*>(field_global);

    int ii,jj;
    ii=f3D_global->dims_[0];
    jj=f3D_global->dims_[1];

    f3D->put_to(0.0);

    iGlobal = 0;
    jGlobal = 0;


    for(int iProcs = 0; iProcs < number_of_procs[0]; iProcs++)
    {
        for(int jProcs = 0; jProcs < number_of_procs[1]; jProcs++)
        {
            procs_rk = iProcs * number_of_procs[1] + jProcs;
            for(int i = 0; i < dims_gather[procs_rk*2]; i++)
            {
                for(int j = 0; j < dims_gather[procs_rk*2 + 1]; j++)
                {
                    iGlobal = iProcs * (dims_gather[0] - 2*oversize[0] -1) + i -oversize[0];
                    jGlobal = jProcs * (dims_gather[1] - 2*oversize[1] -1) + j -oversize[1];
                    if(iProcs == 0 && i < oversize[0] || iProcs == number_of_procs[0] -1 && i > dims_gather[procs_rk*2] - 1 - oversize[0]){
                        iGlobal = abs((int)f3D_global->dims_[0] - abs(iGlobal) - 1);
                    }
                    if(jProcs == 0 && j < oversize[1] || jProcs == number_of_procs[1] -1 && j > dims_gather[procs_rk*2+1] - 1 - oversize[1]){
                        jGlobal = abs((int)f3D_global->dims_[1] - abs(jGlobal) - 1);
                    }
                    iGlobal_gather = send_disp[procs_rk] + i * dims_gather[procs_rk*2+1] + j;
                    if(iGlobal >= ii || jGlobal >= jj) cout<<"error "<<iGlobal<<" "<<iProcs<<" "<<dims_gather[0]<<" "<<oversize[0]<<endl;

                    field_global_gather[iGlobal_gather] = f3D_global->data_3D[iGlobal][jGlobal];

                    //if(f3D_global->data_3D[iGlobal][jGlobal] != 0.0) cout<<"ereeee"; //<<f3D_global->data_3D[iGlobal][jGlobal]<<endl;

                }
            }

        }
    }
    MPI_Scatterv(field_global_gather, &send_cnt[0], &send_disp[0], MPI_DOUBLE, f3D->data_, recv_cnt[smilei_rk], MPI_DOUBLE, 0, SMILEI_COMM_3D);



} // END scatterField
