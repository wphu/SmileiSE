#include "ElectroMagn3D.h"

#include <cmath>

#include <iostream>
#include <sstream>

#include "PicParams.h"
#include "Field3D.h"

#include "SmileiMPI.h"
#include "SmileiMPI_Cart3D.h"

#include "Profile.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Electromagn3D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn3D::ElectroMagn3D(PicParams &params, InputData &input_data, SmileiMPI* smpi) :
ElectroMagn(params, input_data, smpi),
isWestern(smpi->isWestern()),
isEastern(smpi->isEastern()),
isSouthern(smpi->isSouthern()),
isNorthern(smpi->isNorthern())
{
    // local dt to store
    SmileiMPI_Cart3D* smpi3D = static_cast<SmileiMPI_Cart3D*>(smpi);
//    int process_coord_x = smpi3D->getProcCoord(0);


    // --------------------------------------------------
    // Calculate quantities related to the simulation box
    // --------------------------------------------------

    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the x-direction)
    dx       = cell_length[0];
    dt_ov_dx = timestep/dx;
    dx_ov_dt = 1.0/dt_ov_dx;

    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the y-direction)
    dy       = cell_length[1];
    dt_ov_dy = timestep/dy;
    dy_ov_dt = 1.0/dt_ov_dy;

    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the y-direction)
    dz       = cell_length[2];
    dt_ov_dz = timestep/dz;
    dz_ov_dt = 1.0/dt_ov_dz;

    // ----------------------
    // Electromagnetic fields
    // ----------------------
    //! \todo Homogenize 3D/3D dimPrim/dimDual or nx_p/nx_d/ny_p/ny_d

    dimPrim.resize( nDim_field );
    dimDual.resize( nDim_field );
    dim_global.resize( nDim_field );

    // Dimension of the primal and dual grids
    for (size_t i=0 ; i<nDim_field ; i++) {
        // Standard scheme
        dimPrim[i] = n_space[i]+1;
        dimDual[i] = n_space[i]+2;
        // + Ghost domain
        dimPrim[i] += 2*oversize[i];
        dimDual[i] += 2*oversize[i];

        dim_global[i] = n_space_global[i] + 1;
    }
    // number of nodes of the primal and dual grid in the x-direction
    nx_p = n_space[0]+1+2*oversize[0];
    nx_d = n_space[0]+2+2*oversize[0];
    // number of nodes of the primal and dual grid in the y-direction
    ny_p = n_space[1]+1+2*oversize[1];
    ny_d = n_space[1]+2+2*oversize[1];
    // number of nodes of the primal and dual grid in the y-direction
    nz_p = n_space[2]+1+2*oversize[2];
    nz_d = n_space[2]+2+2*oversize[2];

    // Allocation of the EM fields
    Ex_  = new Field3D(dimPrim, "Ex" );
    Ey_  = new Field3D(dimPrim, "Ey" );
    Ez_  = new Field3D(dimPrim, "Ez" );
    Bx_  = new Field3D(dimPrim, "Bx" );
    By_  = new Field3D(dimPrim, "By" );
    Bz_  = new Field3D(dimPrim, "Bz" );
    Bx_m = new Field3D(dimPrim, "Bx_m" );
    By_m = new Field3D(dimPrim, "By_m" );
    Bz_m = new Field3D(dimPrim, "Bz_m" );

    // Total charge currents and densities
    Jx_   = new Field3D(dimPrim, "Jx" );
    Jy_   = new Field3D(dimPrim, "Jy" );
    Jz_   = new Field3D(dimPrim, "Jz" );
    rho_  = new Field3D(dimPrim, "Rho" );
    rho_avg  = new Field3D(dimPrim, "Rho_avg" );

    // Allocation of time-averaged EM fields
    phi_avg = new Field3D(dimPrim, "Phi_avg" );
    Ex_avg  = new Field3D(dimPrim, "Ex_avg");
    Ey_avg  = new Field3D(dimPrim, "Ey_avg");
    Ez_avg  = new Field3D(dimPrim, "Ez_avg");
    Bx_avg  = new Field3D(dimPrim, "Bx_avg");
    By_avg  = new Field3D(dimPrim, "By_avg");
    Bz_avg  = new Field3D(dimPrim, "Bz_avg");



    rho_global = new Field3D(dim_global, "Rho_global");
    phi_global = new Field3D(dim_global, "Phi_global");
    Ex_global = new Field3D(dim_global, "Ex_global");
    Ey_global = new Field3D(dim_global, "Ey_global");
    Ez_global = new Field3D(dim_global, "Ez_global");


    rho_global_avg = new Field3D(dim_global, "Rho_global_avg");
    phi_global_avg = new Field3D(dim_global, "Phi_global_avg");
    Ex_global_avg  = new Field3D(dim_global, "Ex_global_avg");
    Ey_global_avg  = new Field3D(dim_global, "Ey_global_avg");

    rho_global->put_to(0.0);
    phi_global->put_to(0.0);

    Ex_->put_to(0.0);
    Ey_->put_to(0.0);
    Ez_->put_to(0.0);
    Bx_->put_to(params.externB[0]);
    By_->put_to(params.externB[1]);
    Bz_->put_to(params.externB[2]);
    Bx_m->put_to(0.0);
    By_m->put_to(0.0);
    Bz_m->put_to(0.0);
    rho_->put_to(0.0);



    // Allocation of the time-averaged EM fields
    Ex_avg  = new Field3D(dimPrim, "Ex_avg" );
    Ey_avg  = new Field3D(dimPrim, "Ey_avg" );
    Ez_avg  = new Field3D(dimPrim, "Ez_avg" );
    Bx_avg  = new Field3D(dimPrim, "Bx_avg" );
    By_avg  = new Field3D(dimPrim, "By_avg" );
    Bz_avg  = new Field3D(dimPrim, "Bz_avg" );

    // Charge currents currents and density for each species
    for (unsigned int ispec=0; ispec<n_species; ispec++) {
        Jx_s[ispec]  = new Field3D(dimPrim, ("Jx_"+params.species_param[ispec].species_type).c_str());
        Jy_s[ispec]  = new Field3D(dimPrim, ("Jy_"+params.species_param[ispec].species_type).c_str());
        Jz_s[ispec]  = new Field3D(dimPrim, ("Jz_"+params.species_param[ispec].species_type).c_str());
        rho_s[ispec] = new Field3D(dimPrim, ("Rho_"+params.species_param[ispec].species_type).c_str());
        rho_s_avg[ispec]        = new Field3D(dimPrim, ("Rho_"+params.species_param[ispec].species_type+"_avg").c_str());

        rho_s_global[ispec] = new Field3D(dim_global, ("Rho_global_"+params.species_param[ispec].species_type).c_str());
        rho_s_global_avg[ispec] = new Field3D(dim_global, ("Rho_global_"+params.species_param[ispec].species_type+"_avg").c_str());

        Vx_s[ispec]             = new Field3D(dimPrim, ("Vx_"+params.species_param[ispec].species_type).c_str());
        Vx_s_avg[ispec]         = new Field3D(dimPrim, ("Vx_"+params.species_param[ispec].species_type+"_avg").c_str());
        Vx_s_global[ispec]      = new Field3D(dim_global, ("Vx_global_"+params.species_param[ispec].species_type).c_str());
        Vx_s_global_avg[ispec]  = new Field3D(dim_global, ("Vx_global_"+params.species_param[ispec].species_type+"_avg").c_str());

        Vy_s[ispec]             = new Field3D(dimPrim, ("Vy_"+params.species_param[ispec].species_type).c_str());
        Vy_s_avg[ispec]         = new Field3D(dimPrim, ("Vy_"+params.species_param[ispec].species_type+"_avg").c_str());
        Vy_s_global[ispec]      = new Field3D(dim_global, ("Vy_global_"+params.species_param[ispec].species_type).c_str());
        Vy_s_global_avg[ispec]  = new Field3D(dim_global, ("Vy_global_"+params.species_param[ispec].species_type+"_avg").c_str());

        Vz_s[ispec]             = new Field3D(dimPrim, ("Vz_"+params.species_param[ispec].species_type).c_str());
        Vz_s_avg[ispec]         = new Field3D(dimPrim, ("Vz_"+params.species_param[ispec].species_type+"_avg").c_str());
        Vz_s_global[ispec]      = new Field3D(dim_global, ("Vz_global_"+params.species_param[ispec].species_type).c_str());
        Vz_s_global_avg[ispec]  = new Field3D(dim_global, ("Vz_global_"+params.species_param[ispec].species_type+"_avg").c_str());


        Vp_s[ispec]             = new Field3D(dimPrim, ("Vparallel_"+params.species_param[ispec].species_type).c_str());
        Vp_s_avg[ispec]         = new Field3D(dimPrim, ("Vparallel_"+params.species_param[ispec].species_type+"_avg").c_str());
        Vp_s_global[ispec]      = new Field3D(dim_global, ("Vparallel_global_"+params.species_param[ispec].species_type).c_str());
        Vp_s_global_avg[ispec]  = new Field3D(dim_global, ("Vparallel_global_"+params.species_param[ispec].species_type+"_avg").c_str());

        T_s[ispec]              = new Field3D(dimPrim, ("T_"+params.species_param[ispec].species_type).c_str());
        T_s_avg[ispec]          = new Field3D(dimPrim, ("T_"+params.species_param[ispec].species_type+"_avg").c_str());
        T_s_global[ispec]       = new Field3D(dim_global, ("T_global_"+params.species_param[ispec].species_type).c_str());
        T_s_global_avg[ispec]   = new Field3D(dim_global, ("T_global_"+params.species_param[ispec].species_type+"_avg").c_str());


    }


//    ostringstream file_name("");
//    for (unsigned int ispec=0; ispec<n_species; ispec++) {
//        file_name.str("");
//        file_name << "Jx_s" << ispec;
//        Jx_s[ispec]  = new Field3D(dimPrim, 0, false, file_name.str().c_str());
//        file_name.str("");
//        file_name << "Jy_s" << ispec;
//        Jy_s[ispec]  = new Field3D(dimPrim, 1, false, file_name.str().c_str());
//        file_name.str("");
//        file_name << "Jz_s" << ispec;
//        Jz_s[ispec]  = new Field3D(dimPrim, 2, false, file_name.str().c_str());
//        file_name.str("");
//        file_name << "rho_s" << ispec;
//        rho_s[ispec] = new Field3D(dimPrim, file_name.str().c_str());
//    }


    // ----------------------------------------------------------------
    // Definition of the min and max index according to chosen oversize
    // ----------------------------------------------------------------
    index_bc_min.resize( nDim_field, 0 );
    index_bc_max.resize( nDim_field, 0 );
    for (unsigned int i=0 ; i<nDim_field ; i++) {
        index_bc_min[i] = oversize[i];
        index_bc_max[i] = dimDual[i]-oversize[i]-1;
    }
    /*
     MESSAGE("index_bc_min / index_bc_max / nx_p / nx_d" << index_bc_min[0]
            << " " << index_bc_max[0] << " " << nx_p<< " " << nx_d);
     */


    // Define limits of non duplicated elements
    // (by construction 1 (prim) or 2 (dual) elements shared between per MPI process)
    // istart
    for (unsigned int i=0 ; i<3 ; i++)
        for (unsigned int isDual=0 ; isDual<2 ; isDual++)
            istart[i][isDual] = 0;
    for (unsigned int i=0 ; i<nDim_field ; i++) {
        for (unsigned int isDual=0 ; isDual<2 ; isDual++) {
            istart[i][isDual] = oversize[i];
            if (smpi3D->getProcCoord(i)!=0) istart[i][isDual]+=1;
        }
    }

    // bufsize = nelements
    for (unsigned int i=0 ; i<3 ; i++)
        for (unsigned int isDual=0 ; isDual<2 ; isDual++)
            bufsize[i][isDual] = 1;

    for (unsigned int i=0 ; i<nDim_field ; i++) {
        for (int isDual=0 ; isDual<2 ; isDual++)
            bufsize[i][isDual] = n_space[i] + 1;

        for (int isDual=0 ; isDual<2 ; isDual++) {
            bufsize[i][isDual] += isDual;
            if ( smpi3D->getNbrOfProcs(i)!=1 ) {

                if ( ( !isDual ) && (smpi3D->getProcCoord(i)!=0) )
                    bufsize[i][isDual]--;
                else if  (isDual) {
                    bufsize[i][isDual]--;
                    if ( (smpi3D->getProcCoord(i)!=0) && (smpi3D->getProcCoord(i)!=smpi3D->getNbrOfProcs(i)-1) )
                        bufsize[i][isDual]--;
                }

            } // if ( smpi3D->getNbrOfProcs(i)!=1 )
        } // for (int isDual=0 ; isDual
    } // for (unsigned int i=0 ; i<nDim_field

}//END constructor Electromagn3D



// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Electromagn3D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn3D::~ElectroMagn3D()
{
}//END ElectroMagn3D


// ---------------------------------------------------------------------------------------------------------------------
// Save the former Magnetic-Fields (used to center them)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3D::saveMagneticFields()
{
    // Static cast of the fields
    Field3D* Bx3D   = static_cast<Field3D*>(Bx_);
    Field3D* By3D   = static_cast<Field3D*>(By_);
    Field3D* Bz3D   = static_cast<Field3D*>(Bz_);
    Field3D* Bx3D_m = static_cast<Field3D*>(Bx_m);
    Field3D* By3D_m = static_cast<Field3D*>(By_m);
    Field3D* Bz3D_m = static_cast<Field3D*>(Bz_m);

    // Magnetic field Bx^(p,d)
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Bx3D_m)(i,j)=(*Bx3D)(i,j);
        }

    // Magnetic field By^(d,p)
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*By3D_m)(i,j)=(*By3D)(i,j);
        }

    // Magnetic field Bz^(d,d)
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Bz3D_m)(i,j)=(*Bz3D)(i,j);
        }
    }// end for j
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*By3D_m)(nx_p,j)=(*By3D)(nx_p,j);
        }
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Bz3D_m)(nx_p,j)=(*Bz3D)(nx_p,j);
        }

}//END saveMagneticFields



// ---------------------------------------------------------------------------------------------------------------------
// Solve the Maxwell-Ampere equation
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3D::solveMaxwellAmpere()
{
    // Static-cast of the fields
    Field3D* Ex3D = static_cast<Field3D*>(Ex_);
    Field3D* Ey3D = static_cast<Field3D*>(Ey_);
    Field3D* Ez3D = static_cast<Field3D*>(Ez_);
    Field3D* Bx3D = static_cast<Field3D*>(Bx_);
    Field3D* By3D = static_cast<Field3D*>(By_);
    Field3D* Bz3D = static_cast<Field3D*>(Bz_);
    Field3D* Jx3D = static_cast<Field3D*>(Jx_);
    Field3D* Jy3D = static_cast<Field3D*>(Jy_);
    Field3D* Jz3D = static_cast<Field3D*>(Jz_);

    // Electric field Ex^(d,p)
    for (unsigned int i=0 ; i<nx_p ; i++) {
      //    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Ex3D)(i,j) += -timestep*(*Jx3D)(i,j) + dt_ov_dy * ( (*Bz3D)(i,j+1) - (*Bz3D)(i,j) );
        }// end for j
        //    }// end for i

    // Electric field Ey^(p,d)

    //    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Ey3D)(i,j) += -timestep*(*Jy3D)(i,j) - dt_ov_dx * ( (*Bz3D)(i+1,j) - (*Bz3D)(i,j) );
        }// end for j
    //} // end for i

    // Electric field Ez^(p,p)
    //for (unsigned int i=0 ;  i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Ez3D)(i,j) += -timestep*(*Jz3D)(i,j)
            +               dt_ov_dx * ( (*By3D)(i+1,j) - (*By3D)(i,j) )
            -               dt_ov_dy * ( (*Bx3D)(i,j+1) - (*Bx3D)(i,j) );
        } // end for j
    }// end for i



      for (unsigned int j=0 ; j<ny_p ; j++) {
          (*Ex3D)(nx_p,j) += -timestep*(*Jx3D)(nx_p,j) + dt_ov_dy * ( (*Bz3D)(nx_p,j+1) - (*Bz3D)(nx_p,j) );
      }


}//END solveMaxwellAmpere


// ---------------------------------------------------------------------------------------------------------------------
// Center the Magnetic Fields (used to push the particle)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3D::centerMagneticFields()
{
    // Static cast of the fields
    Field3D* Bx3D   = static_cast<Field3D*>(Bx_);
    Field3D* By3D   = static_cast<Field3D*>(By_);
    Field3D* Bz3D   = static_cast<Field3D*>(Bz_);
    Field3D* Bx3D_m = static_cast<Field3D*>(Bx_m);
    Field3D* By3D_m = static_cast<Field3D*>(By_m);
    Field3D* Bz3D_m = static_cast<Field3D*>(Bz_m);

    // Magnetic field Bx^(p,d)
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Bx3D_m)(i,j) = ( (*Bx3D)(i,j) + (*Bx3D_m)(i,j) )*0.5;
        }


    // Magnetic field By^(d,p)

        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*By3D_m)(i,j) = ( (*By3D)(i,j) + (*By3D_m)(i,j) )*0.5;
        }


    // Magnetic field Bz^(d,d)
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Bz3D_m)(i,j) = ( (*Bz3D)(i,j) + (*Bz3D_m)(i,j) )*0.5;
        } // end for j
      } // end for i


      for (unsigned int j=0 ; j<ny_p ; j++) {
          (*By3D_m)(nx_p,j) = ( (*By3D)(nx_p,j) + (*By3D_m)(nx_p,j) )*0.5;
      }
      for (unsigned int j=0 ; j<ny_p ; j++) {
          (*Bz3D_m)(nx_p,j) = ( (*Bz3D)(nx_p,j) + (*Bz3D_m)(nx_p,j) )*0.5;
      } // end for j


}//END centerMagneticFields



// ---------------------------------------------------------------------------------------------------------------------
// Reset/Increment the averaged fields
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3D::incrementAvgFields(unsigned int time_step)
{
    // reset the averaged fields for (time_step-1)%ntime_step_avg == 0
    if ( (time_step-1) % dump_step == 0 ){
        rho_global_avg->put_to(0.0);
        phi_global_avg->put_to(0.0);
        Ex_global_avg->put_to(0.0);
        for (unsigned int ispec=0; ispec<n_species; ispec++) {
            rho_s_avg[ispec]->put_to(0.0);
        }//END loop on species ispec
    }

    // Calculate the sum values for global rho phi Ex and Ey
    if( (time_step % dump_step) > (dump_step - avg_step) || (time_step % dump_step) == 0 )
    {
        for (unsigned int i=0 ; i<dim_global[0]*dim_global[1]*dim_global[2] ; i++) {
            (*rho_global_avg)(i) += (*rho_global)(i);
            (*phi_global_avg)(i) += (*phi_global)(i);
            (*Ex_global_avg)(i)  += (*Ex_global)(i);
        }

        // Calculate the sum values for density of each species
        for (unsigned int ispec=0; ispec<n_species; ispec++) {
            // all fields are defined on the primal grid
            for (unsigned int ix=0 ; ix<dimPrim[0]*dimPrim[1]*dimPrim[2] ; ix++) {
                (*rho_s_avg[ispec])(ix) += (*rho_s[ispec])(ix);
            }
        }//END loop on species ispec
    }

    // calculate the averaged values
    if ( time_step % dump_step == 0 ){
        for (unsigned int i=0 ; i<dim_global[0]*dim_global[1]*dim_global[2]; i++) {
            (*rho_global_avg)(i) /= avg_step;
            (*phi_global_avg)(i) /= avg_step;
            (*Ex_global_avg)(i)  /= avg_step;
        }
        for (unsigned int ispec=0; ispec<n_species; ispec++) {
            for (unsigned int ix=0 ; ix<dimPrim[0]*dimPrim[1] ; ix++) {
                (*rho_s_avg[ispec])(ix) /= avg_step;
            }
        }//END loop on species ispec
    }

}//END incrementAvgFields



// ---------------------------------------------------------------------------------------------------------------------
// Reinitialize the total charge densities and currents
// - save current density as old density (charge conserving scheme)
// - put the new density and currents to 0
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3D::restartRhoJ()
{
    // --------------------------
    // Total currents and density
    // --------------------------

    // static cast of the total currents and densities
    Field3D* Jx3D    = static_cast<Field3D*>(Jx_);
    Field3D* Jy3D    = static_cast<Field3D*>(Jy_);
    Field3D* Jz3D    = static_cast<Field3D*>(Jz_);
    Field3D* rho3D   = static_cast<Field3D*>(rho_);

    // Charge density rho^(p,p) to 0

    rho3D->put_to(0.0);
    Jx3D->put_to(0.0);
    Jy3D->put_to(0.0);
    Jz3D->put_to(0.0);

}//END restartRhoJ


void ElectroMagn3D::restartRhoJs(int ispec, bool currents)
{
    // -----------------------------------
    // Species currents and charge density
    // -----------------------------------
    Field3D* Jx3D_s  = static_cast<Field3D*>(Jx_s[ispec]);
    Field3D* Jy3D_s  = static_cast<Field3D*>(Jy_s[ispec]);
    Field3D* Jz3D_s  = static_cast<Field3D*>(Jz_s[ispec]);
    Field3D* rho3D_s = static_cast<Field3D*>(rho_s[ispec]);

    rho3D_s->put_to(0.0);

    if (currents){
        Jx3D_s->put_to(0.0);
        Jy3D_s->put_to(0.0);
        Jz3D_s->put_to(0.0);
    }
}//END restartRhoJs




// ---------------------------------------------------------------------------------------------------------------------
// Compute the total density and currents from species density and currents
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3D::computeTotalRhoJ()
{

    // static cast of the total currents and densities
    Field3D* Jx3D    = static_cast<Field3D*>(Jx_);
    Field3D* Jy3D    = static_cast<Field3D*>(Jy_);
    Field3D* Jz3D    = static_cast<Field3D*>(Jz_);
    Field3D* rho3D   = static_cast<Field3D*>(rho_);


    // -----------------------------------
    // Species currents and charge density
    // -----------------------------------
    for (unsigned int ispec=0; ispec<n_species; ispec++) {
        Field3D* Jx3D_s  = static_cast<Field3D*>(Jx_s[ispec]);
        Field3D* Jy3D_s  = static_cast<Field3D*>(Jy_s[ispec]);
        Field3D* Jz3D_s  = static_cast<Field3D*>(Jz_s[ispec]);
        Field3D* rho3D_s = static_cast<Field3D*>(rho_s[ispec]);

        // Charge density rho^(p,p) to 0
        for (unsigned int i=0 ; i<nx_p ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
              for (unsigned int k=0 ; k<nz_p ; k++) {
                (*rho3D)(i,j,k) += ( (*rho3D_s)(i,j,k) * species_param[ispec].charge );
              }
            }
        }
/*
        // Current Jx^(d,p) to 0
        for (unsigned int i=0 ; i<nx_p ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*Jx3D)(i,j) += (*Jx3D_s)(i,j);
            }
        }

        // Current Jy^(p,d) to 0
        for (unsigned int i=0 ; i<nx_p ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*Jy3D)(i,j) += (*Jy3D_s)(i,j);
            }
        }

        // Current Jz^(p,p) to 0
        for (unsigned int i=0 ; i<nx_p ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*Jz3D)(i,j) += (*Jz3D_s)(i,j);
            }
        }
*/
    }//END loop on species ispec

}//END computeTotalRhoJ


// ---------------------------------------------------------------------------------------------------------------------
// Compute electromagnetic energy flows vectors on the border of the simulation box
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3D::computePoynting() {

    if (isWestern) {
        unsigned int iEy=istart[0][Ey_->isDual(0)];
        unsigned int iBz=istart[0][Bz_m->isDual(0)];
        unsigned int iEz=istart[0][Ez_->isDual(0)];
        unsigned int iBy=istart[0][By_m->isDual(0)];

        unsigned int jEy=istart[1][Ey_->isDual(1)];
        unsigned int jBz=istart[1][Bz_m->isDual(1)];
        unsigned int jEz=istart[1][Ez_->isDual(1)];
        unsigned int jBy=istart[1][By_m->isDual(1)];

        for (unsigned int j=0; j<=bufsize[1][Ez_->isDual(1)]; j++) {

            double Ey__ = 0.5*((*Ey_)(iEy,jEy+j) + (*Ey_)(iEy, jEy+j+1));
            double Bz__ = 0.25*((*Bz_m)(iBz,jBz+j)+(*Bz_m)(iBz+1,jBz+j)+(*Bz_m)(iBz,jBz+j+1)+(*Bz_m)(iBz+1,jBz+j+1));
            double Ez__ = (*Ez_)(iEz,jEz+j);
            double By__ = 0.5*((*By_m)(iBy,jBy+j) + (*By_m)(iBy+1, jBy+j));

            poynting_inst[0][0] = dy*timestep*(Ey__*Bz__ - Ez__*By__);
            poynting[0][0]+= poynting_inst[0][0];
        }
    }//if Western


    if (isEastern) {
        unsigned int iEy=istart[0][Ey_->isDual(0)]  + bufsize[0][Ey_->isDual(0)] -1;
        unsigned int iBz=istart[0][Bz_m->isDual(0)] + bufsize[0][Bz_m->isDual(0)]-1;
        unsigned int iEz=istart[0][Ez_->isDual(0)]  + bufsize[0][Ez_->isDual(0)] -1;
        unsigned int iBy=istart[0][By_m->isDual(0)] + bufsize[0][By_m->isDual(0)]-1;

        unsigned int jEy=istart[1][Ey_->isDual(1)];
        unsigned int jBz=istart[1][Bz_m->isDual(1)];
        unsigned int jEz=istart[1][Ez_->isDual(1)];
        unsigned int jBy=istart[1][By_m->isDual(1)];

        for (unsigned int j=0; j<=bufsize[1][Ez_->isDual(1)]; j++) {

            double Ey__ = 0.5*((*Ey_)(iEy,jEy+j) + (*Ey_)(iEy, jEy+j+1));
            double Bz__ = 0.25*((*Bz_m)(iBz,jBz+j)+(*Bz_m)(iBz+1,jBz+j)+(*Bz_m)(iBz,jBz+j+1)+(*Bz_m)(iBz+1,jBz+j+1));
            double Ez__ = (*Ez_)(iEz,jEz+j);
            double By__ = 0.5*((*By_m)(iBy,jBy+j) + (*By_m)(iBy+1, jBy+j));

            poynting_inst[1][0] = dy*timestep*(Ey__*Bz__ - Ez__*By__);
            poynting[1][0] -= poynting_inst[1][0];
        }
    }//if Easter

    if (isSouthern) {

        unsigned int iEz=istart[0][Ez_->isDual(0)];
        unsigned int iBx=istart[0][Bx_m->isDual(0)];
        unsigned int iEx=istart[0][Ex_->isDual(0)];
        unsigned int iBz=istart[0][Bz_m->isDual(0)];

        unsigned int jEz=istart[1][Ez_->isDual(1)];
        unsigned int jBx=istart[1][Bx_m->isDual(1)];
        unsigned int jEx=istart[1][Ex_->isDual(1)];
        unsigned int jBz=istart[1][Bz_m->isDual(1)];

        for (unsigned int i=0; i<=bufsize[0][Ez_->isDual(0)]; i++) {
            double Ez__ = (*Ez_)(iEz+i,jEz);
            double Bx__ = 0.5*((*Bx_m)(iBx+i,jBx) + (*Bx_m)(iBx+i, jBx+1));
            double Ex__ = 0.5*((*Ex_)(iEx+i,jEx) + (*Ex_)(iEx+i+1, jEx));
            double Bz__ = 0.25*((*Bz_m)(iBz+i,jBz)+(*Bz_m)(iBz+i+1,jBz)+(*Bz_m)(iBz+i,jBz+1)+(*Bz_m)(iBz+i+1,jBz+1));

            poynting_inst[0][1] = dx*timestep*(Ez__*Bx__ - Ex__*Bz__);
            poynting[0][1] += poynting_inst[0][1];
        }
    }// if South

    if (isNorthern) {
        unsigned int iEz=istart[0][Ez_->isDual(0)];
        unsigned int iBx=istart[0][Bx_m->isDual(0)];
        unsigned int iEx=istart[0][Ex_->isDual(0)];
        unsigned int iBz=istart[0][Bz_m->isDual(0)];

        unsigned int jEz=istart[1][Ez_->isDual(1)]  + bufsize[1][Ez_->isDual(1)] -1;
        unsigned int jBx=istart[1][Bx_m->isDual(1)] + bufsize[1][Bx_m->isDual(1)]-1;
        unsigned int jEx=istart[1][Ex_->isDual(1)]  + bufsize[1][Ex_->isDual(1)] -1;
        unsigned int jBz=istart[1][Bz_m->isDual(1)] + bufsize[1][Bz_m->isDual(1)]-1;

        for (unsigned int i=0; i<=bufsize[0][Ez_->isDual(0)]; i++) {
            double Ez__ = (*Ez_)(iEz+i,jEz);
            double Bx__ = 0.5*((*Bx_m)(iBx+i,jBx) + (*Bx_m)(iBx+i, jBx+1));
            double Ex__ = 0.5*((*Ex_)(iEx+i,jEx) + (*Ex_)(iEx+i+1, jEx));
            double Bz__ = 0.25*((*Bz_m)(iBz+i,jBz)+(*Bz_m)(iBz+i+1,jBz)+(*Bz_m)(iBz+i,jBz+1)+(*Bz_m)(iBz+i+1,jBz+1));

            poynting_inst[1][1] = dx*timestep*(Ez__*Bx__ - Ex__*Bz__);
            poynting[1][1] -= poynting_inst[1][1];
        }
    }//if North

}

void ElectroMagn3D::gatherFields(SmileiMPI *smpi)
{
    SmileiMPI_Cart3D* smpi3D = static_cast<SmileiMPI_Cart3D*>(smpi);
    for(int i = 0; i < rho_s.size(); i++)
    {
        smpi3D->gatherRho( static_cast<Field3D*>(rho_s_global[i]), static_cast<Field3D*>(rho_s[i]) );
    }

}


void ElectroMagn3D::gatherAvgFields(SmileiMPI *smpi)
{
    SmileiMPI_Cart3D* smpi3D = static_cast<SmileiMPI_Cart3D*>(smpi);
    for(int i = 0; i < rho_s.size(); i++)
    {
        smpi3D->gatherRho( static_cast<Field3D*>(rho_s_global_avg[i]), static_cast<Field3D*>(rho_s_avg[i]) );
    }

}
