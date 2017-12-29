#include "Species.h"

#include <cmath>
#include <ctime>
#include <cstdlib>
#include <algorithm>

#include <iostream>

#include <omp.h>

#include "PusherFactory.h"

#include "PartBoundCond.h"
//#include "BoundaryConditionType.h"

#include "ElectroMagn.h"
#include "Interpolator.h"
#include "InterpolatorFactory.h"
#include "Profile.h"
#include "Grid.h"
#include "Grid2D.h"
#include "Projector.h"

#include "SmileiMPI.h"
//#include "SimWindow.h"

// #include "Field.h"
#include "Field1D.h"
#include "Field2D.h"
#include "Field3D.h"
#include "Tools.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Creator for Species
// input: simulation parameters & Species index
// ---------------------------------------------------------------------------------------------------------------------
Species::Species(PicParams& params, int ispec, SmileiMPI* smpi) :
speciesNumber(ispec),
clrw(params.clrw),
oversize(params.oversize),
cell_length(params.cell_length),
species_param(params.species_param[ispec]),
velocityProfile(3,NULL),
temperatureProfile(3,NULL),
ndim(params.nDim_particle),
min_loc(smpi->getDomainLocalMin(0))
{

    particles.species_number = speciesNumber;

    densityProfileType = species_param.density_type;
    chargeProfile         = new Profile(species_param.charge_profile, params.geometry);
    densityProfile        = new Profile(species_param.dens_profile  , params.geometry);

    velocityProfile[0]    = new Profile(species_param.mvel_x_profile, params.geometry);
    velocityProfile[1]    = new Profile(species_param.mvel_y_profile, params.geometry);
    velocityProfile[2]    = new Profile(species_param.mvel_z_profile, params.geometry);

    temperatureProfile[0] = new Profile(species_param.temp_x_profile, params.geometry);
    temperatureProfile[1] = new Profile(species_param.temp_y_profile, params.geometry);
    temperatureProfile[2] = new Profile(species_param.temp_z_profile, params.geometry);

    // ppc means "particles per cell"
    ppcProfile            = new Profile(species_param.ppc_profile   , params.geometry);

    // -------------------
    // Variable definition
    // -------------------
    PI2 = 2.0 * M_PI;

    DEBUG(species_param.species_type);

    electron_species = NULL;

    // atomic number

    // Width of clusters:
    if (params.n_space[0]%clrw != 0)
        ERROR("clrw should divide n_space[0] = " << params.n_space[0] );

    // Arrays of the min and max indices of the particle bins
    bmin.resize(params.n_space[0]/clrw);
    bmax.resize(params.n_space[0]/clrw);
    if (ndim == 3){
        bmin.resize(params.n_space[0]/clrw*params.n_space[1]);
        bmax.resize(params.n_space[0]/clrw*params.n_space[1]);
    }

    //Size in each dimension of the buffers on which each bin are projected
    //In 1D the particles of a given bin can be projected on 6 different nodes at the second order (oversize = 2)

    //Primal dimension of fields.
    f_dim0 =  params.n_space[0] + 2 * oversize[0] +1;
    f_dim1 =  params.n_space[1] + 2 * oversize[1] +1;
    f_dim2 =  params.n_space[2] + 2 * oversize[2] +1;

    if (ndim == 1){
        b_dim0 =  (1 + clrw) + 2 * oversize[0];
        b_dim1 =  1;
        b_dim2 =  1;
        b_lastdim = b_dim0;
    }
    if (ndim == 2){
        b_dim0 =  (1 + clrw) + 2 * oversize[0]; // There is a primal number of bins.
        b_dim1 =  f_dim1;
        b_dim2 =  1;
        b_lastdim = b_dim1;
    }
    if (ndim == 3){
        b_dim0 =  (1 + clrw) + 2 * oversize[0]; // There is a primal number of bins.
        b_dim1 = f_dim1;
        b_dim2 = f_dim2;
        b_lastdim = b_dim2;
    }

    size_proj_buffer = b_dim0*b_dim1*b_dim2;


    if (!params.restart) {

    // unsigned int npart_effective=0;

        // Create particles in a space starting at cell_index
        vector<double> cell_index(3,0);
        for (unsigned int i=0 ; i<params.nDim_field ; i++) {
            if (cell_length[i]!=0)
                cell_index[i] = smpi->getDomainLocalMin(i);
        }

        int starting_bin_idx = 0;
        // does a loop over all cells in the simulation
        // considering a 3d volume with size n_space[0]*n_space[1]*n_space[2]
        //*npart_effective =
        createParticles(params.n_space, cell_index, starting_bin_idx, params );

        //PMESSAGE( 1, smpi->getRank(),"Species "<< speciesNumber <<" # part "<< npart_effective );
    }

    // assign the correct Pusher to Push
    Push = PusherFactory::create( params, ispec );



    // define limits for BC and functions applied and for domain decomposition
    partBoundCond = new PartBoundCond( params, ispec, smpi);

    indexes_of_particles_to_perform_psi.resize(2 * params.nDim_particle);
    //ener_tot = 0.;
    nrj_bc_lost = 0.;
    nrj_mw_lost = 0.;
    nrj_new_particles = 0.;


}//END Species creator


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Species
// ---------------------------------------------------------------------------------------------------------------------
Species::~Species()
{
    delete Push;
    if (partBoundCond) delete partBoundCond;
    if (chargeProfile) delete chargeProfile;
    if (densityProfile) delete densityProfile;
    for (unsigned int i=0; i<velocityProfile.size(); i++)
        delete velocityProfile[i];
    for (unsigned int i=0; i<temperatureProfile.size(); i++)
        delete temperatureProfile[i];

    DEBUG(10,"Species deleted ");
}



// ---------------------------------------------------------------------------------------------------------------------
// For all (np) particles in a mesh initialize its numerical weight (equivalent to a number density)
// ---------------------------------------------------------------------------------------------------------------------
void Species::initWeight(unsigned int nPart, unsigned int ispec, unsigned int iPart, double density)
{
    for (unsigned  p= iPart; p<iPart+nPart; p++) {
        particles.weight(p) = density / nPart;
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// For all (np) particles in a mesh initialize its numerical weight (equivalent to constant)
// ---------------------------------------------------------------------------------------------------------------------
void Species::initWeight_constant(unsigned int nPart, unsigned int ispec, unsigned int iPart, double weight_const)
{
    for (unsigned  p= iPart; p<iPart+nPart; p++) {
        particles.weight(p) = weight_const;
    }
}


// ---------------------------------------------------------------------------------------------------------------------
// For all (np) particles in a mesh initialize its charge state
// ---------------------------------------------------------------------------------------------------------------------
void Species::initCharge(unsigned int nPart, unsigned int ispec, unsigned int iPart, double q)
{
    //short Z = (short)q;
    //double r = q-(double)Z;

    // charges of all particles in a cell are the same
    for (unsigned int p = iPart; p<iPart+nPart; p++)
    {
        particles.charge(p) = q;
    }


    /*
    // if charge is integer, then all particles have the same charge
    if ( r == 0. ) {
        for (unsigned int p = iPart; p<iPart+nPart; p++)
            particles.charge(p) = Z;
    // if charge is not integer, then particles can have two different charges
    } else {
        int tot = 0, Nm, Np;
        double rr=r/(1-r), diff;
        Np = (int)round(r*(double)nPart);
        Nm = (int)nPart - Np;
        for (unsigned int p = iPart; p<iPart+nPart; p++) {
            if (Np > rr*Nm) {
                particles.charge(p) = Z+1;
                Np--;
            } else {
                particles.charge(p) = Z;
                Nm--;
            }
            tot += particles.charge(p);
        }
        diff = ((double)nPart)*q - (double)tot; // missing charge
        if (diff != 0.) {
            WARNING("Could not match exactly charge="<<q<<" for species #"<<ispec<<" (difference of "<<diff<<"). Try to add particles.");
        }
    }
    */
}



// ---------------------------------------------------------------------------------------------------------------------
// For all (np) particles in a mesh initialize their position
//   - either using regular distribution in the mesh (initPosition_type = regular)
//   - or using uniform random distribution (initPosition_type = random)
// ---------------------------------------------------------------------------------------------------------------------
void Species::initPosition(unsigned int nPart, unsigned int iPart, double *indexes, unsigned int ndim,
                           std::vector<double> cell_length, string initPosition_type)
{
    for (unsigned  p= iPart; p<iPart+nPart; p++) {
        for (unsigned  i=0; i<ndim ; i++) {

            // define new position (either regular or random)
            if (initPosition_type == "regular") {
                particles.position(i,p)=indexes[i]+(p-iPart+0.5)*cell_length[i]/nPart;
            } else if (initPosition_type == "random") {
                particles.position(i,p)=indexes[i]+(((double)rand() / RAND_MAX))*cell_length[i];
            }
            particles.position_old(i,p) = particles.position(i,p);
        }// i
    }// p
}



// ---------------------------------------------------------------------------------------------------------------------
// For all (np) particles in a mesh initialize their momentum
//   - at zero (init_momentum_type = cold)
//   - using random distribution (init_momentum_type = maxwell-juettner)
// ---------------------------------------------------------------------------------------------------------------------
void Species::initMomentum(unsigned int nPart, unsigned int iPart, double *temp, double *vel, string initMomentum_type,
                           vector<double>& max_jutt_cumul, PicParams& params)
{

    // average mean-momentum (used to center the distribution)
    double pMean[3]= {0.0,0.0,0.0};

    if (initMomentum_type == "cold") {

        for (unsigned int p= iPart; p<iPart+nPart; p++) {
            for (unsigned int i=0; i<3 ; i++) {
                particles.momentum(i,p) = 0.0;
            }
        }

    } else if (initMomentum_type == "maxwell")
    {
        // initialize using the Maxwell distribution function

        for (unsigned int p= iPart; p<iPart+nPart; p++)
        {
            /*
            double ran;
            do {
                ran = (double)rand() / RAND_MAX;
            }
            while (ran == 0.0);

            double vt = sqrt(2.0 * temp[0] * params.const_e / species_param.mass);
            double psm = vt * sqrt(-log(ran));

            double phi     = 2.0 * M_PI*(double)rand() / RAND_MAX;
            double theta   = M_PI*(double)rand() / RAND_MAX;
            particles.momentum(0,p) = psm*sin(theta)*cos(phi);
            particles.momentum(1,p) = psm*sin(theta)*sin(phi);
            particles.momentum(2,p) = psm*cos(theta);
            */

            /*
            double ran;
            do {
                ran = (double)rand() / RAND_MAX;
            }
            while (ran == 0.0);

            double vt = sqrt(2.0 * temp[0] * params.const_e / species_param.mass);
            double psm = vt * sqrt(-log(ran));


            double phi   = 2.0 * M_PI*(double)rand() / RAND_MAX;
            particles.momentum(0,p) = psm*sin(phi);

            phi   = 2.0 * M_PI*(double)rand() / RAND_MAX;
            particles.momentum(1,p) = psm*sin(phi);

            phi   = 2.0 * M_PI*(double)rand() / RAND_MAX;
            particles.momentum(2,p) = psm*sin(phi);
            */

            double vt = sqrt(2.0 * temp[0] * params.const_e / species_param.mass);
            double x1;
            double x2;

            do {
                x1 = (double)rand() / RAND_MAX;
            }
            while (x1 == 0.0);
            x2 = (double)rand() / RAND_MAX;
            particles.momentum(0,p) = vt * sqrt( -log(x1) ) * sin(2.0 * M_PI * x2);

            do {
                x1 = (double)rand() / RAND_MAX;
            }
            while (x1 == 0.0);
            x2 = (double)rand() / RAND_MAX;
            particles.momentum(1,p) = vt * sqrt( -log(x1) ) * sin(2.0 * M_PI * x2);

            do {
                x1 = (double)rand() / RAND_MAX;
            }
            while (x1 == 0.0);
            x2 = (double)rand() / RAND_MAX;
            particles.momentum(2,p) = vt * sqrt( -log(x1) ) * sin(2.0 * M_PI * x2);

        }


    }else if (initMomentum_type == "maxwell-juettner")
    {
        // initialize using the Maxwell-Juettner distribution function

        for (unsigned int p= iPart; p<iPart+nPart; p++)
        {
            double Renergy=(double)rand() / RAND_MAX;
            double phi=acos(1.0-2.0*(double)rand() / RAND_MAX);
            double theta=2.0*M_PI*(double)rand() / RAND_MAX;

            int il=0;
            int ir=max_jutt_cumul.size();
            while (ir > il+1)  {
                int im=(il+ir)/2;
                if (Renergy > max_jutt_cumul[im]) {
                    il=im;
                } else {
                    ir=im;
                }
            }
            double right_w=(Renergy-max_jutt_cumul[il])/(max_jutt_cumul[il+1]);
            double left_w=1-right_w;

            double Ener=left_w*il*dE +right_w*(il+1)*dE;
            double psm = sqrt(pow(1.0+Ener,2)-1.0);

            particles.momentum(0,p) = psm*cos(theta)*sin(phi);
            particles.momentum(1,p) = psm*sin(theta)*sin(phi);
            particles.momentum(2,p) = psm*cos(phi);
            for (unsigned int i=0; i<3 ; i++)
            {
                pMean[i] += particles.momentum(i,p);
            }
        }//p

        // center the distribution function around pMean
        for (unsigned int p= iPart; p<iPart+nPart; p++)
        {
            for (unsigned int i=0; i<3 ; i++) {
                particles.momentum(i,p) -= pMean[i]/nPart;
            }
        }

        for (unsigned int p= iPart; p<iPart+nPart; p++) {
            particles.momentum(1,p) *= sqrt(temp[1]/temp[0]);
            particles.momentum(2,p) *= sqrt(temp[2]/temp[0]);
        }

    // Rectangular distribution
    } else if (initMomentum_type == "rectangular") {

        for (unsigned int p= iPart; p<iPart+nPart; p++) {
            particles.momentum(0,p) = 2.0 * (2.*(double)rand() / RAND_MAX - 1.) * sqrt(2.0 * temp[0] * params.const_e / species_param.mass);
            particles.momentum(1,p) = (0.00001 * 2.*(double)rand() / RAND_MAX - 1.) * sqrt(2.0 * temp[0] * params.const_e / species_param.mass);
            particles.momentum(2,p) = (0.00001 * 2.*(double)rand() / RAND_MAX - 1.) * sqrt(2.0 * temp[0] * params.const_e / species_param.mass);
            //if( isinf(particles.momentum(0,p)) || isinf(particles.momentum(1,p)) || isinf(particles.momentum(2,p)) ) {
            //    cout<<particles.momentum(0,p)<<" "<<particles.momentum(1,p)<<" "<<particles.momentum(2,p)<<endl;
            //}
        }
    }//END if initMomentum_type

/*
    // Adding the mean velocity (using relativistic composition)
    // ---------------------------------------------------------
    if ( (vel[0]!=0.0) || (vel[1]!=0.0) || (vel[2]!=0.0) ){

        double vx, vy, vz, v2, g, gm1, Lxx, Lyy, Lzz, Lxy, Lxz, Lyz, gp, px, py, pz;

        // mean-velocity
        vx  = -vel[0];
        vy  = -vel[1];
        vz  = -vel[2];
        v2  = vx*vx + vy*vy + vz*vz;
        g   = 1.0/sqrt(1.0-v2);
        gm1 = g - 1.0;

        // compute the different component of the Matrix block of the Lorentz transformation
        Lxx = 1.0 + gm1 * vx*vx/v2;
        Lyy = 1.0 + gm1 * vy*vy/v2;
        Lzz = 1.0 + gm1 * vz*vz/v2;
        Lxy = gm1 * vx*vy/v2;
        Lxz = gm1 * vx*vz/v2;
        Lyz = gm1 * vy*vz/v2;

        // Lorentz transformation of the momentum
        for (unsigned int p=iPart; p<iPart+nPart; p++)
        {
            gp = sqrt(1.0 + pow(particles.momentum(0,p),2) + pow(particles.momentum(1,p),2) + pow(particles.momentum(2,p),2));
            px = -gp*g*vx + Lxx * particles.momentum(0,p) + Lxy * particles.momentum(1,p) + Lxz * particles.momentum(2,p);
            py = -gp*g*vy + Lxy * particles.momentum(0,p) + Lyy * particles.momentum(1,p) + Lyz * particles.momentum(2,p);
            pz = -gp*g*vz + Lxz * particles.momentum(0,p) + Lyz * particles.momentum(1,p) + Lzz * particles.momentum(2,p);
            particles.momentum(0,p) = px;
            particles.momentum(1,p) = py;
            particles.momentum(2,p) = pz;
        }

    }//ENDif vel != 0

*/

    if ( (vel[0]!=0.0) || (vel[1]!=0.0) || (vel[2]!=0.0) ){
        for (unsigned int p=iPart; p<iPart+nPart; p++)
        {
            particles.momentum(0,p) += vel[0];
            particles.momentum(1,p) += vel[1];
            particles.momentum(2,p) += vel[2];
        }

    }//ENDif vel != 0




}//END initMomentum


// ---------------------------------------------------------------------------------------------------------------------
// For all (np) particles in a mesh initialize their momentum
//   - at zero (init_momentum_type = cold)
//   - using random distribution (init_momentum_type = maxwell-juettner)
// ---------------------------------------------------------------------------------------------------------------------
void Species::heat(unsigned int nPart_bin, unsigned int nPart, unsigned int iPart, double temp, PicParams& params)
{
    double momentum_unit[3];
    double momentum_magnitude;
    double momentum_magnitude_add;

    vector<unsigned int> index;
    index.resize(nPart_bin);
    for(int i = 0; i < nPart_bin; i++)
    {
      index[i] = i;
    }
    random_shuffle(index.begin(), index.end());

    for ( int i = 0; i < nPart; i++)
    {
        unsigned int p = iPart + index[i];
        momentum_magnitude = sqrt( particles.momentum(0,p) * particles.momentum(0,p)
                                 + particles.momentum(1,p) * particles.momentum(1,p)
                                 + particles.momentum(2,p) * particles.momentum(2,p) );
        if(momentum_magnitude == 0.0)
        {
          continue;
        }
        momentum_unit[0] = particles.momentum(0,p) / momentum_magnitude;
        momentum_unit[1] = particles.momentum(1,p) / momentum_magnitude;
        momentum_unit[2] = particles.momentum(2,p) / momentum_magnitude;

        momentum_magnitude_add = sqrt(3.0 * temp * params.const_e / species_param.mass + momentum_magnitude * momentum_magnitude) - momentum_magnitude;
        momentum_magnitude += momentum_magnitude_add;

        particles.momentum(0,p) = momentum_unit[0] * momentum_magnitude;
        particles.momentum(1,p) = momentum_unit[1] * momentum_magnitude;
        particles.momentum(2,p) = momentum_unit[2] * momentum_magnitude;
    }


}//END heat





// ---------------------------------------------------------------------------------------------------------------------
// For all (np) particles in a mesh initialize their acceleration for implicit method
// ---------------------------------------------------------------------------------------------------------------------
void Species::initAcceleration_imp(unsigned int nPart, unsigned int iPart)
{
    for (unsigned int p= iPart; p<iPart+nPart; p++) {
        for (unsigned int i=0; i<3 ; i++) {
            particles.al_imp(i,p) = 0.0;
            particles.au_imp(i,p) = 0.0;
        }
    }

}//END initAcceleration



// ---------------------------------------------------------------------------------------------------------------------
// For all particles of the species
//   - interpolate the fields at the particle position
//   - calculate the new velocity
//   - calculate the new position
//   - apply the boundary conditions
//   - increment the currents (projection)
// ---------------------------------------------------------------------------------------------------------------------
void Species::dynamics(double time_dual, unsigned int ispec, ElectroMagn* EMfields, Interpolator* Interp,
                       Projector* Proj, SmileiMPI *smpi, PicParams &params)
{
    Interpolator* LocInterp = InterpolatorFactory::create(params, smpi);

    // Electric field at the particle position
    LocalFields Epart;
    // Magnetic field at the particle position
    LocalFields Bpart;
    // Ionization current
    LocalFields Jion;

    int iloc;
    unsigned int i,j,ibin,iPart;

    //! buffers for currents and charge
    double *b_Jx,*b_Jy,*b_Jz,*b_rho;

    // number of particles for this Species
    unsigned int nParticles = getNbrOfParticles();
    // Reset list of particles to exchange

    int iDirection=-1;

    clearExchList();

    //ener_tot  = 0.;
    //ener_lost = 0.;
    double ener_iPart(0.);
    //bool contribute(true);
    // -------------------------------
    // calculate the particle dynamics
    // -------------------------------
    if (time_dual>species_param.time_frozen) { // moving particle
        double gf = 1.0;

        psi_particles.clear();
        for(int iD=0; iD<indexes_of_particles_to_perform_psi.size(); iD++)
        {
             indexes_of_particles_to_perform_psi[iD].clear();
        }
        for (ibin = 0 ; ibin < (unsigned int)bmin.size() ; ibin++) {
            for (iPart=(unsigned int)bmin[ibin] ; iPart<(unsigned int)bmax[ibin]; iPart++ ) {

                // Interpolate the fields at the particle position
                //(*LocInterp)(EMfields, particles, iPart, &Epart);
                (*LocInterp)(EMfields, particles, iPart, &Epart, &Bpart);

                // Push the particle
                //(*Push)(particles, iPart, Epart);
                (*Push)(particles, iPart, Epart, Bpart);
                //(*Push)(particles, iPart, Epart, Bpart, gf);

                // Apply boundary condition on the particles
                // Boundary Condition may be physical or due to domain decomposition
                // apply returns 0 if iPart is no more in the domain local
                // if omp, create a list per thread
                if ( !partBoundCond->apply( particles, iPart, params.species_param[ispec], ener_iPart, iDirection ) ) {
                    addPartInExchList( iPart );
                    if(iDirection >= 0){
                        addPartInPsiList( iDirection, iPart );
                    }
                }
            }//iPart
        }// ibin

        //for (iPart=0 ; iPart<nParticles; iPart++ ) {
        //    (*Proj)(EMfields->rho_s[ispec], particles, iPart);
        //}

        // copy PSI particles to psi_particles, because after MPi particle exchanging
        // the PSI particles will be erased
        for(int iD=0; iD<indexes_of_particles_to_perform_psi.size(); iD++)
        {
            for(int iPart=0; iPart<indexes_of_particles_to_perform_psi[iD].size(); iPart++)
            {
                int iPart_psi = indexes_of_particles_to_perform_psi[iD][iPart];
                particles.cp_particle(iPart_psi, psi_particles);
            }
        }

    }
    else if (!particles.isTestParticles) { // immobile particle (at the moment only project density)
        for (iPart=0 ; iPart<nParticles; iPart++ ) {
            (*Proj)(EMfields->rho_s[ispec], particles, iPart);
        }
    }//END if time vs. time_frozen

    delete LocInterp;

}//END dynamic




// ---------------------------------------------------------------------------------------------------------------------
// For all particles of the species
//   - implicit method: first push
// ---------------------------------------------------------------------------------------------------------------------
void Species::dynamics_imp_firstPush(double time_dual, unsigned int ispec, ElectroMagn* EMfields, Interpolator* Interp,
                       Projector* Proj, SmileiMPI *smpi, PicParams &params)
{
    Interpolator* LocInterp = InterpolatorFactory::create(params, smpi);

    // Electric field at the particle position
    LocalFields Epart;
    // Magnetic field at the particle position
    LocalFields Bpart;
    // Ionization current
    LocalFields Jion;

    int iloc;
    unsigned int i,j,ibin,iPart;

    //! buffers for currents and charge
    double *b_Jx,*b_Jy,*b_Jz,*b_rho;

    // number of particles for this Species
    unsigned int nParticles = getNbrOfParticles();
    // Reset list of particles to exchange

    int iDirection=-1;

    clearExchList();

    //ener_tot  = 0.;
    //ener_lost = 0.;
    double ener_iPart(0.);
    //bool contribute(true);
    // -------------------------------
    // calculate the particle dynamics
    // -------------------------------
    if (time_dual>species_param.time_frozen) { // moving particle
        double gf = 1.0;

        psi_particles.clear();
        for(int iD=0; iD<indexes_of_particles_to_perform_psi.size(); iD++)
        {
             indexes_of_particles_to_perform_psi[iD].clear();
        }
        for (ibin = 0 ; ibin < (unsigned int)bmin.size() ; ibin++) {
            for (iPart=(unsigned int)bmin[ibin] ; iPart<(unsigned int)bmax[ibin]; iPart++ ) {

                // Interpolate the fields at the particle position
                //(*LocInterp)(EMfields, particles, iPart, &Epart);
                (*LocInterp)(EMfields, particles, iPart, &Epart, &Bpart);

                // Push the particle
                //(*Push)(particles, iPart, Epart);
                Push->firstPush(particles, iPart, Epart);
                //(*Push)(particles, iPart, Epart, Bpart, gf);

                // Apply boundary condition on the particles
                // Boundary Condition may be physical or due to domain decomposition
                // apply returns 0 if iPart is no more in the domain local
                // if omp, create a list per thread
                if ( !partBoundCond->apply( particles, iPart, params.species_param[ispec], ener_iPart, iDirection ) ) {
                    addPartInExchList( iPart );
                    if(iDirection >= 0){
                        addPartInPsiList( iDirection, iPart );
                    }
                }
            }//iPart
        }// ibin

        //for (iPart=0 ; iPart<nParticles; iPart++ ) {
        //    (*Proj)(EMfields->rho_s[ispec], particles, iPart);
        //}

        // copy PSI particles to psi_particles, because after MPi particle exchanging
        // the PSI particles will be erased
        for(int iD=0; iD<indexes_of_particles_to_perform_psi.size(); iD++)
        {
            for(int iPart=0; iPart<indexes_of_particles_to_perform_psi[iD].size(); iPart++)
            {
                int iPart_psi = indexes_of_particles_to_perform_psi[iD][iPart];
                particles.cp_particle(iPart_psi, psi_particles);
            }
        }

    }
    else if (!particles.isTestParticles) { // immobile particle (at the moment only project density)
        for (iPart=0 ; iPart<nParticles; iPart++ ) {
            (*Proj)(EMfields->rho_s[ispec], particles, iPart);
        }
    }//END if time vs. time_frozen

    delete LocInterp;

}//END dynamic


// ---------------------------------------------------------------------------------------------------------------------
// For all particles of the species
//   - implicit method: second Push
// ---------------------------------------------------------------------------------------------------------------------
void Species::dynamics_imp_secondPush(double time_dual, unsigned int ispec, ElectroMagn* EMfields, Interpolator* Interp,
                       Projector* Proj, SmileiMPI *smpi, PicParams &params)
{
    Interpolator* LocInterp = InterpolatorFactory::create(params, smpi);

    // Electric field at the particle position
    LocalFields Epart;
    // Magnetic field at the particle position
    LocalFields Bpart;
    // Ionization current
    LocalFields Jion;

    int iloc;
    unsigned int i,j,ibin,iPart;

    //! buffers for currents and charge
    double *b_Jx,*b_Jy,*b_Jz,*b_rho;

    // number of particles for this Species
    unsigned int nParticles = getNbrOfParticles();
    // Reset list of particles to exchange

    int iDirection=-1;

    clearExchList();

    //ener_tot  = 0.;
    //ener_lost = 0.;
    double ener_iPart(0.);
    //bool contribute(true);
    // -------------------------------
    // calculate the particle dynamics
    // -------------------------------
    if (time_dual>species_param.time_frozen) { // moving particle
        double gf = 1.0;

        for(int iD=0; iD<indexes_of_particles_to_perform_psi.size(); iD++)
        {
             indexes_of_particles_to_perform_psi[iD].clear();
        }
        for (ibin = 0 ; ibin < (unsigned int)bmin.size() ; ibin++) {
            for (iPart=(unsigned int)bmin[ibin] ; iPart<(unsigned int)bmax[ibin]; iPart++ ) {

                // Interpolate the fields at the particle position
                //(*LocInterp)(EMfields, particles, iPart, &Epart);
                (*LocInterp)(EMfields, particles, iPart, &Epart, &Bpart);

                // Push the particle
                //(*Push)(particles, iPart, Epart);
                Push->secondPush(particles, iPart, Epart);
                //(*Push)(particles, iPart, Epart, Bpart, gf);

                // Apply boundary condition on the particles
                // Boundary Condition may be physical or due to domain decomposition
                // apply returns 0 if iPart is no more in the domain local
                // if omp, create a list per thread
                if ( !partBoundCond->apply( particles, iPart, params.species_param[ispec], ener_iPart, iDirection ) ) {
                    addPartInExchList( iPart );
                    if(iDirection >= 0){
                        addPartInPsiList( iDirection, iPart );
                    }
                }
            }//iPart
        }// ibin

        //for (iPart=0 ; iPart<nParticles; iPart++ ) {
        //    (*Proj)(EMfields->rho_s[ispec], particles, iPart);
        //}

        // copy PSI particles to psi_particles, because after MPi particle exchanging
        // the PSI particles will be erased
        for(int iD=0; iD<indexes_of_particles_to_perform_psi.size(); iD++)
        {
            for(int iPart=0; iPart<indexes_of_particles_to_perform_psi[iD].size(); iPart++)
            {
                int iPart_psi = indexes_of_particles_to_perform_psi[iD][iPart];
                particles.cp_particle(iPart_psi, psi_particles);
            }
        }

    }
    else if (!particles.isTestParticles) { // immobile particle (at the moment only project density)
        for (iPart=0 ; iPart<nParticles; iPart++ ) {
            (*Proj)(EMfields->rho_s[ispec], particles, iPart);
        }
    }//END if time vs. time_frozen

    delete LocInterp;

}//END dynamic



void Species::dynamics_EM(double time_dual, unsigned int ispec, ElectroMagn* EMfields, Interpolator* Interp,
                       Projector* Proj, SmileiMPI *smpi, PicParams &params)
{
    Interpolator* LocInterp = InterpolatorFactory::create(params, smpi);

    // Electric field at the particle position
    LocalFields Epart;
    // Magnetic field at the particle position
    LocalFields Bpart;
    // Ionization current
    LocalFields Jion;

    int iloc;
    unsigned int i,j,ibin,iPart;

    //! buffers for currents and charge
    double *b_Jx,*b_Jy,*b_Jz,*b_rho;

    // number of particles for this Species
    unsigned int nParticles = getNbrOfParticles();
    // Reset list of particles to exchange

    int iDirection=-1;

    clearExchList();

    //ener_tot  = 0.;
    //ener_lost = 0.;
    double ener_iPart(0.);
    //bool contribute(true);
    // -------------------------------
    // calculate the particle dynamics
    // -------------------------------
    if (time_dual>species_param.time_frozen) { // moving particle
        double gf = 1.0;
        //Allocate buffer for projection  *****************************
        // *4 accounts for Jy, Jz and rho. * nthds accounts for each thread.
        b_Jx = (double *) malloc(4 * size_proj_buffer * sizeof(double));
        //Point buffers of each thread to the correct position
        b_Jy = b_Jx + size_proj_buffer ;
        b_Jz = b_Jy + size_proj_buffer ;
        b_rho = b_Jz + size_proj_buffer ;

        psi_particles.clear();
        for (ibin = 0 ; ibin < (unsigned int)bmin.size() ; ibin++) {

            // reset all current-buffers
            memset( &(b_Jx[0]), 0, 4*size_proj_buffer*sizeof(double));

            for (iPart=(unsigned int)bmin[ibin] ; iPart<(unsigned int)bmax[ibin]; iPart++ ) {


                //MESSAGE("ipart: "<<iPart);
                //cout<<"ipart: "<<iPart<<endl;
                // Interpolate the fields at the particle position
                //(*LocInterp)(EMfields, particles, iPart, &Epart);
                (*LocInterp)(EMfields, particles, iPart, &Epart, &Bpart);

                //if(Epart.x != 0.0 || Epart.y != 0.0 || Epart.z != 0.0 || Bpart.x != 0.0 || Bpart.y != 0.0 || Bpart.z != 0.0){
                //    cout<<iPart<<"local field not zero  "<<Epart.x<<" "<<Epart.y<<" "<<Epart.z<<" "<<Bpart.x<<" "<<Bpart.y<<" "<<Bpart.z<<endl;
                //}


                // Push the particle
                //(*Push)(particles, iPart, Epart);
                (*Push)(particles, iPart, Epart, Bpart);
                //(*Push)(particles, iPart, Epart, Bpart, gf);

                // Apply boundary condition on the particles
                // Boundary Condition may be physical or due to domain decomposition
                // apply returns 0 if iPart is no more in the domain local
                // if omp, create a list per thread
                if ( !partBoundCond->apply( particles, iPart, params.species_param[ispec], ener_iPart, iDirection ) ) {
                    addPartInExchList( iPart );
                    if(iDirection >= 0){
                        addPartInPsiList( iDirection, iPart );
                    }
                }

                //if (!particles.isTestParticles) {
                //    if (ndim <= 2) {
                //        //(*Proj)(b_Jx, b_Jy, b_Jz, b_rho, particles, iPart, gf, ibin*clrw, b_lastdim);
                //        (*Proj)(EMfields->rho_s[ispec], particles, iPart);
                //    } else {
                //        (*Proj)(EMfields->Jx_s[ispec], EMfields->Jy_s[ispec], EMfields->Jz_s[ispec],
                //                EMfields->rho_s[ispec],particles, iPart, gf);
                //    }
                //}
            }//iPart

/*
            // Copy buffer back to the global array and free buffer****************
            if (!particles.isTestParticles) {
                // this part is dimension dependant !! this is for dim = 1
                if (ndim == 1) {
                    for (i = 0; i < b_dim0 ; i++) {
                        //! \todo Should we care about primal - dual sizes here ?
                        iloc = ibin*clrw + i ;
                        (*EMfields->Jx_s[ispec]) (iloc) +=  b_Jx[i];
                        (*EMfields->Jy_s[ispec]) (iloc) +=  b_Jy[i];
                        (*EMfields->Jz_s[ispec]) (iloc) +=  b_Jz[i];
                        (*EMfields->rho_s[ispec])(iloc) += b_rho[i];
                    }
                } // End if (ndim == 1)
                if (ndim == 2) {
                    for (i = 0; i < 2*oversize[0]+1 ; i++) {
                        iloc = ibin*clrw + i ;
                        //! \todo Here b_dim0 is the dual size. Make sure no problems arise when i == b_dim0-1 for primal arrays.
                        for (j = 0; j < b_dim1 ; j++) {
                            (*EMfields->Jx_s[ispec]) (iloc*(f_dim1  )+j) +=  b_Jx[i*b_dim1+j];   //  primal along y
                            (*EMfields->Jy_s[ispec]) (iloc*(f_dim1+1)+j) +=  b_Jy[i*b_dim1+j];   //+1 because dual along y
                            (*EMfields->Jz_s[ispec]) (iloc*(f_dim1  )+j) +=  b_Jz[i*b_dim1+j];   // primal along y
                            (*EMfields->rho_s[ispec])(iloc*(f_dim1  )+j) += b_rho[i*b_dim1+j];   // primal along y
                        }
                    }
                    for (i = 2*oversize[0]+1; i < clrw ; i++) {
                        iloc = ibin*clrw + i ;
                        //! \todo Here b_dim0 is the dual size. Make sure no problems arise when i == b_dim0-1 for primal arrays.
                        for (j = 0; j < b_dim1 ; j++) {
                            (*EMfields->Jx_s[ispec]) (iloc*(f_dim1  )+j) +=  b_Jx[i*b_dim1+j];   //  primal along y
                            (*EMfields->Jy_s[ispec]) (iloc*(f_dim1+1)+j) +=  b_Jy[i*b_dim1+j];   //+1 because dual along y
                            (*EMfields->Jz_s[ispec]) (iloc*(f_dim1  )+j) +=  b_Jz[i*b_dim1+j];   // primal along y
                            (*EMfields->rho_s[ispec])(iloc*(f_dim1  )+j) += b_rho[i*b_dim1+j];   // primal along y
                        }
                    }
                    for (i = std::max(clrw,2*oversize[0]+1); i < b_dim0 ; i++) {
                        iloc = ibin*clrw + i ;
                        //! \todo Here b_dim0 is the dual size. Make sure no problems arise when i == b_dim0-1 for primal arrays.
                        for (j = 0; j < b_dim1 ; j++) {
                            (*EMfields->Jx_s[ispec]) (iloc*(f_dim1  )+j) +=  b_Jx[i*b_dim1+j];   //  primal along y
                            (*EMfields->Jy_s[ispec]) (iloc*(f_dim1+1)+j) +=  b_Jy[i*b_dim1+j];   //+1 because dual along y
                            (*EMfields->Jz_s[ispec]) (iloc*(f_dim1  )+j) +=  b_Jz[i*b_dim1+j];   // primal along y
                            (*EMfields->rho_s[ispec])(iloc*(f_dim1  )+j) += b_rho[i*b_dim1+j];   // primal along y
                        }
                    }
                } // End if (ndim == 2)
            } // if (!particles.isTestParticles)
*/

        }// ibin

        for (iPart=0 ; iPart<nParticles; iPart++ ) {
            (*Proj)(EMfields->rho_s[ispec], particles, iPart);
        }

        // copy PSI particles to psi_particles, because after MPi particle exchanging
        // the PSI particles will be erased
        for(int iDirection=0; iDirection<indexes_of_particles_to_perform_psi.size(); iDirection++)
        {
            for(int iPart=0; iPart<indexes_of_particles_to_perform_psi[iDirection].size(); iPart++)
            {
                int iPart_psi = indexes_of_particles_to_perform_psi[iDirection][iPart];
                particles.cp_particle(iPart_psi, psi_particles);
            }
        }

        free(b_Jx);

    }
    else if (!particles.isTestParticles) { // immobile particle (at the moment only project density)

        for (iPart=0 ; iPart<nParticles; iPart++ ) {
            (*Proj)(EMfields->rho_s[ispec], particles, iPart);
        }
    }//END if time vs. time_frozen

    delete LocInterp;

}//END dynamic


void Species::Project(double time_dual, unsigned int ispec, ElectroMagn* EMfields, Projector* Proj, SmileiMPI *smpi, PicParams &params)
{

    if (time_dual>species_param.time_frozen) { // moving particle

        for (int ibin = 0 ; ibin < (unsigned int)bmin.size() ; ibin++)
        {
            for (int iPart=(unsigned int)bmin[ibin] ; iPart<(unsigned int)bmax[ibin]; iPart++ )
            {
                (*Proj)(EMfields->rho_s[ispec], particles, iPart, species_param.weight);
                //(*Proj)(EMfields->rho_s[ispec], particles, iPart);
            }//iPart
        }// ibin

    }
    else if (!particles.isTestParticles) {

    }


}//END Project





void Species::absorb2D(double time_dual, unsigned int ispec, Grid* grid, SmileiMPI *smpi, PicParams &params)
{
    double xpn, ypn;
    int ic, jc;
    double dx_inv_, dy_inv_;
    int i_domain_begin, j_domain_begin;

    Grid2D* grid2D = static_cast<Grid2D*>(grid);

    SmileiMPI_Cart2D* smpi2D = static_cast<SmileiMPI_Cart2D*>(smpi);

    dx_inv_   = 1.0/params.cell_length[0];
    dy_inv_   = 1.0/params.cell_length[1];

    i_domain_begin = smpi2D->getCellStartingGlobalIndex(0);
    j_domain_begin = smpi2D->getCellStartingGlobalIndex(1);


    if (time_dual>species_param.time_frozen) { // moving particle

        indexes_of_particles_to_absorb.clear();

        for (int ibin = 0 ; ibin < (unsigned int)bmin.size() ; ibin++) {
            for (int iPart=(unsigned int)bmin[ibin] ; iPart<(unsigned int)bmax[ibin]; iPart++ ) {

                //Locate particle on the primal grid & calculate the projection coefficients
                xpn = particles.position(0, iPart) * dx_inv_;  // normalized distance to the first node
                ic  = floor(xpn);                   // index of the central node

                ypn = particles.position(1, iPart) * dy_inv_;  // normalized distance to the first node
                jc   = floor(ypn);                  // index of the central node

                int i = ic-i_domain_begin; // index of first point for projection in x
                int j = jc-j_domain_begin; // index of first point for projection in y

                if( grid2D->iswall_2D[i][j] == 1 && grid2D->iswall_2D[i+1][j] == 1 && grid2D->iswall_2D[i+1][j+1] == 1
                && grid2D->iswall_2D[i][j+1] == 1 ) {
                    indexes_of_particles_to_absorb.push_back(iPart);
                }

            }//iPart

        }// ibin

        // copy PSI particles to psi_particles, because after MPi particle exchanging
        // the PSI particles will be erased
        for(int iPart=0; iPart<indexes_of_particles_to_absorb.size(); iPart++)
        {
            int iPart_psi = indexes_of_particles_to_absorb[iPart];
            particles.cp_particle(iPart_psi, psi_particles);
        }
        erase_particles_from_bins(indexes_of_particles_to_absorb);
    }
    else if (!particles.isTestParticles) { // immobile particle (at the moment only project density)

    }//END if time vs. time_frozen


}//END absorb2D




// ---------------------------------------------------------------------------------------------------------------------
// Dump for the Species
//! \todo{Check if one wants to keep it like this (MG)}
// ---------------------------------------------------------------------------------------------------------------------
void Species::dump(std::ofstream& ofile)
{
    for (unsigned int i=0; i<(unsigned int)particles.size(); i++ )
    {
        ofile << i ;
        for (unsigned int m=0; m<ndim; m++) ofile << "\t" << particles.position(m,i);
        for (unsigned int m=0; m<3; m++)    ofile << "\t" << particles.momentum(m,i);
        ofile << "\t" << particles.weight(i); //<< "\t" << Push->getMass() << "\t" << Push->getCharge();
        ofile << endl;
    }
    ofile << endl;
}

// ---------------------------------------------------------------------------------------------------------------------
// Sort particles
// This method assume the particle displacement at one timestep is less than one cell_length
// When a particle loss in collision , can not directly move it to bin end and then sort_part
// ---------------------------------------------------------------------------------------------------------------------
void Species::sort_part()
{
    //The width of one bin is cell_length[0] * clrw.

    int p1,p2,bmin_init;
    unsigned int bin;
    double limit;
    int numSort;
    int count;

    count = 0;
    do {
        count++;
        //if(count>3) {
        //    cout<<"count of particles sort: "<<count<<endl;
        //}
        numSort = 0;
        //Backward pass
        for (bin=0; bin<bmin.size()-1; bin++) { //Loop on the bins.
            limit = min_loc + (bin+1)*cell_length[0]*clrw;
            p1 = bmax[bin]-1;
            //If first particles change bin, they do not need to be swapped.
            while (p1 == bmax[bin]-1 && p1 >= bmin[bin]) {
                if (particles.position(0,p1) >= limit ) {
                    bmax[bin]--;
                    numSort++;
                }
                p1--;
            }
            //         Now particles have to be swapped
            for( p2 = p1 ; p2 >= bmin[bin] ; p2-- ) { //Loop on the bin's particles.
                if (particles.position(0,p2) >= limit ) {
                    //This particle goes up one bin.
                    particles.swap_part(p2,bmax[bin]-1);
                    bmax[bin]--;
                    numSort++;
                }
            }
        }
        //Forward pass + Rebracketting
        for (bin=1; bin<bmin.size(); bin++) { //Loop on the bins.
            limit = min_loc + bin*cell_length[0]*clrw;
            bmin_init = bmin[bin];
            p1 = bmin[bin];
            while (p1 == bmin[bin] && p1 < bmax[bin]) {
                if (particles.position(0,p1) < limit ) {
                    bmin[bin]++;
                    numSort++;
                }
                p1++;
            }
            for( p2 = p1 ; p2 < bmax[bin] ; p2++ ) { //Loop on the bin's particles.
                if (particles.position(0,p2) < limit ) {
                    //This particle goes down one bin.
                    particles.swap_part(p2,bmin[bin]);
                    bmin[bin]++;
                    numSort++;
                }
            }

            //Rebracketting
            //Number of particles from bin going down is: bmin[bin]-bmin_init.
            //Number of particles from bin-1 going up is: bmin_init-bmax[bin-1].
            //Total number of particles we need to swap is the min of both.
            p2 = min(bmin[bin]-bmin_init,bmin_init-bmax[bin-1]);
            if (p2 >0) particles.swap_part(bmax[bin-1],bmin[bin]-p2,p2);
            bmax[bin-1] += bmin[bin] - bmin_init;
            bmin[bin] = bmax[bin-1];
        }
    }
    while ( numSort > 0 );

}
/*
void Species::movingWindow_x(unsigned int shift, SmileiMPI *smpi, PicParams& params)
{
    // Update BC positions
    partBoundCond->moveWindow_x( shift*cell_length[0], smpi );
    // Set for bin managment
    min_loc += shift*cell_length[0];

    // Send particles of first bin on process rank-1
    // If no rank-1 -> particles deleted
    clearExchList();

    for (unsigned int ibin = 0 ; ibin < 1 ; ibin++)
        for (int iPart=bmin[ibin] ; iPart<bmax[ibin]; iPart++ ) {
            addPartInExchList( iPart );
            nrj_mw_lost += particles.weight(iPart)*(particles.lor_fac(iPart)-1.0);
        }

    // bin 0 empty
    // Shifts all the bins by 1.
    bmin.erase( bmin.begin() );
    bmax.erase( bmax.begin() );
    // Create new bin at the end
    // Update last values of bmin and bmax to handle correctly received particles
    bmin.push_back( bmax[bmax.size()-1] );
    bmax.push_back( bmax[bmax.size()-1] );
    bmin[0] = 0;

    int iDim(0);
    smpi->exchangeParticles( this, speciesNumber,params, 0, iDim );

    // Create new particles
    if (smpi->isEastern() ) {
        defineNewCells(shift, smpi, params);
    }

}
*/
void Species::defineNewCells(unsigned int shift, SmileiMPI *smpi, PicParams& params)
{
    // does a loop over all cells in the simulation
    // considering a 3d volume with size n_space[0]*n_space[1]*n_space[2]
    vector<double> cell_index(3,0);
    for (unsigned int i=0 ; i<params.nDim_field ; i++) {
        if (cell_length[i]!=0) {
            cell_index[i] = smpi->getDomainLocalMin(i);
        }
    }
    // cell_index[0] goes to end of the domain minus cell to create
    cell_index[0] = smpi->getDomainLocalMax(0) - shift*cell_length[0];

    // Next bin to create
    int new_bin_idx = bmin.size() - 1;

    vector<unsigned int> n_space_created(3,0);
    n_space_created[0] = shift;
    n_space_created[1] = params.n_space[1];
    n_space_created[2] = params.n_space[2];

    //unsigned int npart_effective =
    createParticles(n_space_created, cell_index, new_bin_idx, params );
}


int Species::createParticles(vector<unsigned int> n_space_to_create, vector<double> cell_index, int new_bin_idx, PicParams& params  )
{
    // ---------------------------------------------------------
    // Calculate density and number of particles for the species
    // ---------------------------------------------------------

    // field containing the charge distribution (always 3d)
    Field3D charge(n_space_to_create);
    max_charge = 0.;

    // field containing the density distribution (always 3d)
    Field3D density(n_space_to_create);

    // field containing the temperature distribution along all 3 momentum coordinates (always 3d * 3)
    Field3D temperature[3];

    // field containing the temperature distribution along all 3 momentum coordinates (always 3d * 3)
    Field3D velocity[3];

    // field containing the number of particles in each cell
    Field3D n_part_in_cell(n_space_to_create);

    for (unsigned int i=0; i<3; i++) {
        velocity[i].allocateDims(n_space_to_create);
        temperature[i].allocateDims(n_space_to_create);
    }

    int npart_effective = 0;
    for (unsigned int i=0; i<n_space_to_create[0]; i++) {
        for (unsigned int j=0; j<n_space_to_create[1]; j++) {
            for (unsigned int k=0; k<n_space_to_create[2]; k++) {

                vector<double> x_cell(3,0);
                x_cell[0] = cell_index[0] + (i+0.5)*cell_length[0];
                x_cell[1] = cell_index[1] + (j+0.5)*cell_length[1];
                x_cell[2] = cell_index[2] + (k+0.5)*cell_length[2];

                n_part_in_cell(i,j,k) = round(ppcProfile->valueAt(x_cell));
                if( n_part_in_cell(i,j,k)<=0. ) {
                    n_part_in_cell(i,j,k) = 0.;
                    density(i,j,k) = 0.;
                    continue;
                }

                // assign charge its correct value in the cell
                charge(i,j,k) = chargeProfile->valueAt(x_cell);
                if( charge(i,j,k)>max_charge ) max_charge=charge(i,j,k);
                // assign density its correct value in the cell
                density(i,j,k) = densityProfile->valueAt(x_cell);
                if(density(i,j,k)!=0. && densityProfileType=="charge") {
                    if(charge(i,j,k)==0.) ERROR("Encountered non-zero charge density and zero charge at the same location");
                    density(i,j,k) /= charge(i,j,k);
                }
                density(i,j,k) = abs(density(i,j,k));
                //MESSAGE(1,"density"<<density(i,j,k));
                // for non-zero density define temperature & mean-velocity and increment the nb of particles
                if (density(i,j,k)!=0.0) {

                    // assign the temperature & mean-velocity their correct value in the cell
                    for (unsigned int m=0; m<3; m++) {
                        temperature[m](i,j,k) = temperatureProfile[m]->valueAt(x_cell);
                        //MESSAGE("temp 1 :" <<  temperature[m](i,j,k))
                        velocity[m](i,j,k) = velocityProfile[m]->valueAt(x_cell);
			//MESSAGE(1,"ddd"<<temperature[m](i,j,k)<<velocity[m](i,j,k));
                    }

                    // increment the effective number of particle by n_part_in_cell(i,j,k)
                    // for each cell with as non-zero density
                    npart_effective += n_part_in_cell(i,j,k);
                    //DEBUG(10,"Specie "<< speciesNumber <<" # part "<<npart_effective<<" "<<i<<" "<<j<<" "<<k);

                }//ENDif non-zero density

            }//i
        }//j
    }//k end the loop on all cells

    // defines npart_effective for the Species & create the corresponding particles
    // -----------------------------------------------------------------------

    // if moving_win
    //     particles.create_particles(npart_effective);
    // else {
    //    // reserve included in initialize if particles emty
    //    particles.reserve(round( params->species_param[speciesNumber].c_part_max * npart_effective ), ndim);
    //    particles.initialize(n_existing_particles+npart_effective, params_->nDim_particle);
    // }
    //MESSAGE(1,"density__end");
    int npart_reserve = params.species_param[speciesNumber].c_part_max
                        * params.species_param[speciesNumber].n_part_per_cell_for_weight
                        * n_space_to_create[0] * n_space_to_create[1] * n_space_to_create[2];
    particles.reserve(npart_reserve, ndim);

    int n_existing_particles = particles.size();
    particles.initialize(n_existing_particles+npart_effective, params, speciesNumber);
    psi_particles.initialize(0, params);


    // define Maxwell-Juettner related quantities
    // ------------------------------------------

    // Maxwell-Juettner cumulative function (array)
    std::vector<double> max_jutt_cumul;

    /*
    if (species_param.initMomentum_type=="maxwell-juettner") {
        //! \todo{Pass this parameters in a code constants class (MG)}
        nE     = 20000;
        muEmax = 20.0;

        max_jutt_cumul.resize(nE);
        double mu=species_param.mass/species_param.temperature[0];
        //double mu=species_param.mass/temperature[m](i,j,k);
        double Emax=muEmax/mu;
        dE=Emax/nE;

        double fl=0;
        double fr=0;
        max_jutt_cumul[0]=0.0;
        for (unsigned  i=1; i<nE; i++ ) {
            //! \todo{this is just the isotropic case, generalise to non-isotropic (MG)}
            fr=(1+i*dE)*sqrt(pow(1.0+i*dE,2)-1.0) * exp(-mu*i*dE);
            max_jutt_cumul[i]=max_jutt_cumul[i-1] + 0.5*dE*(fr+fl);
            fl=fr;
        }
        for (unsigned int i=0; i<nE; i++) max_jutt_cumul[i]/=max_jutt_cumul[nE-1];

    }
    */

    // Initialization of the particles properties
    // ------------------------------------------
    unsigned int nPart;
    unsigned int iPart=n_existing_particles;
    double *indexes=new double[params.nDim_particle];
    double *temp=new double[3];
    double *vel=new double[3];

    // start a loop on all cells

    //bmin[bin] point to begining of bin (first particle)
    //bmax[bin] point to end of bin (= bmin[bin+1])
    //if bmax = bmin, bin is empty of particle.
    for (unsigned int i=0; i<n_space_to_create[0]; i++) {
        if (i%clrw == 0) bmin[new_bin_idx+i/clrw] = iPart;
        for (unsigned int j=0; j<n_space_to_create[1]; j++) {
            for (unsigned int k=0; k<n_space_to_create[2]; k++) {
                // initialize particles in meshes where the density is non-zero
                if (density(i,j,k)>0) {

                    if (species_param.initMomentum_type=="maxwell-juettner") {
                        //! \todo{Pass this parameters in a code constants class (MG)}
                        nE     = 20000;
                        muEmax = 20.0;

                        max_jutt_cumul.resize(nE);
                        //double mu=species_param.mass/species_param.temperature[0];
                        double mu=species_param.mass/temperature[0](i,j,k); // For Temperature profile
                        double Emax=muEmax/mu;
                        dE=Emax/nE;

                        double fl=0;
                        double fr=0;
                        max_jutt_cumul[0]=0.0;
                        for (unsigned int l=1; l<nE; l++ ) {
                            //! \todo{this is just the isotropic case, generalise to non-isotropic (MG)}
                            fr=(1.+l*dE)*sqrt(pow(1.0+l*dE,2)-1.0) * exp(-mu*l*dE);
                            max_jutt_cumul[l]=max_jutt_cumul[l-1] + 0.5*dE*(fr+fl);
                            fl=fr;
                        }
                        for (unsigned int l=0; l<nE; l++) max_jutt_cumul[l]/=max_jutt_cumul[nE-1];
                    }

                    temp[0] = temperature[0](i,j,k);
                    vel[0]  = velocity[0](i,j,k);
                    temp[1] = temperature[1](i,j,k);
                    vel[1]  = velocity[1](i,j,k);
                    temp[2] = temperature[2](i,j,k);
                    vel[2]  = velocity[2](i,j,k);
                    nPart = n_part_in_cell(i,j,k);

                    indexes[0]=i*cell_length[0]+cell_index[0];
                    if (ndim > 1) {
                        indexes[1]=j*cell_length[1]+cell_index[1];
                        if (ndim > 2) {
                            indexes[2]=k*cell_length[2]+cell_index[2];
                        }//ndim > 2
                    }//ndim > 1

                    initPosition(nPart, iPart, indexes, params.nDim_particle,
                                 cell_length, species_param.initPosition_type);

                    initMomentum(nPart,iPart, temp, vel,
                                 species_param.initMomentum_type, max_jutt_cumul, params);

                    initAcceleration_imp(nPart, iPart);

                    initWeight(nPart, speciesNumber, iPart, density(i,j,k));
                    initCharge(nPart, speciesNumber, iPart, charge(i,j,k));

                    //calculate new iPart (jump to next cell)
                    iPart+=nPart;
                }//END if density > 0
            }//k end the loop on all cells
        }//j
        if (i%clrw == clrw -1) bmax[new_bin_idx+i/clrw] = iPart;
    }//i

    delete [] indexes;
    delete [] temp;
    delete [] vel;

    // Recalculate former position using the particle velocity
    // (necessary to calculate currents at time t=0 using the Esirkepov projection scheme)
    for (int iPart=n_existing_particles; iPart<n_existing_particles+npart_effective; iPart++) {
        /*897 for (int i=0; i<(int)ndim; i++) {
            particles.position_old(i,iPart) -= particles.momentum(i,iPart)/particles.lor_fac(iPart) * params.timestep;
        }897*/
        nrj_new_particles += particles.weight(iPart)*(particles.lor_fac(iPart)-1.0);
    }

    if (particles.isTestParticles)
	particles.setIds();

    return npart_effective;

}


void Species::insert_particles_to_bins(Particles &insert_Particles, std::vector<int> &count_in_bins)
{
    int n_part_insert;
    int begin_id = 0;
    for(int ibin=0; ibin<bmax.size(); ibin++)
    {
        n_part_insert = count_in_bins[ibin];
        insert_Particles.cp_particles(begin_id, n_part_insert, particles, bmax[ibin]);
        bmax[ibin] += n_part_insert;
        for (int ibin_inLoop=ibin+1 ; ibin_inLoop < bmax.size() ; ibin_inLoop++ ) {
            bmax[ibin_inLoop] += n_part_insert ;
            bmin[ibin_inLoop] = bmax[ibin_inLoop-1] ;
        }
        begin_id += n_part_insert;
    }
}
void Species::erase_particles_from_bins(std::vector<int> &indexs_to_erase)
{
    int ii, iPart;
    for (unsigned int ibin = 0 ; ibin < bmax.size() ; ibin++ ) {
        //        DEBUG(ibin << " bounds " << bmin[ibin] << " " << bmax[ibin]);
        ii = indexs_to_erase.size()-1;
        if (ii >= 0) { // Push lost particles to the end of the bin
            iPart = indexs_to_erase[ii];
            while (iPart >= bmax[ibin] && ii > 0) {
                ii--;
                iPart = indexs_to_erase[ii];
            }
            while (iPart == bmax[ibin]-1 && iPart >= bmin[ibin] && ii > 0) {
                bmax[ibin]--;
                ii--;
                iPart = indexs_to_erase[ii];
            }
            while (iPart >= bmin[ibin] && ii > 0) {
                particles.overwrite_part(bmax[ibin]-1, iPart );
                bmax[ibin]--;
                ii--;
                iPart = indexs_to_erase[ii];
            }
            if (iPart >= bmin[ibin] && iPart < bmax[ibin]) { //On traite la dernière particule (qui peut aussi etre la premiere)
                particles.overwrite_part(bmax[ibin]-1, iPart );
                bmax[ibin]--;
            }
        }
    }
    //Shift the bins in memory
    //Warning: this loop must be executed sequentially. Do not use openMP here.
    for (int unsigned ibin = 1 ; ibin < bmax.size() ; ibin++ ) { //First bin don't need to be shifted
        ii = bmin[ibin]-bmax[ibin-1]; // Shift the bin in memory by ii slots.
        iPart = min(ii,bmax[ibin]-bmin[ibin]); // Number of particles we have to shift = min (Nshift, Nparticle in the bin)
        if(iPart > 0) particles.overwrite_part(bmax[ibin]-iPart,bmax[ibin-1],iPart);
        bmax[ibin] -= ii;
        bmin[ibin] = bmax[ibin-1];
    }


    // Delete useless Particles
    //Theoretically, not even necessary to do anything as long you use bmax as the end of your iterator on particles.
    //Nevertheless, you might want to free memory and have the actual number of particles
    //really equal to the size of the vector. So we do:
    particles.erase_particle_trail(bmax.back());
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
}








/*
void Species::updateMvWinLimits(double x_moved) {
    partBoundCond->updateMvWinLimits(x_moved);
    min_loc += x_moved;
}
//Do we have to project this species ?
 bool Species::isProj(double time_dual, SimWindow* simWindow) {
    //Recompute frozen particles density if
    //moving window is activated, actually moving at this time step, and we are not in a density slope.
    return time_dual > species_param.time_frozen  || (simWindow && simWindow->isMoving(time_dual)) ;
}
*/

void Species::printAvgVelocity()
{
    double vx_avg, vy_avg, vz_avg;
    vx_avg = 0.0;
    vy_avg = 0.0;
    vz_avg = 0.0;
    for (unsigned int p= 0; p<particles.size(); p++)
    {
        vx_avg += abs( particles.momentum(0,p) );
        vx_avg += abs( particles.momentum(0,p) );
        vx_avg += abs( particles.momentum(0,p) );
    }

    vx_avg /= particles.size();
    vy_avg /= particles.size();
    vz_avg /= particles.size();

    MESSAGE("Average velocity of "<<species_param.species_type << "  " << vx_avg << vy_avg << vz_avg);

}


void Species::calDepCharge(ElectroMagn* EMfields, PicParams& params, SmileiMPI* smpi)
{
    unsigned int iPart;
    double depCharge[3][2];
    double depCharge_temp[3][2];
    for ( int iDim = 0 ; iDim<(int)params.nDim_particle ; iDim++ )
    {
        for (int i=0 ; i<indexes_of_particles_to_exchange.size() ; i++) {
            iPart = indexes_of_particles_to_exchange[i];
            if ( particles.position(iDim,iPart) < smpi->getDomainLocalMin(iDim) ) {
                depCharge_temp[iDim][0] += ( particles.charge(iPart) + particles.weight(iPart) );
            }
            else if ( particles.position(iDim,iPart) > smpi->getDomainLocalMax(iDim) ) {
                depCharge_temp[iDim][1] += ( particles.charge(iPart) + particles.weight(iPart) );
            }
            else { // particle will be deleted, if supp_particle particles still in the domain
            }
        }

    }

    if( !smpi->isWestern() ) {
        depCharge_temp[0][0] = 0.0;
    }

    if( !smpi->isEastern() ) {
        depCharge_temp[0][1] = 0.0;
    }

    if( params.nDim_particle >= 2 ) {
        if( !smpi->isSouthern() ) {
            depCharge_temp[1][0] = 0.0;
        }

        if( !smpi->isNorthern() ) {
            depCharge_temp[1][1] = 0.0;
        }
    }

    smpi->reduceDoubleVector( depCharge_temp[0], depCharge[0], 6 );
    for ( int iDim = 0 ; iDim<(int)params.nDim_particle ; iDim++ )
    {
        EMfields->depCharge[iDim][0] += depCharge[iDim][0];
        EMfields->depCharge[iDim][1] += depCharge[iDim][1];
    }

}
