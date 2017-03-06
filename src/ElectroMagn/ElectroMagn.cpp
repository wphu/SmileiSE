#include "ElectroMagn.h"

#include <limits>
#include <iostream>

#include "PicParams.h"
#include "Species.h"
#include "Projector.h"
#include "Field.h"
#include "Profile.h"
#include "SolverFactory.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for the virtual class ElectroMagn
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn::ElectroMagn(PicParams &params, InputData &input_data, SmileiMPI* smpi) :
timestep(params.timestep),
cell_length(params.cell_length),
n_species(params.species_param.size()),
nDim_field(params.nDim_field),
cell_volume(params.cell_volume),
n_space(params.n_space),
n_space_global(params.n_space_global),
oversize(params.oversize)
{
    species_param = params.species_param;

    // initialize poynting vector
    poynting[0].resize(nDim_field,0.0);
    poynting[1].resize(nDim_field,0.0);
    poynting_inst[0].resize(nDim_field,0.0);
    poynting_inst[1].resize(nDim_field,0.0);

    // initialize charge vector: the charge is scalar, not field
    emitCharge[0].resize(nDim_field, 0.0);
    emitCharge[1].resize(nDim_field, 0.0);
    depCharge[0].resize(nDim_field, 0.0);
    depCharge[1].resize(nDim_field, 0.0);
    totCharge[0].resize(nDim_field, 0.0);
    totCharge[1].resize(nDim_field, 0.0);



    // take useful things from params
    for (unsigned int i=0; i<3; i++) {
        DEBUG("____________________ OVERSIZE: " <<i << " " << oversize[i]);
    }

    if (n_space.size() != 3) ERROR("this should not happend");

    Ex_=NULL;
    Ey_=NULL;
    Ez_=NULL;
    Bx_=NULL;
    By_=NULL;
    Bz_=NULL;
    Bx_m=NULL;
    By_m=NULL;
    Bz_m=NULL;
    Jx_=NULL;
    Jy_=NULL;
    Jz_=NULL;
    rho_=NULL;

    Ex_avg=NULL;
    Ey_avg=NULL;
    Ez_avg=NULL;
    Bx_avg=NULL;
    By_avg=NULL;
    Bz_avg=NULL;

    // Species charge currents and density
    Jx_s.resize(n_species);
    Jy_s.resize(n_species);
    Jz_s.resize(n_species);
    rho_s.resize(n_species);
    rho_s_avg.resize(n_species);
    rho_s_global.resize(n_species);
    rho_s_global_avg.resize(n_species);

    Vx_s.resize(n_species);
    Vx_s_avg.resize(n_species);
    Vx_s_global.resize(n_species);
    Vx_s_global_avg.resize(n_species);

    Vy_s.resize(n_species);
    Vy_s_avg.resize(n_species);
    Vy_s_global.resize(n_species);
    Vy_s_global_avg.resize(n_species);

    Vz_s.resize(n_species);
    Vz_s_avg.resize(n_species);
    Vz_s_global.resize(n_species);
    Vz_s_global_avg.resize(n_species);

    T_s.resize(n_species);
    T_s_avg.resize(n_species);
    T_s_global.resize(n_species);
    T_s_global_avg.resize(n_species);
    for (unsigned int ispec=0; ispec<n_species; ispec++) {
        Jx_s[ispec]             = NULL;
        Jy_s[ispec]             = NULL;
        Jz_s[ispec]             = NULL;
        rho_s[ispec]            = NULL;
        rho_s_avg[ispec]        = NULL;
        rho_s_global[ispec]     = NULL;
        rho_s_global_avg[ispec] = NULL;

        Vx_s[ispec]             = NULL;
        Vx_s_avg[ispec]         = NULL;
        Vx_s_global[ispec]      = NULL;
        Vx_s_global_avg[ispec]  = NULL;

        Vy_s[ispec]             = NULL;
        Vy_s_avg[ispec]         = NULL;
        Vy_s_global[ispec]      = NULL;
        Vy_s_global_avg[ispec]  = NULL;

        Vz_s[ispec]             = NULL;
        Vz_s_avg[ispec]         = NULL;
        Vz_s_global[ispec]      = NULL;
        Vz_s_global_avg[ispec]  = NULL;

        Vp_s[ispec]             = NULL;
        Vp_s_avg[ispec]         = NULL;
        Vp_s_global[ispec]      = NULL;
        Vp_s_global_avg[ispec]  = NULL;

        T_s[ispec]              = NULL;
        T_s_avg[ispec]          = NULL;
        T_s_global[ispec]       = NULL;
        T_s_global_avg[ispec]   = NULL;
    }

    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<2; j++) {
            istart[i][j]=0;
            bufsize[i][j]=0;
        }
    }

    //emBoundCond = ElectroMagnBC_Factory::create(params, laser_params);

    //MaxwellFaradaySolver_ = SolverFactory::create(params);

}



// ---------------------------------------------------------------------------------------------------------------------
// Destructor for the virtual class ElectroMagn
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn::~ElectroMagn()
{
    delete Ex_;
    delete Ey_;
    delete Ez_;
    delete Bx_;
    delete By_;
    delete Bz_;
    delete Bx_m;
    delete By_m;
    delete Bz_m;
    delete Jx_;
    delete Jy_;
    delete Jz_;
    delete rho_;

    if (Ex_avg!=NULL) {
        delete Ex_avg;
        delete Ey_avg;
        delete Ez_avg;
        delete Bx_avg;
        delete By_avg;
        delete Bz_avg;
    }

    for (unsigned int ispec=0; ispec<n_species; ispec++) {
      delete Jx_s[ispec];
      delete Jy_s[ispec];
      delete Jz_s[ispec];
      delete rho_s[ispec];
    }


}//END Destructer


// ---------------------------------------------------------------------------------------------------------------------
// Method used to create a dump of the data contained in ElectroMagn
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn::dump()
{
    //!\todo Check for none-cartesian grid & for generic grid (neither all dual or all primal) (MG & JD)

    vector<unsigned int> dimPrim;
    dimPrim.resize(1);
    dimPrim[0] = n_space[0]+2*oversize[0]+1;
    vector<unsigned int> dimDual;
    dimDual.resize(1);
    dimDual[0] = n_space[0]+2*oversize[0]+2;

    // dump of the electromagnetic fields
    Ex_->dump(dimDual);
    Ey_->dump(dimPrim);
    Ez_->dump(dimPrim);
    Bx_->dump(dimPrim);
    By_->dump(dimDual);
    Bz_->dump(dimDual);
    // dump of the total charge density & currents
    rho_->dump(dimPrim);
    Jx_->dump(dimDual);
    Jy_->dump(dimPrim);
    Jz_->dump(dimPrim);
}



// ---------------------------------------------------------------------------------------------------------------------
// Method used to initialize the total charge density
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn::initRhoJ(vector<Species*>& vecSpecies, Projector* Proj)
{
    //! \todo Check that one uses only none-test particles
    // number of (none-test) used in the simulation
    //! \todo fix this: n_species is already a member of electromagn, is it this confusing? what happens if n_species grows (i.e. with ionization)?
    unsigned int n_species = vecSpecies.size();

    //loop on all (none-test) Species
    for (unsigned int iSpec=0 ; iSpec<n_species; iSpec++ ) {
        Particles &cuParticles = vecSpecies[iSpec]->getParticlesList();
        unsigned int n_particles = vecSpecies[iSpec]->getNbrOfParticles();

        DEBUG(n_particles<<" species "<<iSpec);
	if (!cuParticles.isTestParticles) {
	    for (unsigned int iPart=0 ; iPart<n_particles; iPart++ ) {
		// project charge & current densities
		(*Proj)(Jx_s[iSpec], Jy_s[iSpec], Jz_s[iSpec], rho_s[iSpec], cuParticles, iPart,
			cuParticles.lor_fac(iPart));
	    }
	}

    }//iSpec
    DEBUG("before computeTotalRhoJ");
    computeTotalRhoJ();
    DEBUG("projection done for initRhoJ");

}






double ElectroMagn::computeNRJ(unsigned int shift, SmileiMPI *smpi) {
    double nrj(0.);

    if ( smpi->isWestern() ) {
	nrj += Ex_->computeNRJ(shift, istart, bufsize);
	nrj += Ey_->computeNRJ(shift, istart, bufsize);
	nrj += Ez_->computeNRJ(shift, istart, bufsize);

	nrj += Bx_m->computeNRJ(shift, istart, bufsize);
	nrj += By_m->computeNRJ(shift, istart, bufsize);
	nrj += Bz_m->computeNRJ(shift, istart, bufsize);

    }

    return nrj;
}

bool ElectroMagn::isRhoNull(SmileiMPI* smpi)
{
    double norm2(0.);
    double locnorm2(0.);

    // rho_->isDual(i) = 0 for all i
    // istart[i][0] & bufsize[i][0]

    vector<unsigned int> iFieldStart(3,0), iFieldEnd(3,1), iFieldGlobalSize(3,1);
    for (unsigned int i=0 ; i<rho_->isDual_.size() ; i++ ) {
	iFieldStart[i] = istart[i][0];
	iFieldEnd [i] = iFieldStart[i] + bufsize[i][0];
	iFieldGlobalSize [i] = rho_->dims_[i];
    }

    for (unsigned int k=iFieldStart[2]; k<iFieldEnd[2]; k++) {
	for (unsigned int j=iFieldStart[1]; j<iFieldEnd[1]; j++) {
	    for (unsigned int i=iFieldStart[0]; i<iFieldEnd[0]; i++) {
		unsigned int ii=k+ j*iFieldGlobalSize[2] +i*iFieldGlobalSize[1]*iFieldGlobalSize[2];
		locnorm2 += (*rho_)(ii)*(*rho_)(ii);
	    }
	}
    }

    MPI_Allreduce(&locnorm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return (norm2<=0.);

}

string LowerCase(string in){
    string out=in;
    //std::transform(out.begin(), out.end(), out.begin(), ::tolower);
    return out;
}
