#include "PicParams.h"
#include <cmath>
#include "Tools.h"
#include "PyTools.h"
#include "InputData.h"

#include <algorithm>

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// PicParams : open & parse the input data file, test that parameters are coherent
// ---------------------------------------------------------------------------------------------------------------------
PicParams::PicParams(InputData &ifile) {


    // --------------
    // Stop & Restart
    // --------------
    dump_step=1;
    ifile.extract("dump_step", dump_step);

    ntime_step_avg = 0;
    ifile.extract("ntime_step_avg", ntime_step_avg);

    dump_minutes=0.0;
    ifile.extract("dump_minutes", dump_minutes);

    timesteps_collision = 1;
    ifile.extract("timesteps_collision", timesteps_collision);

    timesteps_coulomb = 1;
    ifile.extract("timesteps_coulomb", timesteps_coulomb);

    timesteps_DSMC = 1000000000000;
    ifile.extract("timesteps_DSMC", timesteps_DSMC);

    timesteps_restore = 1000000000;
    ifile.extract("timesteps_restore", timesteps_restore);

    exit_after_dump=true;
    ifile.extract("exit_after_dump", exit_after_dump);

    restart=false;
    ifile.extract("restart", restart);
    if (restart) MESSAGE("Code running from restart"); //! \todo Give info on restart properties

    is_continue = 0;
    ifile.extract("is_continue", is_continue);

    is_calVDF = 0;
    ifile.extract("is_calVDF", is_calVDF);


    //!\todo MG is this still used ?? I cannot find it anywhere
    check_stop_file=false;
    ifile.extract("check_stop_file", check_stop_file);

    dump_file_sequence=2;
    ifile.extract("dump_file_sequence", dump_file_sequence);
    dump_file_sequence=std::max((unsigned int)1,dump_file_sequence);


    // ---------------------
    // Normalisation & units
    // ---------------------

    wavelength_SI = 0.;
    ifile.extract("wavelength_SI",wavelength_SI);


    // method of PIC
    ifile.extract("method", method);
    if (method!="explicit" && method!="implicit") {
        ERROR("Method " << method << " does not exist");
    }

    imp_theta = 0.1;
    ifile.extract("imp_theta", imp_theta);

    imp_iteration_number = 1;
    ifile.extract("imp_iteration_number", imp_iteration_number);

    // -------------------
    // Simulation box info
    // -------------------

    // geometry of the simulation
    ifile.extract("dim", geometry);
    if (geometry!="1d3v" && geometry!="2d3v") {
        ERROR("Geometry " << geometry << " does not exist");
    }
    setDimensions();

    // interpolation order
    ifile.extract("interpolation_order", interpolation_order);

    // projection order
    ifile.extract("projection_order", projection_order);

    /*
    if (interpolation_order!=2 && interpolation_order!=4) {
        ERROR("Interpolation/projection order " << interpolation_order << " not defined");
    }
    if (geometry=="2d3v" && interpolation_order==4) {
        ERROR("Interpolation/projection order " << interpolation_order << " not yet defined in 2D");
    }
    */

    //!\todo (MG to JD) Please check if this parameter should still appear here
    // Disabled, not compatible for now with particles sort
    // if ( !ifile.extract("exchange_particles_each", exchange_particles_each) )
    exchange_particles_each = 1;


    // TIME & SPACE RESOLUTION/TIME-STEPS

    // reads timestep & cell_length
    ifile.extract("timestep", timestep);
    ifile.extract("n_time", n_time);
    ifile.extract("cell_length",cell_length);
    if (cell_length.size()!=nDim_field) {
        ERROR("Dimension of cell_length ("<< cell_length.size() << ") != " << nDim_field << " for geometry " << geometry);
    }
    res_space.resize(nDim_field);
    for (unsigned int i=0;i<nDim_field;i++){
        res_space[i] = 1.0/cell_length[i];
    }

    time_fields_frozen=0.0;
    ifile.extract("time_fields_frozen", time_fields_frozen);

    // testing the CFL condition
    //!\todo (MG) CFL cond. depends on the Maxwell solv. ==> Move this computation to the ElectroMagn Solver
    double res_space2=0;
    for (unsigned int i=0; i<nDim_field; i++) {
        res_space2 += res_space[i]*res_space[i];
    }
    dtCFL=1.0/sqrt(res_space2);
    if ( timestep>dtCFL ) {
        //ERROR("CFL problem: timestep=" << timestep << " should be smaller than " << dtCFL);
    }


    // simulation duration & length
    ifile.extract("sim_time", sim_time);

    ifile.extract("sim_length",sim_length);
    if (sim_length.size()!=nDim_field) {
        ERROR("Dimension of sim_length ("<< sim_length.size() << ") != " << nDim_field << " for geometry " << geometry);
    }

    bcType = "constant";
    ifile.extract("bcType", bcType);

    //! Boundary conditions for ElectroMagnetic Fields
    if ( !ifile.extract("bc_em_type_x", bc_em_type_x)  ) {
        ERROR("Electromagnetic boundary condition type (bc_em_type_x) not defined" );
    }
    if (bc_em_type_x.size()==1) { // if just one type is specified, then take the same bc type in a given dimension
        bc_em_type_x.resize(2); bc_em_type_x[1]=bc_em_type_x[0];
    }

    if ( !ifile.extract("bc_em_value_x", bc_em_value_x)  ) {
        ERROR("Electromagnetic boundary condition type (bc_em_value_x) not defined" );
    }
    if (bc_em_value_x.size()==1) { // if just one type is specified, then take the same bc type in a given dimension
        bc_em_value_x.resize(2); bc_em_value_x[1]=bc_em_value_x[0];
    }


    if ( geometry == "2d3v" || geometry == "3d3v" ) {
        if ( !ifile.extract("bc_em_type_y", bc_em_type_y) )
            ERROR("Electromagnetic boundary condition type (bc_em_type_y) not defined" );
        if (bc_em_type_y.size()==1) { // if just one type is specified, then take the same bc type in a given dimension
            bc_em_type_y.resize(2); bc_em_type_y[1]=bc_em_type_y[0];
        }
    }
    if ( geometry == "3d3v" ) {
        if ( !ifile.extract("bc_em_type_z", bc_em_type_z) )
            ERROR("Electromagnetic boundary condition type (bc_em_type_z) not defined" );
        if (bc_em_type_z.size()==1) { // if just one type is specified, then take the same bc type in a given dimension
            bc_em_type_z.resize(2); bc_em_type_z[1]=bc_em_type_z[0];
        }
    }

    // ------------------------
    // Moving window parameters
    // ------------------------
    if (!ifile.extract("clrw",clrw)) {
        clrw = 1;
    }

    if ( !ifile.extract("externB", externB)  ) {
        ERROR("Extern magnetic field is not defined" );
    }



    // ------------------
    // Species properties
    // ------------------
    readSpecies(ifile);

    global_every=0;

    ifile.extract("every",global_every);

    // --------------------
    // Number of processors
    // --------------------
    if ( !ifile.extract("number_of_procs", number_of_procs) )
        number_of_procs.resize(nDim_field, 0);

    // -------------------------------------------------------
    // Compute usefull quantities and introduce normalizations
    // also defines defaults values for the species lengths
    // -------------------------------------------------------
    computeNormalization();
    compute();
    computeSpecies();

}

void PicParams::readSpecies(InputData &ifile) {
    bool ok;
    for (int ispec = 0; ispec < ifile.nComponents("Species"); ispec++) {
        SpeciesStructure tmpSpec;

        ifile.extract("species_type",tmpSpec.species_type,"Species",ispec);
        if(tmpSpec.species_type.empty()) {
            ERROR("For species #" << ispec << " empty species_type");
        }
        ifile.extract("initPosition_type",tmpSpec.initPosition_type ,"Species",ispec);
        if (tmpSpec.initPosition_type.empty()) {
            ERROR("For species #" << ispec << " empty initPosition_type");
        } else if ( (tmpSpec.initPosition_type!="regular")&&(tmpSpec.initPosition_type!="random") ) {
            ERROR("For species #" << ispec << " bad definition of initPosition_type " << tmpSpec.initPosition_type);
        }

        ifile.extract("initMomentum_type",tmpSpec.initMomentum_type ,"Species",ispec);
        if ( (tmpSpec.initMomentum_type=="mj") || (tmpSpec.initMomentum_type=="maxj") ) {
            tmpSpec.initMomentum_type="maxwell-juettner";
        }
        if (   (tmpSpec.initMomentum_type!="cold")
            && (tmpSpec.initMomentum_type!="maxwell")
            && (tmpSpec.initMomentum_type!="maxwell-juettner")
            && (tmpSpec.initMomentum_type!="rectangular") ) {
            ERROR("For species #" << ispec << " bad definition of initMomentum_type");
        }

        tmpSpec.c_part_max = 1.0;// default value
        ifile.extract("c_part_max",tmpSpec.c_part_max,"Species",ispec);

        if( !ifile.extract("mass",tmpSpec.mass ,"Species",ispec) ) {
            ERROR("For species #" << ispec << ", mass not defined.");
        }

        tmpSpec.dynamics_type = "norm"; // default value
        if (!ifile.extract("dynamics_type",tmpSpec.dynamics_type ,"Species",ispec) )
            WARNING("For species #" << ispec << ", dynamics_type not defined: assumed = 'norm'.");
        if (tmpSpec.dynamics_type!="norm"){
            ERROR("dynamics_type different than norm not yet implemented");
        }

        tmpSpec.time_frozen = 0.0; // default value
        ifile.extract("time_frozen",tmpSpec.time_frozen ,"Species",ispec);
        if (tmpSpec.time_frozen > 0 && \
            tmpSpec.initMomentum_type!="cold") {
            WARNING("For species #" << ispec << " possible conflict between time-frozen & not cold initialization");
        }

        tmpSpec.radiating = false; // default value
        ifile.extract("radiating",tmpSpec.radiating ,"Species",ispec);
        if (tmpSpec.dynamics_type=="rrll" && (!tmpSpec.radiating)) {
            WARNING("For species #" << ispec << ", dynamics_type='rrll' forcing radiating=True");
            tmpSpec.radiating=true;
        }

        if (!ifile.extract("bc_part_type_west",tmpSpec.bc_part_type_west,"Species",ispec) )
            ERROR("For species #" << ispec << ", bc_part_type_west not defined");
        if (!ifile.extract("bc_part_type_east",tmpSpec.bc_part_type_east,"Species",ispec) )
            ERROR("For species #" << ispec << ", bc_part_type_east not defined");

        if (nDim_particle>1) {
            if (!ifile.extract("bc_part_type_south",tmpSpec.bc_part_type_south,"Species",ispec) )
                ERROR("For species #" << ispec << ", bc_part_type_south not defined");
            if (!ifile.extract("bc_part_type_north",tmpSpec.bc_part_type_north,"Species",ispec) )
                ERROR("For species #" << ispec << ", bc_part_type_north not defined");
        }

        // for thermalizing BCs on particles check if thermT is correctly defined
        bool thermTisDefined=false;
        if ( (tmpSpec.bc_part_type_west=="thermalize") || (tmpSpec.bc_part_type_east=="thermalize") ){
            thermTisDefined=ifile.extract("thermT",tmpSpec.thermT,"Species",ispec);
            if (!thermTisDefined) ERROR("thermT needs to be defined for species " <<ispec<< " due to x-BC thermalize");
        }
        if ( (nDim_particle==2) && (!thermTisDefined) &&
             (tmpSpec.bc_part_type_south=="thermalize" || tmpSpec.bc_part_type_north=="thermalize") ) {
            thermTisDefined=ifile.extract("thermT",tmpSpec.thermT,"Species",ispec);
            if (!thermTisDefined) ERROR("thermT needs to be defined for species " <<ispec<< " due to y-BC thermalize");
        }
        if (thermTisDefined) {
            if (tmpSpec.thermT.size()==1) {
                tmpSpec.thermT.resize(3);
                for (unsigned int i=1; i<3;i++)
                    tmpSpec.thermT[i]=tmpSpec.thermT[0];
            }
        } else {
            tmpSpec.thermT.resize(3);
            for (unsigned int i=0; i<3;i++)
                tmpSpec.thermT[i]=0.0;
        }


        tmpSpec.ionization_model = "none"; // default value
        ifile.extract("ionization_model", tmpSpec.ionization_model, "Species",ispec);

        ifile.extract("atomic_number", tmpSpec.atomic_number, "Species",ispec);
        ifile.extract("atomic_mass", tmpSpec.atomic_mass, "Species",ispec);
        ifile.extract("surface_binding_energy", tmpSpec.surface_binding_energy, "Species",ispec);
        ifile.extract("density_solid", tmpSpec.density_solid, "Species",ispec);
        ifile.extract("ne", tmpSpec.ne, "Species",ispec);
        ifile.extract("nz2", tmpSpec.nz2, "Species",ispec);
        ifile.extract("nw", tmpSpec.nw, "Species",ispec);


        tmpSpec.isTest = false; // default value
        ifile.extract("isTest",tmpSpec.isTest ,"Species",ispec);
        if (tmpSpec.ionization_model!="none" && (!tmpSpec.isTest)) {
            ERROR("Disabled for now : test & ionized");
        }

        // Species geometry
        // ----------------

        // Density
        bool ok1, ok2;
        ok1 = extractOneProfile(ifile, "nb_density"    , tmpSpec.dens_profile, ispec);
        ok2 = extractOneProfile(ifile, "charge_density", tmpSpec.dens_profile, ispec);
        if(  ok1 &&  ok2 ) ERROR("For species #" << ispec << ", cannot define both `nb_density` and `charge_density`.");
        if( !ok1 && !ok2 ) ERROR("For species #" << ispec << ", must define `nb_density` or `charge_density`.");
        if( ok1 ) tmpSpec.density_type = "nb";
        if( ok2 ) tmpSpec.density_type = "charge";

        // number density which is a constant
        if( !ifile.extract("nb_density",tmpSpec.density ,"Species",ispec) ) {
            ERROR("For species #" << ispec << ", density not defined.");
        }

        // Number of particles per cell
        if( !extractOneProfile(ifile, "n_part_per_cell", tmpSpec.ppc_profile, ispec) )
            ERROR("For species #" << ispec << ", n_part_per_cell not found or not understood");

        if( !ifile.extract("n_part_per_cell",tmpSpec.n_part_per_cell ,"Species",ispec) ) {
            ERROR("For species #" << ispec << ", n_part_per_cell not defined.");
        }

        // Number of particles per cell for weight
        if( !ifile.extract("n_part_per_cell_for_weight",tmpSpec.n_part_per_cell_for_weight ,"Species",ispec) ) {
            ERROR("For species #" << ispec << ", n_part_per_cell_for_weight not defined.");
        }


        // Charge
        if( !extractOneProfile(ifile, "charge", tmpSpec.charge_profile, ispec) )
            ERROR("For species #" << ispec << ", charge not found or not understood");
        // charge which is a constant
        if( !ifile.extract("charge",tmpSpec.charge ,"Species",ispec) ) {
            ERROR("For species #" << ispec << ", charge not defined.");
        }



        // Mean velocity
        vector<ProfileStructure*> vecMvel;
        extractVectorOfProfiles(ifile, "mean_velocity", vecMvel, ispec);
        tmpSpec.mvel_x_profile = *(vecMvel[0]);
        tmpSpec.mvel_y_profile = *(vecMvel[1]);
        tmpSpec.mvel_z_profile = *(vecMvel[2]);

        // Temperature
        vector<ProfileStructure*> vecTemp;
        extractVectorOfProfiles(ifile, "temperature", vecTemp, ispec);
        tmpSpec.temp_x_profile = *(vecTemp[0]);
        tmpSpec.temp_y_profile = *(vecTemp[1]);
        tmpSpec.temp_z_profile = *(vecTemp[2]);

        //ifile.extract("temperature",tmpSpec.thermT);
        ifile.extract("temperature", tmpSpec.thermT, "Species",ispec);


        // DSMC parameters
        ifile.extract("diameter",tmpSpec.diameter ,"Species",ispec);

        ifile.extract("ref_temperature",tmpSpec.ref_temperature ,"Species",ispec);
        ifile.extract("visc_temp_index",tmpSpec.visc_temp_index ,"Species",ispec);
        ifile.extract("vss_scat_inv",tmpSpec.vss_scat_inv ,"Species",ispec);

        // Save the Species params
        // -----------------------
        species_param.push_back(tmpSpec);
    }
}

bool PicParams::extractProfile(InputData &ifile, PyObject *mypy, ProfileStructure &P)
{
    double val;
    // If the profile is only a double, then convert to a constant function
    if( PyTools::convert(mypy, val) ) {
        // Extract the function "constant"
        PyObject* constantFunction = ifile.extract_py("constant");
        // Create the argument which has the value of the profile
        PyObject* arg = PyTuple_New(1);
        PyTuple_SET_ITEM(arg, 0, PyFloat_FromDouble(val));
        // Create the constant anonymous function
        PyObject * tmp = PyObject_Call(constantFunction, arg, NULL);
        P.py_profile = tmp;
        return true;
    } else if (mypy && PyCallable_Check(mypy)) {
        P.py_profile=mypy;
        return true;
    }
    return false;
}

bool PicParams::extractOneProfile(InputData &ifile, string varname, ProfileStructure &P, int ispec) {
    PyObject *mypy = ifile.extract_py(varname, "Species", ispec);
    if( !extractProfile(ifile, mypy, P) ) return false;
    return true;
}

void PicParams::extractVectorOfProfiles(InputData &ifile, string varname, vector<ProfileStructure*> &Pvec, int ispec)
{
    Pvec.resize(3);
    vector<PyObject*> pvec = ifile.extract_pyVec(varname, "Species", ispec);
    int len = pvec.size();
    if( len==3 ) {
        for(int i=0; i<len; i++) {
            Pvec[i] = new ProfileStructure();
            if( !extractProfile(ifile, pvec[i], *(Pvec[i])) )
                ERROR("For species #" << ispec << ", "<<varname<<"["<<i<<"] not understood");
        }
    } else if ( len==1 ) {
        Pvec[0] = new ProfileStructure();
        if( !extractProfile(ifile, pvec[0], *(Pvec[0])) )
            ERROR("For species #" << ispec << ", "<<varname<<" not understood");
        Pvec[1] = Pvec[0];
        Pvec[2] = Pvec[0];
    } else {
        ERROR("For species #" << ispec << ", "<<varname<<" needs 1 or 3 components.");
    }
}


// ---------------------------------------------------------------------------------------------------------------------
// Compute useful values (normalisation, time/space step, etc...)
// ---------------------------------------------------------------------------------------------------------------------
void PicParams::compute()
{
    // time-related parameters
    // -----------------------
    sim_time = (double)(n_time) * timestep;



    // grid/cell-related parameters
    // ----------------------------
    n_space.resize(3);
    cell_length.resize(3);
    cell_volume=1.0;
    if (nDim_field==res_space.size() && nDim_field==sim_length.size()) {

        // compute number of cells & normalized lengths
        for (unsigned int i=0; i<nDim_field; i++) {
            //!\todo MG Is this really necessary now ?
            cell_length[i] = 1.0/res_space[i];
            n_space[i]     = round(sim_length[i]/cell_length[i]);
            sim_length[i]  = (double)(n_space[i])*cell_length[i]; // ensure that nspace = sim_length/cell_length
            cell_volume   *= cell_length[i];
        }
        // create a 3d equivalent of n_space & cell_length
        for (unsigned int i=nDim_field; i<3; i++) {
            n_space[i]=1;
            cell_length[i]=0.0;
        }
        // compute number of cells per cluster
        n_cell_per_cluster = clrw * n_space[1] * n_space[2];

    } else {
        ERROR("Problem with the definition of nDim_field");
    }

    //!\todo (MG to JD) Are these 2 lines really necessary ? It seems to me it has just been done before
    n_space.resize(3, 1);
    cell_length.resize(3, 0.);	    //! \todo{3 but not real size !!! Pbs in Species::Species}

    n_space_global.resize(3, 1);	//! \todo{3 but not real size !!! Pbs in Species::Species}
    oversize.resize(3, 0);

}


// ---------------------------------------------------------------------------------------------------------------------
// Compute useful values for Species-related quantities
// ---------------------------------------------------------------------------------------------------------------------
void PicParams::computeSpecies()
{
    // Loop on all species
    for (unsigned int ispec=0; ispec< species_param.size(); ispec++) {

        // here I save the dimension of the pb (to use in BoundaryConditionType.h)
        species_param[ispec].nDim_fields = nDim_field;

        // define thermal velocity as \sqrt{T/m}
        species_param[ispec].thermalVelocity.resize(3);
        species_param[ispec].thermalMomentum.resize(3);
        for (unsigned int i=0; i<3; i++) {
            species_param[ispec].thermalVelocity[i]=sqrt(2.0*species_param[ispec].thermT[i]*const_e/species_param[ispec].mass);
            species_param[ispec].thermalMomentum[i]=species_param[ispec].mass * species_param[ispec].thermalVelocity[i];
        }

        species_param[ispec].weight = species_param[ispec].density / species_param[ispec].n_part_per_cell_for_weight;


    }//end loop on all species (ispec)

}

// ---------------------------------------------------------------------------------------------------------------------
// Compute normalization factor of some variables
// ---------------------------------------------------------------------------------------------------------------------
void PicParams::computeNormalization()
{
    const_e = 1.6021766208e-19;
    const_emass = 9.10938356e-31;
    const_c = 299792458.0;
    const_ephi0 = 8.854187817620389e-12;
    const_pi = 3.1415926;
    const_boltz = 1.3806e-23;
    const_h = 6.62606957e-34;

    norm_omiga0 = const_c / wavelength_SI;
    norm_density = const_ephi0*const_emass*norm_omiga0*norm_omiga0 / (const_e*const_e);
    norm_temperature = (const_emass*const_c*const_c) / const_e;
    norm_time = wavelength_SI / const_c;
    norm_length = wavelength_SI;
    norm_efield = const_emass * const_c * norm_omiga0 / const_e;
    norm_voltage = norm_temperature;
    norm_j = const_e / ( norm_time * norm_length * norm_length );

}




// ---------------------------------------------------------------------------------------------------------------------
// Set dimensions according to geometry
// ---------------------------------------------------------------------------------------------------------------------
void PicParams::setDimensions()
{
    if (geometry=="1d3v") {
        nDim_particle=1;
        nDim_field=1;
    } else if (geometry=="2d3v") {
        nDim_particle=2;
        nDim_field=2;
    } else if (geometry=="3d3v") {
        nDim_particle=3;
        nDim_field=3;
    } else if (geometry=="2drz") {
        nDim_particle=3;
        nDim_field=2;
    } else {
        ERROR("Geometry: " << geometry << " not defined");
    }
}



// ---------------------------------------------------------------------------------------------------------------------
// Printing out the data at initialisation
// ---------------------------------------------------------------------------------------------------------------------
void PicParams::print()
{

    // Numerical parameters
    // ---------------------
    MESSAGE("Numerical parameters");
    MESSAGE(1,"Geometry : " << geometry)
    MESSAGE(1,"(nDim_particle, nDim_field) : (" << nDim_particle << ", "<< nDim_field << ")");
    MESSAGE(1,"Interpolation_order : " <<  interpolation_order);
    MESSAGE(1,"(n_time,   timestep) : (" << n_time << ", " << timestep << ")");
    MESSAGE(1,"           timestep  = " << timestep/dtCFL << " * CFL");

    for ( unsigned int i=0 ; i<sim_length.size() ; i++ ){
        MESSAGE(1,"dimension " << i << " - (res_space, sim_length) : (" << res_space[i] << ", " << sim_length[i] << ")");
        MESSAGE(1,"            - (n_space,  cell_length) : " << "(" << n_space[i] << ", " << cell_length[i] << ")");
    }

    // Plasma related parameters
    // -------------------------
    MESSAGE("Plasma related parameters");
    MESSAGE(1,"n_species       : " << species_param.size());
    for ( unsigned int i=0 ; i<species_param.size() ; i++ ) {
        MESSAGE(1,"            species_type : "<< species_param[i].species_type);
    }


}

// Finds requested species in the list of existing species.
// Returns an array of the numbers of the requested species.
// Note that there might be several species that have the same "name" or "type"
//  so that we have to search for all possibilities.
vector<unsigned int> PicParams::FindSpecies( vector<string> requested_species)
{
    bool species_found;
    vector<unsigned int> result;
    unsigned int i;
    vector<string> existing_species;

    // Make an array of the existing species names
    existing_species.resize(0);
    for (unsigned int ispec=0 ; ispec<species_param.size() ; ispec++) {
        existing_species.push_back( species_param[ispec].species_type );
    }

    // Loop over group of requested species
    for (unsigned int rs=0 ; rs<requested_species.size() ; rs++) {
        species_found = false;
        // Loop over existing species
        for (unsigned int es=0 ; es<existing_species.size() ; es++) {
            if (requested_species[rs] == existing_species[es]) { // if found
                species_found = true;
                // Add to the list and sort
                for (i=0 ; i<result.size() ; i++) {
                    if (es == result[i]) break; // skip if duplicate
                    if (es <  result[i]) {
                        result.insert(result.begin()+i,es); // insert at the right place
                        break;
                    }
                }
                // Put at the end if not put earlier
                if (i == result.size()) result.push_back(es);
            }
        }
        if (!species_found)
            ERROR("Species `" << requested_species[rs] << "` was not found.");
    }

    return result;
}
