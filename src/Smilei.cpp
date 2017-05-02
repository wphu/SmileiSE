////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////                                                                                                                ////
////                                                                                                                ////
////                                   PARTICLE-IN-CELL CODE SMILEI                                                 ////
////                    Simulation of Matter Irradiated by Laser at Extreme Intensity                               ////
////                                                                                                                ////
////                          Cooperative OpenSource Object-Oriented Project                                        ////
////                                      from the Plateau de Saclay                                                ////
////                                          started January 2013                                                  ////
////                                                                                                                ////
////                                                                                                                ////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "Smilei.h"

#include <ctime>
#include <cstdlib>
#include <unistd.h>

#include <iostream>
#include <iomanip>

#include "InputData.h"
#include "PicParams.h"

#include "SmileiMPIFactory.h"
#include "GridFactory.h"
#include "SpeciesFactory.h"
#include "PartSourceFactory.h"
#include "CollisionsFactory.h"
#include "ElectroMagnFactory.h"
#include "InterpolatorFactory.h"
#include "ProjectorFactory.h"
#include "SolverFactory.h"
#include "SmileiIOFactory.h"
#include "PSIFactory.h"
#include "ElectroMagnBC_Factory.h"
#include "DiagnosticFactory.h"

#include "Timer.h"
#include <omp.h>

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
//                                                   MAIN CODE
// ---------------------------------------------------------------------------------------------------------------------
int main (int argc, char* argv[])
{
    //cout.setf( ios::fixed,  ios::floatfield ); // floatfield set to fixed

    // Define 2 MPI environments :
    //  - smpiData : to broadcast input data, unknown geometry
    //  - smpi (defined later) : to compute/exchange data, specific to a geometry
    SmileiMPI *smpiData= new SmileiMPI(&argc, &argv );

    // -------------------------
    // Simulation Initialization
    // -------------------------

    // Check for namelists (input files)
    vector<string> namelists(argv + 1, argv + argc);
    if (namelists.size()==0) ERROR("No namelists given!");

    // Send information on current simulation

    MESSAGE("                   _            _");
    MESSAGE(" ___           _  | |        _  \\ \\    ");
    MESSAGE("/ __|  _ __   (_) | |  ___  (_)  | |");
    MESSAGE("\\__ \\ | '  \\   _  | | / -_)  _   | |   Version  : " << __VERSION);
    MESSAGE("|___/ |_|_|_| |_| |_| \\___| |_|  | |   Date     : " << __COMMITDATE);
    MESSAGE("                                /_/    ");

    TITLE("Input data info");
    // Read the namelists file (no check!)
    InputData input_data(smpiData,namelists);

    // Read simulation & diagnostics parameters
    PicParams params(input_data);
    smpiData->init(params);
    smpiData->barrier();
    if ( smpiData->isMaster() ) params.print();
    smpiData->barrier();



    // Geometry known, MPI environment specified
    TITLE("General MPI environement");
    SmileiMPI* smpi = SmileiMPIFactory::create(params, smpiData);
    smpi->barrier();


    //>Initialize Grid
    TITLE("generate grid");
    Grid* grid = NULL;
    grid = GridFactory::create(params, input_data, smpi);
    smpi->barrier();
    smpi->scatterGrid(grid);

    // -------------------------------------------
    // Declaration of the main objects & operators
    // -------------------------------------------

    // ---------------------------
    // Initialize Species & Fields
    // ---------------------------


    TITLE("Initializing particles");
    // Initialize the vecSpecies object containing all information of the different Species
    // ------------------------------------------------------------------------------------
    // vector of Species (virtual)
    vector<Species*> vecSpecies = SpeciesFactory::createVector(params, smpi);
    smpi->barrier();
    for (unsigned int ispec=0 ; ispec<params.species_param.size(); ispec++)
    {
        //vecSpecies[ispec]->printAvgVelocity();
    }

    // Initialize the electromagnetic fields and interpolation-projection operators
    // according to the simulation geometry
    // ----------------------------------------------------------------------------

    TITLE("Initializing ElectroMagn Fields");
    // object containing the electromagnetic fields (virtual)
    ElectroMagn* EMfields = ElectroMagnFactory::create(params, input_data, smpi);
    smpi->barrier();

    TITLE("Initializing Fields Bounary Condition");
    vector<ElectroMagnBC*> vecEmBC = ElectroMagnBCFactory::create(params);
    smpi->barrier();

    TITLE("Creating Diagnostic");
    Diagnostic*  diag  = DiagnosticFactory::create(params, smpi, EMfields);
    smpi->barrier();


    TITLE("Creating Solver");
    Solver* solver = SolverFactory::create(params, input_data, grid, smpi);


    TITLE("Creating PartSource");
    vector<PartSource*> vecPartSource = PartSourceFactory::create(params, input_data, vecSpecies, smpi);
    smpi->barrier();


    // Initialize the collisions (vector of collisions)
    // ------------------------------------------------------------------------------------
    TITLE("Creating Collisions");
    vector<Collisions*> vecCollisions = CollisionsFactory::create(params, input_data, vecSpecies, smpi);
    smpi->barrier();

    TITLE("Creating PSI");
    vector<PSI*> vecPSI = PSIFactory::create(params, input_data, vecSpecies, smpi);
    smpi->barrier();

    TITLE("Creating Interp/Proj");
    // interpolation operator (virtual)
    Interpolator* Interp = InterpolatorFactory::create(params, smpi);

    // projection operator (virtual)
    Projector* Proj = ProjectorFactory::create(params, smpi);
    smpi->barrier();

    //Create mpi i/o environment
    TITLE("Creating IO output environment");
    SmileiIO*  sio  = SmileiIOFactory::create(params, smpi, EMfields, vecSpecies, diag);
    smpi->barrier();

    // ------------------------------------------------------------------------
    // Initialize the simulation times time_prim at n=0 and time_dual at n=+1/2
    // ------------------------------------------------------------------------
    unsigned int stepStart = sio->stepStart, stepStop=params.n_time;
    // time at integer time-steps (primal grid)
    double time_prim = stepStart * params.timestep;
    // time at half-integer time-steps (dual grid)
    double time_dual = (stepStart +0.5) * params.timestep;

    int itime = stepStart;


    TITLE("Solve the field first time before PIC loop");
    (*solver)(EMfields, smpi);
    smpi->barrier();

    // Count timer
    vector<Timer> timer(11);
    timer[0].init(smpi, "Total time");
    timer[1].init(smpi, "EmitLoad");
    timer[2].init(smpi, "Collide");
    timer[3].init(smpi, "Interpolate and Move");
    timer[4].init(smpi, "MPI Exchange Particle");
    timer[5].init(smpi, "Absorb paritcle (2D)");
    timer[6].init(smpi, "Project Particle");
    timer[7].init(smpi, "PSI");
    timer[8].init(smpi, "Diagnostic");
    timer[9].init(smpi, "Fields Solve");
    timer[10].init(smpi,"Write IO");

    // ------------------------------------------------------------------
    //                     HERE STARTS THE PIC LOOP
    // ------------------------------------------------------------------
    TITLE("Time-Loop is started: number of time-steps n_time = " << params.n_time);
    smpi->barrier();
    timer[0].restart();
    while(itime <= stepStop)
    {
        itime++;
        time_prim += params.timestep;
        time_dual += params.timestep;

        // ================== EmitLoad =========================================
        //> add Particle Source: emit from boundary or load in some region
        timer[1].restart();
        for (unsigned int iPS=0 ; iPS<vecPartSource.size(); iPS++)
        {
            vecPartSource[iPS]->emitLoad(params,smpi,vecSpecies,itime, EMfields);
        }
        timer[1].update();


        // ================== Collide =========================================
        timer[2].restart();
        if(itime % params.timesteps_collision == 0)
        {
            for (unsigned int icoll=0 ; icoll<vecCollisions.size(); icoll++)
            {
                vecCollisions[icoll]->collide(params, smpi, EMfields, vecSpecies,itime);
            }
        }
                timer[2].update();

        // ================== Interpolate and Move ===============================
        int tid(0);
        timer[3].restart();
        for (unsigned int ispec=0 ; ispec<params.species_param.size(); ispec++)
        {
            vecSpecies[ispec]->dynamics(time_dual, ispec, EMfields, Interp, Proj, smpi, params);
        }
        timer[3].update();

        // ================== MPI Exchange Particle ============================================
        timer[4].restart();
        for (unsigned int ispec=0 ; ispec<params.species_param.size(); ispec++)
        {
            for ( int iDim = 0 ; iDim<(int)params.nDim_particle ; iDim++ )
            {
                smpi->exchangeParticles(vecSpecies[ispec], ispec, params, tid, iDim);
            }
            vecSpecies[ispec]->sort_part(); // Should we sort test particles ?? (JD)
        }
        timer[4].update();

        // ================== Absorb Particle for 2D ====================================
        timer[5].restart();
        for (unsigned int ispec=0 ; ispec<params.species_param.size(); ispec++)
        {
            if(params.geometry == "2d3v") {
                vecSpecies[ispec]->absorb2D(time_dual, ispec, grid, smpi, params);
            }
        }
        timer[5].update();

        // ================== Project Particle =========================================
        timer[6].restart();
        for (unsigned int ispec=0 ; ispec<params.species_param.size(); ispec++)
        {
            EMfields->restartRhoJs(ispec, 0);
            vecSpecies[ispec]->Project(time_dual, ispec, EMfields, Proj, smpi, params);
        }
        timer[6].update();

        // ================== Plasma Surface Interacton ==================================
        timer[7].restart();
        for (unsigned int ipsi=0 ; ipsi<vecPSI.size(); ipsi++)
        {
            vecPSI[ipsi]->performPSI(params,smpi,vecSpecies,itime, EMfields);
        }
        timer[7].update();

        // ================== Run Diagnostic =============================================
        timer[8].restart();
        diag->run(smpi, vecSpecies, EMfields, itime);
        timer[8].update();

        // ================== Solve Electromagnetic Fields ===============================
        timer[9].restart();
        EMfields->restartRhoJ();
        EMfields->computeTotalRhoJ();
        (*solver)(EMfields, smpi);
        timer[9].update();

        // ================== Write IO ====================================================
        timer[10].restart();
        if(params.ntime_step_avg)
        {
            EMfields->incrementAvgFields(itime, params.ntime_step_avg);
        }
        if(itime % params.dump_step == 0){
            EMfields->gatherAvgFields(smpi);
            MESSAGE("time step = "<<itime);
        }
        sio->write(params, smpi, EMfields, vecSpecies, diag, itime);
        if(itime % params.timesteps_restore == 0)
        {
            sio->storeP(params, smpi, vecSpecies, itime);
        }
        timer[10].update();

    }//END of the time loop
    sio->endStoreP(params, smpi, vecSpecies, itime);
    smpi->barrier();




    // ------------------------------------------------------------------
    //                      HERE ENDS THE PIC LOOP
    // ------------------------------------------------------------------
    timer[0].update();
    for( int i = 0; i < timer.size(); i++)
    {
        TITLE(timer[i].name() << " = " << timer[i].getTime());
    }


    // ------------------------------------------------------------------------
    // check here if we can close the python interpreter
    // ------------------------------------------------------------------------
    TITLE("Cleaning up python runtime environement");
    input_data.cleanup();


    //double timElapsed=smpiData->time_seconds();
    //if ( smpi->isMaster() ) MESSAGE("Time in time loop : " << timElapsed );

    TITLE("Time profiling :");
    //timer[0].update();
    //timer[0].print();


    // ------------------------------------------------------------------
    //                      Temporary validation diagnostics
    // ------------------------------------------------------------------

    // temporary EM fields dump in Fields.h5


    // ------------------------------
    //  Cleanup & End the simulation
    // ------------------------------
    delete Proj;
    delete Interp;
    smpi->barrier();
    delete EMfields;

    for (unsigned int iPS=0 ; iPS<vecPartSource.size(); iPS++) delete vecPartSource[iPS];
    vecPartSource.clear();

    for(unsigned int i=0; i<vecCollisions.size(); i++) delete vecCollisions[i];
    vecCollisions.clear();

    for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) delete vecSpecies[ispec];
    vecSpecies.clear();
    smpi->barrier();
    TITLE("END");
    delete sio;
    delete smpi;
    delete smpiData;
    return 0;

}//END MAIN
