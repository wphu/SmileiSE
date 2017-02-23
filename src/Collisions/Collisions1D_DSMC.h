/*
Collisions1D_DSMC class
!!! The code is mainly from the java code in the website: https://www.particleincell.com/2013/dsmc0-in-java/
!!! For now, the code does not consider the SubCell, if necessary, the SubCell part will be added
*/

#ifndef COLLISIONS1D_DSMC_H
#define COLLISIONS1D_DSMC_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "Collisions1D.h"
#include "H5.h"


using namespace std;

class Collisions1D_DSMC : public Collisions1D
{

public:
    //! Constructor for Collisions1D between two species
    Collisions1D_DSMC(PicParams&,std::vector<Species*>&,SmileiMPI*,unsigned int, vector< vector<unsigned int> >);
    ~Collisions1D_DSMC();

    double cross_section(double ke);

    //! Method called in the main smilei loop to apply collisions at each timestep
    void collide(PicParams&, SmileiMPI* smpi, ElectroMagn* fields, std::vector<Species*>&, int);


    // ================DSMC parameters============================
    vector< vector<unsigned int> > species_group;
    vector< vector<int> > species_count;        //particle number in a Cell of one species
    int NumSpecies;		//MNSP is the maximum number of molecular species
	int NumSpGroups;		//MNSG is the number of species groups for collision sampling
    int NumCells;		//MNC  is the maximum number of sub

	int SampleNumIt;	//NIS is the number of time steps between samples 4

    //sampled_info(N,M,L) sampled information on species_group L in cell M
    //---- N=1 number sum
    //----N=2,3,4 sum of u,v,w
    //----N=5,6,7 sum of u*u,v*v,w*w
    struct SampledInfo
    {
        double num_sum;
        double v_sum[3];
        double v2_sum[3];
    };
    vector< vector<SampledInfo> > sampled_info;

    // [NumCells][NumSpGroups][NumSpGroups]; --CCG(N,M,L,K) is for collisions between species groups L-K in cell M
    struct SpGroupColData
    {
        double sigma_g_max;				//----N=1 is the maximum value of (relative speed)*(coll. cross-section)
        double rem;
    };
    vector< vector< vector<SpGroupColData> > > spg_col_data;

    // IC(N,M,L) information on the molecules of species group L in cell M
    struct CellInfo
    {
        int start;
        int count;
    };
    vector< vector<CellInfo> > cell_info;

    // SPM(N,M,L) information on the interaction between L-M molecules
    struct SpeciesInteraction
    {
        // The molecule and atom diameters are from the following refs:
        // Bird' book: page 410
        // The mathematical theory of non-uniform gases: page238, page237, page228,
        double sigma;			//--N=1  the reference cross-section (diameter in the data)
        double ref_temp;		//--N=2  the reference temperature
        // visc_temp_index is from Bird' book: page411
        double visc_temp_index;				//--N=3  the viscosity-temperature power law
        double vss_scat_inv;				//--N=4  the reciprocal of the VSS scattering parameter
        double reduced_mass;	//--N=5  the reduced mass
        double gamma;			//--N=6  the Gamma function of (5/2 - viscosity-temperature power law)
    };
    vector< vector<SpeciesInteraction> > species_interaction;

    /*comp block*/
	double Spwt;    ////FNUM  is the number of real molecules represented by a simulated mol.
	double DeltaT;     ////DTM is the time step
	int SteadyFlowEstimateCycles;     ////NPS is the estimated number of samples to steady flow

	/*geom block*/
    double CellVolume;
	double CellWidth;      ////CW is the cell width
	int SubCellsPerCell;     ////NSC is the number of sub-cells per cell
	double XMIN;      ////XF is the minimum x coordinate
	double XMAX;      ////XR is the maximum x coordinate

	/*const block*/
	double PI;      ////PI is pi and SPI is the square root of pi
	double SPI;
	double BOLTZ;   ////BOLTZ is the Boltzmann constant

	/*elast block*/
	double VRC[3];
	double VRR;
	double VR;
	int LS,MS;
	int L,M;
    int iSL,iSM;            // Species number of group L, M
    int indexL,indexM;      // particle number in a Cell of species iSL,iSM
    int iL,iM;              // particle number of species iSL,iSM
    Species *sL, *sM;
    Particles *pL, *pM;
	int MSC;
	double CVR;
	//int K,L,M,N;
	int MM,MN,NN;


    double StreamTemp;
    int TotNumSamples;
    int TotNumSelections;
    double CVR_MAX;
    int sample_it;

    void COLLM(vector<Species*>& vecSpecies);
    void SAMPLE_INIT();
    void INDEXM(vector<Species*>& vecSpecies);
    void INIT0(vector<Species*>& vecSpecies);
    void SELECT(int ibin, vector<Species*>& vecSpecies);
    void ELASTIC(vector<Species*>& vecSpecies);
    void SAMPLE0();
    double GAM(double X);



private:
    inline double scatter_particles(Particles* particle1, int iPart1, Particles* particle2, int iPart2);

};


#endif
