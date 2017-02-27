#include "Collisions1D_DSMC.h"
#include "SmileiMPI.h"
#include "Field2D.h"
#include "H5.h"


#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <sstream>

using namespace std;


// Constructor
Collisions1D_DSMC::Collisions1D_DSMC(PicParams& params, vector<Species*>& vecSpecies, SmileiMPI* smpi,
                       unsigned int n_col,
                       vector< vector<unsigned int> > sg)
: Collisions1D(params)
{

    n_collisions    = n_col;
    species_group   = sg;

    DeltaT = params.timesteps_DSMC * params.timestep;


    // Calculate total number of bins In The Current process
    nbins = vecSpecies[0]->bmin.size();
    totbins = nbins;

    SampleNumIt = 4;
    sample_it = 0;
    PI = params.const_pi;
    SPI = sqrt(PI);
    BOLTZ = params.const_boltz;
    CellVolume = pow(params.cell_length[0], 3);
    StreamTemp = 300.0;
    TotNumSelections = 0;
    // Spwt is weight of paritcles
    // !!!Now all collision species should have the same weight
    // the weight has been divided by the cell volume, so equation(11.3) should not contain the cell volume
    Spwt = vecSpecies[0]->species_param.weight;

    // NumSpecies is the number of all species for DSMC, not including electron or ions
    // NumSpGroups is the number of species groups
    NumSpecies = 0;
    NumSpGroups = species_group.size();
    // init SpeciesList[iSpecies]: list of all species participating the collision
    for(int iSG=0; iSG<NumSpGroups; iSG++)
    {
        NumSpecies += species_group[iSG].size();
        for(int iS = 0; iS < species_group[iSG].size(); iS++)
        {
            SpeciesList.push_back( species_group[iSG][iS] );
        }
    }

    // init species_count[totbins][NumSpecies]
    species_count.resize(totbins);
    for(int ibin = 0; ibin < totbins; ibin++)
    {
        species_count[ibin].resize( vecSpecies.size() );
    }

    // init sampled_info[totbins][NumSpGroups], the accumulated info all the time
    // init cell_info[totbins][NumSpGroups], the info at the current time
    // init spg_col_data[totbins][NumSpGroups][NumSpGroups]: max cross section and rem
    sampled_info.resize(totbins);
    cell_info.resize(totbins);
    spg_col_data.resize(totbins);
    for(int ibin = 0; ibin < totbins; ibin++)
    {
        sampled_info[ibin].resize(NumSpGroups);
        cell_info[ibin].resize(NumSpGroups);
        spg_col_data[ibin].resize(NumSpGroups);

        for(int iSG=0; iSG<NumSpGroups; iSG++)
        {
            spg_col_data[ibin][iSG].resize(NumSpGroups);
        }
    }

    species_interaction.resize(NumSpecies);
    for( int iSpec = 0; iSpec < NumSpecies; iSpec++ )
    {
        species_interaction[iSpec].resize(NumSpecies);
    }


    INIT0(vecSpecies);
    SAMPLE_INIT();

}

Collisions1D_DSMC::~Collisions1D_DSMC()
{

}

// Calculates the collisions for a given Collisions1D object
// Now only support 3 groups, each grous has the approximate mass
// like [H, D, T] or [S, Cl, Ar]
void Collisions1D_DSMC::collide(PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies, int itime)
{
    if( (itime % params.timesteps_DSMC) == 0 )
    {
        INDEXM(vecSpecies);
        COLLM(vecSpecies);
        sample_it++;
        while(sample_it == SampleNumIt) {
            sample_it = 0;
            SAMPLE0();
        }
    }


}


void Collisions1D_DSMC::INIT0(vector<Species*>& vecSpecies)
{
    int iSM, iSN;
    int iSL, iSK;
    for(int N=0; N<NumSpecies; N++)
    {
        iSN = SpeciesList[N];
        for(int M=0; M<NumSpecies; M++)
        {
            iSM = SpeciesList[M];
            species_interaction[N][M].sigma = 0.25*PI*pow((vecSpecies[iSN]->species_param.diameter + vecSpecies[iSM]->species_param.diameter), 2);
            species_interaction[N][M].ref_temp = 0.5*(vecSpecies[iSN]->species_param.ref_temperature + vecSpecies[iSM]->species_param.ref_temperature);
            species_interaction[N][M].visc_temp_index = 0.5*(vecSpecies[iSN]->species_param.visc_temp_index + vecSpecies[iSM]->species_param.visc_temp_index);
            species_interaction[N][M].vss_scat_inv = 0.5*(vecSpecies[iSN]->species_param.vss_scat_inv + vecSpecies[iSM]->species_param.vss_scat_inv);
            species_interaction[N][M].reduced_mass = vecSpecies[iSN]->species_param.mass * vecSpecies[iSM]->species_param.mass /
                                                    (vecSpecies[iSN]->species_param.mass + vecSpecies[iSM]->species_param.mass);
            species_interaction[N][M].gamma = GAM(2.5 - species_interaction[N][M].visc_temp_index);
        }
    }

    for(int ibin=0; ibin<totbins; ibin++)
    {
        for (int L=0;L<NumSpGroups;L++)
        {
            for (int K=0;K<NumSpGroups;K++)
            {
                //*--the maximum value of the (rel. speed)*(cross-section) is set to a
				//*--reasonable, but low, initial value and will be increased as necessary
                double VR = 0.5 * sqrt(2.0 * BOLTZ * vecSpecies[SpeciesList[0]]->species_param.ref_temperature /
                            vecSpecies[SpeciesList[0]]->species_param.mass);
                double VRR = VR * VR;
                double sigma_g_max_temp = VR*species_interaction[0][0].sigma*
                				pow( 2.*BOLTZ*species_interaction[0][0].ref_temp/(species_interaction[0][0].reduced_mass*VRR),
                				species_interaction[0][0].visc_temp_index-0.5 ) /
                                species_interaction[0][0].gamma;

                spg_col_data[ibin][L][K].sigma_g_max = sigma_g_max_temp;
                //cout<<"sigma_g_max_temp = "<<sigma_g_max_temp<<endl;
                //spg_col_data[ibin][L][K].sigma_g_max = species_interaction[0][0].sigma * 300. * sqrt(StreamTemp/300.);
                double RF = (double)random() /RAND_MAX;
				spg_col_data[ibin][L][K].rem=RF;
            }
        }
    }
}


void Collisions1D_DSMC::SAMPLE_INIT()
{
    TotNumSamples = 0;
    for(int ibin=0; ibin<totbins; ibin++)
    {
        for(int L=0; L<NumSpGroups; L++)
        {
            sampled_info[ibin][L].num_sum=1.0e-6;
            //for (int i=0;i<3;i++)
            //{
            //    sampled_info[N][L].v_sum[i]=0.0;
            //    sampled_info[N][L].v2_sum[i]=0.0;
            //}
        }
    }
}


void Collisions1D_DSMC::SAMPLE0()
{
    TotNumSamples++;
    for(int MM=0; MM<NumSpGroups; MM++)
    {
        for(int ibin=0; ibin<totbins; ibin++)
        {
            sampled_info[ibin][MM].num_sum += cell_info[ibin][MM].count;

        }
    }

}


// if SubCells need to be considered, spg_cell_info should be added
void Collisions1D_DSMC::INDEXM(vector<Species*>& vecSpecies)
{
    Species *sp;
    // In one species group, different species are treated as the same species
    for(int ibin=0; ibin<totbins; ibin++)
    {
        for(int MM=0; MM<NumSpGroups; MM++)
        {
            cell_info[ibin][MM].count = 0;
            for(int iSpecies=0; iSpecies<species_group[MM].size(); iSpecies++)
            {
                sp = vecSpecies[ species_group[MM][iSpecies] ];
                species_count[ibin][ species_group[MM][iSpecies] ] = ( sp->bmax[ibin] - sp->bmin[ibin] );
                cell_info[ibin][MM].count += species_count[ibin][ species_group[MM][iSpecies] ];
            }
        }
    }
}


// Calculates the collisions for a given Collisions1D object
void Collisions1D_DSMC::COLLM(vector<Species*>& vecSpecies)
{
    double ASEL;
    TotNumCols = 0;
    // Loop on bins(cells)

    for(int ibin = 0; ibin < totbins; ibin++ )
    {
        // Loop on species, NN == MM: the like species; NN != MM unlike species
        // Like species may contain some different species with like mass, like [H, D, T]
        for(NN = 0; NN < NumSpGroups; NN++ )
        {
            for(MM = 0; MM < NumSpGroups; MM++)
            {
                double SN = sampled_info[ibin][MM].num_sum;
                double AVN;
                if(SN > 1.0) {
                    AVN = SN / TotNumSamples;
                }
                else {
                    AVN = cell_info[ibin][MM].count;
                }

                //*--ASEL is the number of pairs to be selected, see eqn (11.5)
                // rem is to handle the problem: the particle number is odd
				ASEL = 0.5 * cell_info[ibin][NN].count * AVN * Spwt * spg_col_data[ibin][NN][MM].sigma_g_max * DeltaT
                              + spg_col_data[ibin][NN][MM].rem;
                int NSEL = ASEL;
                spg_col_data[ibin][NN][MM].rem=ASEL-NSEL;

                if(NSEL > 0) {
                    //*--if there are insufficient molecules to calculate collisions,
					//*--the number NSEL is added to the remainer CCG(2,N,NN,MM)
                    if(  ( (NN!=MM)  &&  (cell_info[ibin][NN].count<1 || cell_info[ibin][MM].count<1) )
					   ||( (NN==MM) && (cell_info[ibin][NN].count<2) )  ) {
                        spg_col_data[ibin][NN][MM].rem=spg_col_data[ibin][NN][MM].rem+NSEL;
                    }
                    else {
                        TotNumSelections = TotNumSelections+NSEL;
                        double CVR_MAX = spg_col_data[ibin][NN][MM].sigma_g_max;
                        for(int ISEL=0;ISEL<NSEL;ISEL++)
						{
							/*note, this sets L and M as well as LS and MS*/
							SELECT(ibin, vecSpecies);

                            //cout<<"ASEL ========== "<<ASEL<<" "<<TotNumSelections<<endl;

							//*--if necessary, the maximum product in CVM is upgraded
							if(CVR > CVR_MAX) CVR_MAX = CVR;

							//*--the collision is accepted with the probability of eqn (11.6)
                            double RF = (double)rand() / RAND_MAX;
							if(RF < CVR / spg_col_data[ibin][NN][MM].sigma_g_max) {
								TotNumCols++;
								//SumColPairSeparations=SumColPairSeparations+Math.abs(part[L].x-part[M].x);
								//COL[LS][MS]++;
								//COL[MS][LS]++;

								ELASTIC(vecSpecies);
							}
						}
						spg_col_data[ibin][NN][MM].sigma_g_max=CVR_MAX;
                    }

                }

            }
        }
    }
    //MESSAGE("TotNumCols = "<<TotNumCols);
    //MESSAGE("TotNumSelections = "<<TotNumSelections);

}


// calculate the the two species, particles and iPart particepating the collision
void Collisions1D_DSMC::SELECT(int ibin, vector<Species*>& vecSpecies)
{
    double RF;
    int K;
    int iSpecies;

    RF = (double)random() / RAND_MAX;
    K = (int)( RF*cell_info[ibin][NN].count );
    indexL = K;
    iSpecies = 0;
    iSL = species_group[NN][iSpecies];
    // iSL and iSM are species number: Species *s = vecSpecies[iSL/iSM]
    // indexL and indexM is the particles indexes
    while (iSpecies < species_group[NN].size() && indexL >= species_count[ibin][iSL]) {
        indexL -= species_count[ibin][iSL];
        iSpecies++;
        iSL = species_group[NN][iSpecies];
    }

    do {
        RF = (double)random() / RAND_MAX;
        K = (int)( RF * cell_info[ibin][MM].count );
        indexM = K;
        iSpecies = 0;
        iSM = species_group[MM][iSpecies];
        while (iSpecies < species_group[MM].size() && indexM >= species_count[ibin][iSM]) {
            indexM -= species_count[ibin][iSM];
            iSpecies++;
            iSM = species_group[MM][iSpecies];
        }
    } while (iSL == iSM && indexL == indexM);

    sL = vecSpecies[iSL];
    sM = vecSpecies[iSM];
    pL = &(sL->particles);
    pM = &(sM->particles);
    iL = sL->bmin[ibin] + indexL;
    iM = sM->bmin[ibin] + indexM;
    for(int i=0; i<3; i++)
    {
        VRC[i] = pL->momentum(i,iL) - pM->momentum(i,iM);
    }
    VRR = VRC[0]*VRC[0] + VRC[1]*VRC[1] + VRC[2]*VRC[2];
    VR = sqrt(VRR);
    // Calculate cross section * relative velocity
    // the collision cross-section is based on eqn (4.63)
    CVR=VR*species_interaction[iSL][iSM].sigma*
    				pow( 2.*BOLTZ*species_interaction[iSL][iSM].ref_temp/(species_interaction[iSL][iSM].reduced_mass*VRR),
    				species_interaction[iSL][iSM].visc_temp_index-0.5 ) /
                    species_interaction[iSL][iSM].gamma;
}


void Collisions1D_DSMC::ELASTIC(vector<Species*>& vecSpecies)
{
    double VRCP[3];		//VRCP(3) are the post-collision components of the relative velocity
    double VCCM[3];		//VCCM(3) are the components of the centre of mass velocity
    double RML = species_interaction[iSL][iSM].reduced_mass / sM->species_param.mass;
    double RMM = species_interaction[iSL][iSM].reduced_mass / sL->species_param.mass;
    double A,B,C,D;
    double OC,SC;

    for (int i=0;i<3;i++)
	{
		VCCM[i] = RML*pL->momentum(i,iL) + RMM*pM->momentum(i,iM);
	}

    //use the VHS logic
    double RF = (double)random() / RAND_MAX;
	B = 2.0 * RF - 1.0;	//B is the cosine of a random elevation angle
	A = sqrt(1.0 - B * B);
	VRCP[0] = B * VR;
	C = 2.0 * PI * RF; //C is a random azimuth angle
	VRCP[1] = A * cos(C) * VR;
	VRCP[2] = A * sin(C) * VR;

    for (int i = 0; i < 3; i++)
	{
        pL->momentum(i,iL) = VCCM[i] + VRCP[i]*RMM;
        pM->momentum(i,iM) = VCCM[i] - VRCP[i]*RML;
	}
}




double Collisions1D_DSMC::GAM(double X)
{
    double A=1.;
    double Y=X;
    if (Y<1.0) { A=A/Y; }
    else
    {
        do {
            Y=Y-1;
            if (Y>=1.) { A=A*Y; }
        } while (Y>=1.);

    }

    double GAM=A*(1.-0.5748646*Y+0.9512363*Y*Y-0.6998588*Y*Y*Y+
            0.4245549*Y*Y*Y*Y-0.1010678*Y*Y*Y*Y*Y);
    return GAM;
}
