#include "Collisions1D_Coulomb.h"
#include "SmileiMPI.h"
#include "Field1D.h"
#include "Field2D.h"
#include "ElectroMagn.h"
#include "H5.h"


#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <sstream>

using namespace std;


// Constructor
Collisions1D_Coulomb::Collisions1D_Coulomb(PicParams& params, vector<Species*>& vecSpecies, SmileiMPI* smpi,
                       unsigned int n_col,
                       vector<unsigned int> sg1,
                       vector<unsigned int> sg2,
                       double coulomb_logarithm,
                       bool intra_col,
                       int debug_every_temp)
: Collisions1D(params)
{
    n_collisions    = n_col;
    species_group1  = sg1;
    species_group2  = sg2;
    coulomb_log     = coulomb_logarithm;
    intra_collisions= intra_col;
    debug_every     = debug_every_temp;
    start           = 0;

    twoPi = 2.0 * const_pi;
    e_ov_ephi0 = const_e / const_ephi0;
    time_coulomb = params.timesteps_coulomb * params.timestep * params.zoom_collision;

    n_Species        = sg1.size();
    n_particle.resize(n_Species);
    density.resize(n_Species);
    n_particle_begin.resize(n_Species);
    index.resize(n_Species);
    bmin.resize(n_Species);

    nbins = vecSpecies[0]->bmin.size();
    debye_length.resize(nbins);
    log_coulomb.resize(nbins);


    cal_collision_pairs(vecSpecies);
}

Collisions1D_Coulomb::~Collisions1D_Coulomb()
{

}


void Collisions1D_Coulomb::collide(PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies, int itime)
{
    unsigned int i1, i2, iS1, iS2;
    Species   *s1, *s2;
    Particles *p1, *p2;
    double gx, gy, gz, g_magnitude, g_square, g_3, g_p, hx, hy, hz, s, A12, cosX, sinX, phi;
    double charge1, charge2, charge;


    if( (itime % params.timesteps_coulomb) == 0 )
    {
        cal_log_coulomb(fields);

        // Loop on bins
        for (unsigned int ibin=0 ; ibin<nbins ; ibin++)
        {
            //cout<<"bins "<<nbins<<"  "<<ibin<<endl;
            // calculate n_particle[], density[], bmin[]
            density_all = 0;
            for(int i = 0; i < n_Species; i++)
            {
              //cout<<"n_Species "<<n_Species<<" "<<i<<endl;
              s1 = vecSpecies[ species_group1[i] ];
              double weight = s1->species_param.weight;
              n_particle[i] = s1->bmax[ibin] - s1->bmin[ibin];
              density[i] = n_particle[i] *weight;
              density_all += density[i];
              bmin[i] = s1->bmin[ibin];
              n_particle_begin[i] = 0;
            }
            if(density_all == 0.0)
            {
              continue;
            }
            density_all_inv = 1.0 / density_all;

            //cout<<"11111"<<endl;
            // calculate index[][]
            for(int i = 0; i < n_collisionPairs - n_Species; i ++)
            {
              iS1 = iSpec1[i];
              iS2 = iSpec2[i];
              n_indexes[i] = n_particle[iS1] * density[iS2] * density_all_inv;
              n_particle_begin[iS1] += n_indexes[i];
              n_particle_begin[iS2] += n_indexes[i];
              isOdd[i] = false;
            }
            //cout<<"22222"<<endl;
            for(int i = n_collisionPairs - n_Species; i < n_collisionPairs; i ++)
            {
              //cout<<"aaaa"<<"i" <<endl;
              int i_index = i - (n_collisionPairs - n_Species);
              iS1 = iSpec1[i];
              n_indexes[i] = n_particle[iS1] - n_particle_begin[iS1];
              if(n_indexes[i] % 2 == 0)
              {
                isOdd[i] = false;
                n_indexes[i] /= 2;

                index[i_index].resize( n_particle[iS1] );
                for (unsigned int j = 0; j < n_particle[iS1]; j++)
                {
                  index[i_index][j] = j;
                }
                random_shuffle(index[i_index].begin(), index[i_index].end());
              }
              else
              {
                isOdd[i] = true;
                n_indexes[i] = (n_indexes[i] + 1) / 2;

                index[i_index].resize( n_particle[iS1] );
                for (unsigned int j = 0; j < n_particle[iS1]; j++)
                {
                  index[i_index][j] = j;
                }
                random_shuffle(index[i_index].begin(), index[i_index].end());
                index[i_index].push_back( n_particle[iS1] - 1 );
              }
            }

            //cout<<"3333"<<endl;
            for(int i = 0; i < n_Species; i++)
            {
              n_particle_begin[i] = 0;
            }

            // loop on collisionPairs
            for(int i = 0; i < n_collisionPairs; i++)
            {
              //cout<<"n_collisionPairs "<<n_collisionPairs<<"  "<<i<<endl;
              iS1 = iSpec1[i];
              iS2 = iSpec2[i];
              s1 = vecSpecies[ species_group1[iS1] ];               s2 = vecSpecies[ species_group1[iS2] ];
              p1 = &(s1->particles);                                p2 = &(s2->particles);

              if(log_coulomb[ibin] == 0.0)
              {
                continue;
              }
              A12 = gamma2[i] * density_all * log_coulomb[ibin];


              // loop on collision pairs of particles
              for (unsigned int j = 0; j < n_indexes[i]; j++)
              {
                  //cout<<"n_indexes[i] "<<n_indexes[i]<<"  "<<j<<endl;
                  if(isIntra[i])
                  {
                      i1 = index[iS1][j + n_particle_begin[iS1]] + bmin[iS1];
                      i2 = index[iS1][j + n_particle_begin[iS2] + n_indexes[i]] + bmin[iS1];
                  }
                  else
                  {
                    i1 = index[iS1][j + n_particle_begin[iS1]] + bmin[iS1];
                    i2 = index[iS2][j + n_particle_begin[iS2]] + bmin[iS2];
                  }

                  // Calculate stuff
                  gx = p1->momentum(0,i1) - p2->momentum(0,i2);
                  gy = p1->momentum(1,i1) - p2->momentum(1,i2);
                  gz = p1->momentum(2,i1) - p2->momentum(2,i2);
                  g_square = gx*gx + gy*gy + gz*gz;
                  g_magnitude = sqrt( g_square );
                  g_3 = g_square * g_magnitude;
                  g_p = sqrt( gy*gy + gz*gz );
                  if(g_p == 0.0 || g_3 == 0.0)
                  {
                      continue;
                  }
                  // the formula below the equation (101)
                  s = A12 * time_coulomb / g_3;
                  if(isOdd[j])
                  {
                    int n_temp = n_indexes[i] * 2 -1;
                    s *= ( n_temp / (n_temp+1) );
                  }

                  // Pick the deflection angles according to Nanbu's theory
                  // ref: improved modeing of relativistic collisions and collisional ionization in paritcle in cell codes
                  cosX = cos_chi(s);
                  if(cosX < -1.0)
                  {
                    cosX = -1.0;
                    cout<<"cosX < -1.0"<<endl;
                  }
                  else if(cosX > 1.0)
                  {
                    cosX = 1.0;
                    cout<<"cosX > 1.0"<<endl;
                  }

                  if(cosX*cosX > 1.0)
                  {
                    cout<<"cosX * cosX > 1.0"<<endl;
                    continue;
                  }
                  sinX = sqrt( 1. - cosX*cosX );


                  //!\todo make a faster rand by preallocating ??
                  phi = twoPi * ((double)rand() / RAND_MAX);
                  hx = g_p * cos(phi);
                  hy = -( gx * gy * cos(phi) + g_magnitude * gz * sin(phi) ) / g_p;
                  hz = -( gx * gz * cos(phi) - g_magnitude * gy * sin(phi) ) / g_p;

                  // Apply the deflection

                  double momentum1[3];
                  double momentum2[3];

                  momentum1[0] = p1->momentum(0,i1) - mr2[i] * ( gx * (1-cosX) + hx * sinX );
                  momentum1[1] = p1->momentum(1,i1) - mr2[i] * ( gy * (1-cosX) + hy * sinX );
                  momentum1[2] = p1->momentum(2,i1) - mr2[i] * ( gz * (1-cosX) + hz * sinX );

                  momentum2[0] = p2->momentum(0,i2) + mr1[i] * ( gx * (1-cosX) + hx * sinX );
                  momentum2[1] = p2->momentum(1,i2) + mr1[i] * ( gy * (1-cosX) + hy * sinX );
                  momentum2[2] = p2->momentum(2,i2) + mr1[i] * ( gz * (1-cosX) + hz * sinX );

                  p1->momentum(0,i1) = momentum1[0];
                  p1->momentum(1,i1) = momentum1[1];
                  p1->momentum(2,i1) = momentum1[2];

                  p2->momentum(0,i2) = momentum2[0];
                  p2->momentum(1,i2) = momentum2[1];
                  p2->momentum(2,i2) = momentum2[2];

              } // end loop on pairs of particles
              n_particle_begin[iS1] += n_indexes[i];
              n_particle_begin[iS2] += n_indexes[i];
            }

        } // end loop on bins
    }

}

void Collisions1D_Coulomb::cal_collision_pairs(vector<Species*>& vecSpecies)
{
  int iS1, iS2;
  Species *s1, *s2;
  double m1, m2;
  double charge1, charge2;
  n_collisionPairs = n_Species + n_Species * (n_Species - 1) / 2;

  n_indexes.resize(n_collisionPairs);
  gamma1.resize(n_collisionPairs);
  gamma2.resize(n_collisionPairs);
  mr1.resize(n_collisionPairs);
  mr2.resize(n_collisionPairs);
  m12.resize(n_collisionPairs);
  m1_inv.resize(n_collisionPairs);
  m2_inv.resize(n_collisionPairs);
  isOdd.resize(n_collisionPairs);


  for(int i = 0; i < n_Species - 1; i++)
  {
    for(int j = i + 1; j < n_Species; j++)
    {
      iSpec1.push_back(i);
      iSpec2.push_back(j);
      isIntra.push_back(false);
    }
  }

  for(int i = 0; i < n_Species; i++)
  {
    iSpec1.push_back(i);
    iSpec2.push_back(i);
    isIntra.push_back(true);
  }

  for(int i = 0; i < n_collisionPairs; i++)
  {
    int iS1 = iSpec1[i];
    int iS2 = iSpec2[i];
    s1 = vecSpecies[ species_group1[iS1] ];             s2 = vecSpecies[ species_group1[iS2] ];
    m1 = s1->species_param.mass;                        m2 = s2->species_param.mass;
    charge1 = s1->species_param.charge;                 charge2 = s2->species_param.charge;
    m1_inv[i] = 1.0 / m1;
    m2_inv[i] = 1.0 / m2;

    // used in equation (104a) and (104b)
    mr1[i] = m1 / (m1 + m2);
    mr2[i] = m2 / (m1 + m2);
    // the reduced mass
    m12[i] = m1 * m2 / (m1 + m2);
    // used in equation (96)
    gamma1[i] = 4.0 * const_pi * const_ephi0 * m12[i] / abs(charge1 * charge2);
    // used in equation (95)
    gamma2[i] = 4.0 * const_pi / ( gamma1[i] * gamma1[i] );
  }

  //cout<<"collisionPairs"<<n_collisionPairs<<" "<<isIntra.size()<<" "<<iSpec1.size()<<endl;


}

void Collisions1D_Coulomb::cal_log_coulomb(ElectroMagn* fields)
{
  double g12_square;
  int iS1, iS2;
  double T1, T2, Vx, Vy, Vz;

  // calculate electron debye_length
  iS1 = 0;
  Field1D* T_electron = static_cast<Field1D*>(fields->T_s[ species_group1[iS1] ]);
  Field1D* density_electron = static_cast<Field1D*>(fields->rho_s[ species_group1[iS1] ]);
  for(int ibin = 0; ibin < nbins; ibin++)
  {
    if( (*density_electron)( ibin + oversize[0] ) == 0.0 )
    {
      debye_length[ibin] = 0.0;
    }
    else
    {
      debye_length[ibin] = sqrt(const_ephi0 * (*T_electron)( ibin + oversize[0] ) / ( (*density_electron)( ibin + oversize[0] ) * const_e));
    }
  }

  for(int i = 0; i < n_collisionPairs; i++)
  {
    iS1 = iSpec1[i];
    iS2 = iSpec2[i];
    Field1D* T1_s = static_cast<Field1D*>(fields->T_s[ species_group1[iS1] ]);
    Field1D* T2_s = static_cast<Field1D*>(fields->T_s[ species_group1[iS2] ]);

    Field1D* Vx1_s = static_cast<Field1D*>(fields->Vx_s[ species_group1[iS1] ]);
    Field1D* Vy1_s = static_cast<Field1D*>(fields->Vy_s[ species_group1[iS1] ]);
    Field1D* Vz1_s = static_cast<Field1D*>(fields->Vz_s[ species_group1[iS1] ]);

    Field1D* Vx2_s = static_cast<Field1D*>(fields->Vx_s[ species_group1[iS2] ]);
    Field1D* Vy2_s = static_cast<Field1D*>(fields->Vy_s[ species_group1[iS2] ]);
    Field1D* Vz2_s = static_cast<Field1D*>(fields->Vz_s[ species_group1[iS2] ]);

    for(int ibin = 0; ibin < nbins; ibin++)
    {
      T1 = (*T1_s)(ibin + oversize[0]) * 3.0 * const_e * m1_inv[i];
      T2 = (*T2_s)(ibin + oversize[0]) * 3.0 * const_e * m2_inv[i];
      Vx = (*Vx1_s)(ibin + oversize[0]) - (*Vx2_s)(ibin + oversize[0]);
      Vy = (*Vy1_s)(ibin + oversize[0]) - (*Vy2_s)(ibin + oversize[0]);
      Vz = (*Vz1_s)(ibin + oversize[0]) - (*Vz2_s)(ibin + oversize[0]);
      // the formula below equation (96)
      g12_square = T1 + T2 + Vx*Vx + Vy*Vy + Vz*Vz;
      if(debye_length[ibin] == 0.0 || g12_square == 0.0)
      {
        log_coulomb[ibin] = 0.0;
      }
      else
      {
        log_coulomb[ibin] = log( gamma1[i] * g12_square * debye_length[ibin] );
      }

    }
  }

}



/*

void Collisions1D_Coulomb::collide(PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies, int itime)
{

    unsigned int nbins = vecSpecies[0]->bmin.size(); // number of bins
    vector<unsigned int> *sg1, *sg2, index1, index2;
    unsigned int bmin1, bmin2;
    unsigned int nspec1, nspec2; // numbers of species in each group
    unsigned int npart1, npart2; // numbers of macro-particles in each group
    unsigned int npairs[3]; // number of pairs of macro-particles
    double n1, n2, n; // densities of particles
    unsigned int i1, i2, iSpec1[3], iSpec2[3], iS1, iS2;
    Species   *s1, *s2;
    Species   *s1_pre, *s2_pre;
    Particles *p1, *p2;
    double m1, m2, m12, mr1, mr2, W1, W2, gamma1, gamma2, gx, gy, gz, g_magnitude, g_square, g_3, g_p,
           g12_square, hx, hy, hz, s, T1, T2, V1, V2, A12,
           cosX, sinX, phi, twoPi, debye_length, debye_squared, time_coulomb;
    double e_ov_ephi0;
    double charge1, charge2, charge, Vx, Vy, Vz, ni, Ti, Vi;
    bool isOdd[3], intra[3];
    int n_particle_begin[2];
    int n_collisionPairs;



    if( (itime % params.timesteps_coulomb) == 0 )
    {
        //!!!=======species_group1 must contain one species like e, D1 and so on;
        // as same for species_group2;
        // species_group1 can be the same as species_group2
        sg1 = &species_group1;
        sg2 = &species_group2;

        iSpec1[0] = (*sg1)[0]; iSpec1[1] = (*sg1)[0]; iSpec1[2] = (*sg2)[0];
        iSpec2[0] = (*sg2)[0]; iSpec2[1] = (*sg1)[0]; iSpec2[2] = (*sg2)[0];

        // =========Calculate some constants in the formulas used below============
        twoPi = 2.0 * const_pi;
        e_ov_ephi0 = const_e / const_ephi0;
        time_coulomb = params.timesteps_coulomb * params.timestep * params.zoom_collision;

        // Loop on bins
        for (unsigned int ibin=0 ; ibin<nbins ; ibin++) {
            s1_pre = vecSpecies[(*sg1)[0]];         s2_pre = vecSpecies[(*sg2)[0]];
            W1 = s1_pre->species_param.weight;      W2 = s2_pre->species_param.weight;
            npart1 = s1_pre->bmax[ibin] - s1_pre->bmin[ibin];
            npart2 = s2_pre->bmax[ibin] - s2_pre->bmin[ibin];
            bmin1 = s1_pre->bmin[ibin];
            bmin2 = s2_pre->bmin[ibin];
            n_particle_begin[0] = 0;
            n_particle_begin[1] = 0;

            n1 = npart1 * W1;
            n2 = npart2 * W2;
            n = n1 + n2;

            if(npart1 + npart2 == 0)
            {
              continue;
            }
            npairs[0] = npart1 * npart2 / (npart1 + npart2);

            npairs[1] = npart1 - npairs[0];
            if(npairs[1] % 2 == 0)
            {
              npairs[1] /= 2;
              isOdd[1] = false;
            }
            else
            {
              if(npairs[1] == 1)
              {
                npairs[1] = 0;
              }
              else
              {
                npairs[1] = (npairs[1] + 1) / 2;
              }
              isOdd[1] = true;

            }

            npairs[2] = npart2 - npairs[0];
            if(npairs[2] % 2 == 0)
            {
              npairs[2] /= 2;
              isOdd[2] = false;
            }
            else
            {
              if(npairs[2] == 1)
              {
                npairs[2] = 0;
              }
              else
              {
                npairs[2] = (npairs[2] + 1) / 2;
              }
              isOdd[2] = true;
            }
            isOdd[0] = false;

            // Shuffle particles to have random pairs
            //    (It does not really exchange them, it is just a temporary re-indexing)


            index1.resize(npart1);
            for (unsigned int i=0; i<npart1; i++)
            {
              index1[i] = i;
            }
            random_shuffle(index1.begin(), index1.end()); // shuffle the index array
            if(isOdd[1])
            {
              unsigned int index_last = index1.size() - 1;
              index1.push_back(index_last);
            }


            index2.resize(npart2);
            for (unsigned int i=0; i<npart2; i++)
            {
              index2[i] = i;
            }
            random_shuffle(index2.begin(), index2.end());
            if(isOdd[2])
            {
              unsigned int index_last = index2.size() - 1;
              index2.push_back(index_last);
            }

            intra[0] = false;
            intra[1] = true;
            intra[2] = true;


            n_collisionPairs = 3;
            for(int n_cP = 0; n_cP < n_collisionPairs; n_cP++)
            {
              iS1 = iSpec1[n_cP];
              iS2 = iSpec2[n_cP];
              s1 = vecSpecies[iS1];            s2 = vecSpecies[iS2];
              p1 = &(s1->particles);              p2 = &(s2->particles);
              m1 = s1->species_param.mass;        m2 = s2->species_param.mass;
              charge1 = s1->species_param.charge; charge2 = s2->species_param.charge;

              // used in equation (104a) and (104b)
              mr1 = m1 / (m1 + m2);
              mr2 = m2 / (m1 + m2);
              // the reduced mass
              m12 = m1 * m2 / (m1 + m2);
              // used in equation (96)
              gamma1 = 4.0 * const_pi * const_ephi0 * m12 / abs(charge1 * charge2);
              // used in equation (95)
              gamma2 = 4.0 * const_pi / ( gamma1 * gamma1 );


              A12 = gamma2 * n * 15.6;
              //cout<<"coulomb: "<<itime<<endl;

              // Now start the real loop on pairs of particles
              // See equations in http://dx.doi.org/10.1063/1.4742167
              // ----------------------------------------------------
              for (unsigned int i=0; i<npairs[n_cP]; i++)
              {
                  if(intra[n_cP])
                  {
                    if(iS1 == 0)
                    {

                      i1 = index1[i + n_particle_begin[iS1]] + bmin1;
                      i2 = index1[i + n_particle_begin[iS2] + npairs[n_cP]] + bmin1;
                    }
                    else if(iS1 == 1)
                    {
                      i1 = index2[i + n_particle_begin[iS1]] + bmin2;
                      i2 = index2[i + n_particle_begin[iS2] + npairs[n_cP]] + bmin2;
                    }

                  }
                  else
                  {
                    i1 = index1[i + n_particle_begin[iS1]] + bmin1;
                    i2 = index2[i + n_particle_begin[iS2]] + bmin2;
                  }


                  // Calculate stuff
                  gx = p1->momentum(0,i1) - p2->momentum(0,i2);
                  gy = p1->momentum(1,i1) - p2->momentum(1,i2);
                  gz = p1->momentum(2,i1) - p2->momentum(2,i2);
                  g_square = gx*gx + gy*gy + gz*gz;
                  g_magnitude = sqrt( g_square );
                  g_3 = g_square * g_magnitude;
                  g_p = sqrt( gy*gy + gz*gz );
                  if(g_p == 0.0 || g_3 == 0.0)
                  {
                      continue;
                  }
                  // the formula below the equation (101)
                  s = A12 * time_coulomb / g_3;
                  if(isOdd[n_cP])
                  {
                    int n_temp = npairs[n_cP] * 2 -1;
                    s *= ( n_temp / (n_temp+1) );
                  }

                  // Pick the deflection angles according to Nanbu's theory
                  // ref: improved modeing of relativistic collisions and collisional ionization in paritcle in cell codes
                  cosX = cos_chi(s);
                  if(cosX < -1.0)
                  {
                    cosX = -1.0;
                    cout<<"cosX < -1.0"<<endl;
                  }
                  else if(cosX > 1.0)
                  {
                    cosX = 1.0;
                    cout<<"cosX > 1.0"<<endl;
                  }

                  if(cosX*cosX > 1.0)
                  {
                    cout<<"cosX * cosX > 1.0"<<endl;
                    continue;
                  }
                  sinX = sqrt( 1. - cosX*cosX );


                  //!\todo make a faster rand by preallocating ??
                  phi = twoPi * ((double)rand() / RAND_MAX);
                  hx = g_p * cos(phi);
                  hy = -( gx * gy * cos(phi) + g_magnitude * gz * sin(phi) ) / g_p;
                  hz = -( gx * gz * cos(phi) - g_magnitude * gy * sin(phi) ) / g_p;

                  // Apply the deflection

                  double momentum1[3];
                  double momentum2[3];

                  momentum1[0] = p1->momentum(0,i1) - mr2 * ( gx * (1-cosX) + hx * sinX );
                  momentum1[1] = p1->momentum(1,i1) - mr2 * ( gy * (1-cosX) + hy * sinX );
                  momentum1[2] = p1->momentum(2,i1) - mr2 * ( gz * (1-cosX) + hz * sinX );

                  momentum2[0] = p2->momentum(0,i2) + mr1 * ( gx * (1-cosX) + hx * sinX );
                  momentum2[1] = p2->momentum(1,i2) + mr1 * ( gy * (1-cosX) + hy * sinX );
                  momentum2[2] = p2->momentum(2,i2) + mr1 * ( gz * (1-cosX) + hz * sinX );

                  p1->momentum(0,i1) = momentum1[0];
                  p1->momentum(1,i1) = momentum1[1];
                  p1->momentum(2,i1) = momentum1[2];

                  p2->momentum(0,i2) = momentum2[0];
                  p2->momentum(1,i2) = momentum2[1];
                  p2->momentum(2,i2) = momentum2[2];

              } // end loop on pairs of particles
              n_particle_begin[iS1] += npairs[n_cP];
              n_particle_begin[iS2] += npairs[n_cP];

              //cout<<n_collisionPairs<<"  "<<n_cP<<endl;
            }
            //cout<<nbins<<"  "<<ibin<<endl;

        } // end loop on bins
    }

}

*/


/*
// non-relativistic, same weight case: Nanbu theory
// the code corresponds to the ref: Probability theory of electron-molecule ion-molecule molecule-molecule and
// coulomb collisions for particle modeling of materials processing plasmas and gases
void Collisions1D_Coulomb::collide(PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies, int itime)
{

    unsigned int nbins = vecSpecies[0]->bmin.size(); // number of bins
    vector<unsigned int> *sg1, *sg2, *sgtmp, index1, index2;
    unsigned int bmin1, bmax1, bmin2, bmax2;
    unsigned int nspec1, nspec2; // numbers of species in each group
    unsigned int npart1, npart2; // numbers of macro-particles in each group
    unsigned int npairs; // number of pairs of macro-particles
    vector<unsigned int> np1, np2; // numbers of macro-particles in each species, in each group
    double n1, n2, n; // densities of particles
    unsigned int i1, i2, ispec1, ispec2, iSpec1, iSpec2;
    Species   *s1, *s2;
    Particles *p1, *p2;
    double m1, m2, m12, mr1, mr2, W1, W2, gamma1, gamma2, gx, gy, gz, g_magnitude, g_square, g_3, g_p,
           g12_square, hx, hy, hz, s, T1, T2, V1, V2, A12,
           cosX, sinX, phi, twoPi, debye_length, debye_squared, time_coulomb;
    double e_ov_ephi0;
    double charge1, charge2, charge, Vx, Vy, Vz, ni, Ti, Vi;


    if( (itime % params.timesteps_coulomb) == 0 )
    {
        //!!!=======species_group1 must contain one species like e, D1 and so on;
        // as same for species_group2;
        // species_group1 can be the same as species_group2
        sg1 = &species_group1;
        sg2 = &species_group2;


        // Loop on bins
        for (unsigned int ibin=0 ; ibin<nbins ; ibin++) {

            // get bin start/end for all necessary species, and number of particles
            for (unsigned int i=0; i<2; i++)
            {
                // try twice to ensure group 1 has more macro-particles
                s1 = vecSpecies[(*sg1)[0]];         s2 = vecSpecies[(*sg2)[0]];
                npart1 = s1->bmax[ibin] - s1->bmin[ibin];
                npart2 = s2->bmax[ibin] - s2->bmin[ibin];
                bmin1 = s1->bmin[ibin];
                bmin2 = s2->bmin[ibin];

                if (npart2 <= npart1)
                {
                    break; // ok if group1 has more macro-particles
                }
                else { // otherwise, we exchange groups and try again
                    sgtmp = sg1; sg1 = sg2; sg2 = sgtmp;
                }
            }
            // now group1 has more macro-particles than group2

            // skip to next bin if no particles
            if (npart1==0 || npart2==0)
            {
                continue;
            }

            iSpec1 = (*sg1)[0];
            iSpec2 = (*sg2)[0];

            p1 = &(s1->particles);              p2 = &(s2->particles);
            m1 = s1->species_param.mass;        m2 = s2->species_param.mass;
            W1 = s1->species_param.weight;      W2 = s2->species_param.weight;
            charge1 = s1->species_param.charge; charge2 = s2->species_param.charge;

            // Shuffle particles to have random pairs
            //    (It does not really exchange them, it is just a temporary re-indexing)
            index1.resize(npart1);
            for (unsigned int i=0; i<npart1; i++) index1[i] = i; // first, we make an ordered array
            //! \todo benchmark and improve the shuffling method ?
            random_shuffle(index1.begin(), index1.end()); // shuffle the index array
            if (intra_collisions) { // In the case of collisions within one species
                n1 = npart1 * W1;
                n2 = 0.0;
                npairs = (int) ceil(((double)npart1)/2.); // half as many pairs as macro-particles
                index2.resize(npairs);
                for (unsigned int i=0; i<npairs; i++) index2[i] = index1[i+npart1-npairs]; // index2 is second half
                index1.resize(npairs); // index1 is first half
            }
            else
            { // In the case of collisions between two species
                n1 = npart1 * W1;
                n2 = npart2 * W2;
                npairs = npart1; // as many pairs as macro-particles in group 1 (most numerous)
                index2.resize(npairs);
                for (unsigned int i=0; i<npart1; i++) index2[i] = i % npart2;
            }
            n = n1 + n2;

            // =========Calculate some constants in the formulas used below============
            twoPi = 2.0 * const_pi;
            e_ov_ephi0 = const_e / const_ephi0;
            time_coulomb = params.timesteps_coulomb * params.timestep * params.zoom_collision;
            // used in equation (104a) and (104b)
            mr1 = m1 / (m1 + m2);
            mr2 = m2 / (m1 + m2);
            // the reduced mass
            m12 = m1 * m2 / (m1 + m2);
            // used in equation (96)
            gamma1 = 4.0 * const_pi * const_ephi0 * m12 / abs(charge1 * charge2);
            // used in equation (95)
            gamma2 = 4.0 * const_pi / ( gamma1 * gamma1 );


            // Pre-calculate some numbers before the big loop
            Field1D* T1_s = static_cast<Field1D*>(fields->T_s[iSpec1]);
            Field1D* T2_s = static_cast<Field1D*>(fields->T_s[iSpec2]);
            Field1D* Vx1_s = static_cast<Field1D*>(fields->Vx_s[iSpec1]);
            Field1D* Vy1_s = static_cast<Field1D*>(fields->Vy_s[iSpec1]);
            Field1D* Vz1_s = static_cast<Field1D*>(fields->Vz_s[iSpec1]);

            Field1D* Vx2_s = static_cast<Field1D*>(fields->Vx_s[iSpec2]);
            Field1D* Vy2_s = static_cast<Field1D*>(fields->Vy_s[iSpec2]);
            Field1D* Vz2_s = static_cast<Field1D*>(fields->Vz_s[iSpec2]);

            if( ibin+oversize[0] < 0 || ibin+oversize[0] >= T1_s->dims_[0] )
            {
                cout<<"error: "<<ibin+oversize[0]<<endl;
            }
            T1 = (*T1_s)(ibin + oversize[0]) * 3.0 * const_e / m1;
            T2 = (*T2_s)(ibin + oversize[0]) * 3.0 * const_e / m2;
            Vx = (*Vx1_s)(ibin + oversize[0]) - (*Vx2_s)(ibin + oversize[0]);
            Vy = (*Vy1_s)(ibin + oversize[0]) - (*Vy2_s)(ibin + oversize[0]);
            Vz = (*Vz1_s)(ibin + oversize[0]) - (*Vz2_s)(ibin + oversize[0]);
            // the formula below equation (96)
            g12_square = T1 + T2 + Vx*Vx + Vy*Vy + Vz*Vz;

            // Calculate debye_length
            debye_squared = 0.0;
            for( int iSpec = 0; iSpec < vecSpecies.size(); iSpec++ )
            {
                Field1D* T1D_s = static_cast<Field1D*>(fields->T_s[iSpec]);
                Field1D* rho1D_s = static_cast<Field1D*>(fields->rho_s[iSpec]);
                charge = vecSpecies[iSpec]->species_param.charge;
                ni = (*rho1D_s)(ibin + oversize[0]);
                Ti =   (*T1D_s)(ibin + oversize[0]);
                if( Ti > 0.0 )
                {
                    debye_squared += ( ni * charge * charge / ( const_ephi0 * const_e * Ti ) );
                }
                else
                {
                    continue;
                }

            }
            if( debye_squared > 0.0 )
            {
                debye_length = sqrt( 1.0 / debye_squared );
            }
            else
            {
                continue;
            }

            // equation (95)
            if( g12_square > 0.0 )
            {
                A12 = gamma2 * n * log( gamma1 * g12_square * debye_length );
                if(A12 < 0.0)
                {
                    cout<<"Coulomb collision: A12 "<<A12<<"  debye_length = "<<debye_length<<endl;
                    continue;
                }
            }
            else
            {
                continue;
            }

            //cout<<"coulomb: "<<itime<<endl;

            // Now start the real loop on pairs of particles
            // See equations in http://dx.doi.org/10.1063/1.4742167
            // ----------------------------------------------------
            for (unsigned int i=0; i<npairs; i++) {

                // find species and index i1 of particle "1"
                i1 = index1[i] + bmin1;
                // find species and index i2 of particle "2"
                i2 = index2[i] + bmin2;

                // Calculate stuff
                gx = p1->momentum(0,i1) - p2->momentum(0,i2);
                gy = p1->momentum(1,i1) - p2->momentum(1,i2);
                gz = p1->momentum(2,i1) - p2->momentum(2,i2);
                g_square = gx*gx + gy*gy + gz*gz;
                g_magnitude = sqrt( g_square );
                g_3 = g_square * g_magnitude;
                g_p = sqrt( gy*gy + gz*gz );
                if(g_p == 0.0 || g_3 == 0.0)
                {
                    continue;
                }
                // the formula below the equation (101)
                s = A12 * time_coulomb / g_3;

                // Pick the deflection angles according to Nanbu's theory
                // ref: improved modeing of relativistic collisions and collisional ionization in paritcle in cell codes
                cosX = cos_chi(s);
                sinX = sqrt( 1. - cosX*cosX );
                while( cosX*cosX > 1.0 || isnan(sinX) )
                {
                    cout<<"Coulomb1D collision: cosX, sinX = "<<cosX<<" "<<sinX<<endl;
                    cosX = cos_chi(s);
                    sinX = sqrt( 1. - cosX*cosX );
                }

                //!\todo make a faster rand by preallocating ??
                phi = twoPi * ((double)rand() / RAND_MAX);
                hx = g_p * cos(phi);
                hy = -( gx * gy * cos(phi) + g_magnitude * gz * sin(phi) ) / g_p;
                hz = -( gx * gz * cos(phi) - g_magnitude * gy * sin(phi) ) / g_p;

                // Apply the deflection

                vector<double> momentum1;
                vector<double> momentum2;
                momentum1.resize(3);
                momentum2.resize(3);

                momentum1[0] = p1->momentum(0,i1) - mr2 * ( gx * (1-cosX) + hx * sinX );
                momentum1[1] = p1->momentum(1,i1) - mr2 * ( gy * (1-cosX) + hy * sinX );
                momentum1[2] = p1->momentum(2,i1) - mr2 * ( gz * (1-cosX) + hz * sinX );

                momentum2[0] = p2->momentum(0,i2) + mr1 * ( gx * (1-cosX) + hx * sinX );
                momentum2[1] = p2->momentum(1,i2) + mr1 * ( gy * (1-cosX) + hy * sinX );
                momentum2[2] = p2->momentum(2,i2) + mr1 * ( gz * (1-cosX) + hz * sinX );

                DEBUGEXEC(
                    int error = 0;
                    for(int idirection = 0; idirection < 3; idirection++)
                    {
                        if( isnan(momentum1[idirection]) || isinf(momentum1[idirection]) )
                        {
                            cout<<"Coulomb Error info: "<<s<<" "<<cosX<<" "<<sinX<<" "<<g_p<<" "<<g_3<<endl;
                            cout<<"Species1: "<<s1->species_param.species_type<<" momentum old "<<p1->momentum(idirection,i1)<<endl;
                            cout<<"momentum new = "<<momentum1[idirection]<<endl;
                            error = 1;
                        }

                        if( isnan(momentum2[idirection]) || isinf(momentum2[idirection]) )
                        {
                            cout<<"Coulomb Error info: "<<s<<" "<<cosX<<" "<<sinX<<" "<<g_p<<" "<<g_3<<endl;
                            cout<<"Species2: "<<s2->species_param.species_type<<" momentum old "<<p2->momentum(idirection,i2)<<endl;
                            cout<<"momentum new = "<<momentum2[idirection]<<endl;
                            error = 1;
                        }
                    }
                    if(error == 1)
                    {
                        continue;
                    }
                );

                p1->momentum(0,i1) = momentum1[0];
                p1->momentum(1,i1) = momentum1[1];
                p1->momentum(2,i1) = momentum1[2];

                p2->momentum(0,i2) = momentum2[0];
                p2->momentum(1,i2) = momentum2[1];
                p2->momentum(2,i2) = momentum2[2];



            } // end loop on pairs of particles

        } // end loop on bins
    }

}

*/




// Calculates the debye length squared in each cluster
// The formula for the inverse debye length squared is sumOverSpecies(density*charge^2/(ephi0*temperature))
void Collisions1D_Coulomb::calculate_debye_length(PicParams& params, vector<Species*>& vecSpecies)
{

    // get info on particle binning
    unsigned int nbins = vecSpecies[0]->bmin.size(); // number of bins
    unsigned int nspec = vecSpecies.size(); // number of species
    unsigned int bmin, bmax;
    double p2, density, density_max, charge, temperature, rmin2;
    Species   * s;
    Particles * p;
    double coeff = params.wavelength_SI/(6.*M_PI*2.8179403267e-15); // normLength/(3*electronRadius) = wavelength/(6*pi*electronRadius)

    debye_length_squared.resize(nbins);

    // Loop on bins
    //! \todo Make OpenMP parallelization (MF & JD)
    for (unsigned int ibin=0 ; ibin<nbins ; ibin++) {

        density_max = 0.;
        debye_length_squared[ibin] = 0.;
        for (unsigned int ispec=0 ; ispec<nspec ; ispec++) { // loop all species
            // Calculation of particles density, mean charge, and temperature
            // Density is the sum of weights
            // Temperature basic definition is the average <v*p> divided by 3
            //    (instead of v*p, we use p^2/gamma)
            s  = vecSpecies[ispec];
            p  = &(s->particles);
            bmin = s->bmin[ibin];
            bmax = s->bmax[ibin];
            density     = 0.;
            charge      = 0.;
            temperature = 0.;
            // loop particles to calculate average quantities
            for (unsigned int iPart=bmin ; iPart<bmax ; iPart++ ) {
                p2 = pow(p->momentum(0,iPart),2)+pow(p->momentum(1,iPart),2)+pow(p->momentum(2,iPart),2);
                density     += p->weight(iPart);
                charge      += p->weight(iPart) * p->charge(iPart);
                temperature += p->weight(iPart) * p2/sqrt(1.+p2);
            }
            if (density <= 0.) continue;
            charge /= density; // average charge
            temperature *= (s->species_param.mass) / (3.*density); // Te in units of eV
            density /= params.n_cell_per_cluster; // density in units of critical density
            // compute inverse debye length squared
            if (temperature>0.) debye_length_squared[ibin] += density*charge*charge/(const_ephi0*temperature);
            // compute maximum density of species
            if (density>density_max) density_max = density;
        }

        // if there were particles,
        if (debye_length_squared[ibin] > 0.) {
            // compute debye length squared in code units
            debye_length_squared[ibin] = 1./(debye_length_squared[ibin]);
            // apply lower limit to the debye length (minimum interatomic distance)
            // I do not know the meanming of below code, so I comment them ......
            //rmin2 = pow(coeff*density_max, -2./3.);
            //if (debye_length_squared[ibin] < rmin2) debye_length_squared[ibin] = rmin2;
        }

      }

  #ifdef  __DEBUG
      // calculate and print average debye length
      double mean_debye_length = 0.;
      for (unsigned int ibin=0 ; ibin<nbins ; ibin++)
          mean_debye_length += sqrt(debye_length_squared[ibin]);
      mean_debye_length /= (double)nbins;
      //DEBUG("Mean Debye length in code length units = " << scientific << setprecision(3) << mean_debye_length);
      mean_debye_length *= params.wavelength_SI/(2.*M_PI); // switch to SI
      DEBUG("Mean Debye length in meters = " << scientific << setprecision(3) << mean_debye_length );
  #endif

}


// relativistic case
// the code corresponds to the ref: improved modeing of relativistic collisions and collisional ionization in paritcle in cell codes
void Collisions1D_Coulomb::collide_relativistic(PicParams& params, SmileiMPI* smpi, vector<Species*>& vecSpecies, int itime)
{

    unsigned int nbins = vecSpecies[0]->bmin.size(); // number of bins
    vector<unsigned int> *sg1, *sg2, *sgtmp, bmin1, bmax1, bmin2, bmax2, index1, index2;
    unsigned int nspec1, nspec2; // numbers of species in each group
    unsigned int npart1, npart2; // numbers of macro-particles in each group
    unsigned int npairs; // number of pairs of macro-particles
    vector<unsigned int> np1, np2; // numbers of macro-particles in each species, in each group
    double n1, n2, n12, n123, n223; // densities of particles
    unsigned int i1, i2, ispec1, ispec2;
    Species   *s1, *s2;
    Particles *p1, *p2;
    double m1, m2, m12, W1, W2, qqm, qqm2, qq2mm, qq12, gamma1, gamma2, gamma_m12, gamma_m12_inv,
           COM_vx, COM_vy, COM_vz, COM_vsquare, COM_gamma,
           term1, term2, term3, term4, term5, term6, coeff1, coeff2, coeff3, coeff4, coeff5, twoPi,
           vcv1, vcv2, px_COM, py_COM, pz_COM, p2_COM, p_COM, gamma1_COM, gamma2_COM,
           logL, bmin, s, vrel, smax,
           cosX, sinX, phi, sinXcosPhi, sinXsinPhi, p_perp, inv_p_perp,
           newpx_COM, newpy_COM, newpz_COM, U, vcp;
    Field2D *smean, *logLmean, *ncol;//, *temperature
    ostringstream name;
    hid_t did;

    sg1 = &species_group1;
    sg2 = &species_group2;




    // Loop on bins
    for (unsigned int ibin=0 ; ibin<nbins ; ibin++) {

        // get bin start/end for all necessary species, and number of particles
        for (unsigned int i=0; i<2; i++) { // try twice to ensure group 1 has more macro-particles
            nspec1 = sg1->size();
            nspec2 = sg2->size();
            bmin1.resize(nspec1); // bin starting point, for each of species group 1
            bmax1.resize(nspec1); // bin  ending  point, for each of species group 1
            np1  .resize(nspec1); // number of particles in that bin
            bmin2.resize(nspec2); // bin starting point, for each of species group 2
            bmax2.resize(nspec2); // bin  ending  point, for each of species group 2
            np2  .resize(nspec2); // number of particles in that bin
            npart1 = 0; npart2 = 0;
            for (ispec1=0 ; ispec1<nspec1 ; ispec1++) {
                s1 = vecSpecies[(*sg1)[ispec1]];
                bmin1[ispec1] = s1->bmin[ibin];
                bmax1[ispec1] = s1->bmax[ibin];
                np1[ispec1] = bmax1[ispec1] - bmin1[ispec1];
                npart1 += np1[ispec1];
            }
            for (ispec2=0 ; ispec2<nspec2 ; ispec2++) {
                s2 = vecSpecies[(*sg2)[ispec2]];
                bmin2[ispec2] = s2->bmin[ibin];
                bmax2[ispec2] = s2->bmax[ibin];
                np2[ispec2] = bmax2[ispec2] - bmin2[ispec2];
                npart2 += np2[ispec2];
            }
            if (npart2 <= npart1) break; // ok if group1 has more macro-particles
            else { // otherwise, we exchange groups and try again
                sgtmp = sg1; sg1 = sg2; sg2 = sgtmp;
            }
        }
        // now group1 has more macro-particles than group2

        // skip to next bin if no particles
        if (npart1==0 or npart2==0) continue;

        // Shuffle particles to have random pairs
        //    (It does not really exchange them, it is just a temporary re-indexing)
        index1.resize(npart1);
        for (unsigned int i=0; i<npart1; i++) index1[i] = i; // first, we make an ordered array
        //! \todo benchmark and improve the shuffling method ?
        random_shuffle(index1.begin(), index1.end()); // shuffle the index array
        if (intra_collisions) { // In the case of collisions within one species
            npairs = (int) ceil(((double)npart1)/2.); // half as many pairs as macro-particles
            index2.resize(npairs);
            for (unsigned int i=0; i<npairs; i++) index2[i] = index1[i+npart1-npairs]; // index2 is second half
            index1.resize(npairs); // index2 is first half
        } else { // In the case of collisions between two species
            npairs = npart1; // as many pairs as macro-particles in group 1 (most numerous)
            index2.resize(npairs);
            for (unsigned int i=0; i<npart1; i++) index2[i] = i % npart2;
        }

        // Calculate density of group 1
        n1 = 0.;
        for (ispec1=0 ; ispec1<nspec1 ; ispec1++)
            for (unsigned int iPart=bmin1[ispec1] ; iPart<bmax1[ispec1] ; iPart++)
                n1 += vecSpecies[(*sg1)[ispec1]]->particles.weight(iPart);
        n1 /= params.n_cell_per_cluster;

        // Calculate density of group 2
        n2 = 0.;
        for (ispec2=0 ; ispec2<nspec2 ; ispec2++)
            for (unsigned int iPart=bmin2[ispec2] ; iPart<bmax2[ispec2] ; iPart++)
                n2 += vecSpecies[(*sg2)[ispec2]]->particles.weight(iPart);
        n2 /= params.n_cell_per_cluster;

        // Calculate the "hybrid" density
        n12 = 0.;
        for (unsigned int i=0; i<npairs; i++) { // for each pair of particles
            // find species and index i1 of particle "1"
            i1 = index1[i];
            for (ispec1=0 ; ispec1<nspec1 ; ispec1++) {
                if (i1 < np1[ispec1]) break;
                i1 -= np1[ispec1];
            }
            i1 += bmin1[ispec1];
            // find species and index i2 of particle "2"
            i2 = index2[i];
            for (ispec2=0 ; ispec2<nspec2 ; ispec2++) {
                if (i2 < np2[ispec2]) break;
                i2 -= np2[ispec2];
            }
            i2 += bmin2[ispec2];
            // sum weights
            n12 += min( vecSpecies[(*sg1)[ispec1]]->particles.weight(i1)
                       ,vecSpecies[(*sg2)[ispec2]]->particles.weight(i2) );
        }
        n12 /= params.n_cell_per_cluster;

        // Pre-calculate some numbers before the big loop
        n123 = pow(n1,2./3.); n223 = pow(n2,2./3.);
        twoPi = 2. * M_PI;
        coeff1 = 0.5*const_h; // h/2
        coeff2 = 1.0 / (4.0*const_pi*const_ephi0*const_c*const_c); //
        coeff3 = params.timestep * n1*n2/n12;
        coeff4 = pow( 4.0*const_pi/3.0 , 1./3. ) * coeff3;
        coeff5 = 1.0 / (4.0 * const_pi * const_ephi0 * pow(const_c,4));     //paraters in equation (17)




        // Now start the real loop on pairs of particles
        // See equations in http://dx.doi.org/10.1063/1.4742167
        // ----------------------------------------------------
        for (unsigned int i=0; i<npairs; i++) {

            // find species and index i1 of particle "1"
            i1 = index1[i];
            for (ispec1=0 ; ispec1<nspec1 ; ispec1++) {
                if (i1 < np1[ispec1]) break;
                i1 -= np1[ispec1];
            }
            i1 += bmin1[ispec1];
            // find species and index i2 of particle "2"
            i2 = index2[i];
            for (ispec2=0 ; ispec2<nspec2 ; ispec2++) {
                if (i2 < np2[ispec2]) break;
                i2 -= np2[ispec2];
            }
            i2 += bmin2[ispec2];

            s1 = vecSpecies[(*sg1)[ispec1]]; s2 = vecSpecies[(*sg2)[ispec2]];
            p1 = &(s1->particles);           p2 = &(s2->particles);
            m1 = s1->species_param.mass;     m2 = s2->species_param.mass;
            W1 = p1->weight(i1);             W2 = p2->weight(i2);

            // Calculate stuff
            m12  = m1 / m2; // mass ratio
            qqm  = p1->charge(i1) * p2->charge(i2) / m1;
            //qqm2 = qqm * qqm;
            qq12 = p1->charge(i1) * p2->charge(i2);
            qq2mm = qq12*qq12 / (m1*m2);


            // Get momenta and calculate gammas
            double const_c2 = const_c * const_c;
            double const_c2_inv = 1.0 / const_c2;
            gamma1 = 1.0 / sqrt(1. - (pow(p1->momentum(0,i1),2) + pow(p1->momentum(1,i1),2) + pow(p1->momentum(2,i1),2))*const_c2_inv);
            gamma2 = 1.0 / sqrt(1. + (pow(p2->momentum(0,i2),2) + pow(p2->momentum(1,i2),2) + pow(p2->momentum(2,i2),2))*const_c2_inv);
            gamma_m12 = gamma1*m1 + gamma2*m2;
            gamma_m12_inv = 1./gamma_m12;

            // Calculate the center-of-mass (COM) frame ==> equation (1)
            // Quantities starting with "COM" are those of the COM itself, expressed in the lab frame.
            // They are NOT quantities relative to the COM.
            // COM_v is Vc in equation (1)
            COM_vx = ( m1 * (p1->momentum(0,i1)) + m2 * (p2->momentum(0,i2)) ) * gamma_m12_inv;
            COM_vy = ( m1 * (p1->momentum(1,i1)) + m2 * (p2->momentum(1,i2)) ) * gamma_m12_inv;
            COM_vz = ( m1 * (p1->momentum(2,i1)) + m2 * (p2->momentum(2,i2)) ) * gamma_m12_inv;
            COM_vsquare = COM_vx*COM_vx + COM_vy*COM_vy + COM_vz*COM_vz;
            COM_gamma = pow( 1.-COM_vsquare*const_c2_inv , -0.5);

            // Change the momentum to the COM frame (we work only on particle 1)
            // Quantities ending with "COM" are quantities of the particle expressed in the COM frame.
            //some terms in equation (2)
            term1 = (COM_gamma - 1.) / COM_vsquare;
            vcv1  = (COM_vx*(p1->momentum(0,i1)) + COM_vy*(p1->momentum(1,i1)) + COM_vz*(p1->momentum(2,i1)))/gamma1;
            vcv2  = (COM_vx*(p2->momentum(0,i2)) + COM_vy*(p2->momentum(1,i2)) + COM_vz*(p2->momentum(2,i2)))/gamma2;
            term2 = (term1*vcv1 - COM_gamma) * m1 * gamma1;
            px_COM = (p1->momentum(0,i1)) + term2*COM_vx;
            py_COM = (p1->momentum(1,i1)) + term2*COM_vy;
            pz_COM = (p1->momentum(2,i1)) + term2*COM_vz;
            p2_COM = px_COM*px_COM + py_COM*py_COM + pz_COM*pz_COM;
            p_COM  = sqrt(p2_COM);
            gamma1_COM = (1.-vcv1*const_c2_inv)*COM_gamma*gamma1;
            gamma2_COM = (1.-vcv2*const_c2_inv)*COM_gamma*gamma2;

            // Calculate some intermediate quantities in equation (17)
            term3 = COM_gamma * gamma_m12_inv;        // before p1*
            term4 = gamma1_COM * gamma2_COM;
            term5 = term4*m1*m2*const_c2/p2_COM + 1.0;

            // Calculate coulomb log if necessary
            // equation (22) and (23)
            logL = coulomb_log;
            if( logL <= 0. ) { // if auto-calculation requested
                bmin = max( coeff1/p_COM , abs(coeff2*qq12*term3*term5) ); // min impact parameter
                logL = 0.5*log(1.+debye_length_squared[ibin]/pow(bmin,2));
                if (logL < 2.) logL = 2.;
            }

            // Calculate the collision parameter s12 (similar to number of real collisions)
            s = coeff3 * coeff5 * logL * qq2mm * term3 * p_COM * term5*term5 / (gamma1*gamma2);

            // Low-temperature correction       // equation (21)
            // !!!not sure vrel is the equation (8), or vrel is the relative velocity in laboratory frame
            vrel = p_COM/ (m1*m2*term3*term4); // relative velocity
            smax = coeff4 * (m1+m2) * vrel / max(m12*n123,n223);
            if (s>smax) s = smax;

            // Pick the deflection angles according to Nanbu's theory
            cosX = cos_chi(s);
            sinX = sqrt( 1. - cosX*cosX );
            //!\todo make a faster rand by preallocating ??
            phi = twoPi * ((double)rand() / RAND_MAX);

            // Calculate combination of angles
            sinXcosPhi = sinX*cos(phi);
            sinXsinPhi = sinX*sin(phi);

            // Apply the deflection
            p_perp = sqrt( px_COM*px_COM + py_COM*py_COM );
            if( p_perp > 1.e-10*p_COM ) { // make sure p_perp is not too small
                inv_p_perp = 1./p_perp;
                newpx_COM = (px_COM * pz_COM * sinXcosPhi - py_COM * p_COM * sinXsinPhi) * inv_p_perp + px_COM * cosX;
                newpy_COM = (py_COM * pz_COM * sinXcosPhi + px_COM * p_COM * sinXsinPhi) * inv_p_perp + py_COM * cosX;
                newpz_COM = -p_perp * sinXcosPhi  +  pz_COM * cosX;
            } else { // if p_perp is too small, we use the limit px->0, py=0
                newpx_COM = p_COM * sinXcosPhi;
                newpy_COM = p_COM * sinXsinPhi;
                newpz_COM = p_COM * cosX;
            }

            // Random number to choose whether deflection actually applies.
            // This is to conserve energy in average when weights are not equal.
            //!\todo make a faster rand by preallocating ??
            U = ((double)rand() / RAND_MAX);

            // Go back to the lab frame and store the results in the particle array
            vcp = COM_vx * newpx_COM + COM_vy * newpy_COM + COM_vz * newpz_COM;
            if( U < W2/W1 ) { // deflect particle 1 only with some probability
                term6 = term1*vcp + gamma1_COM * COM_gamma;
                p1->momentum(0,i1) = newpx_COM + COM_vx * term6;
                p1->momentum(1,i1) = newpy_COM + COM_vy * term6;
                p1->momentum(2,i1) = newpz_COM + COM_vz * term6;
            }
            if( U < W1/W2 ) { // deflect particle 2 only with some probability
                term6 = -m12 * term1*vcp + gamma2_COM * COM_gamma;
                p2->momentum(0,i2) = -m12 * newpx_COM + COM_vx * term6;
                p2->momentum(1,i2) = -m12 * newpy_COM + COM_vy * term6;
                p2->momentum(2,i2) = -m12 * newpz_COM + COM_vz * term6;
            }



        } // end loop on pairs of particles


    } // end loop on bins


}


// Technique given by Nanbu in http://dx.doi.org/10.1103/PhysRevE.55.4642
//   to pick randomly the deflection angle cosine, in the center-of-mass frame.
// It involves the "s" parameter (~ collision frequency * deflection expectation)
//   and a random number "U".
// Technique slightly modified in http://dx.doi.org/10.1063/1.4742167
// ref: improved modeing of relativistic collisions and collisional ionization in paritcle in cell codes
inline double Collisions1D_Coulomb::cos_chi(double s)
{

    double A, invA;
    //!\todo make a faster rand by preallocating ??
    double U = (double)rand() / RAND_MAX;

    if( s < 0.1 ) {
        //if ( U<0.0001 ) U=0.0001; // ensures cos_chi > 0
        while(U == 0 || 1. + s*log(U) <= 0.0)
        {
          U = (double)rand() / RAND_MAX;
        }
        return 1. + s*log(U);
    }
    if( s < 3.  ) {
        // the polynomial has been modified from the article in order to have a better form
        invA = 0.00569578 +(0.95602 + (-0.508139 + (0.479139 + ( -0.12789 + 0.0238957*s )*s )*s )*s )*s;
        A = 1./invA;
        return  invA  * log( exp(-A) + 2.*U*sinh(A) );
    }
    if( s < 6.  ) {
        A = 3.*exp(-s);
        return (1./A) * log( exp(-A) + 2.*U*sinh(A) );
    }
    return 2.*U - 1.;

}
