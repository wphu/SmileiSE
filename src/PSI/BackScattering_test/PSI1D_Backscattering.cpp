#include <cmath>
#include <iostream>
#include "PSI1D_Backscattering.h"

using namespace std;


// Constructor
PSI1D_Backscattering::PSI1D_Backscattering(
    int nz1_in,
    int m1_in,
    int ne_in,
    vector<int> nz2_in,
    vector<int> nw_in )
{
    nz1     = nz1_in;
    m1      = m1_in;
    ne      = ne_in;
    nz2     = nz2_in;
    nw      = nw_in;
}

PSI1D_Backscattering::~PSI1D_Backscattering()
{

}

void PSI1D_Backscattering::scatter(double &rnion, double &reion, double theta, double energy)
{
    double a1, r, th, feps, r0heavy, sklfact, psi, psiden2, psiden1, psinum, fmu, redRN;
    double delta, eps, sb, sa, a2, z2, c, fi0, epsb, beta2, sl, w, zz, z2p, aa, wi, sum;
    double sbeff, saeff, a2eff, z2eff, z1p, z1,eth, eth1, mu, tau;
    double r0, as1, as2;
    if(nz1 == 1)
    {
        if(m1 < 1 || m1 > 3)
        {
            cout<<"BackScattering Error:  1"<<endl;
            return;
        }
        a1 = a1t[ m1-1 ];
    }
    else if(nz1 == 2)
    {
        if(m1 < 3 || m1 > 4)
        {
            cout<<"BackScattering Error:  2"<<endl;
            return;
        }
        a1 = a1t[ m1 ];
    }
    else if(nz1 > 2)
    {
        a1 = a1t[nz1 + 2];
        a2eff = 0.0;
        sum = 0.0;
        for( int i = 0; i < ne; i++ )
        {
            a2 = a2t[ nz2[i] -1 ];
            a2eff += ( nw[i] * a2 );
            sum += nw[i];
        }
        a2 = a2eff / sum;
        delta = a1 / a2;
        if( delta > 4.5 )
        {
            cout<<"BackScattering Message:  3"<<endl;
        }
    }
    else
    {
        cout<<"BackScattering Error:  4"<<endl;
        return;
    }

    z1 = nz1;
    z1p = pow( z1, 2.0/3.0 );
    z2eff = 0.0;
    a2eff = 0.0;
    saeff = 0.0;
    sbeff = 0.0;
    sum = 0.0;

    for( int i = 0; i < ne; i++ )
    {
        z2 = nz2[i];
        a2 = a2t[ nz2[i] -1 ];
        wi = nw[i];
        z2eff = z2eff + wi * z2;
        a2eff = a2eff + wi * a2;
        sum = sum + wi;
        fmu = a1 / a2;
        aa = 1.0 + fmu;
        z2p = pow( z2, 2.0/3.0 );
        zz = sqrt( 1.0 + z1p / z2p );
        eps = fneps( energy,z1,a1,z2,a2 );
        sa = 0.0793 * sqrt(eps) / fmu;
        w = z1 * zz * aa / z2p;
        wi = wi * z1 / w / a2;
        saeff = saeff + wi * sa;
        sl = d[ nz2[i] -1 ] * z1p * pow( aa / zz, 1.5 ) * sa / sqrt(a1);

        if ( z2 - 12.9 <= 0.0 )
        {
            fi0 = 12.0 + 7.0 / z2;
        }
        else
        {
            fi0 = 9.76 + 58.5 / pow( z2, 1.19 );
        }

        tau = eps * pow( z2, 2.0 ) * w / a1 / 3.03055e7;
        beta2 = tau * ( tau + 2.0 ) / pow( tau + 1.0, 2 );
        epsb = 1.022e6 * beta2 / z2 / fi0;

        if ( z1 - 2.9 <= 0.0 )
        {
            c = 100.0 * z1 / z2;
        }
        else
        {
            c = 5.0;
        }

        sb = 61.474 / fmu * w * ( log ( epsb / (1.0 - beta2) + 1.0 + c / epsb ) - beta2 )
             / fi0 / epsb;
        sb = 2.0 * eps / rangen( eps,energy ) + sl * sb / (0.25 * sl + sb );
        sbeff = sbeff + wi * sb;
    }


    z2 = z2eff / sum;
    a2 = a2eff / sum;
    sa = saeff;
    sb = sbeff;
    eps = fneps( energy, z1, a1, z2, a2 );
    if( delta > 0.5 )
    {
        redRN = 1.0 / ( 1 + pow( eps/0.104, 0.577) + pow( eps/0.73, 1.5 ) );
        fmu = a1 / a2;
        psinum = 1.0 + 24.1 * pow( fmu, -3.995 );
        psiden1 = fmu / ( (1 + fmu) * (eps / 1.84 + 1.0) );
        psiden2 = pow( 1.0 + fmu, 3 ) * eps /
                  ( pow(fmu, 3) * (1.0 - 3.0/(2.0*fmu) + 0.9/pow(fmu,2) + 1.0/(2.0*pow(fmu,3)) ) * (eps+13.3) );
        psi = psinum / ( psiden1 + psiden2 );

        // sbeff/saeff corresponds to the rho_a/rho-t in the reference paper ( of course
        // sbeff is not equal to rho_a nor is saeff equal to rho_t).

        sklfact = pow(z1, 0.66666667) / sqrt(a1) * sb/sa * pow(fmu,2) / pow(1.0+fmu,2) / psi;
        r0heavy = redRN / sklfact;

        //write(*,*)'did heavy projectile calculation for delta > 0.5 '

        // =====================Calculate rnion=============================
        r = 1.0 / ( 1.0 + pow(eps/0.133, 0.285) ) + 0.530 / ( 1.0 + pow(85.0/eps, 1.46) );

        //write(*,*)'did heavy proj calculation for rnion'
        if (th != 0.0)
        {
            //write(*,*)'did theta not = 0 calc for rnion'
            r0 = r0heavy;
            as1 = 7.38 * pow(eps, 0.359);
            as2 = 0.836 / pow(eps, 0.087);
            rnion = rne(th, r0, as1, as2);
            //return
        }
        else
        {
            // write(*,*)'did theta = 0 calc. for rnion'
            rnion = r0heavy;
            // write(*,*)r0heavy, rnion, r, reion
            //return
        }

        // ===========Calculate reion================================
        r = 1.0 / ( 1.0 + pow(eps/0.133, 0.285) ) + 0.530 / ( 1.0 + pow(85.0/eps, 1.46) );
        th = theta * 1.745329e-2;
        if (th == 0.0)
        {
            // write(*,*)'did theta=0 calculation of reion'
            reion = r0heavy * r;
            //go to 79
        }
        else
        {
			// write(*,*)'did theta > 0 calculation for reion'
            r0 = r0heavy * r;
            as1 = 17.9 * pow(eps, 0.453);
            as2 = 0.771 / pow(eps, 0.014);
            reion = rne(th, r0, as1, as2);
            //go to 79
        }
    }
    else
    {
        feps = (2.0 * eps - 3.0 ) / ( eps + 1.0 );
        r0 = 0.705 * sa / sb / pow( ( 1.0 + a1/a2 ), feps );
        r0 = r0 / ( 1.0 + pow( (eps/0.047), 0.597 ) + pow( (eps/0.619), 1.5 ) );
        th = theta * 1.745329e-2;


        // ======================Calculate reion=======================
        if (th != 0.0)
        {
            as1 = 17.9 * pow(eps, 0.453);
            as2 = 0.771 / pow(eps, 0.014);
            reion = rne(th, r0, as1, as2);
        }
        else
        {
            reion = r0;
        }


        // ============Calculate rnion===============================
		// write(*,*)'did light projectile calculation for rnion'
        r = 1.0 / ( 1.0 + pow( (eps/0.133), 0.285 ) ) + 0.530 / ( 1.0 + pow(85.0/eps, 1.46) );
        r0 = r0 / r;
        if (th != 0.0)
        {
            as1 = 7.38 * pow(eps, 0.359);
            as2 = 0.836 / pow(eps, 0.087);
            rnion = rne(th, r0, as1, as2);
            //return
        }
        else
        {
            rnion = r0;
        }

        //return


    }

}



double PSI1D_Backscattering::rangen(double epsiln, double energy)
{
    // reduced range of ions for nuclear stopping only
    // W.D.Wilson , L.G.Haggmark  &  J.P.Biersack : PHYS. REV. B 15 , 2458
    // (1977).
    double a = 0.56258;
    double b = 1.1776;
    double c = 0.62680;

    double x = log( b * epsiln );
    double x1 = ( c-1.0 ) * x;
    double x2 = -2.0 * x;
    return ( expint(x1) - expint(x2) ) / a / b;

}


double PSI1D_Backscattering::expint(double x)
{
    // incomplete gamma function : gamma ( 0 , x )
    // gamma( 0 , x ) = - ei ( - x )	for x > 0
    // gamma( 0 , x ) = - eibar( - x )	for x < 0
    double a[4] = { .2677737343, 8.6347608925, 18.059016973, 8.5733287401 };
    double b[4] = { 3.9584969228, 21.0996530827, 25.6329561486, 9.5733223454 };
    double conver = 1e-7;
    double z, p, f, t, q;

    z = x;
    if (x-1.0 < 0.0)
    {
        p = 0.57721566490153;
        if (x < 0.0 && x + 83.5 < 0.0)
        {
            return - 1.e38;
        }
        else if( x == 0.0)
        {
            return - p;
        }
        else if(x > 0.0 && x + 83.5 >= 0.0 || x > 0.0 && x + 83.5 < 0.0 )
        {
            p = p + log(abs(z));
            f = 1.0;
            t = - z;
            q = t + p;
            while ( abs(p-q) - conver > 0.0)
            {
                p = q;
                t = - z * f * t;
                f = f + 1.0;
                t = t / pow(f, 2.0);
                q = t + p;
            }
            return - p;
        }
    }
    else
    {
        if ( x - 84.0 <= 0.0)
        {
            p = (((( z + a[3] ) * z + a[2]) * z + a[1]) * z + a[0] )
                / (((( z + b[3] ) * z + b[2] ) * z + b[1] ) * z + b[0] );
            return exp( -z ) * p / z;
        }
        else
        {
            return 1.e-38;
        }
    }

}
