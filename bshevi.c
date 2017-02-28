/* bshevi.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static doublereal c_b10 = .66666666666666663;
static doublereal c_b15 = 1.5;
static doublereal c_b19 = 1.19;
static doublereal c_b20 = 2.;
static doublereal c_b26 = .577;
static doublereal c_b28 = -3.995;
static doublereal c_b29 = .66666667;
static doublereal c_b30 = .285;
static doublereal c_b31 = 1.46;
static doublereal c_b34 = .453;
static doublereal c_b35 = .014;
static doublereal c_b36 = .597;
static doublereal c_b43 = .359;
static doublereal c_b44 = .087;

/*  a  FORTRAN FUNCTION to evaluate empirical formulas for number and */
/*  energy  BACKSCATTERING  COEFFICIENTS  of Light ions incident on */
/*  Elemental and Compound targets. */

/* *********************************************************************** */

/*  program copied from */

/*  NAYOGA UNIVERSITY REPORT #  IPPJ-AM-41 ( Data On The Backscattering */
/*  Coefficients Of Light Ions From Solids ) - R.Ito, T.Tabata, N.Itoh, */
/*  K.Morita, T.Kato, H.Tawara. */

/*  Modifications : */
/* 	i)  Converted into double precision. */
/* 	ii) Include heavy particle backscattering as described in the */
/* 	    paper " Present Status of Data Compilation on Ion */
/* 	    Backscattering. - Tatsuo Tabata and Rinsuke Ito */
/* 	    Osaka Prefectural Radiation Research Institute, */
/* 	    Shinke-cho, Sakai, Osaka 593 " IPPJ-AM-64, Nagoya */
/*           University Report. */

/* *********************************************************************** */

/*  call RNION for the number backscattering coefficient and */
/*  REION for the energy backscattering coefficient. */

/*  valid upto any mass number for incident ions, provided m1/m2 < 4.8, */
/*  where m1 is the mass of the incident ion and m2 is the mass of the */
/*  target atom. */

doublereal rnion_0_(int n__, doublereal *theta, doublereal *energy, integer *
	nz1, integer *m1, integer *ne, integer *nz2, integer *nw)
{
    /* Initialized data */

    static doublereal a1t[95] = { 1.0072766,2.0135536,3.1055011,3.0149325,
	    4.0015059,6.941,9.01218,10.811,12.011,14.0067,15.9994,18.998403,
	    20.179,22.98977,24.305,26.98154,28.0855,30.97376,32.066,35.453,
	    39.948,39.0893,40.078,44.95591,47.88,50.9415,51.9961,54.938,
	    55.847,58.9332,58.69,63.546,65.39,69.723,72.59,74.9216,78.96,
	    79.904,83.8,85.4678,87.62,88.9059,91.224,92.9064,95.94,99.,101.07,
	    102.9055,106.42,107.8682,112.41,114.82,118.71,121.75,127.6,
	    126.9045,131.29,132.9054,137.33,138.9055,140.12,140.9077,144.24,
	    145.,150.36,151.96,157.25,158.9254,162.5,164.9304,167.26,168.9342,
	    173.04,174.967,178.49,180.9479,183.85,186.207,190.2,192.22,195.08,
	    196.9655,200.59,204.383,207.2,208.9804,210.,210.,222.,223.,226.,
	    227.,232.,231.,238.0289 };
    static doublereal a2t[92] = { 1.00794,4.002602,6.941,9.01218,10.811,
	    12.011,14.0067,15.9994,18.998403,20.179,22.98977,24.305,26.98154,
	    28.0855,30.97376,32.066,35.453,39.948,39.0893,40.078,44.95591,
	    47.88,50.9415,51.9961,54.938,55.847,58.9332,58.69,63.546,65.39,
	    69.723,72.59,74.9216,78.96,79.904,83.8,85.4678,87.62,88.9059,
	    91.224,92.9064,95.94,99.,101.07,102.9055,106.42,107.8682,112.41,
	    114.82,118.71,121.75,127.6,126.9045,131.29,132.9054,137.33,
	    138.9055,140.12,140.9077,144.24,145.,150.36,151.96,157.25,
	    158.9254,162.5,164.9304,167.26,168.9342,173.04,174.967,178.49,
	    180.9479,183.85,186.207,190.2,192.22,195.08,196.9655,200.59,
	    204.383,207.2,208.9804,210.,210.,222.,223.,226.,227.,232.,231.,
	    238.0289 };
    static doublereal d__[92] = { .9341,.6693,.6654,.9712,1.007,1.024,1.111,
	    .9699,.7357,.6842,.8769,1.29,1.395,1.378,1.063,1.123,1.632,1.839,
	    1.642,1.749,1.638,1.523,1.396,1.236,1.072,1.083,.9624,1.085,1.125,
	    1.277,1.525,1.675,1.601,1.762,1.679,1.914,1.696,1.884,1.9,1.993,
	    2.039,1.894,2.001,1.795,1.738,1.534,1.644,1.698,1.816,1.866,2.181,
	    2.027,2.24,2.384,2.108,2.283,2.321,2.159,2.1,2.042,1.986,1.932,
	    1.879,1.931,1.779,1.578,1.492,1.448,1.406,1.365,1.394,1.431,1.348,
	    1.3,1.477,1.439,1.403,1.269,1.376,1.22,1.336,1.504,1.683,1.739,
	    1.751,1.744,1.959,2.115,2.167,2.17,2.084,2.05 };

    /* Format strings */
    static char fmt_600[] = "(\002\tundefined projectile is called\t\002)";

    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2, d__3, d__4;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal), log(
	    doublereal), tan(doublereal);

    /* Local variables */
    static doublereal c__;
    static integer i__, n;
    static doublereal r__, w, a1, a2, r0, z1, z2, aa, sa, sb, th, sl, wi, zz, 
	    fi0, as1, as2, z1p, z2p, fmu, eps, tau, psi, sum, epsb, feps, 
	    a2eff, beta2, z2eff, saeff, sbeff, delta, redrn;
    extern doublereal rangen_(doublereal *, doublereal *);
    static doublereal psinum, psiden1, psiden2, r0heavy, sklfact;

    /* Fortran I/O blocks */
    static cilist io___5 = { 0, 6, 0, fmt_600, 0 };
    static cilist io___12 = { 0, 6, 0, fmt_600, 0 };



/*  number backscattering coefficient of light ions. */

/* meaning of input variables */

/* 	theta 	= angle of incidence in degree. */
/* 	energy	= incident kinetic energy in eV. */
/* 	nz1 	= atomic number of the projectile. */
/* 	m1	= mass number of the projectile. */
/* 	ne	= number of constituent elements in the target. */
/* 	nz2	= array for atomic numbers of the constituents. */
/* 	nw	= array for relative numbers of the constituents. */

/* 	example : for Tritium ions of 1 KeV energy incident on the TiO2 */
/* 		( Titanium Dioxide ) target, Energy=1000 , nz1=1 , */
/* 		m1=3 , ne=2 , nz2(1)=22 , nw(1)=1 , nz2(2)=8 and */
/* 		nw(2)=2 . */

/*  tables prepared */
/* 	a1t	= table of masses of projectiles */
/* 	a2t	= table of atomic weights of elemental targets */
/* 		IUPAC (1983) , C12=12 . */
/* 	d	= table of correction factors for the low energy */
/* 		electronic stopping cross-section. */


/*  data table a1t : * * * * * * * * * * * * * * * * * * * * * * * * * * */

    /* Parameter adjustments */
    --nz2;
    --nw;

    /* Function Body */
    switch(n__) {
	case 1: goto L_reion;
	}

/*  z1 = 1 ( for hydrogen , deuterium and tritium ) */
/*  z1 = 2 ( for helium-3 and helium-4 ) */
/*  z1 = 3-10 */
/*  z1 = 11 - 20 */
/*  z1 = 21 - 30 */
/*  z1 = 31 - 40 */
/*  z1 = 41 - 50 */
/*  z1 = 51 - 60 */
/*  z1 = 61 - 70 */
/*  z1 = 71 - 80 */
/*  z1 = 81 - 90 */
/*  z12 = 91 - 92 */

/*  data table a2t : * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*  z2 = 1 - 10 */
/*  z2 = 11 - 20 */
/*  z2 = 21 - 30 */
/*  z2 = 31 - 40 */
/*  z2 = 41 - 50 */
/*  z2 = 51 - 60 */
/*  z2 = 61 - 70 */
/*  z2 = 71 - 80 */
/*  z2 = 81 - 90 */
/*  z2 = 91 - 92 */

/*  data table d : * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*  z2 = 1 - 10 */
/*  z2 = 11 - 20 */
/*  z2 = 21 - 30 */
/*  z2 = 31 - 40 */
/*  z2 = 41 - 50 */
/*  z2 = 51 - 60 */
/*  z2 = 61 - 70 */
/*  z2 = 71 - 80 */
/*  z2 = 81 - 90 */
/*  z2 = 91 - 92 */
/* L600: */
    n = 1;
    goto L1;

L_reion:
/*  energy backscattering coefficient of light ions */
    n = 0;
L1:
    if (*nz1 == 1) {
	goto L3;
    }
    if (*nz1 == 2) {
	goto L4;
    }
    if (*nz1 > 2) {
	goto L5;
    }
L2:
    s_wsfe(&io___5);
    e_wsfe();
    return ret_val;
L3:
    if (*m1 < 1 || *m1 > 3) {
	goto L2;
    }
    a1 = a1t[*m1 - 1];
    goto L6;
L4:
    if (*m1 < 3 || *m1 > 4) {
	goto L2;
    }
    a1 = a1t[*m1];
    goto L6;
L5:
    a1 = a1t[*nz1 + 2];
/* 	write(*,*)'atomic # of projectile > 2' */
/*  finding effective a2 : */
    a2eff = 0.;
    sum = 0.;
    i__1 = *ne;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a2 = a2t[nz2[i__] - 1];
	a2eff += nw[i__] * a2;
	sum += nw[i__];
    }
    a2 = a2eff / sum;
    delta = a1 / a2;
/* 	write(*,*)'delta = ',delta */
    if (delta > 4.5) {
	s_wsfe(&io___12);
	e_wsfe();
    }
L6:
    z1 = (doublereal) (*nz1);
    z1p = pow_dd(&z1, &c_b10);
    z2eff = 0.;
    a2eff = 0.;
    saeff = 0.;
    sbeff = 0.;
    sum = 0.;
    i__1 = *ne;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z2 = (doublereal) nz2[i__];
	a2 = a2t[nz2[i__] - 1];
	wi = (doublereal) nw[i__];
	z2eff += wi * z2;
	a2eff += wi * a2;
	sum += wi;
	fmu = a1 / a2;
	aa = fmu + 1.;
	z2p = pow_dd(&z2, &c_b10);
	zz = sqrt(z1p / z2p + 1.);
	eps = .032534f * *energy / (z1 * z2 * (1. + a1 / a2) * sqrt(pow_dd(&
		z1, &c_b10) + pow_dd(&z2, &c_b10)));
	sa = sqrt(eps) * .0793f / fmu;
	w = z1 * zz * aa / z2p;
	wi = wi * z1 / w / a2;
	saeff += wi * sa;
	d__1 = aa / zz;
	sl = d__[nz2[i__] - 1] * z1p * pow_dd(&d__1, &c_b15) * sa / sqrt(a1);
	if (z2 - 12.9 <= 0.) {
	    goto L10;
	} else {
	    goto L20;
	}
L10:
	fi0 = 7. / z2 + 12.;
	goto L30;
L20:
	fi0 = 58.5f / pow_dd(&z2, &c_b19) + 9.76f;
L30:
	tau = eps * pow_dd(&z2, &c_b20) * w / a1 / 30305500.;
/* Computing 2nd power */
	d__1 = tau + 1.;
	beta2 = tau * (tau + 2.) / (d__1 * d__1);
	epsb = beta2 * 1.022e6 / z2 / fi0;
	if (z1 - 2.9 <= 0.) {
	    goto L40;
	} else {
	    goto L50;
	}
L40:
	c__ = z1 * 100. / z2;
	goto L60;
L50:
	c__ = 5.;
L60:
	sb = 61.474 / fmu * w * (log(epsb / (1. - beta2) + 1. + c__ / epsb) - 
		beta2) / fi0 / epsb;
	sb = eps * 2. / rangen_(&eps, energy) + sl * sb / (sl * .25 + sb);
	sbeff += wi * sb;
/* L70: */
    }
    z2 = z2eff / sum;
    a2 = a2eff / sum;
    sa = saeff;
    sb = sbeff;
    eps = .032534f * *energy / (z1 * z2 * (1. + a1 / a2) * sqrt(pow_dd(&z1, &
	    c_b10) + pow_dd(&z2, &c_b10)));
    if (delta > .5f) {
	d__1 = eps / .104;
	d__2 = eps / .73;
	redrn = 1. / (pow_dd(&d__1, &c_b26) + 1 + pow_dd(&d__2, &c_b15));
	fmu = a1 / a2;
	psinum = pow_dd(&fmu, &c_b28) * 24.1 + 1.;
	psiden1 = fmu / ((fmu + 1) * (eps / 1.84 + 1.));
/* Computing 3rd power */
	d__1 = fmu + 1.;
/* Computing 3rd power */
	d__2 = fmu;
/* Computing 2nd power */
	d__3 = fmu;
/* Computing 3rd power */
	d__4 = fmu;
	psiden2 = d__1 * (d__1 * d__1) * eps / (d__2 * (d__2 * d__2) * (1. - 
		3. / (fmu * 2.) + .9 / (d__3 * d__3) + 1. / (d__4 * (d__4 * 
		d__4) * 2.)) * (eps + 13.3f));
	psi = psinum / (psiden1 + psiden2);
/*  sbeff/saeff corresponds to the rho_a/rho-t in the reference paper ( of course */
/*  sbeff is not equal to rho_a nor is saeff equal to rho_t). */
/* Computing 2nd power */
	d__1 = fmu;
/* Computing 2nd power */
	d__2 = fmu + 1.;
	sklfact = pow_dd(&z1, &c_b29) / sqrt(a1) * sb / sa * (d__1 * d__1) / (
		d__2 * d__2) / psi;
	r0heavy = redrn / sklfact;
	d__1 = eps / .133;
	d__2 = 85. / eps;
	r__ = 1. / (pow_dd(&d__1, &c_b30) + 1.) + .53f / (pow_dd(&d__2, &
		c_b31) + 1.);
/* 		write(*,*)'did heavy projectile calculation for delta > 0.5 ' */
	if (n == 1) {
	    goto L80;
	}
	th = *theta * .01745329;
	if (th == 0.) {
/* 			write(*,*)'did theta=0 calculation of reion' */
	    ret_val = r0heavy * r__;
	    goto L79;
	} else {
/* 			write(*,*)'did theta > 0 calculation for reion' */
	    r0 = r0heavy * r__;
	    as1 = pow_dd(&eps, &c_b34) * 17.9f;
	    as2 = .771 / pow_dd(&eps, &c_b35);
	    d__1 = tan(th);
	    d__2 = 2. * as2;
	    ret_val = r0 + (1. - r0) / (1. + as1 / pow_dd(&d__1, &d__2));
	    goto L79;
	}
    } else {
    }
/* 	write(*,*)'did light projectile calculation delta < 0.5 ' */
    feps = (eps * 2. - 3.) / (eps + 1.);
    d__1 = a1 / a2 + 1.;
    r0 = sa * .705 / sb / pow_dd(&d__1, &feps);
    d__1 = eps / .047;
    d__2 = eps / .619;
    r0 /= pow_dd(&d__1, &c_b36) + 1. + pow_dd(&d__2, &c_b15);
    th = *theta * .01745329;
    if (n == 1) {
	goto L80;
    }
    if (th == 0.) {
	goto L78;
    }
    as1 = pow_dd(&eps, &c_b34) * 17.9f;
    as2 = .771 / pow_dd(&eps, &c_b35);
    d__1 = tan(th);
    d__2 = 2. * as2;
    ret_val = r0 + (1. - r0) / (1. + as1 / pow_dd(&d__1, &d__2));
    goto L79;
L78:
    ret_val = r0;
L79:
    return ret_val;
L80:
    d__1 = eps / .133;
    d__2 = 85. / eps;
    r__ = 1. / (pow_dd(&d__1, &c_b30) + 1.) + .53f / (pow_dd(&d__2, &c_b31) + 
	    1.);
    if (delta > .5) {
/* 		write(*,*)'did heavy proj calculation for rnion' */
	if (th != 0.) {
/* 			write(*,*)'did theta not = 0 calc for rnion' */
	    r0 = r0heavy;
	    as1 = pow_dd(&eps, &c_b43) * 7.38;
	    as2 = .836 / pow_dd(&eps, &c_b44);
	    d__1 = tan(th);
	    d__2 = 2. * as2;
	    ret_val = r0 + (1. - r0) / (1. + as1 / pow_dd(&d__1, &d__2));
	    return ret_val;
	} else {
/* 			write(*,*)'did theta = 0 calc. for rnion' */
	    ret_val = r0heavy;
/* 			write(*,*)r0heavy, rnion, r, reion */
	    return ret_val;
	}
    } else {
/* 		write(*,*)'did light projectile calculation for rnion' */
	r0 /= r__;
	if (th == 0.) {
	    goto L88;
	}
	as1 = pow_dd(&eps, &c_b43) * 7.38;
	as2 = .836 / pow_dd(&eps, &c_b44);
	d__1 = tan(th);
	d__2 = 2. * as2;
	ret_val = r0 + (1. - r0) / (1. + as1 / pow_dd(&d__1, &d__2));
	return ret_val;
L88:
	ret_val = r0;
	return ret_val;
    }
    return ret_val;
} /* rnion_ */

doublereal rnion_(doublereal *theta, doublereal *energy, integer *nz1, 
	integer *m1, integer *ne, integer *nz2, integer *nw)
{
    return rnion_0_(0, theta, energy, nz1, m1, ne, nz2, nw);
    }

doublereal reion_(doublereal *theta, doublereal *energy, integer *nz1, 
	integer *m1, integer *ne, integer *nz2, integer *nw)
{
    return rnion_0_(1, theta, energy, nz1, m1, ne, nz2, nw);
    }


doublereal rangen_(doublereal *epsiln, doublereal *energy)
{
    /* Initialized data */

    static doublereal a = .56258;
    static doublereal b = 1.1776;
    static doublereal c__ = .6268;

    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static doublereal x, x1, x2;
    extern doublereal expint_(doublereal *);

/*  reduced range of ions for nuclear stopping only */
/*  W.D.Wilson , L.G.Haggmark  &  J.P.Biersack : PHYS. REV. B 15 , 2458 */
/*  (1977). */
    x = log(b * *epsiln);
    x1 = (c__ - 1.) * x;
    x2 = x * -2.;
    ret_val = (expint_(&x1) - expint_(&x2)) / a / b;
    return ret_val;
} /* rangen_ */


doublereal expint_(doublereal *x)
{
    /* Initialized data */

    static doublereal a[4] = { .2677737343,8.6347608925,18.059016973,
	    8.5733287401 };
    static doublereal b[4] = { 3.9584969228,21.0996530827,25.6329561486,
	    9.5733223454 };
    static doublereal conver = 1e-7;

    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double log(doublereal), pow_dd(doublereal *, doublereal *), exp(
	    doublereal);

    /* Local variables */
    static doublereal f, p, q, t, z__;

/*  incomplete gamma function : gamma ( 0 , x ) */
/*  gamma( 0 , x ) = - ei ( - x )	for x > 0 */
/*  gamma( 0 , x ) = - eibar( - x )	for x < 0 */
    z__ = *x;
    if (*x - 1. >= 0.) {
	goto L40;
    } else {
	goto L10;
    }
L10:
    p = .57721566490153f;
    if (*x < 0.) {
	goto L12;
    } else if (*x == 0) {
	goto L30;
    } else {
	goto L15;
    }
L12:
    if (*x + 83.5 >= 0.) {
	goto L15;
    } else {
	goto L13;
    }
L13:
    ret_val = -1e38;
    return ret_val;
L15:
    p += log((abs(z__)));
    f = 1.;
    t = -z__;
L20:
    q = t + p;
    if ((d__1 = p - q, abs(d__1)) - conver <= 0.) {
	goto L30;
    } else {
	goto L22;
    }
L22:
    p = q;
    t = -z__ * f * t;
    f += 1.;
    t /= pow_dd(&f, &c_b20);
    goto L20;
L30:
    ret_val = -p;
    return ret_val;
L40:
    if (*x - 84. <= 0.) {
	goto L43;
    } else {
	goto L41;
    }
L41:
    ret_val = 1e-38;
    return ret_val;
L43:
    p = ((((z__ + a[3]) * z__ + a[2]) * z__ + a[1]) * z__ + a[0]) / ((((z__ + 
	    b[3]) * z__ + b[2]) * z__ + b[1]) * z__ + b[0]);
    ret_val = exp(-z__) * p / z__;
    return ret_val;
} /* expint_ */

