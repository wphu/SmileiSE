c  a  FORTRAN FUNCTION to evaluate empirical formulas for number and
c  energy  BACKSCATTERING  COEFFICIENTS  of Light ions incident on
c  Elemental and Compound targets.
c
c***********************************************************************
c
c  program copied from
c
c  NAYOGA UNIVERSITY REPORT #  IPPJ-AM-41 ( Data On The Backscattering
c  Coefficients Of Light Ions From Solids ) - R.Ito, T.Tabata, N.Itoh,
c  K.Morita, T.Kato, H.Tawara.
c
c  Modifications :

c	i)  Converted into double precision.

c	ii) Include heavy particle backscattering as described in the
c	    paper " Present Status of Data Compilation on Ion
c	    Backscattering. - Tatsuo Tabata and Rinsuke Ito
c	    Osaka Prefectural Radiation Research Institute,
c	    Shinke-cho, Sakai, Osaka 593 " IPPJ-AM-64, Nagoya
c           University Report.
c
c***********************************************************************
c
c  call RNION for the number backscattering coefficient and
c  REION for the energy backscattering coefficient.
c
c  valid upto any mass number for incident ions, provided m1/m2 < 4.8,
c  where m1 is the mass of the incident ion and m2 is the mass of the
c  target atom.
c
      function rnion( theta,energy,nz1,m1,ne,nz2,nw )
c
c  number backscattering coefficient of light ions.
c
c meaning of input variables
c
c	theta 	= angle of incidence in degree.
c	energy	= incident kinetic energy in eV.
c	nz1 	= atomic number of the projectile.
c	m1	= mass number of the projectile.
c	ne	= number of constituent elements in the target.
c	nz2	= array for atomic numbers of the constituents.
c	nw	= array for relative numbers of the constituents.
c
c	example : for Tritium ions of 1 KeV energy incident on the TiO2 
c		( Titanium Dioxide ) target, Energy=1000 , nz1=1 ,
c		m1=3 , ne=2 , nz2(1)=22 , nw(1)=1 , nz2(2)=8 and 
c		nw(2)=2 .
c
c  tables prepared
c	a1t	= table of masses of projectiles
c	a2t	= table of atomic weights of elemental targets
c		IUPAC (1983) , C12=12 .
c	d	= table of correction factors for the low energy
c		electronic stopping cross-section.
c
      implicit real*8 (a-h,o-z)
      real*8 a1t(95) , a2t(92) , d(92)
      dimension nz2(ne) ,nw(ne)
c
c  data table a1t : * * * * * * * * * * * * * * * * * * * * * * * * * * 
c
      data a1t /
c  z1 = 1 ( for hydrogen , deuterium and tritium )
     # 1.0072766,2.0135536,3.1055011,
c  z1 = 2 ( for helium-3 and helium-4 )
     # 3.0149325,4.0015059, 
c  z1 = 3-10
     # 6.941,9.01218,10.811,12.011,14.0067,15.9994,18.998403,20.179,
c  z1 = 11 - 20
     # 22.98977,24.305,26.98154,28.0855,30.97376,32.066,35.453,39.948,
     # 39.0893,40.078,
c  z1 = 21 - 30
     # 44.95591,47.88,50.9415,51.9961,54.9380,55.847,58.9332,58.69,
     # 63.546,65.39,
c  z1 = 31 - 40
     # 69.723,72.59,74.9216,78.96,79.904,83.80,85.4678,87.62,88.9059,
     # 91.224,
c  z1 = 41 - 50
     # 92.9064,95.94,99.0,101.07,102.9055,106.42,107.8682,112.41,114.82,
     # 118.710,
c  z1 = 51 - 60
     # 121.75,127.60,126.9045,131.29,132.9054,137.33,138.9055,140.12,
     # 140.9077,144.24,
c  z1 = 61 - 70
     # 145.0,150.36,151.96,157.25,158.9254,162.50,164.9304,167.26,
     # 168.9342,173.04,
c  z1 = 71 - 80
     # 174.967,178.49,180.9479,183.85,186.207,190.2,192.22,195.08,
     # 196.9655,200.59,
c  z1 = 81 - 90
     # 204.383,207.2,208.9804,210.0,210.0,222.0,223.0,226.0,227.0,232.0,
c  z12 = 91 - 92
     # 231.0,238.0289 /
c
c  data table a2t : * * * * * * * * * * * * * * * * * * * * * * * * * * 
c
      data a2t /
c  z2 = 1 - 10
     # 1.00794,4.002602,6.941,9.01218,10.811,12.011,14.0067,15.9994,
     # 18.998403,20.179,
c  z2 = 11 - 20
     # 22.98977,24.305,26.98154,28.0855,30.97376,32.066,35.453,39.948,
     # 39.0893,40.078,
c  z2 = 21 - 30
     # 44.95591,47.88,50.9415,51.9961,54.9380,55.847,58.9332,58.69,
     # 63.546,65.39,
c  z2 = 31 - 40
     # 69.723,72.59,74.9216,78.96,79.904,83.80,85.4678,87.62,88.9059,
     # 91.224,
c  z2 = 41 - 50
     # 92.9064,95.94,99.0,101.07,102.9055,106.42,107.8682,112.41,114.82,
     # 118.710,
c  z2 = 51 - 60
     # 121.75,127.60,126.9045,131.29,132.9054,137.33,138.9055,140.12,
     # 140.9077,144.24,
c  z2 = 61 - 70
     # 145.0,150.36,151.96,157.25,158.9254,162.50,164.9304,167.26,
     # 168.9342,173.04,
c  z2 = 71 - 80
     # 174.967,178.49,180.9479,183.85,186.207,190.2,192.22,195.08,
     # 196.9655,200.59,
c  z2 = 81 - 90
     # 204.383,207.2,208.9804,210.0,210.0,222.0,223.0,226.0,227.0,232.0,
c  z2 = 91 - 92
     # 231.0,238.0289 /
c
c  data table d : * * * * * * * * * * * * * * * * * * * * * * * * * * 
c
      data d / 
c  z2 = 1 - 10
     # 0.9341d0,0.6693d0,0.6654d0,0.9712d0,1.007d0,1.024d0,1.111d0,
     # 0.9699d0,0.7357d0,0.6842d0,
c  z2 = 11 - 20
     # 0.8769d0,1.290d0,1.395d0,1.378d0,1.063d0,1.123d0,1.632d0,1.839d0,
     # 1.642d0,1.749d0,
c  z2 = 21 - 30
     # 1.638d0,1.523d0,1.396d0,1.236d0,1.072d0,1.083d0,0.9624d0,1.085d0,
     # 1.125d0,1.277d0,
c  z2 = 31 - 40
     # 1.525d0,1.675d0,1.601d0,1.762d0,1.679d0,1.914d0,1.696d0,1.884d0,
     # 1.9d0,1.993d0,
c  z2 = 41 - 50
     # 2.039d0,1.894d0,2.001d0,1.795d0,1.738d0,1.534d0,1.644d0,1.698d0,
     # 1.816d0,1.866d0,
c  z2 = 51 - 60
     # 2.181d0,2.027d0,2.240d0,2.384d0,2.108d0,2.283d0,2.321d0,2.159d0,
     # 2.1d0,2.042d0,
c  z2 = 61 - 70
     # 1.986d0,1.932d0,1.879d0,1.931d0,1.779d0,1.578d0,1.492d0,1.448d0,
     # 1.406d0,1.365d0,
c  z2 = 71 - 80
     # 1.394d0,1.431d0,1.348d0,1.3d0,1.477d0,1.439d0,1.403d0,1.269d0,
     # 1.376d0,1.22d0,
c  z2 = 81 - 90
     # 1.336d0,1.504d0,1.683d0,1.739d0,1.751d0,1.744d0,1.959d0,2.115d0,
     # 2.167d0,2.170d0,
c  z2 = 91 - 92
     # 2.084d0,2.050d0 /

      f(z1,z2) = dsqrt ( z1 ** (2.d0/3.d0) + z2 ** (2.d0/3.d0) )
      fneps(e,z1,a1,z2,a2) = 0.032534 * e / (z1 * z2 * ( 1.d0 + a1/a2 )
     #	* f(z1,z2))
      rne(th) = r0 + ( 1.d0 - r0 ) / ( 1.d0 + as1 / dtan(th) **
     #	( 2.d0 * as2 ) )
 600  format('	undefined projectile is called	')
      n = 1
      go to 1

      entry reion(theta,energy,nz1,m1,ne,nz2,nw)

c  energy backscattering coefficient of light ions

      n = 0
 1    if (nz1.eq.1) go to 3
      if (nz1.eq.2) go to 4
      if (nz1.gt.2) go to 5
 2    write(6,600)
      return
 3    if (m1.lt.1.or.m1.gt.3) go to 2
      a1 = a1t(m1)
      go to 6
 4    if (m1.lt.3.or.m1.gt.4) go to 2
      a1 = a1t( m1+1 )
      go to 6
 5    a1 = a1t( nz1+3 )
c	write(*,*)'atomic # of projectile > 2'
c  finding effective a2 :
	a2eff = 0.0d0
	sum = 0.0d0
	do i = 1,ne
		a2 = a2t(nz2(i))
		a2eff = a2eff + nw(i) * a2
		sum = sum + nw(i)
	enddo
	a2 = a2eff / sum
	delta = a1 / a2
c	write(*,*)'delta = ',delta
	if (delta.gt.4.5d0) write(6,600)
 6    z1 = nz1
      z1p = z1 ** (2.d0/3.d0)
      z2eff = 0.d0
      a2eff = 0.d0
      saeff = 0.d0
      sbeff = 0.d0
      sum = 0.d0

      do 70 i=1,ne
      z2 = nz2(i)
      a2 = a2t(nz2(i))
      wi = nw(i)
      z2eff = z2eff + wi * z2
      a2eff = a2eff + wi * a2
      sum = sum + wi
      fmu = a1/a2
      aa = 1.d0 + fmu
      z2p = z2 ** (2.d0/3.d0)
      zz = dsqrt ( 1.d0 + z1p/z2p )
      eps = fneps(energy,z1,a1,z2,a2)
      sa = 0.0793 * dsqrt (eps) / fmu
      w = z1 * zz * aa / z2p
      wi = wi * z1 / w / a2
      saeff = saeff + wi * sa
      sl = d(nz2(i)) * z1p * ( aa/zz ) ** 1.5d0 * sa / dsqrt(a1)
      if ( z2 - 12.9d0 ) 10 , 10 , 20
 10   fi0 = 12.d0 + 7.d0 / z2
      go to 30
 20   fi0 = 9.76 + 58.5 / z2 ** 1.19
 30   tau = eps * z2 ** 2.d0 * w / a1 / 3.03055d7
      beta2 = tau * ( tau + 2.d0 ) / ( tau + 1.d0 ) ** 2
      epsb = 1.022d6 * beta2 / z2 / fi0
      if ( z1 - 2.9d0 ) 40,40,50
 40   c = 100.d0 * z1 / z2
      go to 60
 50   c = 5.d0
 60   sb = 61.474d0 / fmu * w * ( dlog ( epsb/(1.d0-beta2) + 1.d0 +
     #	c/epsb ) - beta2 ) / fi0 / epsb
      sb = 2.d0 * eps / rangen(eps,energy) + sl * sb / (0.25d0 * sl +
     #  sb )
      sbeff = sbeff + wi * sb
 70   continue
      z2 = z2eff / sum
      a2 = a2eff / sum
      sa = saeff
      sb = sbeff
      eps = fneps(energy,z1,a1,z2,a2)
	if (delta.gt.0.5) then
		redRN = 1.0d0 /
     @			( 1+(eps/0.104d0)**0.577d0 + (eps/0.73d0)**1.5 )
		fmu = a1/a2
		psinum = 1.0d0+24.1d0*fmu**(-3.995d0)
		psiden1 = fmu / ( (1+fmu) * (eps/1.84d0+1.0d0) )
		psiden2 = (1.0d0+fmu)**3 * eps / ( fmu**3 * (1.0d0 -
     @			3.0d0/(2.0d0*fmu) + 0.9d0/fmu**2 +
     @			1.0d0/(2.0d0*fmu**3) ) * (eps+13.3) )
		psi = psinum / ( psiden1 + psiden2 )

c  sbeff/saeff corresponds to the rho_a/rho-t in the reference paper ( of course
c  sbeff is not equal to rho_a nor is saeff equal to rho_t).

		sklfact = z1**0.66666667d0 / dsqrt(a1) * sb/sa *
     @			fmu**2 / (1.0d0+fmu)**2 / psi
		r0heavy = redRN / sklfact
		r = 1.d0 / ( 1.d0 + (eps/0.133d0 ) ** 0.285 ) + 0.530 /
     #				 ( 1.d0 + (85.0d0/eps) **1.46d0 )
c		write(*,*)'did heavy projectile calculation for delta > 0.5 '
      		if ( n.eq.1 ) go to 80
		th = theta * 1.745329d-2
		if (th.eq.0.0d0) then
c			write(*,*)'did theta=0 calculation of reion'
			reion = r0heavy * r
			go to 79
		else
c			write(*,*)'did theta > 0 calculation for reion'
			r0 = r0heavy * r
      			as1 = 17.9 * eps ** 0.453d0
      			as2 = 0.771d0 / eps ** 0.014d0
      			reion = rne(th)
			go to 79
		endif
	else
		continue
	endif

c	write(*,*)'did light projectile calculation delta < 0.5 '

      feps = (2.d0 * eps - 3.d0 ) / ( eps + 1.d0 )
      r0 = 0.705d0 * sa / sb / ( 1.0d0 + a1/a2 ) ** feps
      r0 = r0 / (1.d0 + (eps/0.047d0) ** 0.597d0 + ( eps / 0.619d0 ) ** 
     #	1.5d0 )
      th = theta * 1.745329d-2
      if (n.eq.1) go to 80
      if (th.eq.0) go to 78
      as1 = 17.9 * eps ** 0.453d0
      as2 = 0.771d0 / eps ** 0.014d0
      reion = rne(th)
	go to 79
 78   reion = r0
 79   return
 80   r = 1.d0 / ( 1.d0 + (eps/0.133d0 ) ** 0.285 ) + 0.530 / ( 1.d0 +
     #	(85.d0/eps) ** 1.46d0 )
	if (delta.gt.0.5d0) then
c		write(*,*)'did heavy proj calculation for rnion'
		if (th.ne.0.0d0) then
c			write(*,*)'did theta not = 0 calc for rnion'
			r0 = r0heavy
      			as1 = 7.38d0 * eps ** 0.359d0
      			as2 = 0.836d0 /eps ** 0.087d0
			rnion = rne(th)
			return
		else
c			write(*,*)'did theta = 0 calc. for rnion'
			rnion = r0heavy
c			write(*,*)r0heavy, rnion, r, reion
			return
		endif
	else
c		write(*,*)'did light projectile calculation for rnion'
		r0 = r0 / r
		if (th.eq.0.d0) go to 88
		as1 = 7.38d0 * eps ** 0.359d0
		as2 = 0.836d0 /eps ** 0.087d0
		rnion = rne(th)
		return
  88   		rnion = r0
		return
	endif

	end
c
      function rangen(epsiln,energy)
c  reduced range of ions for nuclear stopping only
c  W.D.Wilson , L.G.Haggmark  &  J.P.Biersack : PHYS. REV. B 15 , 2458
c  (1977).
      implicit real*8 (a-h,o-z)
      data a,b,c / 0.56258d0 , 1.1776d0 , 0.62680d0 /
      x = dlog ( b * epsiln )
      x1 = ( c-1.d0 ) * x
      x2 = -2.d0 * x
      rangen = ( expint(x1) - expint(x2) ) / a / b
      return
      end
c
      function expint(x)
      implicit real*8(a-h,o-z)
c  incomplete gamma function : gamma ( 0 , x )
c  gamma( 0 , x ) = - ei ( - x )	for x > 0
c  gamma( 0 , x ) = - eibar( - x )	for x < 0
      real*8 a(4),b(4),conver,f,p,q,t,z
      data a / 0.2677737343d0 , 8.6347608925 , 18.059016973d0 ,
     #	8.5733287401 /
      data b / 3.9584969228 , 21.0996530827 , 25.6329561486 ,
     #	9.5733223454 /
      data conver / 1.d-7 /
      z = dble(x)
      if (x-1.d0) 10,40,40
 10   p = 0.57721566490153
      if (x) 12,30,15
 12   if ( x + 83.5d0 ) 13,15,15
 13   expint = - 1.d38
      return
 15   p = p + dlog(dabs(z))
      f = 1.d0
      t = - z
 20   q = t + p
      if ( dabs(p-q) - conver ) 30,30,22
 22   p = q
      t = - z * f * t
      f = f + 1.d0
      t = t / f ** 2.d0
      go to 20
 30   expint = - p
      return
 40   if ( x - 84.d0 ) 43,43,41
 41   expint = 1.d-38
      return
 43   p = (((( z + a(4) ) * z + a(3)) * z + a(2)) * z + a(1) ) /
     #	(((( z + b(4) ) * z + b(3)) * z + b(2)) * z + b(1) )
      expint = dexp ( -z ) * p / z
      return
      end
