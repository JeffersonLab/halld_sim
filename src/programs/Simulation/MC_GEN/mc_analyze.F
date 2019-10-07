c*******************************************************************************
c	MC_ANALYZE.  Part of MC_GEN
c	Analyze a single event.   Histogram any and all variables and
c	correlations of interest.  The parent routine, MC_GEN is strictly
c	and event generator.  This code analyzes these events one by one.
c
c	R. A. Schumacher,   Carnegie Mellon University
c	Original code was developed ~1991
c
c	Modifications:
c		10-28-91	- create this code from some older version
c		11-12-92	- reconfigure for f0 production
c		7-18-95		- bring up to date for HBOOK and PAW
c				- reconfigure back to CLAS hyperons
c				- more weak decay stuff
c		2-27-96		- additions to look at quasi-free Xi production
c		3-8-96		- additions for dE/dx calculations
c		8-12-97		- additions for Cx' and Cz' polarizations
c				- additions for CLAS spectrometer
c               8-18-97         - UNIX version
c                               - fine tuning to simulate CLAS triggers
c               5-21-98         - fix up for Linux running
c               5-30-00         - stuff for better hyperon polarization analysis
c               1-20-01         - add MC_CLAS routine to export events in TXT2BOS
c                                 compatible format for comparison with GSIM
c               2-5-01          - henceforth the instrumental smearing of the
c                                 tracks will have been done in mc_generator.
c                                 Various calls to PSMEAR have been removed.
c               2-26-01         - bug fix: definition of PHI was rotated 180deg.
c               5-23-01         - add stuff for pi+ photoproduction
c               6-14-01         - add stuff for pi0 photoproduction
c               6-28-01         - wire maps for drift chambers and DC knockout
c               7-6-01          - add trigger mask to simplify switching reactions
c               9-1-01          - turn off the TOF counter off-by-one blurring
c               9-8-01          - improve TOF counter matching algorithm
c               9-10-01         - introduce swimming of particles thru field
c               4-2-02          - fix polarization info: fix array bounds bug
c               7-27-02         - introduce E. Pasyuk eloss code into K+ reconstruction
c               10-07-04        - turn off tracking in CLAS
c               6-06-05         - add K+p missing mass off proton to see pi-
c               11-16-05        - mods to look at Lambda(1520)
c               1-18-06         - mods to look at Lambda(1405)
c               8-10-06         - turned off close-particle zapping
c               8-21-06         - modify run control such that we can 
c                                 skip all analysis by setting IDOCUTS to 0
c               10-17-06        - improve/expand x(1280) analysis
c               12-06-06        - fix bug in event weighting when combinatorics are present
c                                 array 'acceptlist' introduced
c               1-22-07         - add flag IWRITEIT to argument list to permit writing
c                                 to disk of selected events we want in the BOS file
c               12-17-07        - more code to analyze Lambda(1520)
c               1-28-08         - adapt code to compile with gfortran
c               2-02-08         - restructure to build in a GlueX acceptance function 
c               2-15-08         - building up code for analysis of some charm channels
c                               - disable obsolete wire-hit checking in CLAS
c               11-19-08        - histograms for Lambda(1520) tests changed 
c               3-18-09         - add histograms to look at Lambda(1405) lineshape
c               7-8-09          - histogram to examine eta' --> rho gamma
c               5-20-10         - kinematic tricks with the eta pi pi Dalitz plot
c               11-8-11         - fix acceptance function to handle reversed torus field
c               3-20-12         - new code to look at K0 Lambda in g14
c               5-25-12         - change return flag iaccept_CLAS to include momentum selection of track
c               3-24-14         - add code to look at N Delta dibaryon states in g13 and g14
c               9-02-14         - start recovery from disastrous loss of ALL /monte code for unknown 
c                                 reasons.  Last backup of source files was from about 4-10-14, and
c                                 quite a lot of coding was lost, as were some non-code ancillary files.
c               12-17-14        - add code to look at gamma d --> K0 Lambda with boost to c.m. frame
c                2-27-15        - look at g + n --> pi- p missing momentum resolution
c                4-29-15        - look at g + p --> K+ Lambda --> K+ p lepton neutrino  semi-leptonic decays
c                7-14-15        - write out event by event data for x-->eta pi pi
c                9-30-15        - code tweak to satisfy gfortran on the Mac
c               11-06-15        - revisit 4-body decays of the f(1285)      
c                1-27-16        - trace bug that is changing PLIST entries
c                3-16-16        - tweaks for f1(1285) analysis
c                6-20-16        - set up for p-pbar p detection
c                8-28-16        - more detailed calculations of p p-bar IM in the CM frame
c                6-26-17        - start of setup for Lambda-anti-Lambda analysis
c	         2-10-18        - start setup of Lambda p scattering analysis
c                6-14-18        - expand code for Lambda-anti-Lambda analysis
c                8-24-18	- resurrect the K0Y0 code
c                2-01-19        - adapt the code to treat Sigma0-anti-Sigma0, too
c	         7-01-19        - morph to keep up with developing MC_GEN code
c                7-02-19        - huge amounts of analysis code for old projects was removed       
c
c	Coding for the list of particles generated in the decay chain:
c	All variables in this list are measured in the laboratory frame.
c	
c	PLIST(I,MARK)	
c	I = index into list of a given particle's parameters
c	MARK = counts the particles produced
c	I = 
c	1	PX
c	2	PY
c	3	PZ
c	4	MASS
c	5	TOTAL ENERGY
c	6	LUND ID FOR THIS PARTICLE
c	7	X OF CREATION
c	8	Y OF CREATION
c	9	Z OF CREATION
c	10	X COMPONENT OF POLARIZATION
c	11	Y COMPONENT OF POLARIZATION
c	12	Z COMPONENT OF POLARIZATION
c	13	PX AT DECAY POINT
c	14	PY AT DECAY POINT
c	15	PZ AT DECAY POINT
c       16      MASS (same as (4))
c	17	TOTAL ENERGY AT DECAY POINT
c	18	X OF DECAY
c	19	Y OF DECAY
c	20	Z OF DECAY
c	21	CTAU OF PARTICLE
c       22      Charge of the particle
c       23      spare
c       24      spare
c
c       Each particle must have an associated particle identification (ID) code.
c       The very original version of this code was based on the LUND numerology.
c       However, this coding is totally unrelated to the GEANT (PDG) coding.
c       Hence, this code effectively has its own particle ID suite.
c       The user can make up new particle ID ad libitum in the definition files.
c       The ordering of particles is somewhat haphazard, reflecting the history
c       of how the code was developed.
c       
c	Old LUND particle ID's 
c	1	gamma
c	2	neutrino
c	7	electron
c	9	muon
c	17	pi+
c	18	K+
c	19	K0
c	23	pi0
c	24	eta
c	27	rho+
c	33	rho0
c	34	omega
c	35	phi
c	38	K_long
c	57	Lambda
c	41	proton
c	42	neutron
c	43	Sigma+
c	44	Sigma0
c	45	Sigma-
c	46	Cascade-
c	47	Cascade0
c	61	Delta(1232)++
c       Additional values used by various definition files used here.
c	20	K*0(892)
c	21	K*+(892)
c	48	Sigma0(1385)
c	49	Cascade-(1530)
c	50	3He
c	51	L3H
c	52	Deuteron
c	53      Carbon_12
c	54      Carbon_11
c	55	Boron 11
c	56	Sigma+(1385)
c	58	Lambda(1405)
c	59	Lambda(1520)
c	69	f0
c	70      D0
c       71      D+
c	74	x(1285)		
c	75	a0+(980)
c   	76	N(1535) S11
c	77	a2(1320)
c       80      Lambda_c +
c       90      J/psi
c       91      Baryonium
c       92      Diproton
c       93      Lambda anti-Lambda state
c       94      Lambda anti-Lambda state
c       95      proton anti-Lambda state
c       98      "Fat Meson"
c       99      "Fat Baryon"
c      100      "Box-like" mass distribution
c
c
c	Top of active code
c  
	subroutine mc_analyze(plist,imark)
	parameter (nwpawc = 20000000)
	parameter (maxnumlis=50)  !Maximum number of particles in the decay list
	parameter (nkbinpol = 10)    !for defining kaon angle binning FOR polarizations
	parameter (nebinpol = 34)    !for defining photon energy scale
	parameter (ebinminpol = 0.8) !GeV
	parameter (ebinmaxpol = 2.5) !GeV
c	parameter (nkbin = 20)    !for defining kaon angle binning   FOR g1C
c	parameter (nebin = 80)    !for defining photon energy scale
c	parameter (ebinmin = 0.4875) !GeV
c	parameter (ebinmax = 2.4875) !GeV
c	parameter (nkbin = 18)     !for defining kaon angle binning   FOR g1C
c	parameter (nebin = 100)    !for defining photon energy scale
c	parameter (ebinmin = 0.4875) !GeV
c	parameter (ebinmax = 2.9875) !GeV
c	parameter (cbinmin = -0.85)  !GeV
c	parameter (cbinmax = +0.95)  !GeV
	parameter (nkbin = 10)     !for defining kaon angle binning   FOR g11
c	parameter (nebin = 15)     !for defining photon energy scale
	parameter (nebin = 60)     !for defining photon energy scale
	parameter (ebinmin = 1.200) !GeV
c	parameter (ebinmin = 0.350) !GeV
	parameter (ebinmax = 4.000) !GeV
c	parameter (ebinmax = 2.500) !GeV
	parameter (cbinmin = -0.85)  !GeV
	parameter (cbinmax = +0.95)  !GeV
	parameter (rmneut = 0.93956536)   !GeV neutron mass
	parameter (rmmuon = 0.105658)     !GeV muon mass
	parameter (rmelectron = 0.000511) !GeV electron mass
	character*256 prg_name, histnm, filenm, prefix
	real plist(24,maxnumlis)
	real nhatpol(4),nbeam(4)
	real zhatprime(4),xhatprime(4)
	real zhprime(4),xhprime(4)
	real beta(3),ptemp(5),ptemp2(5),ptemp3(5),pcorout(3),pcorin(3)
	real betapair(3),betatophyp(3),betatoahyp(3)
	real psigp(4),psigm(4)
	real psystem(5)
	real pimp1ap(5), pimp2ap(5), pimp1apcm(5), pimp2apcm(5), pimp1p2(5)
	real ptemp5(5),ptemp5b(5),ptemp5c(5),pxmeson(5),pxcmpip(5),pxcmpim(5)
	real ptemp6(5),ptemp7(5),pK01(5),pK02(5),pK01cm(5),pK02cm(5)
	real plambda1(5),plambdas(5),plambda2(5),plambda1cm(5),plambda2cm(5)
	real pspect(5),pinitial(5)
c	real ppione(5),ppitwo(5)
	real plam(5),plam2(5),pkaon(5),pik(5),pproton(5),pneutrino(5)
	real pkaoncm(5),pprotoncm(5),pleptoncm(5)
	real pkaonhyp(5),pprotonhyp(5),pleptonhyp(5)
	real ppiminuscm(5),ppiminushyp(5),ppiminushyp2(5)
	real pnothing(5),pnothingcm(5),pnothinghyp(5),pnothinghyp2(5)
	real pnothings(5),pnothingcms(5),pnothinghyps(5),pnothinghyp2s(5)
	real pkaonhyp2(5),pprotonhyp2(5),pleptonhyp2(5)
	real plambda(5),plambdacm(5),plambdacms(5),pneutrinocm(5)
	real plambdahyp(5),pneutrinohyp(5),pneutrinos(5)
	real plambdahyp2(5),pneutrinohyp2(5)
	real phyperon(5),betahyp(3),betahyp2(3),betahyps(3),pneut(5),pneuts(5)
	real plambdahyps(5),pprotonhyps(5),ppiminushyps(5)
	real pkgen(5),plgen(5)
	real vertex(3),rvertex(3),dvertex(3),avertex(3),vv(3)
	real alampath(3),plampath(3),dlampath(3)
	real precoil(5),pa0p(5),pa0m(5)
	real precval(5),pxcmeta(5),ppipi(5),ppipipip(5),ppipipim(5)
	real pbeamlab(5),pxcmbeam(5)
	real acceptlist(maxnumlis)
	real dirmeson(5)
	real ptotrest(5),pollamlam(5),nhatdircm(5)
	real p1cm(5),p1bb(5),p2cm(5),p2bb(5),apcm(5),apbb(5),pbeamcm(5),pbeambb(5)
	real aplab(5)
	real plamcm(5),alamcm(5),pimcm(5),pipcm(5),plamalam(5)
	real plamlab(5),alamlab(5),pimlab(5),piplab(5),plamalamlab(5)
	real ppplam(5),ppalam(5),plampair(5),alampair(5),p2hyp(5)
	real pphyp(5),aphyp(5)
	real p1pair(5),p2pair(5)
	real xhatphyp(5),yhatphyp(5),zhatphyp(5)  !need 5 entries for 3 vectors to avoid cross product issue
	real xhatahyp(5),yhatahyp(5),zhatahyp(5)
	real corrsum(3,5),correrr(3,5)
	
	logical ifirst,ifirstprocess,iacc_piplus,iacc_piminus,iacc_kplus,
     1      iacc_kminus,
     1      iacc_kstar,iacc_phi,iacc_l3h,iacc_proton,iacc_firstk,
     1      istop_cas,iacc_photon,iacc_photon2,
     1      iacc_electron,iacc_positron,
     1      ihackflag,iacc_kplus2,iwinner,iacc,
     1      iacc_firstp,iacc_firstpip,iacc_firstpim,iacc_firstkm,
     1      iacc_firstphoton, 
     1      iacc_piplus2,iacc_proton2,iacc_piminus2,iacc_kminus2,
     1      iacc_neutral,iacc_neutron,iacc_deuteron,iacc_antiproton
c
	logical iaccept_GlueX
	logical iexpt_trigger,ifiddle
c	data twopi   /6.283185308/
	data piover2 /1.570796327/
	data rmK0  /0.497614/   !GeV/c^2
	data rmkaon/0.493677/   !GeV/c^2
	data rmpion/0.13957/  !GeV/c^2
	data rmasslambda,rmasssigma/1.11568,1.19255/
	data rmneut,rmprot/0.93956536,0.9382703/
	data clight/29.9792458/ !speed of light in cm/nsec
	data rade /57.29577951/
	data ifirst/.true./
	data ifirstprocess/.true./
c	data rmanot  / 0.9847/ !mass of the a_0 +- (from PDG)
	data rmanot  / 0.960 / !mass of the a_0 +-
c	data gammaa0 / 0.1000/ !average width of the a_0 +-
c	data gammaa0 / 0.0750/ !average width of the a_0 +-
	data gammaa0 / 0.0500/ !average width of the a_0 +-
	common /pawc/hmemor(nwpawc)
	common /stuff/poldir(4),xx,icount,ikeep,ihackflag,iattempt
        common /fermi/prec,qfmass,pfermi,ndynamics
	common /beam/pbeam,ppbeam,rmbeam,rmtarg,jcode,xsigma,ysigma,ncode
        common /analysis/dpoverp_in,dpoverp_out,ptotrest
	common /run_control/n_arg,prg_name,filenm,histnm,prefix,
     1                      irestart,lfl,lprefix,maxbosoutput
	common /detector/acceptweight
	common /misc/cthcm_gen,itoken !for tests and hacks
	common /monte/ iseed,xflat,xnormal,xbw
	common /correlations/corrsum,correrr
c     1     polsumplamX,polsumplamY,polsumplamZ,
c     1     polsumalamX,polsumalamY,polsumalamZ,
c     1	   corrsumXX,corrsumYY,corrsumZZ,corrsumZY,corrsumYZ,
c     1	   corrsumXY,corrsumYX,corrsumXZ,corrsumZX
c
c	write(6,*)'Entering analyzer',icount,ikeep
c	if(ifirst.eqv..true.)then
	if(ifirst)then
	   ifirst = .false.
	   write(6,*)'First call to event analyzer'
	   write(6,*)'     Version 10.01'
	   call flush(6)
c       
	   call hlimit(nwpawc)
c       
	   if(irestart.eq.1)then
	      write(6,*)'Opening existing histogram file for update.'
	      write(6,*)'Histograms from file:', histnm(1:lfl)
	      call hropen(1,'MCDIR',histnm(1:lfl),'U',1024,istat)
	      call hcdir('//PAWC',' ') !Memory
	      call hmdir('RAWEVENT',' ')	
	      call hmdir('ACCEPTED',' ')
	      call hmdir('WEIGHTED',' ')
	      call hcdir('//PAWC/RAWEVENT',' ')
	      call hcdir('//MCDIR/RAWEVENT',' ')
	      call hrin(0,99999,0)
	      call hcdir('//PAWC/ACCEPTED',' ')
	      call hcdir('//MCDIR/ACCEPTED',' ')
	      call hrin(0,0,0)
	      call hcdir('//PAWC/WEIGHTED',' ')
	      call hcdir('//MCDIR/WEIGHTED',' ')
	      call hrin(0,0,0)
c	      call hldir('//PAWC','T')   !prints out all the titles...
c	      call hldir('//MCDIR','T')
	   else   
           do 20 iloop=1,3
	      if(iloop.eq.1)then
		 write(6,*)'Opening new or overwritten histogram file.'
		 write(6,*)'Histograms to file:', histnm(1:lfl)
		 call hropen(1,'MCDIR',histnm(1:lfl),
     1  	      'N',1024,istat)
		 write(6,*)'Booking RAW histograms...'
		 call hmdir('RAWEVENT',' ') !Disk
		 call hmdir('ACCEPTED',' ')
		 call hmdir('WEIGHTED','S')
		 call hrend('MCDIR') 
		 call hcdir('//PAWC',' ') !Memory
		 call hmdir('RAWEVENT',' ')	
		 call hmdir('ACCEPTED',' ')
		 call hmdir('WEIGHTED','S')
		 call hcdir('//PAWC/RAWEVENT',' ')
	      elseif(iloop.eq.2)then
		 write(6,*)'Booking ACCEPTED histograms...'
		 call hcdir('//PAWC/ACCEPTED',' ')
	      else
		 write(6,*)'Booking WEIGHTED histograms...'
		 call hcdir('//PAWC/WEIGHTED',' ')
	      end if
c
c--------------------------------------------------------------------------
c             A huge number of histograms were deleted from this place.
c             Streamlined the program for GlueX era.
c--------------------------------------------------------------------------
c
	      call hbook1(28,'Beam Momentum (GeV/c)' , 100,  3.5, 12.0,0.)	      
c
c------ Histograms for GlueX baryonium study (p pbar)--------------
c
	      call hbook1(131,'P1 Lab Momentum (GeV/c)'     ,100,0.,8.0,0.)
	      call hbook1(132,'P2 Lab Momentum (GeV/c)'     ,100,0.,8.0,0.)
	      call hbook1(133,'p&^-! Lab Momentum (GeV/c)'  ,100,0.,8.0,0.)
	      call hbook2(134,'P2 vs. P1 LAB Momenta (GeV/c)',
     1	                            50,0.,8.0,50,0.,8.0,0.)
	      call hbook1(135,'P1 Lab Angle (deg)'     ,90,0.,90.0,0.)
	      call hbook1(136,'P2 Lab Angle (deg)'     ,90,0.,90.0,0.)
	      call hbook1(137,'p&^-!  Lab Angle (deg)'     ,90,0.,90.0,0.)
c       

	      call hbook2(209,'PROTON  [Q] (deg) vs. p (GeV/c)'       
     1  					,100,0.,8.,100,0.,90.,0.)
	      call hbook2(230,'ANTI-PROTON  [Q] (deg) vs. p (GeV/c)'       
     1  					,100,0.,8.,100,0.,90.,0.)
              
	      call hbook2(8850,'IM (p2 p&^-!) vs IM (p1 p2)'
     1                                ,100,1.6,4.0,100,1.6,4.0,0)
	      call hbook1(8851, 'IM (p1 p&^-!) (GeV/c^2!) '
     1  						 ,100,1.8,4.0,0.)
	      call hbook1(8852, 'IM (p2 p&^-!) (GeV/c^2!) '
     1  						 ,100,1.8,4.0,0.)
	      call hbook1(8871, 'IM (p1 p2) (GeV/c^2!) '
     1  						 ,100,1.8,4.0,0.)
	      call hbook2(8853,'IM (p2 p&^-!) vs IM (p1 p&^-!)'
     1                                ,100,1.6,4.0,100,1.6,4.0,0)
	      call hbook2(8854,'IM^2! (p2 p&^-!) vs IM^2! (p1 p&^-!)' !Dalitz
     1                                ,100,3.0,11.0,100,3.0,11.0,0)
	      call hbook2(8855,'Longitudinal Momentum (van Hove) Plot'
     1                                ,100,-2.5,2.5,100,-2.5,2.5,0)
	      call hbook2(8856,'Proton Momenta P1 vs. P2 in c.m. (GeV/c)'
     1                                ,100,0.0,2.0,100,0.0,2.0,0)
	      call hbook1(8857, 'p?p1! (c.m.)  (GeV/c)'
     1  						 ,100,0.0,2.0,0.)
	      call hbook1(8858, 'p?p2! (c.m.)  (GeV/c)'
     1  						 ,100,0.0,2.0,0.)
	      call hbook2(8859,'IM^2! (p2 p&^-!) vs IM^2! (p1 p&^-!)' !Munged Dalitz
     1                                ,50,3.0,11.0,50,3.0,11.0,0)
	      call hbook1(8860, 'p?p1! (c.m.)  (RESCALED) (GeV/c)'
     1  						 ,100,0.0,2.0,0.)
	      call hbook1(8861, 'p?p1! (c.m.)  (GeV/c) (RESCALED)'
     1  						 ,100,0.0,2.0,0.)
	      call hbook1(8862, 'p?p2! (c.m.)  (GeV/c) (RESCALED)'
     1  						 ,100,0.0,2.0,0.)
	      call hbook2(8863,'IM^2! (p2 p&^-!) vs p?p1! RESCALED ' !Munged Combo plot
     1                                ,50,0.0,2.0,50,3.0,11.0,0)
	      call hbook1(8864, 'p?pbar! (c.m.)  (GeV/c)'
     1  						 ,100,0.0,2.0,0.)
	      call hbook1(8865,'p1 cos [Q]?cm! ',200,-1.,1.0,0.)
	      call hbook1(8866,'p2 cos [Q]?cm! ',200,-1.,1.0,0.)
	      call hbook1(8867,'p&^-! cos [Q]?cm! ',200,-1.,1.0,0.)
	      call hbook2(8868,'cos [Q]?cm! (p1 vs. p2)',
	1    	                       50,-1.,1.0,50,-1.,1.0,0.)
	      call hbook2(8869,'cos [Q]?cm! (pbar vs. p1)',
	1    	                       50,-1.,1.0,50,-1.,1.0,0.)
	      call hbook2(8870,'cos [Q]?cm! (pbar vs. p2)',
	1    	                       50,-1.,1.0,50,-1.,1.0,0.)
	      call hbook1(8875,'P1 [Q]?cm! ',90,0.,180.0,0.)
	      call hbook1(8876,'P2 [Q]?cm! ',90,0.,180.0,0.)
	      call hbook1(8877,'p&^-! [Q]?cm! ',90,0.,180.0,0.)
	      call hbook2(8878,'[Q]?cm! (p1 vs. p2)',
	1    	                       90,0.,180.0,90,0.,180.0,0.)
	      call hbook2(8879,'[Q]?cm! (pbar vs. p1)',
	1	                       90,0.,180.0,90,0.,180.0,0.)
	      call hbook2(8880,'[Q]?cm! (pbar vs. p2)',
	1    	                       90,0.,180.0,90,0.,180.0,0.)

	      call hbook2(8881,'[Q]?cm! (deg) vs. p?cm!   Proton1',
	1    	                       100,0.0,2.0,90,0,180.,0.)
	      call hbook2(8882,'[Q]?cm! (deg) vs. p?cm!   Proton2',
	1    	                       100,0.0,2.0,90,0,180.,0.)
	      call hbook2(8883,'[Q]?cm! (deg) vs. p?cm!  Antiproton',
	1    	                       100,0.0,2.0,90,0,180.,0.)
	      call hbook1(8884,'[Q]?cm! IM(p1 pbar)',90,0.,180.0,0.)
	      call hbook1(8885,'[Q]?cm! IM(p2 pbar)',90,0.,180.0,0.)
	      call hbook1(8886,'cos[Q]?bb! Anti-proton (Adair)' ,100,-1.0,1.0,0.)
	      call hbook1(8887,'cos[Q]?bb! Anti-proton (Hel.)',100,-1.0,1.0,0.)
	      call hbook1(8888,'cos[Q]?bb! Anti-proton (GJ) ',100,-1.0,1.0,0.)
	      call hbook1(8889,'cos [Q]?cm! IM(p1 pbar)',100,-1.0,1.0,0.)
	      call hbook1(8890,'cos [Q]?cm! IM(p2 pbar)',100,-1.0,1.0,0.)
	      call hbook1(8891,'cos [Q]?cm! p1  ',100,-1.0,+1.0,0.0)
	      call hbook1(8892,'cos [Q]?cm! p2  ',100,-1.0,+1.0,0.0)
	      call hbook1(8893,'cos [Q]?cm! pbar',100,-1.0,+1.0,0.0)
	      call hbook2(8894,'cos[Q]?cm! (p2 vs. p1) BEFORE SWAPPING',
	1	                           50,-1.0,1.0,50,-1.0,1.0,0.)
	      call hbook1(8895, '-(t) (GeV^2!)  p&^-! p',
	1	                         100,0.0,3.0,0.)
	      call hbook1(8896, '-(t - tmin) (GeV^2!)  p&^-! p',
	1	                         100,0.0,3.0,0.)
	      call hbook2(8897,'Beam Energy (GeV) vs. IM(p?2! p&^-!)',
	1	   220,1.8,4.0,90,3.0,12.0,0.)
	      call hbook2(8898,'-(t - tmin) (GeV^2!)  p&^-! p vs. IM(p?2! p&^-!)',
	1	   110,1.8,4.0,50,0.0,3.0,0.)
	      call hbook2(8899,'-t (GeV^2!)  p&^-! p vs. IM(p?2! p&^-!)',
	1	   110,1.8,4.0,50,0.0,3.0,0.)
	      call hbook2(8900,'-t (GeV^2!) vs. p&^-! cos [Q]?cm!',
	1	   50,-1.,+1.0,50,0.0,3.0,0.)
	      call hbook2(8901,'Beam Energy (GeV) vs. p&^-! cos [Q]?cm!',
	1	   50,-1.,+1.0,90,3.0,12.0,0.)
c 
c------Histograms for GlueX baryonium study (Lambda anti-Lambda)--------------
c

	      call hbook1(601,'X Vertex Distribution',100,-3.,+3.,0.)
	      call hbook1(602,'Y Vertex Distribution',100,-3.,+3.,0.)
	      call hbook1(603,'Z Vertex Distribution',700,-50.,+20.,0.)

	      call hbook1(604,'X Primary Distribution (cm)',100,-3.,+3.,0.)
	      call hbook1(605,'Y Primary Distribution',100,-3.,+3.,0.)
	      call hbook1(606,'Z Primary Distribution',100, 0.0,200.,0.)
	      call hbook1(607,'X [L] Distribution (cm)',100,-3.,+3.,0.)
	      call hbook1(608,'Y [L] Distribution',100,-3.,+3.,0.)
	      call hbook1(609,'Z [L] Distribution',100,40.,200.,0.)
	      call hbook1(610,'X [L]&^-! Distribution (cm)',100,-3.,+3.,0.)
	      call hbook1(611,'Y [L]&^-! Distribution',100,-3.,+3.,0.)
	      call hbook1(612,'Z [L]&^-! Distribution',100,40.,200.,0.)
	      call hbook1(613,'[L] Lab Pathlength (cm)'     ,100,0.,100.,0.)
	      call hbook1(614,'[L]&^-! Lab Pathlength (cm)',100,0.,100.,0.)
	      call hbook1(615,'[L] to [L]&^-! Lab Distance (cm)',100,0.,100.,0.)
	      call hbook1(616,'[L] Rest Lifetime (ns)'    ,200,0.,2.0,0.)
	      call hbook1(617,'[L]&^-! Rest Lifetime (ns)',200,0.,2.0,0.)
	      call hbook2(618,'[L] vs. [L]&^-! Rest Lifetimes (ns)',
     1		                 30,0.,1.0,30,0.,1.0,0.)

	      call hbook1(9601,'cos[Q]?Y! of p in [L] frm',
     1                                           50, -1.0, 1.0, 0.)
	      call hbook1(9602,'cos[Q]?Y! of p&^-! in [L]&^-! frm.',
     1                                           50, -1.0, 1.0, 0.)
	      call hbook1(9603,'cos[Q] of [L] in pair frm.',
     1	                                         50, -1.0, 1.0, 0.)
              call hbook1(9604,'cos[Q] of [L]&^-! in pair frm.',
     1                                           50, -1.0, 1.0, 0.)
              call hbook1(9605,'cos[Q]?X! of p in [L] frm.',
     1                                           50, -1.0, 1.0, 0.)
              call hbook1(9606,'cos[Q]?Y! of p in [L] frm.',
     1                                           50, -1.0, 1.0, 0.)
              call hbook1(9607,'cos[Q]?Z! of p in [L] frm.',
     1                                           50, -1.0, 1.0, 0.)
              call hbook1(9608,'cos[Q]?X! of p&^-! in [L]&^-! frm.',
     1                                           50, -1.0, 1.0, 0.)
              call hbook1(9609,'cos[Q]?Y! of p&^-! in [L]&^-! frm.',
     1                                           50, -1.0, 1.0, 0.)
              call hbook1(9610,'cos[Q]?Z! of p&^-! in [L]&^-! frm.',
     1                                           50, -1.0, 1.0, 0.)
              call hbook1(9611,'Momentum in [L][L]&^-! pair rest frame (GeV/c)',
     1                                           100, 0.0, 1.5, 0.)
              call hbook1(9612,'p [L]&^-! in pair rest frame (GeV/c)',
     1                                           100, 0.0, 1.5, 0.)
	      call hbook1(9613,'cos[Q] of [L] in pair frm.',
     1	                                         50, -1.0, 1.0, 0.)
              call hbook1(9614,'cos[Q] of [L]&^-! in pair frm.',
     1                                           50, -1.0, 1.0, 0.)
              call hbook1(9615,'Momentum in [L][L]&^-! pair rest frame (GeV/c)',
     1                                           100, 0.0, 1.5, 0.)
              call hbook1(9616,'Mom. [L][L]&^-! (GeV/c) (q-weighted)',
     1                                           100, 0.0, 1.5, 0.)
              call hbook1(9617,
     1         'Mom. [L][L]&^-! (GeV/c) (q and k weighted)',
     1                                           100, 0.0, 1.5, 0.)
	      
	      call hbook1(9801,'Lambda IM (p [p]^-!) ',
	1	                         100,1.05,1.15,0)
	      call hbook1(9802,'Anti-Lambda IM (p&^-! [p]^+!) ',
	1	                         100,1.05,1.15,0)
c	      call hbook1(9801,'Lambda IM (p [p]^-!) ',
c	1	                         100,1.05,1.50,0)
c	      call hbook1(9802,'Anti-Lambda IM (p&^-! [p]^+!) ',
c	1	                         100,1.05,1.50,0)

	      call hbook2(9850,'IM (p [L]&^-!) vs IM (p [L])'
	1	                            ,100,2.0,4.0,100,2.0,4.0,0)
	      call hbook2(9855,'IM (p [L]&^-!) vs IM ([L]&^-! [L])'
	1	                            ,100,2.0,4.0,100,2.0,4.0,0)
	      call hbook2(9854,'IM^2! (p [L]&^-!) vs IM^2! (p [L]) (Dalitz Plot)' !Dalitz
	1       	                    ,100,3.0,11.0,100,3.0,11.0,0)

              call hbook1(9851, 'IM ([L]&^-![L]) (GeV/c^2!) ',
	1	                         100,2.0,4.0,0.)
	      call hbook1(9852, 'IM ([L]&^-!p) (GeV/c^2!) ',
	1	                         100,2.0,4.0,0.)
              call hbook1(9853, 'IM ([L]&^-![L]) (GeV/c^2!) (q-weighted) ',
	1	                         100,2.0,4.0,0.)
              call hbook1(9858, 'IM ([L]&^-![L]) (GeV/c^2!) (q and k weighted)',
	1	                         100,2.0,4.0,0.)
	      call hbook1(9856, '-(t - tmin) (GeV^2!)  [L]&^-! [L]',
	1	                         100,0.0,3.0,0.)
	      call hbook1(9857, '-(t - tmin) (GeV^2!)  [L]&^-! p',
	1	                         100,0.0,3.0,0.)

	      call hbook2(9869,'cos [Q]?cm! anti-Lambda vs. Lambda',
	1    	                       50,-1.,1.0,50,-1.,1.0,0.)
              
	      call hbook2(9870,'Long. Mom. Standard (van Hove Plot)',
	1    	                         100,-2.5,2.5,100,-2.5,2.5,0)
	      call hbook2(9871,'Long. Mom. Inv. Mass (Schu Plot)',
	1    	                         100,-3.5,3.5,100,-3.5,3.5,0)
	      call hbook2(9872,'Long. Mom. Kin. Enrgy. (Schu2 Plot)',
	1    	                         100,-2.5,2.5,100,-2.5,2.5,0)
	      call hbook2(9873,'Long. Mom. Trans. Mom. (Schu3 Plot)',
	1    	                         100,-2.5,2.5,100,-2.5,2.5,0)
	      call hbook2(9874,'I.M. (by sector pair) vs. [w] (degrees)',
	1    	                         100,0.0,360,50,2.0,3.5,0)

              
	      call hbook2(9880,'p?cm! (GeV/c) anti-Lambda vs. Lambda',
	1    	                         100,0.0,2.0,100,0,2.,0.)

	      call hbook2(9881,'cos [Q]?cm!  vs. p?cm!   Lambda',
	1    	                         100,0.0,2.0,100,-1.0,+1.0,0.)
	      call hbook2(9883,'cos [Q]?cm!  vs. p?cm!   anti-Lambda',
	1    	                         100,0.0,2.0,100,-1.0,+1.0,0.)
	      call hbook2(9884,'cos [Q]?cm! vs. p?cm! ([L]&^-! [L] Pair)',
	1    	                         100,0.0,2.0,100,-1.0,+1.0,0.)

	      
		
	      call hbook1(9885,'cos [Q]?cm! Lambda'     ,100 ,-1.0,1.0,0.)
	      call hbook1(9886,'cos [Q]?cm! Anti-Lambda',100 ,-1.0,1.0,0.)
	      call hbook1(9887,'Overall Missing Mass (GeV)',100 ,-1.E-1,+5.E-2,0.)
	      call hbook1(9888,'cos [Q]?cm! Between Planes',100 ,-1.0,1.0,0.)
	      call hbook1(9889,'[Q]?cm! Between Planes',18 ,0.0,180.0,0.)
	      
	      call hbook1(9620,'Overall Missing Energy (GeV)',100,-0.5,0.5,0.)
	      call hbook2(9621,'Missing Energy vs. Missing Mass',
     1	                           50,-1.E-1,+5.E-2, 50,-0.5, 0.5,0)
c	      
c	      For Lambda anti-Lambda computationns
c
	      write(6,*)'Initializing correlation sums...'
c	      polsumplamX = 0
c	      polsumplamY = 0
c	      polsumplamZ = 0
c	      polsumalamX = 0
c	      polsumalamY = 0
c	      polsumalamZ = 0
c	      corrsumXX  = 0
c	      corrsumYY  = 0
c	      corrsumZZ  = 0
c	      corrsumZX  = 0
c	      corrsumXZ  = 0
c	      corrsumXY  = 0
c	      corrsumYX  = 0
c	      corrsumZY  = 0
c	      corrsumYZ  = 0
	      do ii=1,3
		 do jj=1,5
		    corrsum(ii,jj) = 0
		    correrr(ii,jj) = 0
		 enddo
	      enddo
	      call hbook2(9890,'Correlation Matrix [L]&^-! vs [L]', 
     1	                      5,0.5,5.5 ,3,0.5,3.5,0.)
	      
c             Lab momentum and angle for Lambda-anti-Lambda
	      call hbook1(9901,'Mom p1 (GeV) (lab)'    ,100,0.0,2.0,0.)
	      call hbook1(9902,'Mom p2 (GeV) (lab)'    ,100,0.0,8.0,0.)
	      call hbook1(9903,'Mom p-bar (GeV) (lab)' ,100,0.0,8.0,0.)
	      call hbook1(9904,'Mom [p]^+! (GeV) (lab)',100,0.0,2.0,0.)
	      call hbook1(9905,'Mom [p]^-! (GeV) (lab)',100,0.0,2.0,0.)
	      call hbook1(9906,'cos[Q] (deg) p1 (lab)          ',100,0.5,1.0,0.)
	      call hbook1(9907,'cos[Q] (deg) p2 (lab)          ',100,0.5,1.0,0.)
	      call hbook1(9908,'cos[Q] (deg) p-bar (lab)       ',100,0.5,1.0,0.)
	      call hbook1(9909,'cos[Q] (deg) [p]^+! (lab)      ',100,0.5,1.0,0.)
	      call hbook1(9910,'cos[Q] (deg) [p]^-! (lab)      ',100,0.5,1.0,0.)
c             CM momentum and angle for Lambda-anti-Lambda
	      call hbook1(9911,'Mom p1 (GeV) (cm)'    ,100,0.0,2.0,0.)
	      call hbook1(9912,'Mom p2 (GeV) (cm)'    ,100,0.0,2.0,0.)
	      call hbook1(9913,'Mom p-bar (GeV) (cm)' ,100,0.0,2.0,0.)
	      call hbook1(9914,'Mom [p]^+! (GeV) (cm)',100,0.0,1.0,0.)
	      call hbook1(9915,'Mom [p]^-! (GeV) (cm)',100,0.0,1.0,0.)
	      call hbook1(9916,'cos[Q] (deg) p1 (cm)      ',100,-1.,1.0,0.)
	      call hbook1(9917,'cos[Q] (deg) p2 (cm)      ',100,-1.,1.0,0.)
	      call hbook1(9918,'cos[Q] (deg) p-bar (cm)   ',100,-1.,1.0,0.)
	      call hbook1(9919,'cos[Q] (deg) [p]^+! (cm)  ',100,-1.,1.0,0.)
	      call hbook1(9920,'cos[Q] (deg) [p]^-! (cm)  ',100,-1.,1.0,0.)
c       
c------Histograms for GlueX Lambda scattering study -------------
c
c	      call hbook1(10887,'Missing Mass^2! using K^+! mass (GeV)',
cc     1	                        100,-0.02,0.02,0.)
c     1	                        100,-3.0,3.0,0.)
c	      call hbook1(10888,'Missing Mass^2! using [p]^+! mass (GeV)',
cc     1	                        100,-0.02,0.02,0.)
c     1	                        100,-3.0,3.0,0.)
c	      call hbook2(10889,
c     1	            'Missing Mass^2! using [p]^+! vs. K^+! mass (GeV)',
c     1	                        100,-3.0,3.0,100,-3.0,3.0,0.)
c	      call hbook1(10890,'Missing Mass using K^+! mass (GeV)',
c     1	                        100,-0.5,0.5,0.)
c	      call hbook1(10100,'IM(p pi-) ([L])'     ,100 ,1.0,1.2,0.)
c	      call hbook1(10100,'IM(p pi-) ([L])'     ,200 ,1.0,1.6,0.)

 20	   continue
	end if    !histogram file opening
c
	end if    !first call to analyzer
c
c	Top of loop to decide whether this event should be histogrammed at all
c	Do the experiment specific cuts, or skip to line 60
c	for full phase-space acceptance
c
	icount = icount + 1   !Counts how many events are examined
c	write(6,*)'Call to analyzer',icount
c	call prevent(plist,imark+1,ikeep)   !temporary	
c
	iacc_piplus      = .false.
	iacc_piplus2     = .false.
	iacc_piminus     = .false.
	iacc_piminus2    = .false.
	iacc_kplus       = .false.
	iacc_kplus2      = .false.
	iacc_kminus      = .false.
	iacc_kminus2     = .false.
	iacc_kstar       = .false.
	iacc_phi         = .false.
	iacc_l3h         = .false.
	iacc_proton      = .false.
	iacc_proton2     = .false.
	iacc_antiproton  = .false.
	iacc_firstk      = .false.
	iacc_firstkm     = .false.
	iacc_firstpip    = .false.
	iacc_firstpim    = .false.
	iacc_firstp      = .false.
	iacc_firstphoton = .false.
	iacc_electron    = .false.
	iacc_positron    = .false.
	iacc_neutral     = .false.
	iacc_neutron     = .false.
	iacc_deuteron    = .false.
	istop_cas        = .false.
	iexpt_trigger    = .false.
	acceptweight_deuteron   = 0.0
	acceptweight_electron   = 0.0
	acceptweight_eta        = 0.0
	acceptweight_kaon       = 0.0
	acceptweight_kaon2      = 0.0
	acceptweight_neutron    = 0.0
	acceptweight_neutral    = 0.0
	acceptweight_muon       = 0.0
	acceptweight_piplus     = 0.0
	acceptweight_piplus2    = 0.0
	acceptweight_piminus    = 0.0
	acceptweight_piminus2   = 0.0
	acceptweight_proton     = 0.0
	acceptweight_proton2    = 0.0
	acceptweight_antiproton = 0.0
	acceptweight_photon     = 0.0
	acceptweight_photon1    = 0.0
	acceptweight_photon2    = 0.0
	weight       = 0.0
c
	if(ihackflag)then
	   write(6,*)'> Hack masses for tests....'
c          Turn two pions into kaons for events that have only pions...
c	   ipipcount = 0
c	   ipimcount = 0
c	   iskipit = 0
c       ... lots of code cut here
	endif
c
	pxtot = 0.
	pytot = 0.
	pztot = 0.
	acceptlist(1) = -1 !just to initialize value; no meaning intended here
c	call prevent(plist,imark+1,ikeep)
c
c       Do acceptance checks for each of the tracks and for the whole event
c
	do 50 i = 2,imark
	lund = plist(6,i)
	icharge = plist(22,i)
c
c	Get lab angles, in degrees, and other particle parameters
c
c       Take values from formation time
c	pmom = sqrt(plist(5,i)**2 - plist(4,i)**2) * 1000.	!MeV/c
c	px = plist(1,i) * 1000.	!MeV/c
c	py = plist(2,i) * 1000.
c	pz = plist(3,i) * 1000.
c
c       Some instrumental smearing was also done in mc_gen
c
	px = plist(13,i) * 1000. !MeV/c
	py = plist(14,i) * 1000.
	pz = plist(15,i) * 1000.
	pmom = sqrt(px*px +py*py +pz*pz) !MeV/c
	if(pmom.gt.0)then
		arg=pz/pmom
		if(arg.gt. 1.0)arg= 1.0
		if(arg.lt.-1.0)arg=-1.0
		theta_xz = atan2(px,pz) * rade
		theta_yz = atan2(py,pz) * rade
		theta_polar = acos(arg) * rade
		phi = atan2(py,px) * rade
		if(phi.lt.0.0)phi = phi + 360.  !azimuthal angle from 0 to 360.
	else
		theta_polar = 1000.
		theta_xz    = 1000.
		theta_yz    = 1000.
		phi         = 1000.
	end if	
	xx   = plist(18,i) - plist(7,i)
	yy   = plist(19,i) - plist(8,i)
	zz   = plist(20,i) - plist(9,i)
	path = 1.e15
	vertex(1) = plist(7,i)  
	vertex(2) = plist(8,i)
	vertex(3) = plist(9,i)
	if(xx.lt.1.e15 .and. yy.lt.1.e15 .and. zz.lt.1.e15)then
	   path = xx*xx + yy*yy + zz*zz
	   path = sqrt(path)	!pathlength between creation and decay
	end if
	betamag = sqrt(plist(13,i)**2. + plist(14,i)**2. 
     1                                 + plist(15,i)**2.)/plist(17,i)
	gamma = 1. - betamag**2.
	if(gamma.gt.0)then
	   gamma = 1./sqrt(gamma)
	else
	   gamma = 1.e32
c	   write(6,*)'Bad scene...',betamag
c	   write(6,*)plist(13,i),plist(14,i),plist(15,i),plist(17,i)
	end if
	betagammactau = betamag*gamma*plist(21,i)   !in centimeters
cc	write(6,*)betamag,gamma,plist(21,i),betagammactau,pmom,plist(17,i)
c
c       Single track acceptance function call
c       Apply the acceptance function to each track in the event
c	the logical variable IACC is .true. when the track is accepted
c	with any degree of probability, i.e. in the geometric and momentum range.
c	The common-block variable ACCEPTWEIGHT is more nuanced, giving a real value
c	that is in proportion to the acceptance probability.   Lately, this has
c	been used for decaying particles like pions and kaons to include their
c	decay probabilities in the acceptance.
c
c	write(6,*)'>Calling acceptance',lund,theta_polar,phi,pmom,icharge
	iacc = iaccept_GlueX(theta_polar,phi,pmom,icharge,vertex,path,
     1	                                          betagammactau,lund)
c	write(6,*)'Single track>',icount,iacc,lund,acceptweight,
c    1        theta_polar,phi,pmom,icharge

	iexpt_trigger = .true.                     !always triggers
c
	if(iacc)then
	   acceptlist(i) = acceptweight	!for later weighting of i_th track
	else
	   acceptlist(i) = 0.0
	endif
	if(iacc .and. lund.eq.17)then 
	   if(iacc_firstpip.eqv..true.)then
	      iacc_piplus2 = .true.   !flags the second pi+ in the event
	      acceptweight_piplus2 = acceptweight
	   else
	      iacc_firstpip  = .true.
	      iacc_piplus    = .true.
	      acceptweight_piplus = acceptweight
c	      write(6,*)'First pi plus...', acceptweight
	   endif
	   
        else if(iacc .and. lund.eq.-17)then
	   if(iacc_firstpim.eqv..true.)then
	      iacc_piminus2 = .true.   !flags the second pi- in the event
	      acceptweight_piminus2 = acceptweight
	   else
   	      iacc_firstpim   = .true.
	      iacc_piminus    = .true.
	      acceptweight_piminus = acceptweight
c	      write(6,*)'First pi minus...', acceptweight
	   endif

        else if(iacc .and. lund.eq.41)then
	   if(iacc_firstp.eqv..true.)then
	      iacc_proton2 = .true.   !flags the second proton in the event
	      acceptweight_proton2 = acceptweight
c	      write(6,*)'Second proton...', acceptweight
	   else
   	      iacc_firstp   = .true.
	      iacc_proton   = .true.
	      acceptweight_proton = acceptweight
c	      write(6,*)'First proton...', acceptweight
	   endif

        else if(iacc .and. lund.eq.-41)then   !anti-proton
	   acceptweight_antiproton = acceptweight
	   iacc_antiproton = .true.

        else if(iacc .and. lund.eq.42)then   !neutron
	   iacc_neutral = .true. 
	   acceptweight_neutral = acceptweight
	   iacc_neutron = .true. 
	   acceptweight_neutron = acceptweight

        else if(iacc .and. lund.eq.52)then   !deuteron
	   acceptweight_deuteron = acceptweight
	   iacc_deuteron = .true.

        else if(iacc .and. abs(lund).eq.9)then   !muon
	   acceptweight_muon   = acceptweight
	   acceptweight_lepton = acceptweight
	   iacc_muon = .true.

        else if(iacc .and. abs(lund).eq.7)then   !electron
	   acceptweight_electron = acceptweight
	   acceptweight_lepton   = acceptweight
	   iacc_electron = .true.

       	else if(iacc .and. lund.eq.18)then
	   if(iacc_firstk.eqv..true.)then
	      iacc_kplus2 = .true. !flags that two kaons were seen
	   end if
	   iacc_kplus    = .true.
	   iacc_firstk   = .true.
	   acceptweight_kaon = acceptweight
c	   write(6,*)'First kaon...', acceptweight

	else if(iacc .and. lund.eq.-18)then
	   iacc_kminus   = .true.

	else if(iacc .and. lund.eq.51)then
	   iacc_l3h      = .true.

        else if(iacc .and. lund.eq.1)then     
	   if(iacc_firstphoton.eqv..true.)then
	      iacc_photon2 = .true. !flags that a second photon was seen
	      acceptweight_photon2 = acceptweight
	   end if
	   iacc_photon      = .true.
	   iacc_firstphoton = .true.
	   iacc_neutral     = .true. 
	   acceptweight_photon  = acceptweight
	   acceptweight_photon1 = acceptweight
	   acceptweight_neutral = acceptweight

        else if(iacc .and. lund.eq.7)then
	   iacc_electron = .true.
	   acceptweight_electron = acceptweight

        else if(iacc .and. lund.eq.-7)then
	   iacc_positron = .true.
	   acceptweight_positron = acceptweight
	end if

 50	continue
c
c	Now we have looped through all the particles in the event and
c	must decide whether to accept the event as a whole, depending
c	upon how restrictive our trigger is.  So far we have only got 
c	the acceptance flags for each of the separate particles in the 
c	event.
c       
c       Since we are doing a Monte Carlo calculation we can "peek"
c       at the identity of the fiducially-accepted particles and make
c       our event selection based on our knowledge of the particles.
c       In the analysis of real data, we have to make particle identification
c       (PID) plots first, then apply some cuts, etc.  For present 
c       purposes this is not necessary: assume perfect particle ID.
c
c       Software Trigger:  (quite a hack...)
c
c        write(6,*)'>>> ',imark
	iprocess  = 0  ! flag to say whether this event is to be processed
	itrigmask = 1           ! event is accepted in any case 
        if(imark.eq.4)then
           itrigmask = 2  ! event is accepted if we have p, pbar, p   -  hard code this for now
        elseif(imark.eq.8)then
           itrigmask = 3        ! event is accepted if we have p, pbar, p, pi+, pi-   -  hard code this for now
        else
           write(6,*)'Analyzer needs to be told what to accept...'
        endif
	if(iexpt_trigger)then
	   if(itrigmask.eq.1)then
	      iprocess = 1
	   elseif((iacc_proton.and.iacc_proton2.and.iacc_antiproton) !for GlueX p pbar reaction
     1 	      .and. itrigmask.eq.2)then
	      iprocess = 1
	   elseif((iacc_proton.and.iacc_proton2.and.iacc_antiproton.and.
     1	      iacc_piplus.and.iacc_piminus) !for GlueX Lambda anti-Lambda reaction
     1	      .and. itrigmask.eq.3)then
	      iprocess = 1
c	      write(6,*)'Triggermask got to here....'
	   endif
	endif
c
	if(iprocess.eq.1)then
	   if(ifirstprocess)then
	      write(6,*)'The first event(s) to pass the trigger selection:'
	      call prevent(plist,imark+1,ikeep)
	      write(6,*)' '
	      ifirstprocess = .false.
	   endif
	   iwinner = .true.
	else
	   iwinner = .false.
	end if
c
	nbeam(1) = 0     !defines the beam direction for various computations
	nbeam(2) = 0
	nbeam(3) = 1
	nbeam(4) = 0
	pbeamlab(1) = 0.0
	pbeamlab(2) = 0.0
	pbeamlab(3) = ppbeam
	pbeamlab(4) = 0.0
	pbeamlab(5) = ppbeam
c
c	Top of loop to look at particles in this event
c
c 60	continue
 	continue
c       
	do 1000 iloop=1,3
	if(iloop.eq.1)then
	   call hcdir('//PAWC/RAWEVENT',' ')
	   weight  = 1.0
	   ifiddle = .true.	!flag for printing diagnostics later
c	   write(6,*)'> Filling raw histograms'
	elseif(iloop.eq.2)then
	   if(iwinner.eqv..true.)then
	      call hcdir('//PAWC/ACCEPTED',' ')
	   else
c	      write(6,*)'> Event not accepted; returning to calling program.'
	      return
	   end if
	   weight = 1.0
	   ifiddle = .false.  
c	   write(6,*)'> Filling accepted histograms'
	else
	   call hcdir('//PAWC/WEIGHTED',' ')
	   if(itrigmask.eq.10)then
	      weight = acceptweight_proton *
     1	      acceptweight_proton2 *
     1        acceptweight_antiproton 
	   elseif(itrigmask.eq.11)then
	      weight = acceptweight_proton *
     1	      acceptweight_proton2 *
     1        acceptweight_antiproton * 
     1        acceptweight_piplus * 
     1        acceptweight_piminus
	      ifiddle = .false.  
	   end if
	endif
c	write(6,*)'>>>Top of looping',iloop,itrigmask,weight

c**************************************************************************************
c**************************************************************************************
c       A major amount of code was cut away here, corresponding to the study of many
c       reaction channels investigated with mc_generator.  The code was removed
c       in July 2019 to streamline everything for the GlueX baryon-antibaryon studies.
c**************************************************************************************
c**************************************************************************************
c
c	Reconstruct  p p pbar events (for summer 2016+ baryonium exploration for GlueX)
cc 
cc 	Apply a hack up front to analyze this reaction only if exactly 3 tracks 
cc 	are present in the event (+1 for initial state +1 for intermediate state)
cc
	if(imark.gt.5)then
c		write(6,*)'Skipping event with >3 final state tracks', imark
	   goto 910
	endif
c
c	write(6,*)'Starting p - pbar code'
c	call prevent(plist,imark+1,ikeep)   !temporary	

 	do i=2,imark
	lundi = plist(6,i)	!proton
	if(lundi.eq.41)then
c	  write(6,*)'Proton 1    :',(plist(ix,i),  ix=1,4)
	   acceptweight_proton1 = acceptlist(i)
	   do j=i+1,imark
	      lundj = plist(6,j) !proton
	      if(lundj.eq.41)then
c	 	   write(6,*)'Proton 2       :',(plist(ix,j),  ix=1,4)
		 acceptweight_proton2 = acceptlist(j)
		 do k=2,imark
		    lundk = plist(6,k) !antiproton
		    if(lundk.eq.-41)then
c                   write(6,*)'anti-proton:',(plist(ix,k),  ix=1,4)
		       acceptweight_antipro = acceptlist(k)
c
c	       write(6,*)i,j,k
c	       write(6,*)'Proton 1    :',(plist(ix,i),ix=1,4),acceptweight_proton1
c	       write(6,*)'Proton 2    :',(plist(ix,j),ix=1,4),acceptweight_proton2
c	       write(6,*)'Anti-proton :',(plist(ix,k),ix=1,4),acceptweight_antipro
c
		 if(iloop.eq.1 .and. iprocess.eq.1) ikeep = ikeep + 1 ! This event has "passed all cuts"
c       
c                We want the overall c.m quantities
c               
		 e0 = etot(plist(1,1)) !defines the beam plus target
		 ww = plist(4,1) !W of system 
		 do n=1,3
		    beta(n) = plist(n,1)/e0
c                 write(6,*)'beta',beta
		 end do
ccc                 call getbeta(plist(1,1),beta)  !not tested yet...
		 call boost(plist(1,i),beta,p1cm) !proton 1
		 call boost(plist(1,j),beta,p2cm) !proton 2
		 call boost(plist(1,k),beta,apcm) !anti-proton 
c
c                Now do the tricky bit in trying to decide whether to swap the identity of the two protons
c                We do this because of the fact that the real data always has two protons that we have to
c		 associate with the antiproton somehow
c       
c 		 Calculations for Proton 1
c		     
		 pmomp1cm  = absp(p1cm)
		 plongp1cm = p1cm(3)
		 costhp1cm = dot(p1cm,nbeam)/sqrt(psq(p1cm))
		 thp1cm    = acos(costhp1cm)*rade
c
c                Calculations for Proton 2
c		     
		 pmomp2cm  = absp(p2cm)
		 plongp2cm = p2cm(3)
		 costhp2cm = dot(p2cm,nbeam)/sqrt(psq(p2cm))
		 thp2cm    = acos(costhp2cm)*rade
c
		 call hf2(8894,costhp1cm,costhp2cm,weight)
		 if(thp2cm.gt.thp1cm)then !this is what we do with the real data
c		 if(thp2cm.lt.thp1cm)then  !hack to try reverse assignment
		    do n=1,5
		       ptemp(n) = p1cm(n)
		       p1cm(n)  = p2cm(n)
		       p2cm(n)  = ptemp(n)
		       ptemp(n) = plist(n,i)
		       plist(n,i) = plist(n,j)
		       plist(n,j) = ptemp(n)
		    enddo
		 endif
c       
c                Now proceed to do all the invariant mass calculations
c		       
		 call invarmass(plist(1,i),plist(1,j),pimp1p2, rmp1p2) !IM of proton1 and proton2
		 call invarmass(plist(1,i),plist(1,k),pimp1ap, rmp1ap) !IM of proton1 and antiproton
		 call invarmass(plist(1,j),plist(1,k),pimp2ap, rmp2ap) !IM of proton2 and antiproton
		 call hf1(8851,rmp1ap,weight) !look at p1 pbar mass distribution 
		 call hf1(8852,rmp2ap,weight) !look at p2 pbar mass distribution
		 call hf1(8871,rmp1p2,weight) !look at p1 p2   mass distribution
		 call hf2(8850,rmp1p2,rmp2ap,weight) !correlation (not Dalitz plot)
		 call hf2(8853,rmp1ap,rmp2ap,weight) !correlation (not Dalitz plot)
		 rmp1ap2 = rmp1ap**2.
		 rmp2ap2 = rmp2ap**2.
		 call hf2(8854,rmp1ap2,rmp2ap2,weight) ! Dalitz plot
		 call hf2(8897,rmp2ap,ppbeam,weight) !See how invariant mass range varies with beam energy
c       
c                Fill a bunch of "basic" momentum and angle histograms for comparison to real data
c
		 call hf1(28,ppbeam,weight) !beam momentum		
		 call getangles(plist(1,i),pval1,thp1lab,costhp1lab,phip1lab)
		 call getangles(plist(1,j),pval2,thp2lab,costhp2lab,phip2lab)
		 call getangles(plist(1,k),pvala,thpalab,costhpalab,phipalab)
		       
c                pval1 = absp(plist(1,i))   !proton1
		 call hf1(131,pval1,weight) !proton 1
c		 pval2 = absp(plist(1,j))   !proton2
		 call hf1(132,pval2,weight) !proton 2
c		 pvala = absp(plist(1,k))    !anti-proton
		 call hf1(133,pvala,weight) !anti-proton 
		 call hf2(134,pval1,pval2,weight)
		 call hf1(135,thp1lab,weight) !proton 1
		 call hf1(136,thp2lab,weight) !proton 2
		 call hf1(137,thpalab,weight) !anti-proton
                 call hf2(209,pval1,thp1lab,weight)
                 call hf2(230,pvala,thpalab,weight)
c       
c                Calculations for Proton 1 (repeat these in case protons were swapped)
c		     
		 pmomp1cm  = absp(p1cm)
		 plongp1cm = p1cm(3)
		 costhp1cm = dot(p1cm,nbeam)/sqrt(psq(p1cm))
		 thp1cm    = acos(costhp1cm)*rade
c
c                Calculations for Proton 2
c		     
		 pmomp2cm  = absp(p2cm)
		 plongp2cm = p2cm(3)
		 costhp2cm = dot(p2cm,nbeam)/sqrt(psq(p2cm))
		 thp2cm    = acos(costhp2cm)*rade
c
c                Calculations for the antiproton
c		     
		 pmompbarcm  = absp(apcm)
		 plongpbarcm = apcm(3)
		 costhpbarcm = dot(apcm,nbeam)/sqrt(psq(apcm))
		 thpbarcm    = acos(costhpbarcm)*rade

		 call boost(pbeamlab,beta,pbeamcm) !beam photon

		 call boost(pimp1ap,beta,pimp1apcm) !p1 anti-proton combination
		 call boost(pimp2ap,beta,pimp2apcm) !p2 anti-proton combination
		 call getangles
     1 		    (pimp1apcm,pmomimp1ap,thimp1apcm,costhimp1apcm,phiimp1apcm)
		 call getangles
     1		    (pimp2apcm,pmomimp2ap,thimp2apcm,costhimp2apcm,phiimp2apcm)
		 call gett(pbeamcm,pimp2apcm,ttp2ap,ttprimep2ap) !Mandelstam t for beam and proton/anti-proton pair
		 call hf1(8895,-ttp2ap,weight)
		 call hf1(8896,       -ttprimep2ap,weight)
		 call hf2(8898,rmp2ap,-ttprimep2ap,weight) !See how invariant mass range varies with t
		 call hf2(8899,rmp2ap,-ttp2ap     ,weight) !See how invariant mass range varies with t
		 call hf2(8900,costhpbarcm,-ttp2ap,weight) !See how anti-proton angle    varies with t
		 call hf2(8901,costhpbarcm,ppbeam,weight)  !See how anti-proton angle    varies with beam energy
c
		 call hf2(8856,pmomp2cm,pmomp1cm,weight)
		 call hf1(8857,pmomp1cm,weight)
		 call hf1(8858,pmomp2cm,weight)
		 call hf1(8864,pmompbarcm,weight)
c
c		 Angle correlations in CM frame
c
		 call hf1(8865,costhp1cm,weight)
		 call hf1(8866,costhp2cm,weight)
		 call hf1(8867,costhpbarcm,weight)
		 call hf2(8868,costhp2cm,costhp1cm,weight)
		 call hf2(8869,costhp1cm,costhpbarcm,weight)
		 call hf2(8870,costhp2cm,costhpbarcm,weight)

		 call hf1(8875,thp1cm,weight)
		 call hf1(8876,thp2cm,weight)
		 call hf1(8877,thpbarcm,weight)
		 call hf2(8878,thp2cm,thp1cm,weight)
		 call hf2(8879,thp1cm,thpbarcm,weight)
		 call hf2(8880,thp2cm,thpbarcm,weight)

		 call hf2(8881,pmomp1cm,   thp1cm,weight)
		 call hf2(8882,pmomp2cm,   thp2cm,weight)
		 call hf2(8883,pmompbarcm, thpbarcm,weight)

		 call hf1(8884,thimp1apcm,weight)
		 call hf1(8885,thimp2apcm,weight)
		 call hf1(8889,costhimp1apcm,weight)
		 call hf1(8890,costhimp2apcm,weight)
		 call hf1(8891,costhp1cm,weight)
		 call hf1(8892,costhp2cm,weight)
		 call hf1(8893,costhpbarcm,weight)
c
c                Rescale the Dalitz plot... not strictly kinematically kosher
c
c                Egamma = 8.5 GeV --> W = 4.08 GeV
c
		 wtarget = 4.08 !GeV
		 twomp2  = (2.*rmprot)**2.
		 scale = ((wtarget-rmprot)**2.-twomp2)/((ww-rmprot)**2.-twomp2)
		 rmp1ap2hack = (rmp1ap2-twomp2)*scale + twomp2
		 rmp2ap2hack = (rmp2ap2-twomp2)*scale + twomp2
		 call hf2(8859,rmp1ap2hack,rmp2ap2hack,weight) ! Dalitz plot with hack		       
c
		 pcm     = (1./(2.*ww))*sqrt((ww**2. - 0.0)*(ww**2 - twomp2)) !nominal cm two-body momentum
		 ptarget = (1./(2.*wtarget))*
	1		   sqrt((wtarget**2. - 0.0)*(wtarget**2 - twomp2)) !scaled cm two-body momentum
		 sratio = ptarget/pcm

		 pmomp1cmscale = pmomp1cm * sratio
		 pmomp2cmscale = pmomp2cm * sratio
		 call hf1(8861,pmomp1cmscale,weight)
		 call hf1(8862,pmomp2cmscale,weight)

		 call hf2(8863,pmomp1cmscale,rmp2ap2hack,weight)
c
c                Move the "created" baryons to their mutual rest frame, the "baryon-baryon" or bb frame
c
c                Also boost the recoil (target) proton and the beam photon
c       
c		 e0 = etot(pimp2apcm) !use the created-baryon pair
c		 do n=1,3
c		    beta(n) = pimp2apcm(n)/e0
c	 	 end do
cc                write(6,*)'beta',beta
                 call getbeta(pimp2apcm,beta)  !not tested yet...
		 call boost(p1cm,beta,p1bb) !proton 1
		 call boost(p2cm,beta,p2bb) !proton 2
		 call boost(apcm,beta,apbb) !antiproton
		 call boost(pbeamcm,beta,pbeambb) !beam photon (for GJ system)
		 call getangles(p2bb,p2mombb,thp2bb,costhp2bb,phip2bb)
		 call getangles(apbb,apmombb,thapbb,costhapbb,phiapbb)
		 costhhelicity = dot(apbb,beta)/ !s-channel helicity frame polar angle of pbar
     1  			       sqrt(psq(apbb)*psq(beta)) 
		 costhgj       = dot(apbb,pbeambb)/ !Gottfried Jackson polar angle of pbar
     1  			       sqrt(psq(apbb)*psq(pbeambb))
		 costhadair    = costhapbb !Adair system angle of pbar

		 call hf1(8886,costhadair,weight)
		 call hf1(8887,costhhelicity,weight)
		 call hf1(8888,costhgj,weight)

c                write(6,*)wtarget,ww,rmprot,twomp2,scale,pmomp1rescale,pmomp2cmscale
c       
c                Using the longitudinal center-of-mass momenta for the three particles
c                proceed to compute the VanHove plot values
c       
cc                rrmax = (abs(plongp1cm)+abs(plongp2cm)+abs(plongpbarcm))/sqrt(3.)
c		 plongsum = plongp1cm + plongp2cm + plongpbarcm	!should be zero
c                 q1 = plongp1cm
c		  q2 = plongp2cm
c		  q3 = plongpbarcm
c		 q1 = plongp1cm
c		 q2 = plongp2cm
c		 q3 = plongpbarcm
c		 rr = 0.81650*sqrt(q1**2. + q2**2. + q3**2.) !sqrt(2/3) = 0.8165
c		 rad1   = q1/rr
c		 rad2   = q2/rr
c		 rad3   = q3/rr
c		 omega1 = asin(rad1)
c		       omega2 = asin(rad2) - 2.09440 !2pi/3 !not debugged
c		       omega3 = asin(rad3) - 4.18879 !4pi/3
c		       omega  = (omega1 + omega2 + omega3)/3. !use average angle derived from all three particle
c		 omega  = omega1
c       
c              Determine the sector of plot occupied by this event numbered ccw from q1=0 axis
c
c		       if    (q1.ge.0 .and. q2.le.0 .and. q3.ge.0)then
c			     isector = 1
c		       elseif(q1.ge.0 .and. q2.le.0 .and. q3.le.0)then
c			     isector = 2
c		       elseif(q1.ge.0 .and. q2.ge.0 .and. q3.le.0)then
c			     isector = 3
c		       elseif(q1.le.0 .and. q2.ge.0 .and. q3.le.0)then
c			     isector = 4
c		       elseif(q1.le.0 .and. q2.ge.0 .and. q3.ge.0)then
c			     isector = 5
c		       elseif(q1.le.0 .and. q2.le.0 .and. q3.ge.0)then
c			     isector = 6
c		       else
c			     write(6,*)'Bloody Murder'
c			     call exit()
c		       endif
c		       if(isector.eq.3 .or.(isector.eq.2 .and. (abs(q3).gt.abs(q2))))then
c			     omega = 3.14159265 - omega
c		       endif
c		       if(isector.eq.4 .or.(isector.eq.5 .and. (abs(q3).lt.abs(q2))))then
c			     omega = 3.14159265 - omega
c		       endif
c		       xvalue = abs(rad1)*cos(omega) 
c		       yvalue = abs(rad1)*sin(omega) 
c		       xvalue = abs(q1)*cos(omega) 
c		       yvalue = abs(q1)*sin(omega) 
c		       if(iloop.eq.3)then
c			  write(6,*)' '
c			  write(6,*) '>',plongsum, plongp1cm, plongp2cm, plongpbarcm
c			  write(6,*)e0,rr,rad1,rad2,rad3,omega,isector
c			  write(6,*)xvalue,yvalue
c		       endif
c		       call hf2(8855,xvalue,yvalue,weight)
	      endif
	   enddo
	endif
	enddo
	endif
	enddo
c
 910	continue  !used this to skip events not consistent with 3 tracks in final state

c       
c	Reconstruct  Lambda anti-Lambda proton events (for summer 2017+)
c       The same code is supposed to handle Sigma0 anti-Sigma0 events.
c
c 	Skip events that don't have the right number of particles in the tree
c       This is really a hack to speed up the code; could be skipped
c	6    - pure 5 body phase space
c	9    - intermediate resonance decaying to L Lbar
c	13   - intermediate resonance decaying to S0 S0 bar
c
c	write(6,*)'Got to Lambda-anti-Lambda code:', imark
c	if(.not.(imark.eq.9 .or. imark.eq.13 .or. imark.eq.6 .or. imark.eq.8))then
c		write(6,*)'Skipping event with wrong number of tracks'
c		goto 920
c	endif

c     
c       Begin 5-track selection test
c
c	write(6,*)'> Analyze New Event:'  !Hack to see which proton goes with Lambda and which is recoil
c	do i=2,imark
c	   if(iloop.eq.1)then
c	      write(6,*)i, plist(6,i)
c	   endif
c	enddo
c	write(6,*)'Just before 5-track search'
c	call prevent(plist,imark+1,ikeep)   !hack	
c
c       Hacky cheat to see whether we are simulating Mechanism 1, 2, 3 or 4
c
	if(iloop.eq.1)then
	   ipart2 = plist(6,2)  !the LundID of the second particle in the decay table for this event
	   if(ipart2.eq.41)then	!proton in slot 2 of the particle table
	      imechanism = 1
	   elseif(ipart2.eq.57)then !Lambda in slot 2 of the particle table
	      imechanism = 2
	   elseif(ipart2.eq.96)then !anti-Lambda p in slot 2 of the particle table
	      imechanism = 3
	   elseif(ipart2.eq.94)then !"wide" anti-Lambda Lambda in slot 2 of the particle table
	      imechanism = 4
	   endif
c          write(6,*)'Reaction mechanism:',imechanism
	endif
c
 	do i=2,imark
	   lundi = plist(6,i)	!proton
c          write(6,*)'i',i,lundi
 	   if(lundi.eq.41)then
c 	   write(6,*)'Proton 1    :',(plist(ix,i),  ix=1,5)
	      acceptweight_proton1 = acceptlist(i)
	      do j=i+1,imark
		 lundj = plist(6,j) !proton
c 	         write(6,*)'j',j,lundj
		 if(lundj.eq.41)then
c 	 	 write(6,*)'Proton 2    :',(plist(ix,j),  ix=1,5)
		    acceptweight_proton2 = acceptlist(j)
		    do k=2,imark
		      lundk = plist(6,k) !antiproton
c                      write(6,*)'k',k,lundk
		       if(lundk.eq.-41)then
c 		       write(6,*)'anti-proton:',(plist(ix,k),  ix=1,5)
			  acceptweight_antipro = acceptlist(k)
			  do l=2,imark
			     lundl = plist(6,l) !piplus
c 			     write(6,*)'l',l,lundl
			     if(lundl.eq.17)then
c 			     write(6,*)'pi-plus    :',(plist(ix,l),  ix=1,5)
				acceptweight_piplus = acceptlist(l)
				do m=2,imark
				   lundm = plist(6,m) !piminus
c 				   write(6,*)'m',m,lundm
				   if(lundm.eq.-17)then
c                         write(6,*)'pi-minus   :',(plist(ix,m),  ix=1,5)
c				      acceptweight_piminus = acceptlist(m)
c  	       if(iloop.eq.1)then
c		  write(6,*)'Analyze next event:'
c		  write(6,*)i,j,k,l,m
c		  write(6,*)'Proton 1    :',(plist(ix,i),ix=1,5),acceptweight_proton1
c		  write(6,*)'Proton 2    :',(plist(ix,j),ix=1,5),acceptweight_proton2
c		  write(6,*)'Anti-proton :',(plist(ix,k),ix=1,5),acceptweight_antipro
c		  write(6,*)'Pion plus   :',(plist(ix,l),ix=1,5),acceptweight_piplus
c		  write(6,*)'Pion minus  :',(plist(ix,m),ix=1,5),acceptweight_piminus
c	       endif
c
	        if(iloop.eq.1 .and. iprocess.eq.1) ikeep = ikeep + 1	! This event has "passed all cuts"
c       
c               First let us look at the raw lab-frame momentum and angle distributions
c
		call getangles(plist(1,i),p1momlab,thp1lab,costhp1lab,phip1lab) !proton 1
		call getangles(plist(1,j),p2momlab,thp2lab,costhp2lab,phip2lab) !proton 2
		call getangles(plist(1,k),apmomlab,thaplab,costhaplab,phiaplab) !antiproton 
		call getangles(plist(1,l),pipmomlab,thpiplab,costhpiplab,phipiplab) !pi+
		call getangles(plist(1,m),pimmomlab,thpimlab,costhpimlab,phip2pimlab) !pi-
		call hf1(9901,p1momlab,weight)
		call hf1(9902,p2momlab,weight)
		call hf1(9903,apmomlab,weight)
		call hf1(9904,pipmomlab,weight)
		call hf1(9905,pimmomlab,weight)
		call hf1(9906,costhp1lab,weight)
		call hf1(9907,costhp2lab,weight)
		call hf1(9908,costhaplab,weight)
		call hf1(9909,costhpiplab,weight)
		call hf1(9910,costhpimlab,weight)
c       
c               We want the overall c.m quantities
c               
c		e0 = etot(plist(1,1)) !defines the beam plus target
		ww = plist(4,1) !W of system
                call getbeta(plist(1,1),beta)
c                write(6,*)'Beta>',beta
		if(iloop.eq.1) then !as generated 
		   call boost(plist(1,i),beta,p1cm) !proton 1
		   call boost(plist(1,j),beta,p2cm) !proton 2
		   call boost(plist(1,k),beta,apcm) !anti-proton 
		   call boost(plist(1,l),beta,pipcm) !pi plus 
		   call boost(plist(1,m),beta,pimcm) !pi minus 
		   call boost(pbeamlab,  beta,pbeamcm) !beam photon
		   do icopy=1,5
		      aplab(icopy)  = plist(icopy,k)
		      piplab(icopy) = plist(icopy,l)
		   enddo
		else		!as reconstructed
		   call boost(plist(13,i),beta,p1cm) !proton 1
		   call boost(plist(13,j),beta,p2cm) !proton 2
		   call boost(plist(13,k),beta,apcm) !anti-proton 
		   call boost(plist(13,l),beta,pipcm) !pi plus 
		   call boost(plist(13,m),beta,pimcm) !pi minus 
		   call boost(pbeamlab,   beta,pbeamcm) !beam photon
		   do icopy=1,5
		      aplab(icopy)  = plist(icopy+12,k)
		      piplab(icopy) = plist(icopy+12,l)
		   enddo
		endif
		call boost(plist(1,1),beta,psystem)
c
c               Next let us look at the raw center-of-mass frame momentum and angle distributions
c
		call getangles(p1cm, p1momcm,thp1cm,costhp1cm,phip1cm) !proton1
		call getangles(p2cm, p2momcm,thp2cm,costhp2cm,phip2cm) !proton2
		call getangles(apcm, apmomcm,thapcm,costhapcm,phiapcm) !antiproton 
		call getangles(pipcm,pipmomcm,thpipcm,costhpipcm,phipipcm) !pi+
		call getangles(pimcm,pimmomcm,thpimcm,costhpimcm,phip2pimcm) !pi-
		call hf1(9911,p1momcm,weight)
		call hf1(9912,p2momcm,weight)
		call hf1(9913,apmomcm,weight)
		call hf1(9914,pipmomcm,weight)
		call hf1(9915,pimmomcm,weight)
		call hf1(9916,costhp1cm,weight)
		call hf1(9917,costhp2cm,weight)
		call hf1(9918,costhapcm,weight)
		call hf1(9919,costhpipcm,weight)
		call hf1(9920,costhpimcm,weight)
c       
c               Now proceed to do all the invariant mass calculations using CM variables and also LAB variables (as check)
c               Must be careful to combine the correct particles together using 'truth' info from the event generation step;
c               this next section of code is quite a hack:  to work correctly the input definition file must create particles
c               in a fixed order to match the code below.
c
		call invarmass(apcm, pipcm, alamcm,  rmalam) !IM of anti-proton and piplus --> anti-Lambda  CM
		call invarmass(aplab,piplab,alamlab, rmalamlab) !IM of anti-proton and piplus --> anti-Lambda  LAB
		if(imechanism.eq.1)then !Mechanism 1 (Lam-anti-Lam) created at upper vertex
		   rvertex(1) = plist(7,i)  !recoil proton vertex  
		   rvertex(2) = plist(8,i)
		   rvertex(3) = plist(9,i)
		   dvertex(1) = plist(7,j)  !decay proton production vertex
		   dvertex(2) = plist(8,j)
		   dvertex(3) = plist(9,j)
		   call invarmass(plist(1,j),plist(1,m),plamlab,rmplamlab) ! --> Lambda (lab variables)
		   call invarmass(p2cm,pimcm,plamcm, rmplam) !IM of proton2 and piminus    --> Lambda 
		   call invarmass(p1cm,plamcm,ppplam, rmpplam) !IM of proton1 (recoil) and Lambda 
		   call invarmass(p1cm,alamcm,ppalam, rmpalam) !IM of proton1 (recoil) and anti-Lambda 
c		   write(6,*)'One'
		elseif(imechanism.eq.2)then ! Mechanism 2 (p-anti-Lam) created at upper vertex
		   rvertex(1) = plist(7,j)  !recoil proton (primary) vertex
		   rvertex(2) = plist(8,j)
		   rvertex(3) = plist(9,j)
		   dvertex(1) = plist(7,i)  !decay proton production vertex
		   dvertex(2) = plist(8,i)
		   dvertex(3) = plist(9,i)
		   call invarmass(plist(1,i),plist(1,m),plamlab,rmplamlab) ! --> Lambda (lab variables)
		   call invarmass(p1cm,pimcm,plamcm, rmplam) !IM of proton1 and piminus --> Lambda (bachelor) 
		   call invarmass(p2cm,plamcm,ppplam, rmpplam) !IM of proton2 (created) and Lambda 
		   call invarmass(p2cm,alamcm,ppalam, rmpalam) !IM of proton2 (created) and anti-Lambda
c		   write(6,*)'Two'
		elseif(imechanism.eq.3)then !Mechanism 3 (anti-Lam p) created
		   rvertex(1) = plist(7,i)  !recoil proton vertex   - check for correctness!
		   rvertex(2) = plist(8,i)
		   rvertex(3) = plist(9,i)
		   dvertex(1) = plist(7,j)  !decay proton production vertex
		   dvertex(2) = plist(8,j)
		   dvertex(3) = plist(9,j)
		   call invarmass(plist(1,j),plist(1,m),plamlab,rmplamlab)
		   call invarmass(p2cm,pimcm,plamcm, rmplam) !IM of proton2 and piminus    --> Lambda 
		   call invarmass(p1cm,plamcm,ppplam, rmpplam) !IM of proton1 (created) and Lambda (bachelor) 
		   call invarmass(p1cm,alamcm,ppalam, rmpalam) !IM of proton1 (created) and anti-Lambda 
c		   write(6,*)'Three'
		elseif(imechanism.eq.4)then !Mechanism 4 (anti-Lam Lambda) created at lower vertex
		   rvertex(1) = plist(7,i)  !recoil proton vertex   - check for correctness!
		   rvertex(2) = plist(8,i)
		   rvertex(3) = plist(9,i)
		   dvertex(1) = plist(7,j)  !decay proton production vertex
		   dvertex(2) = plist(8,j)
		   dvertex(3) = plist(9,j)
		   call invarmass(plist(1,j),plist(1,m),plamlab,rmplamlab)
		   call invarmass(p2cm,pimcm,plamcm, rmplam) !IM of proton2 and piminus    --> Lambda 
		   call invarmass(p1cm,plamcm,ppplam, rmpplam) !IM of proton1 (recoil) and Lambda 
		   call invarmass(p1cm,alamcm,ppalam, rmpalam) !IM of proton1 (recoil) and anti-Lambda 
c		   write(6,*)'Four'
		else
		   write(6,*)'This hack is messed up at top...',i
		   call exit()
		endif
		call invarmass(plamcm,alamcm,plamalam, rmplamalam)     !IM of Lambda and anti-Lambda 
		call gett(pbeamcm,plamalam,ttplamalam,ttprimeplamalam) !Mandelstam t for beam and Lambda anti-Lambda pair
		call gett(pbeamcm,ppalam,ttp1alam,ttprimep1alam)      !Mandelstam t for beam and proton anti-Lambda pair
		if(imechanism.eq.1)then  !dreadful hack to make all reaction sequences work
		   call invarmass(plamalam,p1cm,ptemp, rmtemp) !IM of Lambda, anti-Lambda, and remaining proton 
                   call invarmass(p1cm,alamcm,ppalam, rmpalam) !IM of proton and anti-Lambda 
		elseif(imechanism.eq.2)then
		   call invarmass(plamalam,p2cm,ptemp, rmtemp) !IM of Lambda, anti-Lambda, and remaining proton 
                   call invarmass(p2cm,alamcm,ppalam, rmpalam) !IM of proton and anti-Lambda 
		elseif(imechanism.eq.3)then 
		   call invarmass(plamalam,p1cm,ptemp, rmtemp) !IM of Lambda, anti-Lambda, and remaining proton 
                   call invarmass(p1cm,alamcm,ppalam, rmpalam) !IM of proton and anti-Lambda 
		elseif(imechanism.eq.4)then 
		   call invarmass(plamalam,p1cm,ptemp, rmtemp) !IM of Lambda, anti-Lambda, and remaining proton 
                   call invarmass(p1cm,alamcm,ppalam, rmpalam) !IM of proton and anti-Lambda 
		else
		   write(6,*)'This hack is messed up at bottom...',i
		   call exit()
		endif
            	call missingmass(psystem,ptemp,rmissplamalam,rmissplamalam2)
c                write(6,*)'Overall missing mass^2 :',rmissplamalam2
		call hf1(9887,rmissplamalam2,weight) !look at the overall missing mass distribution
		rmissenergy = psystem(5) - ptemp(5)  
         	call hf1(9620,rmissenergy,weight)    !look at the overall missing energy distribution 
         	call hf2(9621,rmissplamalam2,rmissenergy,weight) !correlation of missing energy and missing mass
c       
		call getangles(plamcm,pmomplamcm,thplamcm,costhplamcm,phiplamcm) !Lambda
		call getangles(alamcm,pmomalamcm,thalamcm,costhalamcm,phialamcm) !anti-Lambda

		call getangles(plamlab,pmomplamlab,thplamlab,costhplamlab,phiplamlab) !Lambda
		call getangles(alamlab,pmomalamlab,thalamlab,costhalamlab,phialamlab) !anti-Lambda

		call getangles
     1		(plamalam,pmomplamalamcm,thplamalamcm,costhplamalamcm,phiplamalamcm) !Combined
c       
c               Plots of invariants
c
		call hf1(9801,rmplam,weight) !look at the Lambda  mass distribution 
		call hf1(9802,rmalam,weight) !look at the anti-Lambda  mass distribution 
c               if(costhplamcm.gt.0.0)then
		   call hf1(9851,rmplamalam,weight) !look at the anti-Lambda Lambda mass distribution
		   call hf1(9856,-ttprimeplamalam,weight) !look at the t distribution for the produced pair 
c		else
		   call hf1(9852,rmpalam,weight)   !look at the anti-Lambda proton mass distribution
		   call hf1(9857,-ttprimep1alam,weight)   !look at the t distribution for the produced pair 
c		endif
                call hf2(9850,rmpplam,rmpalam,weight) !correlation (not Dalitz plot)
		call hf2(9855,rmplamalam,rmpalam,weight) !correlation (not Dalitz plot)
		rmpplam2 = rmpplam**2.
		rmpalam2 = rmpalam**2.
		call hf2(9854,rmpplam2,rmpalam2,weight) ! Dalitz plot
c
c               Look at quantities in the LAB frame
c               
		call hf1(28,ppbeam,weight) !beam momentum		
		avertex(1) = plist(7,k)	   !anti-proton production vertex (assumes known particle index in PLIST, won't work for Phase Space, etc)
		avertex(2) = plist(8,k)
		avertex(3) = plist(9,k)
		plampath(1) = dvertex(1) - rvertex(1)  !Lambda pathlength
		plampath(2) = dvertex(2) - rvertex(2)
		plampath(3) = dvertex(3) - rvertex(3)
		alampath(1) = avertex(1) - rvertex(1)  !anti-Lambda pathlength
		alampath(2) = avertex(2) - rvertex(2)
		alampath(3) = avertex(3) - rvertex(3)
		dlampath(1) = avertex(1) - dvertex(1)  !separation of the two...
		dlampath(2) = avertex(2) - dvertex(2)
		dlampath(3) = avertex(3) - dvertex(3)
		plamlength = sqrt(plampath(1)**2.+plampath(2)**2.+plampath(3)**2.)
		alamlength = sqrt(alampath(1)**2.+alampath(2)**2.+alampath(3)**2.)
		alamlamlen = sqrt(dlampath(1)**2.+dlampath(2)**2.+dlampath(3)**2.)
		tauplam = plamlength*rmasslambda/(clight*pmomplamlab)  !use speed of light in cm/nsec
		taualam = alamlength*rmasslambda/(clight*pmomalamlab)  !
		call hf1(604,rvertex(1),weight)
		call hf1(605,rvertex(2),weight)
		call hf1(606,rvertex(3),weight)
		call hf1(607,dvertex(1),weight)
		call hf1(608,dvertex(2),weight)
		call hf1(609,dvertex(3),weight)
		call hf1(610,avertex(1),weight)
		call hf1(611,avertex(2),weight)
		call hf1(612,avertex(3),weight)
		call hf1(613,plamlength,weight)
		call hf1(614,alamlength,weight)
		call hf1(615,alamlamlen,weight)
		call hf1(616,tauplam,weight)
		call hf1(617,taualam,weight)
		call hf2(618,taualam,tauplam,weight)
c
c               Look at momentum and angular correlations in the overall CM frame
c
		call hf2(9881,pmomplamcm, costhplamcm,weight)
		call hf2(9883,pmomalamcm, costhalamcm,weight)
		call hf2(9884,pmomplamalamcm, costhplamalamcm,weight)
		call hf2(9869,costhplamcm,costhalamcm,weight)
		call hf2(9880,pmomplamcm,pmomalamcm,weight)
c                
c               Define the normal direction to the CM production plane
c
		call crossprod(pbeamcm,plamalam,nhatdircm,temp)
c       
c		Define the direction out of the Lambda - anti-Lambda production plane in overall CM
c
c		if(costhplamcm.gt.0.0)then    !this will bias the PAIR frame distributions
		if(.true.)then                !do this in any case
		   call crossprod(plamcm,alamcm,pollamlam,openingangle)
		   costh_nhatcm_nlamlamcm = dot(nhatdircm,pollamlam)
	1		/(absp(nhatdircm)*absp(pollamlam))
		   planeangle = acos(costh_nhatcm_nlamlamcm)*rade
c	  	    write(6,*)'Lambda      in overall CM:',plamcm
c		    write(6,*)'Anti-Lambda in overall CM:',alamcm
c		    write(6,*)'Polarization axis, angle :',pollamlam, openingangle
c		    write(6,*)'Production normal, angle :',nhatdircm, temp
c		    write(6,*)'Cosine between planes    :',costh_nhatcm_nlamlamcm 
c		    write(6,*)'Angle  between planes    :',planeangle
c
c                   Various projections
c
                   call hf1(9885,costhplamcm,weight)
                   call hf1(9886,costhalamcm,weight)
                   call hf1(9888,costh_nhatcm_nlamlamcm,weight)
                   call hf1(9889,planeangle,weight)
c
c                  Boost the Lambda anti-Lambda pair to their mutual rest frame and examine their momentum
c                  and angular correlations there
c
		   call getbeta(plamalam,betapair)
		   call boost(plamcm, betapair, plampair)
		   call boost(alamcm, betapair, alampair)
		   call getangles
	1	     (plampair,pmomplampair,thplampair,costhplampair,phiplpair) !Lambda XYZ
		   call getangles
	1	     (alampair,pmomalampair,thalampair,costhalampair,phialpair) !anti-Lambda XYZ
		   call getdecayangle(plampair,betapair,costhplamdecay,thetaplamdecay) !     Lambda wrt pair axis
		   call getdecayangle(alampair,betapair,costhalamdecay,thetaalamdecay) !anti-Lambda wrt pair axis
		   call hf1(9603,costhplamdecay,weight)
		   call hf1(9604,costhalamdecay,weight)
		   call hf1(9611,pmomplampair,weight) !momentum in pair rest frame
		   call hf1(9612,pmomalampair,weight) !momentum in pair rest frame (better be the same...)
                   weightq = pmomplampair    ! momentum of decaying pair in its own rest frame
                   weightk = pmomplamalamcm  ! momentum of pair in the CM frame
                   phasespaceq = weightq !/rmplamalam !phase space factor for two-body decay of pair
                   phasespacek = weightk !/ww         !phase space factor for CM system of two bodies
                   pscale = 0.5   !arbitrary scaling
                   reweight  = weight/phasespaceq * pscale
                   reweight2 = reweight/phasespacek
                   call hf1(9616,pmomplampair,reweight) !momentum in pair rest frame with phase space divided out
                   call hf1(9617,pmomplampair,reweight2) !momentum in pair rest frame with phase space divided out
		   call hf1(9853,rmplamalam,reweight)  !look at the anti-Lambda Lambda mass distribution
		   call hf1(9858,rmplamalam,reweight2) !look at the anti-Lambda Lambda mass distribution
c
                   if(costhplamcm.gt.0.0)then !Lambda distribution in pair frame for forward Lambdas in CM frame.
		      call hf1(9613,costhplamdecay,weight)
		      call hf1(9614,costhalamdecay,weight)
                      call hf1(9615,pmomplampair,weight) !momentum in pair rest frame
		   endif
		   
c		   if(iloop.eq.1)then
c		      write(6,*)'Boost Lambda-anti-Lambda to their mutual rest frame'
c		      write(6,*)' Boosting         : ',plamalam
c		      write(6,*)' Pair beta        : ',betapair
c		      write(6,*)' Lambda CM        : ',plamcm 
c		      write(6,*)' anti-Lambda CM   : ',alamcm 
c		      write(6,*)' Lambda PAIR      : ',plampair !check back-to-backness
c		      write(6,*)' anti-Lambda PAIR : ',alampair
c		      write(6,*)' costh XYZ   : ',costhplampair ,costhalampair
c		      write(6,*)' costh Decay : ',costhplamdecay,costhalamdecay
c                     endif
c

c
c               Using the longitudinal center-of-mass momenta for the three particles,
c               proceed to compute the Van Hove plot values, plus some alternatives
c       
                p1long   = plamcm(3)
                p2long   = alamcm(3)
                p1trans  = sqrt(plamcm(1)**2. + plamcm(2)**2)
                p2trans  = sqrt(alamcm(1)**2. + alamcm(2)**2)
                if(imechanism.eq.1)then
                   p3long   = p1cm(3)
                   p3trans  = sqrt(p1cm(1)**2. + p1cm(2)**2)
                elseif(imechanism.eq.2)then
                   p3long   = p2cm(3)
                   p3trans  = sqrt(p2cm(1)**2. + p2cm(2)**2)
                elseif(imechanism.eq.3)then
                   p3long   = p1cm(3)
                   p3trans  = sqrt(p1cm(1)**2. + p1cm(2)**2)
                elseif(imechanism.eq.4)then
                   p3long   = p1cm(3)
                   p3trans  = sqrt(p1cm(1)**2. + p1cm(2)**2)
                endif
                plongsum = p1long + p2long + p3long !should be zero
                q1 = p1long
                q2 = p2long
                q3 = p3long
                scale0 = 0.81650 !sqrt(2/3) = 0.81650
c                scale1 = 1.22474 !sqrt(3/2) = 1.22474
                qq = scale0 * sqrt(q1**2. + q2**2. + q3**2.) 
c
c               The pi term is to cope with the principle value issue.
c               
                omega = atan2(-sqrt(3.)*q1,(q1+2.*q2)) + 3.14159265
c       
c               Hackery to determine the sector of this event, numbered ccw from q1=0 axis
c
                isector = 0
                if    (q1.ge.0 .and. q2.le.0 .and. q3.ge.0)then
                   isector = 1
                   pairmass = rmpplam
                elseif(q1.ge.0 .and. q2.le.0 .and. q3.le.0)then
                   isector = 2
                   pairmass = rmpalam
                elseif(q1.ge.0 .and. q2.ge.0 .and. q3.le.0)then
                   isector = 3
                   pairmass = rmplamalam
                elseif(q1.le.0 .and. q2.ge.0 .and. q3.le.0)then
                   isector = 4
                   pairmass = rmpplam
                elseif(q1.le.0 .and. q2.ge.0 .and. q3.ge.0)then
                   isector = 5
                   pairmass = rmpalam
                elseif(q1.le.0 .and. q2.le.0 .and. q3.ge.0)then
                   isector = 6
                   pairmass = rmplamalam
                else
                   if(iloop.eq.1)
     1              write(6,*)'Cannot sort right...',q1,q2,q3
c                   call exit()
                endif
c
                if(isector.gt.0)then
                   cosvalue = cos(omega)
                   sinvalue = sin(omega)
c
c                  Standard van Hove plot
c               
                   xvalue = qq*cosvalue
                   yvalue = qq*sinvalue
                   call hf2(9870,xvalue,yvalue,weight)
c     
c     
c                  Alternative longitudinal momentum plots:
c                  Invariant mass of pair by each sector with phase space weighting
c                  (SchuPlot)
c                   
                   xvalue = pairmass*cosvalue
                   yvalue = pairmass*sinvalue
c                   reweight3 = reweight2 / pairmass  !should make plot flater for pure phase space
                   call hf2(9871,xvalue,yvalue,weight)
c
c                  Alternative longitudinal momentum plots:
c                  Available kinetic energy
c                   
c                   ww = plist(4,1) !W of system
                   rkinetic = ww - (2.0 * rmasslambda + rmprot)
                   if(rkinetic.lt.0)then
                      write(6,*)'WTF',rkinetic
                      call exit()
                   endif
                   xvalue = rkinetic*cosvalue
                   yvalue = rkinetic*sinvalue
                   call hf2(9872,xvalue,yvalue,weight)
c     
c                  Alternative longitudinal momentum plots:
c                  Transverse momentum
c
                   ptrans = sqrt(p1trans**2.+p2trans**2.+p3trans**2.)
                   xvalue = ptrans*cosvalue
                   yvalue = ptrans*sinvalue
                   call hf2(9873,xvalue,yvalue,weight)
                endif
                angle = omega*rade
                call hf2(9874,angle,pairmass,weight)

c       
c               Construct the proton and anti-proton 4-vectors in the rest frame of each respective
c               hyperon.   This requires boosting them first from c.m. to the mutual pair rest frame
c               and then to the rest frame of each hyperon separately.
c
		   call boost(p1cm, betapair, p1pair) !boost both protons to hyperon pair rest frame
		   call boost(p2cm, betapair, p2pair) 
		   call boost(apcm, betapair, appair) !the anti-proton that is part of the anti-lambda
		   call getbeta(plampair,betatophyp)
		   if(imechanism.eq.2)then  !hack to resolve combinatoric issue
		      call boost(p1pair, betatophyp, pphyp) !the Lambda's proton in the Lambda rest frame
		   elseif(imechanism.eq.1 .or. imechanism.eq.3 .or. imechanism.eq.4)then
		      call boost(p2pair, betatophyp, pphyp) !the Lambda's proton in the Lambda rest frame
		   else
		      write(6,*)'very confused...'
		      call exit()
		   endif
		   call getbeta(alampair,betatoahyp)
		   call boost(appair, betatoahyp, aphyp) !the anti_Lambda's anti-proton in anti-Lambda rest frame
c	   	   write(6,*)'Proton in Hyperon rest frame          :',pphyp
c		   write(6,*)'Anti-proton in Anti-Hyperon rest frame:',aphyp
c
c                  Project the (anti) proton direction onto the polarization axis and also onto the
c                  coordinate axes to get the direction cosines.   In each case, we need unit vectors for
c                  the respective hyperon rest frame axes.  Use the direction of the hyperon as the Z axis.
c
		   pmagpphyp = absp(pphyp)
		   pmagaphyp = absp(aphyp)
		   pmagpollamlam = absp(pollamlam)  ! this may already be a unit vector, but make sure...
c
		   do ii=1,3
		      zhatphyp(ii) = betatophyp(ii)
		      yhatphyp(ii) = pollamlam(ii)
		      zhatahyp(ii) = betatoahyp(ii)
		      yhatahyp(ii) = pollamlam(ii)
		   enddo
		   call crossprod(yhatphyp,zhatphyp,xhatphyp,pangle)
		   call crossprod(yhatahyp,zhatahyp,xhatahyp,aangle)
		   xnormphyp = absp(xhatphyp)
		   ynormphyp = absp(yhatphyp)
		   znormphyp = absp(zhatphyp)
		   xnormahyp = absp(xhatahyp)
		   ynormahyp = absp(yhatahyp)
		   znormahyp = absp(zhatahyp)
		   do ii=1,3
		      xhatphyp(ii) = xhatphyp(ii)/xnormphyp
		      yhatphyp(ii) = yhatphyp(ii)/ynormphyp
		      zhatphyp(ii) = zhatphyp(ii)/znormphyp
		      xhatahyp(ii) = xhatahyp(ii)/xnormahyp
		      yhatahyp(ii) = yhatahyp(ii)/ynormahyp
		      zhatahyp(ii) = zhatahyp(ii)/znormahyp
		   enddo
c		   if(iloop.eq.1)then
c   		      write(6,*)'Lambda frame basis vectors'
c		      write(6,*)'x',xhatphyp, xnormphyp, pangle*rade
c		      write(6,*)'y',yhatphyp, ynormphyp
c		      write(6,*)'z',zhatphyp, znormphyp
c		      write(6,*)'anti Lambda frame basis vectors'
c		      write(6,*)'x',xhatahyp, xnormahyp, aangle*rade
c		      write(6,*)'y',yhatahyp, ynormahyp
c		      write(6,*)'z',zhatahyp, znormahyp
c                  endif
		   costhx_pphyp = dot(pphyp,xhatphyp)/pmagpphyp
		   costhy_pphyp = dot(pphyp,yhatphyp)/pmagpphyp
		   costhz_pphyp = dot(pphyp,zhatphyp)/pmagpphyp
		   costhx_aphyp = dot(aphyp,xhatahyp)/pmagaphyp
		   costhy_aphyp = dot(aphyp,yhatahyp)/pmagaphyp
		   costhz_aphyp = dot(aphyp,zhatahyp)/pmagaphyp
c       
c               Plot polarization-axis event distributions and direction cosines
c
c		   if(iloop.eq.1)then
c		      write(6,*)'Analyzer:'
c		      write(6,*)'costhy_pphyp:',costhy_pphyp
c		      write(6,*)'costhy_aphyp:',costhy_aphyp
c		      write(6,*)'pphyp       :',pphyp
c		      write(6,*)'aphyp       :',aphyp
c		      write(6,*)'p1pair      :',p1pair
c		      write(6,*)'p2pair      :',p2pair
c		      write(6,*)'p1 lab      :',(plist(ix,i),ix=1,5),acceptweight_proton1
c		      write(6,*)'p2 lab      :',(plist(ix,j),ix=1,5),acceptweight_proton2
c                   endif

		   call hf1(9601,costhy_pphyp,weight)
		   call hf1(9602,costhy_aphyp,weight)
		   call hf1(9605,costhx_pphyp,weight)

		   call hf1(9608,costhx_aphyp,weight)
		   call hf1(9607,costhz_aphyp,weight)
		   call hf1(9610,costhz_aphyp,weight)
c
c                  Sum this event with the rest to compute (at the end of the analysis)
c                  the various polarizations and spin correlations
c
		   if(iloop.eq.1)then !do this for raw data, i.e. events before GlueX acceptance applied
c		      polsumplamX = polsumplamX + costhx_pphyp
c		      polsumplamY = polsumplamY + costhy_pphyp
c		      polsumplamZ = polsumplamZ + costhz_pphyp
c		      polsumalamX = polsumalamX + costhx_aphyp
c		      polsumalamY = polsumalamY + costhy_aphyp
c		      polsumalamZ = polsumalamZ + costhz_aphyp

		      corrsum(1,4) = corrsum(1,4) + costhx_pphyp
		      corrsum(2,4) = corrsum(2,4) + costhy_pphyp
		      corrsum(3,4) = corrsum(3,4) + costhz_pphyp
		      corrsum(1,5) = corrsum(1,5) + costhx_aphyp
		      corrsum(2,5) = corrsum(2,5) + costhy_aphyp
		      corrsum(3,5) = corrsum(3,5) + costhz_aphyp
		      correrr(1,4) = correrr(1,4) + (costhx_pphyp)**2.
		      correrr(2,4) = correrr(2,4) + (costhy_pphyp)**2.
		      correrr(3,4) = correrr(3,4) + (costhz_pphyp)**2.
		      correrr(1,5) = correrr(1,5) + (costhx_aphyp)**2.
		      correrr(2,5) = correrr(2,5) + (costhy_aphyp)**2.
		      correrr(3,5) = correrr(3,5) + (costhz_aphyp)**2.

c		      corrsumXX  = corrsumXX + costhx_aphyp*costhx_pphyp
c		      corrsumYY  = corrsumYY + costhy_aphyp*costhy_pphyp
c		      corrsumZZ  = corrsumZZ + costhz_aphyp*costhz_pphyp
c		      corrsumXY  = corrsumXY + costhx_aphyp*costhy_pphyp
c		      corrsumYX  = corrsumYX + costhy_aphyp*costhx_pphyp ! same as previous
c		      corrsumXZ  = corrsumXZ + costhx_aphyp*costhz_pphyp
c		      corrsumZX  = corrsumZX + costhz_aphyp*costhx_pphyp ! same as previous
c		      corrsumYZ  = corrsumYZ + costhy_aphyp*costhz_pphyp
c		      corrsumZY  = corrsumZY + costhz_aphyp*costhy_pphyp ! same as previous

		      corrsum(1,1)  = corrsum(1,1) +  costhx_aphyp*costhx_pphyp
		      correrr(1,1)  = correrr(1,1) + (costhx_aphyp*costhx_pphyp)**2.
		      corrsum(2,2)  = corrsum(2,2) +  costhy_aphyp*costhy_pphyp
		      correrr(2,2)  = correrr(2,2) + (costhy_aphyp*costhy_pphyp)**2.
		      corrsum(3,3)  = corrsum(3,3) +  costhz_aphyp*costhz_pphyp
		      correrr(3,3)  = correrr(3,3) + (costhz_aphyp*costhz_pphyp)**2.
		      corrsum(1,2)  = corrsum(1,2) +  costhx_aphyp*costhy_pphyp
		      correrr(1,2)  = correrr(1,2) + (costhx_aphyp*costhy_pphyp)**2.
		      corrsum(2,1)  = corrsum(2,1) +  costhy_aphyp*costhx_pphyp 
		      correrr(2,1)  = correrr(2,1) + (costhy_aphyp*costhx_pphyp)**2. 
		      corrsum(1,3)  = corrsum(1,3) +  costhx_aphyp*costhz_pphyp
		      correrr(1,3)  = correrr(1,3) + (costhx_aphyp*costhz_pphyp)**2.
		      corrsum(3,1)  = corrsum(3,1) +  costhz_aphyp*costhx_pphyp 
		      correrr(3,1)  = correrr(3,1) + (costhz_aphyp*costhx_pphyp)**2. 
		      corrsum(2,3)  = corrsum(2,3) +  costhy_aphyp*costhz_pphyp
		      correrr(2,3)  = correrr(2,3) + (costhy_aphyp*costhz_pphyp)**2.
		      corrsum(3,2)  = corrsum(3,2) +  costhz_aphyp*costhy_pphyp 
		      correrr(3,2)  = correrr(3,2) + (costhz_aphyp*costhy_pphyp)**2.
		   endif
		   call hf2(9890,4.,1.,costhx_pphyp)
		   call hf2(9890,4.,2.,costhy_pphyp)
		   call hf2(9890,4.,3.,costhz_pphyp)
		   call hf2(9890,5.,1.,costhx_aphyp)
		   call hf2(9890,5.,2.,costhy_aphyp)
		   call hf2(9890,5.,3.,costhz_aphyp)
		      
		   call hf2(9890,1.,1.,costhx_aphyp*costhx_pphyp)
		   call hf2(9890,2.,2.,costhy_aphyp*costhy_pphyp)
		   call hf2(9890,3.,3.,costhz_aphyp*costhz_pphyp)

		   call hf2(9890,2.,1.,costhy_aphyp*costhx_pphyp)
		   call hf2(9890,1.,2.,costhx_aphyp*costhy_pphyp)
		   call hf2(9890,1.,3.,costhx_aphyp*costhz_pphyp)
		   call hf2(9890,3.,1.,costhz_aphyp*costhx_pphyp)
		   call hf2(9890,2.,3.,costhy_aphyp*costhz_pphyp)
		   call hf2(9890,3.,2.,costhz_aphyp*costhy_pphyp)
		   goto 920	!we  have successfully processed this event... bail out
	     endif
c       
c
c 					Step out of the 5-level deep loop of track selection
c
	                                endif
					enddo
				endif
				enddo
			endif
			enddo
		endif
		enddo
	endif
	enddo
 920	continue
c
c       End of event-examination loop
c       
 1000   continue
c
c        write(6,*)'Returning from analyzer'
 2000	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c	This version is for the GlueX spectrometer at Jefferson Lab
c       This routine is mainly to estimate the geometrical acceptance
c       of the detector.
c
c	We use the usual coordinates for the experiment
c	1 = X = to the left as viewed in beam direction
c	2 = Y = up
c	3 = Z = beam direction
c	angles are in degrees
c
c       The cut-off values are pure fiction at this point
c       
		
	function iaccept_GlueX(theta,phi,pmom,icharge,vertex,path,
     1	                                      betagammactau,lund)
	logical iaccept_GlueX,iacc
	real vertex(3)
	common /gluex_local/ipcount
	common /detector/acceptweight
	data ipcount/0/
c
	ipcount = ipcount+1
	if(ipcount.eq.1)then
	   write(6,*)'Using GlueX acceptance function',ipcount
	endif
	iaccept_GlueX = .false.
	iacc          = .true.
	acceptweight  = 1.00   !this parameter is here in case we're allowing for decay lifetimes...
	pmomgev = pmom/1000.   !particle momentum in GeV/c 
        pminvalcharged = .15   !minimum accepted momentum in GeV
        pminvalpion    = .15   !minimum accepted momentum in GeV
        pminvalproton  = .35   !minimum accepted momentum in GeV
        pminvalaproton = .35   !minimum accepted momentum in GeV
        pminvalneutral = 0.03  !GeV
c
	if(icharge.ne.0)then  !charged particles
	   if(abs(lund).eq.17)then !pions
	      if(pmomgev.lt.pminvalpion)then
		 iacc = .false.
		 acceptweight = 0.0
	      end if 
	   elseif(abs(lund).eq.41)then !protons
	      if(pmomgev.lt.pminvalproton)then
		 iacc = .false.
		 acceptweight = 0.0
	      end if
	   elseif(abs(lund).eq.-41)then !anti-protons
	      if(pmomgev.lt.pminvalaproton)then
		 iacc = .false.
		 acceptweight = 0.0
	      end if
	   elseif(pmomgev.lt.pminvalcharged)then
	      iacc = .false.
	      acceptweight = 0.0
	   end if
	   if(theta.gt.150. .or. theta.lt.2.0)then  !RASch 8-31-16; cut on lab polar angle of charged tracks
	      iacc = .false.
	      acceptweight = 0.0
	   end if
c	   write(6,*)'>',lund,icharge,pmomgev,pminvalcharged,iacc
	else    !neutral particles
	   if(lund.eq.1)then  !this is the ID for a photon
c	      write(6,*)'>Got a photon',theta,pmomgev
	      if(theta.lt.2.0 .or. theta.gt.125.)then
		 iacc = .false.
		 acceptweight = 0.0
c		 write(6,*)'Failed on angle...'
	      endif	      
	      if(pmomgev.lt.pminvalneutral)then
		 iacc = .false.
		 acceptweight = 0.0
c		 write(6,*)'Failed on energy...'
	      endif
	   else               !some other neutral particle
	      continue
	   endif
	endif

c
	if(iacc.eqv..true.)then
		iaccept_GlueX = .true.
	end if
	return
	x = vertex(1)  !fool the gfortran compiler
	x = path 
	x = phi
	x = theta
	x = betagammactau
	x = lund
	end


c
c	Compute missing mass
c	1 = px; 2 = py; 3=pz; 4=rmass; 5=total energy
c	Input is in GeV, output is in GeV/c2
c
	subroutine missingmass(p1,p2,rmiss,rmiss2)
	real p1(5),p2(5)
	px = p1(1) - p2(1)	
	py = p1(2) - p2(2)	
	pz = p1(3) - p2(3)	
	e1 = p1(5)
	e2 = p2(5)
	ee = e1 - e2
	rmiss2= ee*ee - (px*px + py*py + pz*pz)
	if(rmiss2.gt.0)then
		rmiss = sqrt(rmiss2)
	else
		rmiss = -0.000001 !hack for later plotting
	end if
	return
	end

c
c	Compute missing mass and return the 4-vector as well.
c	1 = px; 2 = py; 3=pz; 4=rmass; 5=total energy 
c	Input is in GeV, output is in GeV/c2
c
c       (Compute the total energy as element 5 as of 5-4-15)
c
	subroutine missmass(p1,p2,p3,rmiss,rmiss2)
	real p1(5),p2(5),p3(5)
	p3(1) = p1(1) - p2(1)	
	p3(2) = p1(2) - p2(2)	
	p3(3) = p1(3) - p2(3)	
	e1 = p1(5)
	e2 = p2(5)
	ee = e1 - e2
	rmiss2= ee*ee - (p3(1)**2. + p3(2)**2. + p3(3)**2.)
	if(rmiss2.gt.0)then
		rmiss = sqrt(rmiss2)
		p3(5) = ee
	else
		rmiss = -0.00001
		p3(5) = 0
	end if
	p3(4) = rmiss
c	write(6,*)'Missmass     p1>',p1
c	write(6,*)'             p2>',p2
c	write(6,*)'             p3>',p3
c	write(6,*)'e1,e2,ee,rm,rm2>',e1,e2,ee,rmiss,rmiss2
	return
	end

c****************************************************************************
c       Convert a particle 4-vector to angles in standard coordinates
c       Z - beam axis, when given a vector in lab coordinates
c       Y - up in lab frame
c       X - to left at seen along beam
c
c       R. A. Sch. CMU
c       5-1-2015 - first iteration.  I have no idea why I didn't write this
c                  code a very long time ago.
c
c       Input:
c       particle - input giving px,py,pz,mass,energy
c
c       Output:
c       ppart - scalar momentum
c       theta - polar angle, in degrees
c       costheta - cosine of polar angle
c       phi - azimuthal angle
c
        subroutine getangles(particle,ppart,theta,costheta,phi)
	real particle(5) !px,py,px,mass,total energy
	data rade /57.29577951/
c
	ppart = absp(particle)
	px = particle(1) 
	py = particle(2) 
	pz = particle(3) 
	if(ppart.gt.0)then
	   costheta = pz/ppart
	   if(costheta.gt. 1.0)costheta = 1.0
	   if(costheta.lt.-1.0)costheta =-1.0
c	   theta_xz = atan2(px,pz) * rade
c	   theta_yz = atan2(py,pz) * rade
	   theta = acos(costheta) * rade
	   phi = atan2(py,px) * rade
	   if(phi.lt.0.0)phi = phi + 360. !azimuthal angle from 0 to 360.
	else
	   costheta = 2.0
	   theta = 1000.
c	   theta_xz = 1000.
c	   theta_yz = 1000.
c	   phi = 1000.
	end if	
	return
	end


c*******************************************************************************
c       In a particular frame, take a particle 4-vector and a direction vector
c       called 'BETA" and compute the angle of one with respect to the other
c
c       R.A. Sch,  CMU,  8-2017
c
c       PVECTOR (input) is defined as PX, PY, PZ, MASS, ENERGY
c       BETA (input) is the boost 3-vector
c       COSTHDECAY - cosine of angle between a boost direction and the direction of a
c	             decay particle
c       THETADECAY - same angle in degrees
c 
        subroutine getdecayangle(pvector,beta,costhdecay,thetadecay)
c
	real pvector(5),beta(3)
	data rade /57.29577951/
c
	costhdecay = dot(pvector,beta)/(absp(pvector)*absp(beta))
	thetadecay = acos(costhdecay) * rade
c	write(6,*)'Getdecayangle: decay angle calculation'
c	write(6,*)' pvector   ',pvector
c	write(6,*)' beta      ',beta
c	write(6,*)' COSTHdecay',costhdecay
c	write(6,*)' Thetadecay',thetadecay
	return
	end

c*******************************************************************************
c       Compute Mandelstam t for two given 4-vectors
c       P1 is the incident particle
c       P3 is the "produced" particle that couples at a vertex with P1
c
c       tt is the Mandelstam invariant t;  it will be < 0
c       ttprime is the reduced value of t, namely t - tmin
c
c       R.A.Sch. 8-2017
c
	subroutine gett(p1,p3,tt,ttprime)
	real p1(5),p3(5),p0(5)
	tt = (p3(1)-p1(1))**2. + (p3(2)-p1(2))**2. + (p3(3)-p1(3))**2.
	tt = (p3(5)-p1(5))**2. - tt
c
	pmagp3 = sqrt(p3(1)**2.+p3(2)**2.+p3(3)**2.)
	p0(1) = 0.0
	p0(2) = 0.0
	p0(3) = pmagp3
	p0(4) = p3(4)
	p0(5) = p3(5)
	ttprime = (p0(1)-p1(1))**2. + (p0(2)-p1(2))**2. + (p0(3)-p1(3))**2.
	ttprime = (p0(5)-p1(5))**2. - ttprime
	ttprime = tt - ttprime
c	write(6,*)'Computing t',tt,ttprime
c	write(6,*)p1
c	write(6,*)p3
	return
	end
	
cccccccccccccccccccccc  END of MC_ANALYZE CODE  ccccccccccccccccccccccccccccccccccc
