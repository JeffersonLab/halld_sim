c*******************************************************************************
c	MC_ANALYZE.  Part of MC_GENERATOR
c	Do checks on the data
c
c	Original:
c		10-28-91	
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
c                                 Various calls in this subroutine to PSMEAR 
c                                 have been removed.
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
c       
	subroutine mc_analyze(plist,imark,iwriteit)
c
c	Coding for the list of particles generated in the decay chain:
c	All variables in this list are measured in the lab.
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
c	A list of particle ID codes for reference:
c	LUND particle ID's  (partial list)
c	1	gamma
c	2	neutrino
c	7	electron
c	9	muon -
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
c	57	lambda
c	41	proton
c	42	neutron
c	43	sigma+
c	44	sigma0
c	45	sigma-
c	46	cascade-
c	47	cascade0
c	61	Delta(1232)++
cThe following ones are not standard (I made 'em up)
c	20	K*0(892)
c	21	K*+(892)
c	48	sigma0(1385)
c	49	cascade(1530)
c	50	3He
c	51	L3H
c	52	Deuteron
c	53      Carbon_12
c	54      Carbon_11
c	55	Boron 11
c	56	sigma+(1385)
c	58	Lambda(1405)
c	59	Lambda(1520)
c	69	f0
c	74	x(1285)		
c	75	a0+(980)
c   	76	N(1535) S11
c	77	a2(1320) 
c  
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
	real pp1plam(5),pp1alam(5),plampair(5),alampair(5),p2hyp(5)
	real pphyp(5),aphyp(5)
	real p1pair(5),p2pair(5)
	real xhatphyp(3),yhatphyp(3),zhatphyp(3)
	real xhatahyp(3),yhatahyp(3),zhatahyp(3)
	real corrsum(3,5),correrr(3,5)
	
	logical ifirst,iacc_piplus,iacc_piminus,iacc_kplus,iwriteit,
     1      iacc_kminus,
     1      iacc_kstar,iacc_phi,iacc_l3h,iacc_proton,iacc_firstk,
     1      istop_cas,iacc_photon,iacc_photon2,
     1      iacc_electron,iacc_positron,
     1      ihackflag,iacc_kplus2,iwinner,iacc,
     1      idogluex,
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
c	data rmanot  / 0.9847/ !mass of the a_0 +- (from PDG)
	data rmanot  / 0.960 / !mass of the a_0 +-
c	data gammaa0 / 0.1000/ !average width of the a_0 +-
c	data gammaa0 / 0.0750/ !average width of the a_0 +-
	data gammaa0 / 0.0500/ !average width of the a_0 +-
	common /pawc/hmemor(nwpawc)
	common /stuff/poldir(4),xx,idocuts,icount,ikeep,ihackflag,
     1                iattempt
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
	
c	write(6,*)'Entering analyzer',icount,ikeep
	if(idocuts.eq.1 .and. ifirst.eqv..true.)then
	   ifirst = .false.
	   write(6,*)'First call to event analyzer'
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
	      call hbook1(1,'PHOTON KINETIC ENERGY (MeV)'     ,80,0.,1000.,0.)
	      call hbook1(2,'NEUTRINO KINETIC ENERGY (MeV)'   ,80,0.,1000.,0.)
c              call hbook1(3,'ELECTRON KINETIC ENERGY (MeV)'   ,80,0.,100.,0.)
	      call hbook1(4,'MUON KINETIC ENERGY (MeV)'       ,80,0.,1500.,0.)
	      call hbook1(5,'[p]^+! Kinetic Energy (MeV)'     ,80,0.,1500.,0.)
	      call hbook1(6,'PI0  KINETIC ENERGY (MeV)'       ,50,0.,2500.,0.)
	      call hbook1(1006,'ETA KINETIC ENERGY (MeV)'     ,50,0.,2500.,0.)
	      call hbook1(7,'K^+! KINETIC ENERGY (MeV)'       ,80,0.,1500.,0.)
	      call hbook1(8,'KAON0  KINETIC ENERGY (MeV)'     ,80,0.,1500.,0.)
	      call hbook1(9,'PROTON KINETIC ENERGY (MeV)'     ,50,0.,2000.,0.)
	      call hbook1(10,'NEUTRON KINETIC ENERGY (MeV)'   ,80,0.,1500.,0.)
	      call hbook1(11,'LAMBDA KINETIC ENERGY (MeV)'    ,80,0.,1500.,0.)
c	      call hbook1(12,'CASCADE- KINETIC ENERGY (MeV)'  ,80,0.,500.,0.)
c	      call hbook1(13,'CASCADE0 KINETIC ENERGY (MeV)'  ,80,0.,500.,0.)
	      call hbook1(14,'[p]^-! Kinetic Energy (MeV)'    ,80,0.,1500.,0.)
	      call hbook1(15,'KAON -  KINETIC ENERGY (MeV)'   ,80,0.,1500.,0.)
	      call hbook1(16,'SIGMA +  KINETIC ENERGY (MeV)'  ,80,0.,1500.,0.)
	      call hbook1(17,'SIGMA 0  KINETIC ENERGY (MeV)'  ,80,0.,1500.,0.)
	      call hbook1(18,'SIGMA -  KINETIC ENERGY (MeV)'  ,80,0.,1500.,0.)
c              call hbook1(19,'f0 KINETIC ENERGY (MeV)'        ,80,0.,100.,0.)
c              call hbook1(20,'PHI  KINETIC ENERGY (MeV)'      ,80,0.,100.,0.)
	      call hbook1(21,'L3H  KINETIC ENERGY (MeV)'      ,80,0.,100.,0.)
	      call hbook1(22,'Deuteron Kinetic Energy (MeV)'  ,80,0.,500.,0.)
c	      call hbook1(24,'CAS(1530) KINETIC ENERGY (MeV)'	,80,0.,1000.,0.)
	      call hbook1(25,'K*(892) KINETIC ENERGY (MeV)'	,70,0.,700.,0.)
	      call hbook1(26,'BEAM ENERGY (GeV)'   ,nebin,ebinmin,ebinmax,0.)
	      call hbook1(27,'cos([Q]?c.m.!) K^+!'            ,100,-1.,1.,0.)
c	      call hbook1(28,'Beam Momentum (GeV/c)' ,nebin,ebinmin,ebinmax,0.)
c	      call hbook1(28,'Beam Momentum (GeV/c)' , 100,  0.4, 5.0,0.)
	      call hbook1(28,'Beam Momentum (GeV/c)' , 100,  3.5, 12.0,0.)
	      call hbook1(29,'cos([Q]?c.m.!) [p]^+!'            ,20,-1.,1.,0.)
	      call hbook1(30,'cos([Q]?c.m.!^meson!)'           ,20,-1.,1.,0.)
	      call hbook1(1030,'cos([Q]?c.m.!) [c]'             ,20,-1.,1.,0.)
c
	      call hbook1(101,'PHOTON Momentum (GeV/c)'     ,100,0.,1.0,0.)
	      call hbook1(102,'NEUTRINO MOMENTUM (MeV/c)'   ,180,0.,1800.,0.)
	      call hbook1(103,'Electron Momentum (GeV/c)'   ,100,0.,1.5,0.)
c	      call hbook1(104,'Muon Momentum (GeV/c)'       ,100,0.,1.0,0.)
	      call hbook1(104,'Lepton Momentum (GeV/c)'     ,100,0.,1.0,0.)
	      call hbook1(105,'[p]^+!  Momentum (GeV/c)'    ,100,0.,1.5,0.)
	      call hbook1(106,'[p]^0! Momentum (GeV/c)'     ,100,0.,1.5,0.)
	      call hbook1(1106,'ETA MOMENTUM (MeV/c)'       ,125,0.,2500.,0.)
	      call hbook1(107,'K^+! Momentum (GeV/c)'       ,100,0.,6.0,0.)
	      call hbook1(108,'KAON0  MOMENTUM (MeV/c)'     ,180,0.,1800.,0.)
c	      call hbook1(109,'PROTON Momentum (GeV/c)'     ,100,0.,2.0,0.)
	      call hbook1(109,'PROTON Momentum (GeV/c)'     ,100,0.,6.0,0.)
	      call hbook1(110,'NEUTRON Momentum (GeV/c)'    ,100,0.,2.0,0.)
	      call hbook1(111,'LAMBDA MOMENTUM (MeV/c)'     ,180,0.,2000.,0.)
	      call hbook1(112,'CASCADE - MOMENTUM (MeV/c)'  ,180,0.,2500.,0.)
	      call hbook1(113,'CASCADE0  MOMENTUM (MeV/c)'  ,180,0.,2500.,0.)
c	      call hbook1(114,'[p]^-!  Momentum (GeV/c)'    ,100,0.,2.,0.)
	      call hbook1(114,'[p]^-!  Momentum (GeV/c)'    ,100,0.,1.5,0.)
	      call hbook1(115,'K^-!  Momentum (GeV/c)'      ,100,0.,6.0,0.)
	      call hbook1(116,'SIGMA +  MOMENTUM (MeV/c)'   ,180,0.,1800.,0.)
	      call hbook1(117,'SIGMA 0  MOMENTUM (MeV/c)'   ,180,0.,1800.,0.)
	      call hbook1(118,'SIGMA -  MOMENTUM (MeV/c)'   ,180,0.,2000.,0.)
c		call hbook1(119,'f0 MOMENTUM (MeV/c)'         ,180,0.,1800.,0.)
c		call hbook1(120,'PHI  MOMENTUM (MeV/c)'       ,180,0.,1800.,0.)
	      call hbook1(121,'L3H  MOMENTUM (MeV/c)'       ,180,0.,1800.,0.)
	      call hbook1(122,'Deuteron Momentum (MeV/c)'   ,180,0.,2000.,0.)
	      call hbook1(124,'CAS(1530) MOMENTUM (MeV/c)'  ,180,0.,2500.,0.)
	      call hbook1(125,'K*(892) MOMENTUM (MeV/c)'    ,180,0.,1800.,0.)
	      call hbook1(126,'K*(892) MASS (MeV)'        ,150,700.,1000.,0.)
	      call hbook1(127,'11C RECOIL MOMENTUM (MeV/c)' ,50, 0., 500.,0.)
	      call hbook1(130,'ANTI-PROTON Momentum (GeV/c)'  ,100,0.,8.0,0.)
	      call hbook1(131,'P1 Lab Momentum (GeV/c)'     ,100,0.,8.0,0.)
	      call hbook1(132,'P2 Lab Momentum (GeV/c)'     ,100,0.,8.0,0.)
	      call hbook1(133,'p&^-! Lab Momentum (GeV/c)'  ,100,0.,8.0,0.)
	      call hbook2(134,'P2 vs. P1 LAB Momenta (GeV/c)',
     1	                            50,0.,8.0,50,0.,8.0,0.)
	      call hbook1(135,'P1 Lab Angle (deg)'     ,90,0.,90.0,0.)
	      call hbook1(136,'P2 Lab Angle (deg)'     ,90,0.,90.0,0.)
	      call hbook1(137,'p&^-!  Lab Angle (deg)'     ,90,0.,90.0,0.)
c       
	      call hbook2(201,'PHOTON [Q] (deg) vs. p (GeV/c)'     
     1  					,100,0.,1.,90,0.,180.,0.)
c		call hbook2(202,'NEUTRINO ANGLE (deg) vs. P (MeV/c)'   
c     1  					,40,0.,1800.,40,0.,180.,0.)
	      call hbook2(203,'Electron [Q] (deg) vs. p (GeV/c)'   
     1  					,100,0.,1.,90,0.,180.,0.)
	      call hbook2(204,'Muon [Q] (deg) vs. p (GeV/c)'       
     1  					,40,0.,1.,40,0.,180.,0.)
	      call hbook2(205,'[p]^+! [Q] (deg) vs. p (GeV/c)'       
     1  					,100,0.,1.5,90,0.,180.,0.)
	      call hbook2(206,'[p]^0! [Q] (deg) vs. p (GeV/c)'       
     1  					,100,0.,1.5,90,0.,180.,0.)
	      call hbook2(1206,'ETA ANGLE (deg) vs. P (MeV/c)'       
     1  					,50,0.,2500.,40,0.,180.,0.)
	      call hbook2(207,'K^+!  [Q] (deg) vs. p?K! (GeV/c)'       
     1  					,100,0.,6.,100,0.,90.,0.)
	      call hbook2(2070,'K^+!  ANGLE (deg) vs. P (MeV/c)'    
     1  					,36,0.,180.,20,0.,2000.,0.)
	      call hbook2(2071,'K^+! (GeV/c)  vs. E?[g]! (GeV)',    
     1                              60,0.,5.0, 60, 0.0, 4.0,0.)
	      call hbook2(208,'KAON0  ANGLE (deg) vs. P (MeV/c)'     
     1  					,40,0.,1800.,40,0.,180.,0.)
c	      call hbook2(209,'PROTON  [Q] (deg) vs. p (GeV/c)'       
c     1  					,100,0.,2.,100,0.,90.,0.)
	      call hbook2(209,'PROTON  [Q] (deg) vs. p (GeV/c)'       
     1  					,100,0.,8.,100,0.,90.,0.)
	      call hbook2(2090,'ANTI-PROTON  [Q] (deg) vs. p (GeV/c)'       
     1  					,100,0.,8.,100,0.,90.,0.)
              call hbook2(210,'NEUTRON [Q] (deg) vs. P (GeV/c)'   
     1  					,100,0.,2.,100,0.,90.,0.)
	      call hbook2(211,'[L] ANGLE (deg) vs. P (MeV/c)'    
     1  					,40,0.,2000.,40,0., 90.,0.)
	      call hbook2(212,'CASCADE - ANGLE (deg) vs. P (MeV/c)'
     1  					,40,0.,2500.,40,0., 90.,0.)
	      call hbook2(213,'CASCADE0  ANGLE (deg) vs. P (MeV/c)'
     1  					,40,0.,2500.,40,0., 90.,0.)
	      call hbook2(214,'[p]^-!  [Q] (deg) vs. p! (GeV/c)'       
     1  					,100,0.,1.5,100,0.,180.,0.)
	      call hbook2(215,'K^-!  [Q] (deg) vs. p (GeV/c)'       
     1  					,100,0.,6.,100,0.,90.,0.) 
	      call hbook2(216,'{S}^+! ANGLE (deg) vs. P (MeV/c)'  
     1  					,40,0.,1800.,40,0., 90.,0.)
	      call hbook2(217,'SIGMA 0  ANGLE (deg) vs. P (MeV/c)'  
     1  					,40,0.,1800.,40,0., 90.,0.)
	      call hbook2(218,'SIGMA -  ANGLE (deg) vs. P (MeV/c)'  
     1  					,40,0.,2000.,40,0., 90.,0.)
c		call hbook2(219,'f0 ANGLE (deg) vs. P (MeV/c)'        
c     1  					,40,0.,1800.,40,0., 90.,0.)
c		call hbook2(220,'PHI  ANGLE (deg) vs. P (MeV/c)'      
c     1  					,40,0.,1800.,40,0., 90.,0.)
	      call hbook2(221,'L3H  ANGLE (deg) vs. P (MeV/c)'      
     1  					,40,0.,1800.,40,0., 90.,0.)
	      call hbook2(222,'Deuteron Angle (deg) vs. P (MeV/c)'      
     1  					,40,0.,2000.,40,0., 90.,0.)
	      call hbook2(223,'PROTON ANGLE (deg) vs. P (MeV/c)'     
     1  					,80,0.,1800.,90,0., 180.,0.)
c	      call hbook2(224,'CAS(1530) ANGLE (deg) vs. P (MeV/c)'
c     1  					,40,0.,2500.,40,0., 90.,0.)
	      call hbook2(225,'K*(892) ANGLE (deg) vs. P (MeV/c)'
     1  					,40,0.,1800.,40,0., 90.,0.)
	      call hbook2(226,'cos([Q]?c.m.!) vs. E?[g]!) (GeV)',    
     1  	   	    nebin,ebinmin,ebinmax,nkbin,cbinmin,cbinmax,0.)
	      call hbook2(2226,'cos([Q]?c.m.!) vs. E?[g]! (GeV)',    
     1  	   	    nebin,ebinmin,ebinmax,nkbin,cbinmin,cbinmax,0.)
              call hbook2(3226,'cos([Q]?c.m.!^[p]^o!) vs. E?[g]! (GeV)',    
     1  	   	    nebin,ebinmin,ebinmax,nkbin,cbinmin,cbinmax,0.)
	      call hbook2(1226,'cos([Q]?c.m.!^[c]!) vs. E?[g]! (GeV)',    
     1  	   	    nebin,ebinmin,ebinmax,nkbin,cbinmin,cbinmax,0.)
	      call hbook2(4226,'cos([Q]?c.m.!^[c]!) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,10,-1.0,1.0,0.)

	      call hbook2(4236,'cos([Q]?c.m.!^[c]!) vs. W (GeV) [c]^,!-->[rg] ',    
     1  	   10,1.8,2.8,10,-1.0,1.0,0.)
	      call hbook2(4246,'cos([Q]?c.m.!^[c]!) vs. W (GeV) [c]^,!-->[cpp] ',    
     1  	   10,1.8,2.8,10,-1.0,1.0,0.)
	      call hbook2(4227,'cos([Q]?c.m.!^[c]!) vs. W (GeV) x1280',    
     1  	   10,1.8,2.8,10,-1.0,1.0,0.)

	      call hbook2(4601,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)
	      call hbook2(4602,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)
	      call hbook2(4603,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)
	      call hbook2(4604,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)
	      call hbook2(4605,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)
	      call hbook2(4606,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)
	      call hbook2(4607,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)
	      call hbook2(4608,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)
	      call hbook2(4609,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)
	      call hbook2(4610,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)

	      call hbook2(4611,'MM^2! (GeV/c)^2! vs. W (GeV) [c]^,! ',    
     1  	   10,1.8,2.8,100,-.2,0.5,0.)
	      call hbook2(4612,'MM^2! (GeV/c)^2! vs. W (GeV) [c]^,! ',    
     1  	   10,1.8,2.8,100,-.2,0.5,0.)
	      call hbook2(4613,'MM^2! (GeV/c)^2! vs. W (GeV) [c]^,! ',    
     1  	   10,1.8,2.8,100,-.2,0.5,0.)
	      call hbook2(4614,'MM^2! (GeV/c)^2! vs. W (GeV) [c]^,! ',    
     1  	   10,1.8,2.8,100,-.2,0.5,0.)
	      call hbook2(4615,'MM^2! (GeV/c)^2! vs. W (GeV) [c]^,! ',    
     1  	   10,1.8,2.8,100,-.2,0.5,0.)
	      call hbook2(4616,'MM^2! (GeV/c)^2! vs. W (GeV) [c]^,! ',    
     1  	   10,1.8,2.8,100,-.2,0.5,0.)
	      call hbook2(4617,'MM^2! (GeV/c)^2! vs. W (GeV) [c]^,! ',    
     1  	   10,1.8,2.8,100,-.2,0.5,0.)
	      call hbook2(4618,'MM^2! (GeV/c)^2! vs. W (GeV) [c]^,! ',    
     1  	   10,1.8,2.8,100,-.2,0.5,0.)
	      call hbook2(4619,'MM^2! (GeV/c)^2! vs. W (GeV) [c]^,! ',    
     1  	   10,1.8,2.8,100,-.2,0.5,0.)
	      call hbook2(4620,'MM^2! (GeV/c)^2! vs. W (GeV) [c]^,! ',    
     1  	   10,1.8,2.8,100,-.2,0.5,0.)


	      call hbook2(4701,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)
	      call hbook2(4702,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)
	      call hbook2(4703,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)
	      call hbook2(4704,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)
	      call hbook2(4705,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)
	      call hbook2(4706,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)
	      call hbook2(4707,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)
	      call hbook2(4708,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)
	      call hbook2(4709,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)
	      call hbook2(4710,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)

	      call hbook2(4711,'MM^2! (GeV/c)^2! vs. W (GeV) [c]^,! ',    
     1  	   10,1.8,2.8,100,-.2,0.5,0.)
	      call hbook2(4712,'MM^2! (GeV/c)^2! vs. W (GeV) [c]^,! ',    
     1  	   10,1.8,2.8,100,-.2,0.5,0.)
	      call hbook2(4713,'MM^2! (GeV/c)^2! vs. W (GeV) [c]^,! ',    
     1  	   10,1.8,2.8,100,-.2,0.5,0.)
	      call hbook2(4714,'MM^2! (GeV/c)^2! vs. W (GeV) [c]^,! ',    
     1  	   10,1.8,2.8,100,-.2,0.5,0.)
	      call hbook2(4715,'MM^2! (GeV/c)^2! vs. W (GeV) [c]^,! ',    
     1  	   10,1.8,2.8,100,-.2,0.5,0.)
	      call hbook2(4716,'MM^2! (GeV/c)^2! vs. W (GeV) [c]^,! ',    
     1  	   10,1.8,2.8,100,-.2,0.5,0.)
	      call hbook2(4717,'MM^2! (GeV/c)^2! vs. W (GeV) [c]^,! ',    
     1  	   10,1.8,2.8,100,-.2,0.5,0.)
	      call hbook2(4718,'MM^2! (GeV/c)^2! vs. W (GeV) [c]^,! ',    
     1  	   10,1.8,2.8,100,-.2,0.5,0.)
	      call hbook2(4719,'MM^2! (GeV/c)^2! vs. W (GeV) [c]^,! ',    
     1  	   10,1.8,2.8,100,-.2,0.5,0.)
	      call hbook2(4720,'MM^2! (GeV/c)^2! vs. W (GeV) [c]^,! ',    
     1  	   10,1.8,2.8,100,-.2,0.5,0.)

	      call hbook2(5601,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)
	      call hbook2(5602,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)
	      call hbook2(5603,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)
	      call hbook2(5604,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)
	      call hbook2(5605,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)
	      call hbook2(5606,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)
	      call hbook2(5607,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)
	      call hbook2(5608,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)
	      call hbook2(5609,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)
	      call hbook2(5610,'p?[p+]! (GeV/c) vs. W (GeV) [c]^,!',    
     1  	   10,1.8,2.8,20,0.0,2.5,0.)

              call hbook2(13260,'cos([Q]?c.m.!^[c]!) vs. E?[g]! (GeV)',    
     1  	   	    nebin,ebinmin,ebinmax,nkbin,cbinmin,cbinmax,0.)
	      call hbook2(1326,'E?[g]! vs. (GeV)cos([Q]?c.m.!^[c]!)',    
     1  	                               10,-1.,1.,50,.750,2.000,0.)
	      call hbook2(228,'[Q]?cm! vs. Phi?lab! '    
     1  					,50,0.,180.,50,0.,360.,0.)
	      call hbook2(229,'cos([Q]?c.m.!) vs. P?[g]! (GeV/c)'    
     1  					,100,1.0,4.2, 20,-1.,1.,0.)
	      call hbook2(230,'ANTI-PROTON  [Q] (deg) vs. p (GeV/c)'       
     1  					,100,0.,8.,100,0.,90.,0.)
	      call hbook2(2290,'P?[g]! (GeV/c) vs. cos([Q]?c.m.!)'    
     1  	        	,110,-1.1,1.1,nebin,ebinmin,ebinmax,0.)
	      call hbook2(1070,'K+  Reconstructed vs. Thrown (MeV/c)' 
     1		                               ,40,0.,2000.,40,0.,2000.,0.)
              call hbook2(1071,'K+ p(generated) - p(recon) (MeV/c)' 
     1		                            ,200,0.,2000.,40,-100.,100.,0.)
	      call hbook1(1072,'Missing Mass off Kaon'  ,80,1.0,1.3,0.)
	      call hbook2(1073,'Kaon Lab Angle vs. Missing Mass' 
     1		                                ,100,1.0,1.3,50,0.,100.,0.)
c
c		call hbook1(31,'PHOTON DECAY DISTANCE (cm)'  ,80,0.,200.,0.)
c		call hbook1(32,'NEUTRINO DECAY DISTANCE (cm)',80,0.,100.,0.)
c		call hbook1(33,'ELECTRON DECAY DISTANCE (cm)',80,0.,100.,0.)
		call hbook1(34,'MUON DECAY DISTANCE (cm)'    ,80,0.,1e5,0.)
		call hbook1(35,'PI+- DECAY DISTANCE (cm)'    ,80,0.,6e3,0.)
		call hbook1(36,'PI0  DECAY DISTANCE (cm)'    ,80,0.,10.,0.)
		call hbook1(1036,'ETA  DECAY DISTANCE (cm)'  ,80,0.,10.,0.)
		call hbook1(37,'K^+-! DECAY DISTANCE (cm)'   ,80,0.,1200.,0.)
		call hbook1(38,'KAON0  DECAY DISTANCE (cm)'  ,80,0.,30.,0.)
c		call hbook1(39,'PROTON DECAY DISTANCE (cm)'  ,80,0.,100.,0.)
c		call hbook1(40,'NEUTRON DECAY DISTANCE (cm)' ,80,0.,50.,0.)
		call hbook1(41,'LAMBDA DECAY DISTANCE (cm)'  ,80,0.,30.,0.)
		call hbook1(42,'CASCADE- DECAY DISTANCE (cm)',80,0.,7.,0.)
		call hbook1(43,'CASCADE0 DECAY DISTANCE (cm)',80,0.,5.,0.)
c
c		call hbook1(70,'COS([Q]_LAM_LAM)' ,40,-1.05,1.05,0.)
c		call hbook1(71,'COS([Q]_P_P)'     ,40,-1.05,1.05,0.)
c		call hbook1(72,'COS([Q]_PI-_PI-)',40,-1.05,1.05,0.)
c		call hbook1(73,'COS([Q]_GAM_GAM) ONE PI0',40,-1.05,1.05,0.)
c		call hbook1(74,'COS([Q]_GAM_GAM) ANY PI0',40,-1.05,1.05,0.)
c		call hbook1(75,'[Q]_K_K (lab)'            ,90,0.0,90.0,0.)
c		call hbook1(76,'[Q]_pi_pi (lab)'          ,180,0.0,180.0,0.)
		call hbook1(77,'COS([Q]p wrt NHAT) in [L]cm',40,-1.05,1.05,0.)
		call hbook1(78,'COS([Q]K wrt BEAM), s-frame',40,-1.05,1.05,0.)
		call hbook2(79,'COS p (wrt Nhat in L cm) vs COS K (CM)'    
     1  				,80,-1.00,1.00,20,-1.00, 1.00,0.)
		call hbook1(80,'COS([Q]p wrt Xprime) in [L]cm',40,-1.05,1.05,0.)
		call hbook1(81,'COS([Q]p wrt Zprime) in [L]cm',40,-1.05,1.05,0.)
		call hbook2(82,'COS p (wrt xprime in L cm) vs COS K (CM)'  
     1  				,80,-1.00,1.00,20,-1.00, 1.00,0.)
		call hbook2(83,'COS p (wrt zprime in L cm) vs COS K (CM)'  
     1  				,80,-1.00,1.00,20,-1.00, 1.00,0.)
		call hbook2(84,'Mass(K+pi-)**2 vs. Mass(pi-p)**2'  
     1  				,60,1.00,4.00,60,0.20, 2.00,0.)
		call hbook2(85,'Index of Kinematic Bins vs. Proton Angle'  
     1         ,10,-1.00,1.00,nebinpol*nkbinpol,0.5, nebinpol*nkbinpol+0.5,0.)
		call hbook2(850,'COS p (wrt Nhat in L cm) vs COS K (CM)'    
     1  				,80,-1.00,1.00,20,-1.00, 1.00,0.)
		call hbook2(86,'Index of Kinematic Bins vs. Proton Angle'  
     1         ,10,-1.00,1.00,nebinpol*nkbinpol,0.5, nebinpol*nkbinpol+0.5,0.)
		call hbook2(860,'COS p (wrt Nhat in L cm) vs COS K (CM)'    
     1  				,80,-1.00,1.00,20,-1.00, 1.00,0.)
		call hbook2(87,'Index of Kinematic Bins vs. Proton Angle'  
     1  	,10,-1.00,1.00,nebinpol*nkbinpol,0.5, nebinpol*nkbinpol+0.5,0.)
		call hbook2(870,'COS p (wrt Nhat in L cm) vs COS K (CM)'    
     1  				,80,-1.00,1.00,20,-1.00, 1.00,0.)
		call hbook2(990,'Kinematic Bins vs. Proton Angle: Cx, K p'  
     1         ,10,-1.00,1.00,nebinpol*nkbinpol,0.5, nebinpol*nkbinpol+0.5,0.)
		call hbook2(9900,'COS p (wrt Cx in L cm) vs COS K (CM)'    
     1  				,80,-1.00,1.00,20,-1.00, 1.00,0.)
		call hbook2(991,'Kinematic Bins vs. Proton Angle: Cz, K p'  
     1         ,10,-1.00,1.00,nebinpol*nkbinpol,0.5, nebinpol*nkbinpol+0.5,0.)
c		call hbook2(9910,'COS p (wrt Cz in L cm) vs COS K (CM)'    
c     1  				,80,-1.00,1.00,20,-1.00, 1.00,0.)
		call hbook2(992,'Kinematic Bins vs. Proton Angle: Cx, K pi'  
     1         ,10,-1.00,1.00,nebinpol*nkbinpol,0.5, nebinpol*nkbinpol+0.5,0.)
c		call hbook2(9920,'COS p (wrt Cx in L cm) vs COS K (CM)'    
c     1  				,80,-1.00,1.00,20,-1.00, 1.00,0.)
		call hbook2(993,'Kinematic Bins vs. Proton Angle: Cz, K pi'  
     1         ,10,-1.00,1.00,nebinpol*nkbinpol,0.5, nebinpol*nkbinpol+0.5,0.)
		call hbook2(9930,'COS p (wrt Cx in L cm) vs COS K (CM)'    
     1  				,80,-1.00,1.00,20,-1.00, 1.00,0.)
		call hbook2(994,'Kinematic Bins vs. Proton Angle: Cx, K p pi'  
     1  		  ,10,-1.00,1.00,nebin*10,0.5, nebin*10.+0.5,0.)
c		call hbook2(9940,'COS p (wrt Cx in L cm) vs COS K (CM)'    
c     1  				,80,-1.00,1.00,20,-1.00, 1.00,0.)
c		call hbook2(995,'Kinematic Bins vs. Proton Angle: Cz, K p pi'  
c     1  		  ,10,-1.00,1.00,nebin*10,0.5, nebin*10.+0.5,0.)
c		call hbook2(9950,'COS p (wrt Cx in L cm) vs COS K (CM)'    
c     1  				,80,-1.00,1.00,20,-1.00, 1.00,0.)
c		call hbook2(9960,'cos(th pz) vs cos(th px) in L cm'    
c     1  				,80,-1.00,1.00,80,-1.00, 1.00,0.)
c
		call hbook1(51,'PHOTON cos([Q])'      ,40,-1.05,1.05,0.)
c		call hbook1(52,'NEUTRINO cos([Q])'    ,40,-1.05,1.05,0.)
c		call hbook1(53,'ELECTRON cos([Q])'    ,40,-1.05,1.05,0.)
c		call hbook1(54,'MUON cos([Q])'        ,40,-1.05,1.05,0.)
		call hbook1(55,'[p]^+-! cos([Q]?lab!)' ,40,-1.05,1.05,0.)
		call hbook1(56,'[p]^o! cos([Q]?lab!)'  ,40,-1.05,1.05,0.)
		call hbook1(1056,'ETA cos([Q]?lab!)'   ,40,-1.05,1.05,0.)
		call hbook1(57,'K^+-! cos([Q]?lab!)'  ,40,-1.05,1.05,0.)
		call hbook1(58,'KAON0  cos([Q])'      ,40,-1.05,1.05,0.)
		call hbook1(59,'PROTON cos([Q]?cm!)'  ,40,-1.05,1.05,0.)
		call hbook1(60,'NEUTRON cos([Q]?cm!)' ,40,-1.05,1.05,0.)
		call hbook1(61,'LAMBDA cos([Q])'      ,40,-1.05,1.05,0.)
c		call hbook1(62,'CASCADE - cos([Q])'   ,40,-1.05,1.05,0.)
c		call hbook1(63,'CASCADE0  cos([Q])'   ,40,-1.05,1.05,0.)
c
		call hbook1(99,'COS([Q] P) C.M.'    ,40,-1.05,1.05,0.)
c
		call hbook2(450,'MM2 ([g],p[p]^-!) vs. MM2 ([g],p[p]^+!) ',
     1                                   125,0.4,1.4,125,0.4,1.4,0.)
		call hbook2(451,'MM^2! ([g],p[p]^-!) vs. MM^2! ([g],p[p]^+!) (smeared)',
     1                         	         125,0.4,1.4,125,0.4,1.4,0.)
		call hbook2(452,'p([p]^+!) vs. MM^2! ([g],[p]^+![p]^-!) (GeV/c)^2!'
     1  				 ,100,-0.2,0.5,40,0.0,2.,0.)
		call hbook2(453,'p([p]^-!) vs. MM^2! ([g],[p]^+![p]^-!) (GeV/c)^2!'
     1  				 ,100,-0.2,0.5,40,0.0,2.,0.)
		call hbook2(454,'MM^2!([g],p[p]^-!) vs MM^2!([g],p[p]^+!)'
     1  				 ,80,0.0,0.8,80,0.0,0.8,0.)
		call hbook2(455,'MM^2!([g],p[p]^-!) vs MM^2!([g],p[p]^+!)'
     1  				 ,80,0.0,0.8,80,0.0,0.8,0.)
		call hbook2(456,'M^2! ([c][p]^-!) vs. M^2! ([c][p]^+!)',
c     1                         	         125,0.4,1.4,125,0.4,1.4,0.)
     1            	                          50,0.4,1.4, 50,0.4,1.4,0.)

		call hbook1(457,'M^2! ([c][p])'   ,100,0.4,1.4,0.)
		call hbook1(458,'M^2! ([c][p])'   ,100,0.4,1.4,0.)
		call hbook1(459,'IM2 ([c][p])'   ,100,0.4,1.4,0.)
		call hbook2(461,'IM^2! ([c][p]^-!) vs. IM^2! ([c][p]^+!)',
     1                         	         125,0.4,1.4,125,0.4,1.4,0.)
		call hbook2(462,'IM^2! ([c][p]^-!) vs. IM^2! ([c][p]^+!)',
c     1                         	         125,0.4,1.4,125,0.4,1.4,0.)
     1            	                          50,0.4,1.4, 50,0.4,1.4,0.)

c
c		call hbook1(501,'MASS (KK) (MeV/c^2!)'     ,250,850.,1100.,0.)
c		call hbook1(502,'MASS (pi pi) (GeV/c^2!)'  ,250,0.0,2.0,0.)
		call hbook1(503,'MASS (pi-p) (GeV/c^2!)'   ,100,1.0,1.4,0.)
		call hbook1(504,'MM ([g]K^+!) (GeV/c^2!) (smeared)'
     1  						 ,320,0.9,2.5,0.)
		call hbook1(505,'MISSING MASS ([g]K^+!)(GeV/c^2!) K* DATA'
     1  						 ,300,1.2,1.8,0.)
		call hbook1(507,'MISSING MASS ([g]K^+!) (GeV/c^2!) PHI DATA'
     1  						 ,150,1.2,1.8,0.)
c		call hbook1(508,'MISSING MASS (quasi-free) (GeV/c^2!)'
c     1  						 ,150,1.2,1.8,0.)
c		call hbook1(509,'MISSING MASS (quasi-free) (GeV/c^2!)'
c     1  						 ,150,1.2,1.8,0.)
c		call hbook1(510,'MISSING MASS (K+p) (GeV/c^2!)'
c     1  						 ,1000,0.9,1.8,0.)
c		call hbook1(511,'MISSING MASS ((pipi)p) (GeV/c^2!)'
c     1  						 ,1000,0.9,1.8,0.)
		call hbook1(512,'EXCITATION ENERGY (MeV)'
     1  					 ,100,-50.,50.,0.)
		call hbook1(513,'EXCITATION ENERGY (MeV)'
     1  					 ,100,-50.,500.,0.)
		call hbook1(514,'InvMass (p pi-) (GeV/c^2!)'  ,400,1.0,2.0,0.)
		call hbook2(515,'InvMass(p pi-) vs. MissMass(gamma K)'  
     1  				,100,0.5,2.0,100,1.0,2.0,0.)
		call hbook1(516,'InvMass (p pi+) (GeV/c^2!)'  ,100,1.0,2.0,0.)
		call hbook1(517,'([g],[p]^+!) MISSING MASS (neutron) (GeV/c^2!)'
     1  						 ,300,0.6,1.4,0.)
		call hbook1(518,'([g],p) MISSING MASS ([p]^o!) (GeV/c^2!)'
     1  						 ,50,0.0,0.5,0.)
		call hbook1(1518,'([g],p) MISSING MASS ([c]) (GeV/c^2!)'
     1  						 ,50,0.0,0.8,0.)
		call hbook1(2518,'([g],p) MISSING MASS ([c]) (GeV/c^2!)'
     1  						 ,80,0.7,1.5,0.)
		call hbook1(519,'([g],p) MISSING MASS SQUARED ([p]^o!) (GeV/c^2!)^2!'
     1  						 ,50,-0.1,0.15,0.)
		call hbook1(1519,'([g],p) MISSING MASS SQUARED ([c]) (GeV/c^2!)^2!'
     1  						 ,50,-0.1,0.15,0.)
		call hbook1(2519,'([g],p) MISSING MASS SQUARED ([c]) (GeV/c^2!)^2!'
     1  						 ,50,-0.1,3.0,0.)
		call hbook1(520,'([g],pK^+!) MISSING MASS  (GeV/c^2!)'
     1  						 ,100,0.0,1.0,0.)
		call hbook1(521,'([g],pK^+!) MM SQUARED ([p]^-!) (GeV/c^2!)^2!'
     1  						 ,100,-0.1,1.0,0.)
		call hbook1(530,'MM ([g],pK^+![p]^-!) (GeV/c^2!)'
     1  						 ,100,0.0,0.4,0.)
		call hbook1(531,'MM ([g],pK^+![p]^-!) (GeV/c^2!) (smeared)'
     1  						 ,200,0.0,0.8,0.)
		call hbook1(532,'MM ([g]K^+!) (GeV/c^2!)'
     1  						 ,800,1.0,1.2,0.)
c     1  						 ,320,0.9,2.5,0.)
c
c------Histograms for detail MC study of L(1520) decay in Case 3 --------------

	      call hbook1(7601,'K^+! Momentum (GeV/c)'       ,100,0.,1.0,0.)
	      call hbook1(7602,'[p]^+! Momentum (GeV/c)'     ,100,0.,1.0,0.)
	      call hbook1(7603,'[p]^-! Momentum (GeV/c)'     ,100,0.,1.0,0.)

c	      call hbook1(7611,'K^+! c.m. angle'             ,90,0.,180.0,0.)
c	      call hbook1(7612,'[p]^+! c.m. angle'           ,90,0.,180.0,0.)
c	      call hbook1(7613,'[p]^-! c.m. angle'           ,90,0.,180.0,0.)
	      call hbook1(7611,'K^+! c.m. angle'             ,90,-1.0,1.0,0.)
	      call hbook1(7612,'[p]^+! c.m. angle'           ,90,-1.0,1.0,0.)
	      call hbook1(7613,'[p]^-! c.m. angle'           ,90,-1.0,1.0,0.)


c----Histograms for g14 K0 Lambda study---------------------------

		call hbook1(2531,'IM ([p]^+![p]^-! 1) (GeV/c^2!) '
     1  						 ,250,0.25,0.75,0.)
		call hbook1(2532,'IM ([p]^+![p]^-! 2) (GeV/c^2!) '
     1  						 ,250,0.25,0.75,0.)
		call hbook1(2533,'IM ([p]^+![p]^-! 2) (GeV/c^2!) '
     1  						 ,250,0.25,0.75,0.)
		call hbook1(3531,'IM ([p]^+![p]^-! 1) (GeV/c^2!) (Smeared)'
     1  						 ,250,0.25,0.75,0.)
		call hbook1(3532,'IM ([p]^+![p]^-! 2) (GeV/c^2!) (Smeared)'
     1  						 ,250,0.25,0.75,0.)
		call hbook1(3533,'IM ([p]^+![p]^-! both) (GeV/c^2!) (Smeared)'
     1  						 ,250,0.25,0.75,0.)
		call hbook1(4531,'IM ([p]^+![p]^-! 1) (GeV/c^2!) '
     1  						 ,250,0.25,0.75,0.)
		call hbook1(5531,'IM ([p]^+![p]^-! 1) (GeV/c^2!) (Selected)'
     1  						 ,250,0.25,0.75,0.)
		call hbook1(5532,'IM ([p]^+![p]^-! 2) (GeV/c^2!) (Selected)'
     1  						 ,250,0.25,0.75,0.)
		call hbook1(5533,'IM ([p]^+![p]^-! both) (GeV/c^2!) (Selected)'
     1  						 ,250,0.25,0.75,0.)

		call hbook1(2541,'IM (p ![p]^-! 1) (GeV/c^2!) '
     1  						 ,200,1.05,1.35,0.)
		call hbook1(2542,'IM (p ![p]^-! 2) (GeV/c^2!) '
     1  						 ,200,1.05,1.35,0.)
		call hbook1(2543,'IM (p ![p]^-! both) (GeV/c^2!) '
     1  						 ,200,1.05,1.35,0.)
		call hbook1(2544,'MM ([g]n, K^0!)X (GeV/c^2!) '
     1  						 ,200,0.9,1.4,0.)
		call hbook2(2545,'MM([g],K^0![p]^-!p) vs. MM([g],K^0!0)'
     1  						 ,64,0.9,1.3,64,0.9,1.1,0.)

		call hbook1(3541,'IM (p ![p]^-! 1) (GeV/c^2!) (Smeared)'
     1  						 ,200,1.05,1.35,0.)
		call hbook1(3542,'IM (p ![p]^-! 2) (GeV/c^2!) (Smeared)'
     1  						 ,200,1.05,1.35,0.)
		call hbook1(3543,'IM (p ![p]^-! both) (GeV/c^2!) (Smeared)'
     1  						 ,200,1.05,1.35,0.)
		call hbook1(3544,'MM ([g]n, K^0!)X (GeV/c^2!) (Smeared)'
     1  						 ,200,0.9,1.4,0.)
		call hbook1(5541,'IM (p ![p]^-! 1) (GeV/c^2!) (Selected)'
     1  						 ,200,1.05,1.35,0.)
		call hbook1(5542,'IM (p ![p]^-! 2) (GeV/c^2!) (Selected)'
     1  						 ,200,1.05,1.35,0.)
		call hbook1(5543,'IM (p ![p]^-! both) (GeV/c^2!) (Selected)'
     1  						 ,200,1.05,1.35,0.)
		call hbook1(5544,'MM ([g]n, K^0!)X (GeV/c^2!) (Selected)'
     1  						 ,200,0.9,1.4,0.)

		call hbook1(2550,'MM Spectator (GeV/c^2!) '
     1  						 ,200,0.8,1.2,0.)
		call hbook1(2560,'p K^0!1 (GeV/c) (LAB) '
     1  						 ,200,0.0,2.0,0.)
		call hbook1(2561,'p K^0!2 (GeV/c) (LAB) '
     1  						 ,200,0.0,2.0,0.)
		call hbook1(2562,'p [L]!1 (GeV/c) (LAB) '
     1  						 ,200,0.0,2.0,0.)
		call hbook1(2563,'p [L]!2 (GeV/c) (LAB) '
     1  						 ,200,0.0,2.0,0.)
		call hbook1(2564,'p K^0!1 (GeV/c) (CM) '
     1  						 ,200,0.0,2.0,0.)
		call hbook1(2565,'p K^0!2 (GeV/c) (CM) '
     1  						 ,200,0.0,2.0,0.)
		call hbook1(2566,'p [L]!1 (GeV/c) (CM) '
     1  						 ,200,0.0,2.0,0.)
		call hbook1(2567,'p [L]!2 (GeV/c) (CM) '
     1  						 ,200,0.0,2.0,0.)
		call hbook1(2568,'p x (GeV/c) (CM) '
     1  						 ,200,-2.0,2.0,0.)
		call hbook1(2569,'p y (GeV/c) (CM) '
     1  						 ,200,-2.0,2.0,0.)
		call hbook1(2570,'p z (GeV/c) (CM) '
     1  						 ,200,-2.0,2.0,0.)
		call hbook1(2571,'p tot (GeV/c) (CM) '
     1  						 ,200,0.0,2.0,0.)

c		call hbook1(3550,'MM Spectator (GeV/c^2!) (Smeared)'
c     1  						 ,200,0.8,1.2,0.)
		call hbook1(3550,'Missing Mass m?x! (GeV/c^2!)'
     1  						 ,200,0.8,1.4,0.)
		call hbook1(2551,'Spectator Momentum (GeV/c) '
     1  						 ,200,0.0,1.0,0.)
		call hbook2(3545,'MM([g],K^0![p]^-!p) vs. MM([g],K^0!) (Smeared)'
     1  						 ,64,0.9,1.3,64,0.9,1.1,0.)
		call hbook1(3551,'Spectator Momentum (GeV/c^2!) (Smeared)'
     1  						 ,200,0.0,1.0,0.)
		call hbook2(3552,'[p]^-! 2(K0) vs. [p]^-! 1(Lambda) (GeV/c!)'
     1  						 ,64,0.0,1.0,64,0.0,1.0,0.)
		call hbook2(3553,'[p]^-! 2((K0) vs. [p]^-! 1(L) (Smeared)'
     1  						 ,64,0.0,1.0,64,0.0,1.0,0.)





c----Histograms for g14 pi Delta study---------------------------

		call hbook1(6531,'IM ([p]^+! n) (GeV/c^2!) '
     1  						 ,250,0.8,3.0,0.)
		call hbook1(6532,'IM ([p]^-! n) (GeV/c^2!) '
     1  						 ,250,0.8,3.0,0.)
		call hbook2(6631,'IM^2! ([p]^-! n) vs. ([p]^+! n)(GeV)^2! '
     1  				       ,100,0.6,4.0,100,0.6,4.0,0.)
		call hbook1(6533,'Neutron Momentum (GeV/c) '
     1  						 ,200,0.0,2.0,0.)
		call hbook1(6534,'Neutral beta '
     1  						 ,200,0.0,1.1,0.)
		call hbook1(6535,'IM ([p]^-! [p]^+!) (GeV/c^2!) '
     1  						 ,250,0.25,1.25,0.)



c-----------------------------------------------------------------
		call hbook1(5321,'MM ([g]K^+!) (GeV/c^2!)'
     1  						 ,800,0.9,2.5,0.)
		call hbook1(5322,'MM ([g]K^+!) (GeV/c^2!)'
     1  					 ,800,0.9,2.5,0.)
		call hbook1(5323,'MM ([g]K^+!) (GeV/c^2!)'
     1  						 ,800,0.9,2.5,0.)
		call hbook1(5324,'MM ([g]K^+!) (GeV/c^2!)'
     1  						 ,800,0.9,2.5,0.)
		call hbook1(5325,'MM ([g]K^+!) (GeV/c^2!)'
     1  						 ,800,0.9,2.5,0.)
		call hbook1(5326,'MM ([g]K^+!) (GeV/c^2!)'
     1  						 ,800,0.9,2.5,0.)
		call hbook1(5327,'MM ([g]K^+!) (GeV/c^2!)'
     1  						 ,800,0.9,2.5,0.)
		call hbook1(5328,'MM ([g]K^+!) (GeV/c^2!)'
     1  						 ,800,0.9,2.5,0.)
		call hbook1(5329,'MM ([g]K^+!) (GeV/c^2!)'
     1  						 ,800,0.9,2.5,0.)



		call hbook1(533,'MM ([g]K^+!) (GeV/c^2!) (smeared)'
     1  						 ,320,1.0,1.2,0.)
		call hbook1(534,'MM ([g],K^+![p]^-!) (GeV/c^2!)'
     1  						 ,70,0.9,1.6,0.)
		call hbook1(535,'MM ([g],K^+![p]^-!) (GeV/c^2!) (smeared)'
     1  						 ,70,0.9,1.6,0.)
		call hbook1(536,'MM ([g],pK^+!) (GeV/c^2!)'
     1  						 ,100,0.0,1.0,0.)
		call hbook1(537,'MM ([g],pK^+!) (GeV/c^2!) (smeared)'
     1  						 ,100,0.0,1.0,0.)
		call hbook1(538,'IM (K^+![p]^-!) (GeV/c^2!)'
     1  						 ,100,0.6,1.4,0.)
		call hbook1(5381,'IM (K^+![p]^-!) (GeV/c^2!)'
     1  						 ,100,0.6,1.4,0.)
		call hbook1(5382,'IM (K^+![p]^-!) (GeV/c^2!)'
     1  						 ,100,0.6,1.4,0.)
		call hbook1(5383,'IM (K^+![p]^-!) (GeV/c^2!)'
     1  						 ,100,0.6,1.4,0.)
		call hbook1(5384,'IM (K^+![p]^-!) (GeV/c^2!)'
     1  						 ,100,0.6,1.4,0.)
		call hbook1(5385,'IM (K^+![p]^-!) (GeV/c^2!)'
     1  						 ,100,0.6,1.4,0.)
		call hbook1(5386,'IM (K^+![p]^-!) (GeV/c^2!)'
     1  						 ,100,0.6,1.4,0.)
		call hbook1(5387,'IM (K^+![p]^-!) (GeV/c^2!)'
     1  						 ,100,0.6,1.4,0.)
		call hbook1(5388,'IM (K^+![p]^-!) (GeV/c^2!)'
     1  						 ,100,0.6,1.4,0.)
		call hbook1(5389,'IM (K^+![p]^-!) (GeV/c^2!)'
     1  						 ,100,0.6,1.4,0.)
		call hbook1(539,'IM (K^+![p]^-!) (GeV/c^2!) (smeared)'
     1  						 ,100,0.6,1.4,0.)
		call hbook1(540,'IM (p[p]^-!) (GeV/c^2!)'
     1  						 ,1000,1.0,2.0,0.)
		call hbook1(541,'IM (p[p]^-!) (GeV/c^2!) (smeared)'
     1                                           ,1000,1.0,2.0,0.)
		call hbook1(542,'IM (p[p]^-!) (GeV/c^2!) (smeared)'
     1                                           ,1000,1.0,2.0,0.)
		call hbook2(550,'MM [p]- vs. IM [L]'
     1  				 ,100,1.05,1.25,100,0.0,0.4,0.)


		call hbook1(570,'MM ([g]p) (GeV/c^2!)'
     1  						 ,250,0.6,1.8,0.)
		call hbook1(571,'MM ([g]p) (GeV/c^2!) (smeared)'
     1  						 ,240,0.6,1.8,0.)
		call hbook1(572,'MM^2! ([g]p[p][p]) (GeV/c^2!)^2!'
     1  						 ,400,-0.05,0.5,0.)
		call hbook1(573,'MM^2! ([g]p[p][p]) (GeV/c^2!)^2! (smeared)'
     1  						 ,400,-0.05,0.5,0.)
		call hbook1(574,'Case 1, cos([Q]^c.m.!?x!)'   ,20,-1.,1.,0.)
		call hbook1(5740,'Case 1, cos([Q]^c.m.!?x!)'   ,20,0.8,1.,0.)
		call hbook1(575,'Case 2, cos([Q]^c.m.!?x!)'   ,20,-1.,1.,0.)
		call hbook1(576,'MM ([g],p) (GeV/c^2!) (smeared)'
     1  						 ,1600,0.8,1.6,0.)
		call hbook1(5760,'MM ([g]p) (GeV/c^2!)'
     1  						 ,350,1.10,1.45,0.)

		call hbook2(580,'Case 1, MM^2!([g],p[p]^+![p]^-!) vs MM([g],p)'
     1  				 ,100,0.2,1.6,100,-0.05,0.6,0.)
		call hbook2(581, 
     1                'MM^2!([g],p[p]^+![p]^-!) vs MM([g],p) (smeared)'
     1  				 ,200,0.8,1.6,160,-0.05,0.6,0.)
c     1  				 ,100,0.2,1.6,100,-0.05,0.6,0.)
		call hbook1(582,'IM ([p]^+![p]^-!) (GeV/c)^2!'
     1  						 ,200,0.0,2.0,0.)
		call hbook2(583, 
     1                'MM([g],p[p]^+![p]^-!) vs MM([g],p) (smeared)'
     1  				 ,200,0.8,1.6,160,0.0,1.0,0.)
c     1  				 ,100,0.2,1.6,100,-0.05,0.6,0.)
		call hbook2(591,'Case 1, [c][p]^+![p]^-!,  cos([Q]?x!) vs W bin',
     1  	                                   20,1.8,2.8,9,-0.8,1.0,0.)
		call hbook2(592,'Case 2, K^+![p]^-!(K^0!),  cos([Q]?x!) vs W bin',
     1  	                                   20,1.8,2.8,9,-0.8,1.0,0.)
		call hbook1(593,'MM^2! ([g]p[p][p]) (GeV/c^2!)^2! (smeared) for [c]'
     1  						 ,400,-.05,0.5,0.)
		call hbook1(594,'MM^2! ([g]p[p][p]) (GeV/c^2!)^2! (smeared) for x'
     1  						 ,400,-.05,0.5,0.)
		call hbook1(595,'IM ([p][p]) (GeV/c^2!) '
     1  						 ,400,0.0,1.0,0.)
		call hbook1(596,'IM ([p][p]) (GeV/c^2!) for [g] recoil'
     1  						 ,400,0.0,1.0,0.)
		call hbook1(597,'MM^2!([g],p[p]^+![p]^-!)(GeV/c^2!)^2! (smeared) for [c]'
     1  						 ,650,-0.05,0.6,0.)
		call hbook1(598,'MM^2!([g],p[p]^+![p]^-!)(GeV/c^2!)^2! (smeared) for [c]'
     1  						 ,650,-0.05,0.6,0.)
		call hbook1(5970,'MM([g],p[p]^+![p]^-!)(GeV/c^2!) (smeared) for [c]'
     1  						 ,600,.20,0.8,0.)
		call hbook1(5980,'MM([g],p[p]^+![p]^-!)(GeV/c^2!) (smeared) for [c]'
     1  						 ,600,.20,0.8,0.)
c
c       Block of histograms for background correction calculation in f1 and eta' cross sections.
c       For each of 8 W bins there is a 2D histogram of cos(angle) vs eta mass.
c
		call hbook2(10001,'MM([g],p[p]^+![p]^-!)(GeV/c^2!) (smeared) for [c]'
     1  					      ,600,.20,0.8,  9,-0.8, +1.0,0.)
		call hbook2(10002,'MM([g],p[p]^+![p]^-!)(GeV/c^2!) (smeared) for [c]'
     1  					      ,600,.20,0.8,  9,-0.8, +1.0,0.)
		call hbook2(10003,'MM([g],p[p]^+![p]^-!)(GeV/c^2!) (smeared) for [c]'
     1  					      ,600,.20,0.8,  9,-0.8, +1.0,0.)
		call hbook2(10004,'MM([g],p[p]^+![p]^-!)(GeV/c^2!) (smeared) for [c]'
     1  					      ,600,.20,0.8,  9,-0.8, +1.0,0.)
		call hbook2(10005,'MM([g],p[p]^+![p]^-!)(GeV/c^2!) (smeared) for [c]'
     1  					      ,600,.20,0.8,  9,-0.8, +1.0,0.)
		call hbook2(10006,'MM([g],p[p]^+![p]^-!)(GeV/c^2!) (smeared) for [c]'
     1  					      ,600,.20,0.8,  9,-0.8, +1.0,0.)
		call hbook2(10007,'MM([g],p[p]^+![p]^-!)(GeV/c^2!) (smeared) for [c]'
     1  					      ,600,.20,0.8,  9,-0.8, +1.0,0.)
		call hbook2(10008,'MM([g],p[p]^+![p]^-!)(GeV/c^2!) (smeared) for [c]'
     1  					      ,600,.20,0.8,  9,-0.8, +1.0,0.)
c
		call hbook2(20001,'MM([g],p[p]^+![p]^-!)(GeV/c^2!) (smeared) for [c]'
     1  					      ,600,.20,0.8,  9,-0.8, +1.0,0.)
		call hbook2(20002,'MM([g],p[p]^+![p]^-!)(GeV/c^2!) (smeared) for [c]'
     1  					      ,600,.20,0.8,  9,-0.8, +1.0,0.)
		call hbook2(20003,'MM([g],p[p]^+![p]^-!)(GeV/c^2!) (smeared) for [c]'
     1  					      ,600,.20,0.8,  9,-0.8, +1.0,0.)
		call hbook2(20004,'MM([g],p[p]^+![p]^-!)(GeV/c^2!) (smeared) for [c]'
     1  					      ,600,.20,0.8,  9,-0.8, +1.0,0.)
		call hbook2(20005,'MM([g],p[p]^+![p]^-!)(GeV/c^2!) (smeared) for [c]'
     1  					      ,600,.20,0.8,  9,-0.8, +1.0,0.)
		call hbook2(20006,'MM([g],p[p]^+![p]^-!)(GeV/c^2!) (smeared) for [c]'
     1  					      ,600,.20,0.8,  9,-0.8, +1.0,0.)
		call hbook2(20007,'MM([g],p[p]^+![p]^-!)(GeV/c^2!) (smeared) for [c]'
     1  					      ,600,.20,0.8,  9,-0.8, +1.0,0.)
		call hbook2(20008,'MM([g],p[p]^+![p]^-!)(GeV/c^2!) (smeared) for [c]'
     1  					      ,600,.20,0.8,  9,-0.8, +1.0,0.)
c
		call hbook1(599,'IM ([p][p]) (GeV/c^2!) '
     1  						 ,100,0.0,1.0,0.)
		call hbook1(601,'X Vertex Distribution',100,-3.,+3.,0.)
		call hbook1(602,'Y Vertex Distribution',100,-3.,+3.,0.)
		call hbook1(603,'Z Vertex Distribution',700,-50.,+20.,0.)

		call hbook1(604,'X Primary Distribution (cm)',100,-3.,+3.,0.)
		call hbook1(605,'Y Primary Distribution',100,-3.,+3.,0.)
		call hbook1(606,'Z Primary Distribution',100, 40.,200.,0.)
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

		call hbook1(621,'W (GeV)  [c]^,!',    100,1.8,2.8,0.)
		call hbook1(622,'W (GeV)  x(1280)',      100,1.8,2.8,0.)
		call hbook1(623,'W (GeV)  [c]^,!',    10 ,1.8,2.8,0.)
		call hbook1(624,'W (GeV)  x(1280)',      10 ,1.8,2.8,0.)
		call hbook1(625,'E?[g]! (GeV)  [c]^,! "5# [cpp]',
     1         		                                 100,1.0,4.0,0.)
		call hbook1(626,'E?[g]! (GeV)  [c]^,! "5# [rg]', 
     1		                                         100,1.0,4.0,0.)
                call hbook1(627,'cos [q]^cm! [c]^,!',  10 ,-1.0,1.0,0.)
		call hbook1(628,'cos [q]^cm! x(1280)', 10 ,-1.0,1.0,0.)
		call hbook1(6260,'E?[g]! (GeV)', 
     1		                                         100,0.0,3.0,0.)



		call hbook2(700,'IM (K^+![p]^-!) vs. MM([g],K)',
     1                                     100,0.9,2.5,100,0.6,1.4,0.)	

                call hbook1(705,'[p]^+!  Momentum (MeV/c)',
     1  					      180,0.,1500.,0.)
	        call hbook1(709,'PROTON Momentum (MeV/c)' ,
     1		                                       60,0.,3000.,0.)
                call hbook1(714,'[p]^-!  Momentum (MeV/c)',
     1                                                180,0.,1500.,0.)
                call hbook1(715,'[p]^+!  Momentum (MeV/c)',
     1                                                180,0.,1500.,0.)
                call hbook1(716,'[g] pperp (GeV/c)',
     1                                                   100,0.,1.,0.)
	        call hbook1(719,'PROTON Momentum (MeV/c)' ,
     1		                                       60,0.,3000.,0.)
                call hbook1(724,'[p]^-!  Momentum (MeV/c)',
     1                                                180,0.,1500.,0.)
		call hbook2(770,'cos [Q]^c.m.!?&&&&&meson! vs. MM([g],p)',
     1                      240,0.6,1.8, 10,-1.0,1.0, 0.)                      
		call hbook2(771,'cos [Q]^c.m.!?&&&&&meson! vs. MM([g],p)',
     1                      240,0.6,1.8, 10,-1.0,1.0, 0.)                      
		call hbook2(772,'cos [Q]^c.m.!?&&&&&meson! vs. MM([g],p)',
     1                      240,0.6,1.8, 10,-1.0,1.0, 0.)                      
		call hbook2(773,'cos [Q]^c.m.!?&&&&&meson! vs. MM([g],p)',
     1                      240,0.6,1.8, 10,-1.0,1.0, 0.)                      
		call hbook2(774,'cos [Q]^c.m.!?&&&&&meson! vs. MM([g],p)',
     1                      240,0.6,1.8, 10,-1.0,1.0, 0.)                      
		call hbook2(775,'cos [Q]^c.m.!?&&&&&meson! vs. MM([g],p)',
     1                      240,0.6,1.8, 10,-1.0,1.0, 0.)                      
		call hbook2(776,'cos [Q]^c.m.!?&&&&&meson! vs. MM([g],p)',
     1                      240,0.6,1.8, 10,-1.0,1.0, 0.)                      
		call hbook2(777,'cos [Q]^c.m.!?&&&&&meson! vs. MM([g],p)',
     1                      240,0.6,1.8, 10,-1.0,1.0, 0.)                      
		call hbook2(778,'cos [Q]^c.m.!?&&&&&meson! vs. MM([g],p)',
     1                      240,0.6,1.8, 10,-1.0,1.0, 0.)                      
		call hbook2(779,'cos [Q]^c.m.!?&&&&&meson! vs. MM([g],p)',
     1                      240,0.6,1.8, 10,-1.0,1.0, 0.)                      
		call hbook2(780,'cos [Q]^c.m.!?&&&&&meson! vs. MM([g],p)',
     1                      240,0.6,1.8, 10,-1.0,1.0, 0.)                      


		call hbook1(800,'MM ([g]K^+![p]^+![p]^-!) (GeV/c^2!) '
     1  					      ,120,0.7,1.3,0.)
		call hbook1(801,'MM ([g]K^+![p]^+![p]^-!) (GeV/c^2!) (smeared)'
     1  						 ,120,0.7,1.3,0.)
		call hbook1(802,'IM (n)[p]^+! (GeV/c^2!) (smeared)'
     1  						 ,1000,1.0,1.5,0.)
		call hbook1(803,'IM (n)[p]^-! (GeV/c^2!) (smeared)'
     1  						 ,1000,1.0,1.5,0.)
		call hbook2(804,'IM (n)[p]^+! vs.(n)[p]^-!) (GeV/c^2!) (smeared)',
     1  	                                    100,1.0, 1.5,100,1.0,1.5,0.)
		call hbook1(805,'MM [g],K^+![p]^-! (GeV/c^2!) '
     1  						 ,1000,1.0,1.5,0.)
		call hbook1(806,'MM [g],K^+![p]^-! (GeV/c^2!) (smeared)'
     1  						 ,1000,1.0,1.5,0.)
		call hbook2(807,
     1		   'IM^2! (n)[p]^+! vs.(n)[p]^-! (GeV/c^2!)^2! (smeared)',
     1  	                                 100,1.0, 2.25,100,1.0,2.25,0.)
		call hbook2(808,
     1		   'IM [S]^+![p]^-! vs.K^+![p]^+! (GeV/c^2!)^2!',
     1  	                           50,0.55, 1.40,50,1.30,2.20,0.)
		call hbook2(809,
     1		   'IM [S]^-![p]^+! vs.K^+![p]^-! (GeV/c^2!)^2!',
     1  	                           50,0.55, 1.40,50,1.30,2.20,0.)
		call hbook2(810,
     1		   'IM^2! [S]^+![p]^-! vs.K^+![p]^+! (GeV/c^2!)^2!',
     1  	                           50,0.30, 2.90,50,1.7,5.8,0.)
		call hbook2(811,
     1		   'IM^2! [S]^-![p]^+! vs.K^+![p]^-! (GeV/c^2!)^2!',
     1  	                           50,0.30, 2.90,50,1.7,5.8,0.)
		call hbook2(81200,
     1		   'IM^2! [S]^+![p]^-! vs.[S]^+!K^+! (GeV/c^2!)^2!',
     1  	                           50,2.50, 7.50,50,1.50,6.00,0.)
		call hbook2(81300,
     1		   'IM^2! [S]^-![p]^+! vs.[S]^-!K^+! (GeV/c^2!)^2!',
     1  	                           50,2.50, 7.50,50,1.50,6.00,0.)
		call hbook1(8020,'IM ([L])[p]^+! (GeV/c^2!) (smeared)'
     1  						 ,1000,1.0,1.5,0.)
		call hbook1(8030,'IM ([L])[p]^-! (GeV/c^2!) (smeared)'
     1  						 ,1000,1.0,1.5,0.)
		call hbook2(8040,
     1		   'IM^2! ([L])[p]^+! vs.([L])[p]^-! (GeV/c^2!)^2! (smeared)',
     1  	                                    100,1.5, 2.1,100,1.5,2.1,0.)
c		call hbook2(8040,'IM ([L])[p]^+! vs.([L])[p]^-!) (GeV/c^2!) (smeared)',
c     1  	                                    100,1.0, 1.5,100,1.0,1.5,0.)
		call hbook1(8050,'MM [g],K^+![p]^-! (GeV/c^2!) '
     1  						 ,1000,1.0,1.5,0.)
		call hbook1(8060,'MM [g],K^+![p]^-! (GeV/c^2!) (smeared)'
     1  						 ,1000,1.0,1.5,0.)
		call hbook2(8070,
     1		   'IM^2! ([L])[p]^-! vs.([L])[p]^+! (GeV/c^2!)^2! (smeared)',
     1  	                                    100,1.5, 2.0,100,1.5,2.0,0.)
		call hbook2(8071,
     1		   'IM^2! ([L])[p]^-! vs.([L])[p]^+! (GeV/c^2!)^2! (LOW mass)',
     1  	                                    100,1.5, 2.0,100,1.5,2.0,0.)
		call hbook2(8072,
     1		   'IM^2! ([L])[p]^-! vs.([L])[p]^+! (GeV/c^2!)^2! (HIGH mass)',
     1  	                                    100,1.5, 2.0,100,1.5,2.0,0.)
		call hbook2(8073,
     1		   'IM (K^+![p]^-!) vs.MM([g],K^+!)',
     1  	                                    100,1.2, 2.5,100,0.5,2.0,0.)
		call hbook2(8074,
     1		   'IM (K^+![p]^-!) vs.MM(K^+! [p]^-!)',
     1  	                                    100,1.2, 1.8,100,0.6,1.3,0.)
		call hbook2(8075,
     1		   'IM (K^+![p]^-!) vs.MM(K^+! [p]^+!)',
     1  	                                    100,1.2, 1.8,100,0.6,1.3,0.)

		call hbook2(812,'All Events,  W bin vs. cos([Q]?K!)',
     1  	                                  20,-1.0,1.0,10,1.85,2.85,0.)
		call hbook2(8120,'All Events,  W bin vs. cos([Q]?K!)',
     1  	                                  20,-1.0,1.0,10,1.80,2.80,0.)
c     1  	                                  10,-1.0,1.0,14,1.8,3.2,0.)
		call hbook2(8121,'All Events, E[g] bin vs. cos([Q]?K0!) lab',
     1  	                                  20,-1.0,1.0,10,0.800,2.280,0.)
		call hbook2(8122,'All Events, E[g] bin vs. cos([Q]?K0!) lab',
     1  	                                  20,-1.0,1.0,10,0.440,2.280,0.)
		call hbook2(8123,'All Events, E[g] bin vs. cos([Q]?K0!) cm',
     1  	                                  20,-1.0,1.0,10,0.440,2.280,0.)
		call hbook2(813,'All Events,  W bin vs. MM([g],K^+!)',
     1  	                                  80,1.3,1.7, 14,1.8,3.2,0.)
		call hbook2(814,'All Events,  cos([Q]?K!) vs. MM([g],K^+!)',
     1                                            80,1.3,1.7, 20,-1.0,1.0,0.)
c		call hbook2(839,'Case 2, [L] vs. [S]^+!',
c     1  	                                  50,1.0,1.5,50,1.0,1.5,0.)

c		call hbook1(900,'FERMI RECOIL MOMENTUM'  ,100,0.0,500.,0.)
c		call hbook1(901,'EFFECTIVE RECOIL MASS'  ,600,000.,3000.,0.)
c		call hbook1(902,'RANDOM NUMBER'  ,400, -.1,1.1,0.)
c		call hbook1(903,'NORMAL RANDOM NUMBER'  ,200,-4.,4.,0.)
c		call hbook1(904,'B.W.  RANDOM NUMBER'  ,500,-12.,12.,0.)
c		call hbook1(9040,'B.W.  RANDOM NUMBER'  ,500,-12.,12.,0.)

c
c               Stuff for Lambda(1520)
c

		call hbook1(1530,'MM ([g],pK^+!K^-!) (GeV/c^2!)'
     1  						 ,100,-0.01,0.01,0.)
		call hbook1(1531,'MM ([g],pK^+!K^-!) (GeV/c^2!) (smeared)'
     1  						 ,200,-0.01,0.01,0.)
		call hbook1(1532,'MM ([g]K^+!) (GeV/c^2!)'
     1  						 ,320,0.9,2.5,0.)
		call hbook1(1533,'MM ([g]K^+!) (GeV/c^2!) (smeared)'
     1  						 ,320,0.9,2.5,0.)
		call hbook1(1534,'W (GeV) for [L](1520)', 100 ,1.8,3.2,0.)
c
c               Histograms related to charm search
c
                call hbook1(7001,'IM K^+! [p]^-! (D^0!bar)'    ,100,1.2,2.2,0)
                call hbook1(7003,'IM K^+! [p]^-! (D^0!bar)'    ,100,1.2,2.2,0)

                call hbook1(7002,'IM p [K]^-! [p]^+! ([L]?c^+!)',100,2.0,2.5,0)
                call hbook1(7004,'IM p [K]^-! [p]^+! ([L]?c^+!)',100,2.0,2.5,0)

                call hbook1(7005,'IM [K]^+! [p]^-! [p]^-! (D^-!)',100,1.2,2.2,0)
                call hbook1(7006,'IM [K]^+! [p]^-! [p]^-! (D^-!)',100,1.2,2.2,0)
                call hbook1(7007,'IM [K]^+! [p]^-! [p]^-! (D^-!)',100,0.0,2.3,0)

                call hbook1(7008,'Energy Conservation',100,-1.0,1.0,0)

                call hbook1(7009,'IM p [K]^-! [p]^+! ([L]?c^+!)',100,2.0,2.5,0)
                call hbook1(7010,'IM [K]^+! [p]^-! [p]^-! (D^-!)',100,1.2,2.2,0)

                call hbook1(7011,'IM e^+! [p]^-! e^-! (J/psi)',100,2.7,3.5,0)

                call hbook1(7023,'IM K^+! [p]^-! (D^0!bar)'    ,100,1.2,2.2,0)
                call hbook1(7024,'IM p [K]^-! [p]^+! ([L]?c^+!)',100,2.0,2.5,0)

                call hbook1(7033,'IM K^+! [p]^-! (D^0!bar)'    ,100,1.2,2.2,0)
                call hbook1(7034,'IM p [K]^-! [p]^+! ([L]?c^+!)',100,2.0,2.5,0)
                call hbook1(7035,'IM [g][g]([p]?0!)',100,0.0,0.2,0)
                call hbook1(7036,'IM p [K]^-! [p]^+! ([L]?c^+! [p]^0!)'
     1                                                          ,100,2.1,2.6,0)

                call hbook1(7043,'IM K^+! [p]^-! [p]^-! (D^-!)'    
     1                                                          ,100,1.2,2.2,0)
                call hbook1(7044,'IM p [K]^-! [p]^+! ([L]?c^+!)',100,2.0,2.5,0)
                call hbook1(7046,'IM p [K]^-! 2[p]^+! ([L]?c^+! [p]^+!)'
     1                                                          ,100,2.1,2.6,0)


                call hbook2(7100,
     1                 'IM p K^-! [p]^+! vs. IM [K]^+! [p]^-! [p]^-!',
     1                                 100, 2.0, 2.5, 100, 1.5, 2.3,0)
                call hbook2(7101,
     1                 'IM p K^-! [p]^+! vs. IM [K]^+! [p]^-!',
     1                                 100, 1.2, 2.2, 100, 2.0, 2.5,0)
                call hbook2(7102,
     1                 'IM p K^-! [p]^+! vs. IM [K]^+! [p]^-!',
     1                                 100, 1.2, 2.2, 100, 2.0, 2.5,0)
                call hbook2(7103,
     1                 'IM p K^-! [p]^+! vs. IM [K]^+! [p]^-! [p]^-!',
     1                                 100, 1.2, 2.2, 100, 2.0, 2.5,0)
		call hbook2(9000,'[p][p] relative momentum vs. cos [p][p] in x frame'
     1                                ,100,-1.0,1.0,60,0.0,1.5,0)
		call hbook2(9001,
     1       '[p][p] relative momentum vs. cos [p][p] in [p][p] frame'
     1                                ,100,-1.0,1.0,60,0.0,1.5,0)
		
c
c------Histograms for g13 and g14 dibaryon study --------------

		call hbook1(8801, 'IM (d[p]^+!) (GeV/c^2!) '
     1  						 ,100,2.0,3.0,0.)
		call hbook1(8802, 'IM (d[p]^-!) (GeV/c^2!) '
     1  						 ,100,2.0,3.0,0.)
		call hbook2(8803,'IM (d[p]^+!) vs IM (d[p]^-!)'
     1                                ,50,2.0,3.0,50,2.0,3.0,0)
		call hbook1(8804, 'IM ([p]^+![p]^-!) (GeV/c^2!) '
     1  						 ,100,0.25,1.25,0.)


	      call hbook1(8601,'Deuteron Momentum (GeV/c)'       ,100,0.,1.0,0.)
	      call hbook1(8611,'d c.m. angle'           ,40,-1.0,+1.0,0.)
	      call hbook1(8612,'[p]^+! c.m. angle'      ,40,-1.0,+1.0,0.)
	      call hbook1(8613,'[p]^-! c.m. angle'      ,40,-1.0,+1.0,0.)
	      call hbook2(8614,'c.m. cos [Q] [p]^-! vs cos [Q] [p]^+!'
     1                                ,50,-1.0,1.0,50,-1.0,1.0,0)
	      call hbook2(8615,'Momentum c.m. [p]^-! vs [p]^+!'
     1                                ,50,0.0,1.0,50,0.0,1.0,0)
	      call hbook2(8616,'c.m. cos [Q] [p]^+-! vs cos [Q] d'
     1                                ,50,-1.0,1.0,50,-1.0,1.0,0)
	      call hbook2(8617,'E[g] vs. [p][p] mass'
     1                                ,50,0.25,1.25,50,0.4,2.5,0)
	      call hbook2(8618,'Deuteron mtm, lab vs. cm'
     1                                ,50,0.0,1.0,50,0.0,2000.0,0)

c
c------Histograms for pi- p in g14 study --------------

	      call hbook1(2626,'Beam Momentum (GeV)'   ,50, 0.5, 3.0, 0.)
	      call hbook1(2650,'Spectator Mass (GeV/c^2!) '
     1  						 ,200,0.8,1.2,0.)
	      call hbook1(3650,'Spectator Mass (GeV/c^2!) (Smeared)'
     1  						 ,200,0.8,1.2,0.)
	      call hbook1(2651,'Spectator Momentum (GeV/c) '
     1  						 ,200,0.0,1.0,0.)
	      call hbook1(3651,'Spectator Momentum (GeV/c^2!) (Smeared)'
     1  						 ,200,0.0,1.0,0.)
	      call hbook2(3122,'All Events, E[g] bin vs. cos([Q]?[p]-!)',
     1  	                                  20,-1.0,1.0,10,0.440,2.280,0.)
c
c --------Histograms for semi-leptonic decay study----------

	      call hbook1(632,'MM^2! ([g]K^+!p l) (GeV/c^2!)^2!'
     1  						 ,200,-0.004,0.004,0.)
	      call hbook1(633,'MM^2! ([g]K^+!p l) (GeV/c^2!)^2! (smeared)'
     1  						 ,200,-0.004,0.004,0.)
	      call hbook1(634,'Neutrino momentum (GeV)'
     1  						 ,100,0.0,0.2,0.)
	      call hbook1(6340,'Neutrino momentum (GeV) (smeared)'
     1  						 ,100,0.0,0.2,0.)
	      call hbook1(6341,'Nothing  momentum (GeV) (smeared)'
     1  						 ,100,0.0,0.2,0.)
	      call hbook1(6342,'cos[q] proton lepton ([L] frame)'
     1  						 ,100,-1.0,1.0,0.)
	      call hbook1(6343,'cos[q] proton pion ([L] frame)'
     1  						 ,100,-1.0,1.0,0.)
	      call hbook1(6320,'MM^2! ([g]K^+!p[p]) (GeV/c^2!)^2!'
     1  						 ,200,-0.01,0.01,0.)
	      call hbook1(6330,'MM^2! ([g]K^+!p[p]) (GeV/c^2!)^2! (smeared)'
     1  						 ,200,-0.01,0.01,0.)
	      call hbook1(635,'Neutrino momentum (GeV) (smeared)'
     1  						 ,100,0.0,0.2,0.)
	      call hbook2(636,'Momentum vs. mass^2! of missing part. (smeared)',
     1  	                          50,-0.01,0.01,50,0.0,0.5,0.)
	      call hbook1(637,'Nothing momentum (GeV)'
     1  						 ,100,0.0,0.2,0.)
	      call hbook1(638,'Nothing momentum (GeV) (smeared)'
     1  						 ,100,0.0,0.2,0.)
	      call hbook2(639,'Momentum vs. mass^2! of Nothing (smeared)',
     1  	                          50,-0.01,0.01,50,0.0,0.5,0.)
	      call hbook1(640,'MM^2!([g], K^+!p) (smeared)',
     1  	                          50,-0.05,0.05,0.)
	      call hbook1(641,'MM^2!([g], K^+!p) (smeared)',
     1  	                          50,-0.05,0.05,0.)
	      call hbook2(642,'leptonic mass vs. 0  mass (smeared)',
     1  	                          50,-.01,.01, 50,-0.05,0.05,0.)
	      call hbook2(643,'[p] mass vs. 0  mass (smeared)',
     1  	                          50,-.01,.01, 50,-0.05,0.05,0.)
	      call hbook1(1038,'Proton Momentum (GeV/c)'  ,100,0.,0.25,0.)
	      call hbook1(1039,'Proton Momentum (GeV/c)'  ,100,0.,0.25,0.)
	      call hbook1(1040,'Lepton Momentum (GeV/c)'  ,100,0.,0.2,0.)
	      call hbook1(1041,'[p]^-!  Momentum (GeV/c)' ,100,0.,0.2,0.)
	      call hbook1(1042,'HYP-frame Opening angle'  ,100,-1.0,1.0,0.)
	      call hbook2(2040,'Lepton cos[Q] vs. p (GeV/c)'       
     1  					,50,0.,0.2,50,-1.0,1.0,0.)
	      call hbook2(2041,'[p]^-! cos[Q] vs. p (GeV/c)'       
     1  					,50,0.,0.2,50,-1.0,1.0,0.)
	      call hbook2(2042,'CM coplanarity vs. p (GeV/c)'       
     1  					,50,0.,0.2,50,-1.0,1.0,0.)
	      call hbook2(2050,'HYP frame plepton vs costh l[m]'       
     1  					,50,-1.0,1.0,50,0.0,0.15,0.)
	      call hbook2(2043,'CM coplanarity vs. p (GeV/c)'       
     1  					,50,0.,0.2,50,-1.0,1.0,0.)
	      call hbook2(6300,'[D] TOF vs. p (GeV/c) for muons'       
     1  			            ,100,0.0,0.6,120,-5.0,5.0,0.)
	      call hbook2(6301,'[D] TOF vs. p (GeV/c) for electrons'       
     1  			            ,100,0.0,0.6,120,-5.0,5.0,0.)
	      call hbook1(6302,'Invariant Mass p l (GeV/c)^2!' ,
     1                                              100,1.0,1.15,0.)
	      call hbook1(6303,'Invariant Mass p [p] (fake) (GeV/c)^2!' ,
     1                                              100,1.0,1.15,0.)
	      call hbook1(6304,'Invariant Mass p [p] (GeV/c)^2!' ,
     1                                              100,1.0,1.15,0.)
	      call hbook1(6305,'Invariant Mass p l (fake) (GeV/c)^2!' ,
     1                                              100,1.0,1.15,0.)
c
c------Histograms for GlueX baryonium study (p pbar)--------------
c
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
	      call hbook2(8855,'Van Hove Plot'
     1                                ,100,-1.5,1.5,100,-1.5,1.5,0)
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
	      
c	      call hbook1(9801,'Lambda IM (p [p]^-!) ',
c	1	                         100,1.05,1.15,0)
c	      call hbook1(9802,'Anti-Lambda IM (p&^-! [p]^+!) ',
c	1	                         100,1.05,1.15,0)
	      call hbook1(9801,'Lambda IM (p [p]^-!) ',
	1	                         100,1.05,1.50,0)
	      call hbook1(9802,'Anti-Lambda IM (p&^-! [p]^+!) ',
	1	                         100,1.05,1.50,0)
	      call hbook1(9851, 'IM ([L]&^-![L]) (GeV/c^2!) ',
	1	                         100,2.0,4.0,0.)
	      call hbook1(9852, 'IM ([L]&^-!p) (GeV/c^2!) ',
	1	                         100,2.0,4.0,0.)
	      call hbook1(9856, '-(t - tmin) (GeV^2!)  [L]&^-! [L]',
	1	                         100,0.0,3.0,0.)
	      call hbook1(9857, '-(t - tmin) (GeV^2!)  [L]&^-! p',
	1	                         100,0.0,3.0,0.)

	      call hbook2(9880,'p?cm! (GeV/c) anti-Lambda vs. Lambda',
	1    	                         100,0.0,2.0,100,0,2.,0.)

	      call hbook2(9881,'cos [Q]?cm!  vs. p?cm!   Lambda',
	1    	                         100,0.0,2.0,100,-1.0,+1.0,0.)
	      call hbook2(9883,'cos [Q]?cm!  vs. p?cm!   anti-Lambda',
	1    	                         100,0.0,2.0,100,-1.0,+1.0,0.)
	      call hbook2(9884,'cos [Q]?cm! vs. p?cm! ([L]&^-! [L] Pair)',
	1    	                         100,0.0,2.0,100,-1.0,+1.0,0.)

	      
	      call hbook2(9869,'cos [Q]?cm! anti-Lambda vs. Lambda',
	1    	                       50,-1.,1.0,50,-1.,1.0,0.)

	      call hbook2(9850,'IM (p [L]&^-!) vs IM (p [L])'
	1	                            ,100,2.0,4.0,100,2.0,4.0,0)
	      call hbook2(9855,'IM (p [L]&^-!) vs IM ([L]&^-! [L])'
	1	                            ,100,2.0,4.0,100,2.0,4.0,0)
	      call hbook2(9854,'IM^2! (p [L]&^-!) vs IM^2! (p [L])' !Dalitz
	1       	                    ,100,3.0,11.0,100,3.0,11.0,0)
		
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
	      call hbook1(10887,'Missing Mass^2! using K^+! mass (GeV)',
c     1	                        100,-0.02,0.02,0.)
     1	                        100,-3.0,3.0,0.)
	      call hbook1(10888,'Missing Mass^2! using [p]^+! mass (GeV)',
c     1	                        100,-0.02,0.02,0.)
     1	                        100,-3.0,3.0,0.)
	      call hbook2(10889,
     1	            'Missing Mass^2! using [p]^+! vs. K^+! mass (GeV)',
     1	                        100,-3.0,3.0,100,-3.0,3.0,0.)
	      call hbook1(10890,'Missing Mass using K^+! mass (GeV)',
     1	                        100,-0.5,0.5,0.)
c	      call hbook1(10100,'IM(p pi-) ([L])'     ,100 ,1.0,1.2,0.)
	      call hbook1(10100,'IM(p pi-) ([L])'     ,200 ,1.0,1.6,0.)
	      

	      
 20	   continue
	end if    !histogram file opening
c
c       In Fall 2015, we wanted to write out events for the eta pi+ pi- dalitz plot event by event
c
c	write(6,*)'Opening file for saving Dalitz plot data for x(1280)'
c	open(unit=4, file= 'mc_dalitz_data.txt',status='unknown') !for x(1280) study
c
	end if    !first call to analyzer

c
c	Top of loop to decide whether this event should be histogrammed at all
c	Do the experiment specific cuts, or skip to line 60
c	for full phase-space acceptance
c
	iwriteit = .false.    !Don't write out this event upon return to main program
	icount = icount + 1   !Counts how many events are examined
c	write(6,*)'Call to analyzer',icount
c	call prevent(plist,imark+1,ikeep)   !temporary	
	if(idocuts.eq.0)then
	   ikeep = ikeep + 1
	   goto 2000		!Bail out of doing data analysis here
	endif
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
c	   do i=1,imark
c	      lund = plist(6,i)
c	      if(lund.eq.+18 .or. lund.eq.-18)then
c		 iskipit = 1
c	      endif
c	   enddo
c	   if(iskipit.eq.0)then
c	      x = ran(iseed)
c	      if(x.gt.0.5)then
cc		 write(6,*)'> Stepping up'
c		 imin  = 1
c		 imax  = imark
c		 istep = +1
c	      else
cc		 write(6,*)'> Stepping down'
c		 imin  = imark
c		 imax  = 1
c		 istep = -1
c	      endif
c	      do i=imin,imax,istep
c		 lund = plist(6,i)
c		 if(lund.eq.+17 .and. ipipcount.eq.0)then
c		    lund = +18
c		    plist(6,i)  = lund
c		    plist(4,i)  = 0.49368 !Hack to test PID errors
c		    plist(5,i)  = plist(5,i)
cc		    plist(5,i)  = etot(plist(1,i))
c		    plist(16,i) = plist(4,i) ! copy of mass; NO smearing done
c		    plist(17,i) = plist(17,i) !total energy if wrong mass assumed 
cc		    plist(17,i) = etot(plist(13,i)) !total energy if wrong mass assumed 
c		    ipipcount = ipipcount + 1
cc                    write(6,*)'> Flipped pi+ to K+'
c c                   write(6,*)lund,(plist(iii,i),iii=1,22)
c  c                  call exit()
c		 endif
c		 if(lund.eq.-17 .and. ipimcount.eq.0)then
c		    lund = -18
c		    plist(6,i) = lund
c		    plist(4,i)  = 0.49368 !Hack to test PID errors
c		    plist(5,i)  = plist(5,i)
cc		    plist(5,i)  = etot(plist(1,i))
c		    plist(16,i) = plist(4,i) ! copy of mass; NO smearing done
c		    plist(17,i) = plist(17,i) 
c		    plist(17,i) = etot(plist(13,i)) 
c		    ipimcount = ipimcount + 1
cc                    write(6,*)'> Flipped pi- to K-'
c		 end if
c	      enddo
c	   endif
c	endif

c	if(idocuts.ne.1)goto 60	!Bail out of doing data cuts here
c
	pxtot = 0.
	pytot = 0.
	pztot = 0.
	acceptlist(1) = -1 !just to initialize value; no meaning intended here
cx	call prevent(plist,imark+1,ikeep)
	do 50 i = 2,imark
	lund = plist(6,i)
c
c	Get lab angles, in degrees, and other particle parameters
c
c       Take values from formation time
c	pmom = sqrt(plist(5,i)**2 - plist(4,i)**2) * 1000.	!MeV/c
c	px = plist(1,i) * 1000.	!MeV/c
c	py = plist(2,i) * 1000.
c	pz = plist(3,i) * 1000.
c
c       Take values from AFTER energy losses were computed (in mc_generator)
c       Instrumental smearing was also done in mc_generator
c
	icharge = plist(22,i)
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
		theta_xz = 1000.
		theta_yz = 1000.
		phi = 1000.
	end if	
	xx   = plist(18,i) - plist(7,i)
	yy   = plist(19,i) - plist(8,i)
	zz   = plist(20,i) - plist(9,i)
	path = 1.e15
c	vertex(1) = 0 ! plist(7,i)  !for track swimming, X and Y must be zero...
c	vertex(2) = 0 ! plist(8,i)
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
c       
c	write(6,*)'>Calling acceptance',lund,theta_polar,phi,pmom,icharge
	iacc = iaccept_GlueX(theta_polar,phi,pmom,icharge,vertex,path,
     1	                                          betagammactau,lund)
c	write(6,*)'Single track>',icount,iacc,lund,acceptweight,
c     1        theta_polar,phi,pmom,icharge

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
c
c	if(lund.eq.18)then              !kaons
c	   isectork     = isector  	!save for later plotting
c	   numtofbark   = numtofbar     !ditto
c	   pmomk        = pmom
c	   theta_polark = theta_polar
c	   phik         = phi
c	   pxk          = px
c	   pyk          = py
c	   pzk          = pz
c	   else if(lund.eq.-18)then	!kaon minuses
c	   isectorkm    = isector  	!save for later plotting
c	   numtofbarkm  = numtofbar     !ditto
c	   pmomkm       = pmom
c	   theta_polarkm= theta_polar
c	   phikm        = phi
c	   pxkm         = px
c	   pykm         = py
c	   pzkm         = pz
c	else if(lund.eq.41)then         !protons
c	   isectorp     = isector 	!save for later plotting
c	   numtofbarp   = numtofbar     !ditto
c	   pmomp        = pmom
c	   theta_polarp = theta_polar
c	   phip         = phi
c	else if(lund.eq.-17)then        !pions
c	   isectorpim    = isector 	!save for later plotting
c	   numtofbarpim  = numtofbar    !ditto
c	   pmompim       = pmom
c	   theta_polarpim= theta_polar
c	   phipim        = phi
c	else if(lund.eq.+17)then        !pions
c	   isectorpi     = isector 	!save for later plotting
c	   numtofbarpi   = numtofbar    !ditto
c	   pmompip       = pmom
c	   theta_polarpi = theta_polar
c	   phipip        = phi
c	else if(lund.eq.1)then          !photons
c	   isectorph     = isector 	!save for later plotting
c	   pmomph        = pmom
c	   theta_polarph = theta_polar
c	   phiph         = phi
c	endif
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
c       Software Trigger:
c
	iprocess  = 0		!flag to say whether this event is to be processed
	itrigmask = 1  ! event is accepted in any case  -  hard code this now now
	itrigmask = 10 ! event is accepted is we have p, pbar, p   -  hard code this for now
	itrigmask = 11 ! event is accepted is we have p, pbar, p, pi+, pi-   -  hard code this for now

	if(iexpt_trigger)then
	   if(itrigmask.eq.1)then
	      iprocess = 1
	   elseif((iacc_proton.and.iacc_proton2.and.iacc_antiproton) !for GlueX p pbar reaction
     1 	      .and. itrigmask.eq.10)then
	      iprocess = 1
	   elseif((iacc_proton.and.iacc_proton2.and.iacc_antiproton.and.
     1	      iacc_piplus.and.iacc_piminus) !for GlueX Lambda anti-Lambda reaction
     1	      .and. itrigmask.eq.11)then
	      iprocess = 1
	      write(6,*)'Triggermask got to here....'
	   endif
	endif
c
	if(iprocess.eq.1)then
	   ikeep = ikeep + 1
	   if(ikeep.le.1)then
	      write(6,*)'The first event(s) to pass cuts:'
	      call prevent(plist,imark+1,ikeep)
	      write(6,*)' '
	   end if
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
c	write(6,*)'>>>Top of looping',iloop,itrigmask,weight
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
     1	      accept_weight_proton2 *
     1        acceptweight_antiproton 
	   elseif(itrigmask.eq.11)then
	      weight = acceptweight_proton *
     1	      accept_weight_proton2 *
     1        acceptweight_antiproton * 
     1        acceptweight_piplus * 
     1        acceptweight_piminus
	      ifiddle = .false.  
	   end if
	endif
	
c	do 100 i=1,imark
	do 100 i=2,imark  !skip the parent particle
	lund = plist(6,i)
c
c	Formation-time variables
c
	rmas = plist(4,i) * 1000.                               !MeV
	eeee = plist(5,i) * 1000.				!MeV
	tkin = eeee - rmas					!MeV
	pmom = sqrt(eeee**2 - rmas**2)				!MeV/c
	pmomgev = pmom/1000.
	px = plist(1,i) * 1000.	!MeV/c
	py = plist(2,i) * 1000.
	pz = plist(3,i) * 1000.
c
c	Decay-time variables
c
	edec = plist(17,i) * 1000.				!MeV
	tdec = edec - rmas					!MeV
	pdec = sqrt(edec**2 - rmas**2)				!MeV/c
	pxdec = plist(13,i) * 1000.	!MeV/c
	pydec = plist(14,i) * 1000.
	pzdec = plist(15,i) * 1000.
c
c       Force rest of this code to use smeared variables unless we're
c       on the first loop through when we want to use "generated" 
c       values of the variables.
c
	if(iloop.gt.1)then
	   eeee = edec
	   tkin = tdec
	   pmom = pdec
	   px   = pxdec
	   py   = pydec
	   pz   = pzdec
	endif
c
	if(pmom.ne.0)then
		arg=pz/pmom
		if(arg.gt. 1.0)arg= 1.0
		if(arg.lt.-1.0)arg=-1.0
		theta_polar = acos(arg) * rade	!lab angles, in degrees
		theta_xz = atan2(px,pz) * rade
		theta_yz = atan2(py,pz) * rade
		phi = atan2(py,px) * rade
		if(phi.lt.0.0)phi = phi + 360.  !azimuthal angle from 0 to 360.
		costhetapolar = arg
	else
		theta_polar = 1000.
		theta_xz = 1000.
		theta_yz = 1000.
		phi = 1000.
		costhetapolar = -999.
	end if
c
c	Distance from center of event at which particle decays
c
	rx = plist(18,i)
	ry = plist(19,i)
	rz = plist(20,i)
	big = 1.e15
	if(abs(rx).lt.big .and. abs(ry).lt.big .and. 
     1                                        abs(rz).lt.big)then
		rr = sqrt(rx*rx + ry*ry + rz*rz)
	else
		rr = 1.e32
	end if
	vertex(1) = plist(7,i)  !was not used for swimming track
	vertex(2) = plist(8,i)
	vertex(3) = plist(9,i)  !WAS used for swimming track
c
c
c	Do all the routine histograms
c
c	call hf1(28,ppbeam,weight)  !beam momentum
	if(         lund.eq.18)then			!K+
		call hf1(7,tkin,weight)
		call hf1(107,pmomgev,weight)
c		call hf1(107,pdec,weight)
		call hf1(57,arg,weight)
		call hf1(37,rr,weight)
		call hf2(207,pmomgev,theta_polar,weight)
		call hf2(2070,theta_polar,pmom,weight)
		call hf2(2071,ppbeam,pmomgev,weight) 
	     else if(lund.eq.41)then                       !proton
		call hf1(9,tkin,weight)
		call hf1(109,pmomgev,weight)
c		call hf1(39,rr,weight)
		call hf1(59,costhetapolar,weight)
		call hf2(209,pmomgev,theta_polar,weight)
		call hf2(2090,pmomgev,theta_polar,weight)
c		call hf2(1209,pdec,theta_polar,weight)
		call hf2(1209,theta_xz,theta_yz,weight)
		call hf2( 223,pmom,theta_polar,weight)
		if(abs(xx).le.1.0)call hf1(99,xx,weight)
	     else if(lund.eq.-41)then                     !anti-proton
		call hf1(130,pmomgev,weight)
		call hf2(230,pmomgev,theta_polar,weight)
	     else if(abs(lund).eq.42)then                !neutron
		call hf1(10,tkin,weight)
		call hf1(110,pmomgev,weight)
c		call hf1(40,rr,weight)
		call hf1(60,costhetapolar,weight)
		call hf2(210,pmomgev,theta_polar,weight)
	     else if(abs(lund).eq.57)then                !Lambda
		call hf1(11,tkin,weight)
		call hf1(111,pmom,weight)
		call hf1(41,rr,weight)
		call hf1(61,costhetapolar,weight)
		call hf2( 211,pmom,theta_polar,weight)
		call hf2(1211,pdec,theta_polar,weight)
c               Save generated lambda data for comparison with GSIM
		plgen(1) = plist(5,i)
		do k=1,4
		   plgen(k+1) = plist(k,i)
		enddo
	     else if(    lund .eq.17)then                 !pi+
		call hf1(5,tkin,weight)
		call hf1(105,pmomgev,weight)
		call hf1(35,rr,weight)
		call hf1(55,arg,weight)
		call hf2(205,pmomgev,theta_polar,weight)
	     else if(    lund .eq.-17)then
c		call hf1(14,tkin,weight)
		call hf1(114,pmomgev,weight)
		call hf1(35,rr,weight)
		call hf1(55,arg,weight)
		call hf2(214,pmomgev,theta_polar,weight)
		pmomgev_pion = pmomgev  !save for later cross check
	     else if(abs(lund).eq.23)then                !pi_0
		call hf1(6,tkin,weight)
		call hf1(106,pmomgev,weight)
		call hf1(36,rr,weight)
		call hf1(56,arg,weight)
		call hf2(206,pmomgev,theta_polar,weight)
	     else if(abs(lund).eq. 1)then
		call hf1(1,tkin,weight)
		call hf1(101,pmomgev,weight)
		call hf1(51,arg,weight)
		call hf2(201,pmomgev,theta_polar,weight)
	     else if(abs(lund).eq. 2)then
		call hf1(2,tkin,weight)
		call hf1(102,pmom,weight)
c		call hf1(32,rr,weight)
	     else if(abs(lund).eq. 7)then
		call hf1(3,tkin,weight)
		call hf1(103,pmomgev,weight)
		call hf1(33,rr,weight)
		call hf2(203,pmomgev,theta_polar,weight)
c	else if(abs(lund).eq. 9)then
c		call hf1(4,tkin,weight)
c		call hf1(104,pmom,weight)
c		call hf1(34,rr,weight)
c		call hf1(54,arg,weight)
	     else if(abs(lund).eq.24)then               !eta
		call hf1(1006,tkin,weight)
		call hf1(1106,pmom,weight)
		call hf1(1036,rr,weight)
		call hf1(1056,arg,weight)
		call hf2(1206,pmom,theta_polar,weight)
	     else if(    lund .eq.19)then               !K0
		call hf1(8,tkin,weight)
		call hf1(108,pmom,weight)
		call hf1(38,rr,weight)
		call hf2(208,pmom,theta_polar,weight)
	     else if(    lund .eq.20)then !K*(892)
		call hf1(25,tkin,weight)
		call hf1(125,pmom,weight)
		call hf1(126,rmas,weight)
		call hf2(225,pmom,theta_polar,weight)
	     else if(    lund .eq.-18)then               !K-
		call hf1(15,tkin,weight)
		call hf1(115,pmomgev,weight)
		call hf1(37,rr,weight)
		call hf2(215,pmomgev,theta_polar,weight)
c	else if(abs(lund).eq.19)then
c		call hf1(8,tkin,weight)
c		call hf1(108,pmom,weight)
c		call hf1(38,rr,weight)
c	else if(abs(lund).eq.46)then			!Cascade -
c		call hf1(  12,tkin,weight)
c		call hf1(1012,tdec,weight)
c		call hf1( 112,pmom,weight)
c		call hf1(1112,pdec,weight)
c		call hf1(42,rr,weight)
c		call hf2(1042,rx,ry,weight)
c		call hf2(2042,rz,ry,weight)
c		if(tdec.eq.0)then
c			call hf2(3042,rx,ry,weight)
c			call hf2(4042,rz,ry,weight)
c		end if
c		call hf2( 212,pmom,theta_polar,weight)
c		call hf2(1212,pdec,theta_polar,weight)
c	else if(abs(lund).eq.47)then			!Cascade 0
c		call hf1(13,tkin,weight)
c		call hf1(113,pmom,weight)
c		call hf1(43,rr,weight)
c		call hf2(213,pmom,theta_polar,weight)
c	else if(abs(lund).eq.49)then			!Cascade(1530)
c		call hf1(24,tkin,weight)
c		call hf1(124,pmom,weight)
c		call hf2(224,pmom,theta_polar,weight)
c	else if(abs(lund).eq.69)then			!f0
c		call hf1(19,tkin,weight)
c		call hf1(119,pmom,weight)
c		call hf2(219,pmom,theta_polar,weight)
c	else if(abs(lund).eq.35)then			!phi
c		call hf1(20,tkin,weight)
c		call hf1(120,pmom,weight)
c		call hf2(220,pmom,theta_polar,weight)
	     else if(abs(lund).eq.43)then !Sigma +
		call hf1(16,tkin,weight)
		call hf1(116,pmom,weight)
		call hf2(216,pmom,theta_polar,weight)
	     else if(abs(lund).eq.44)then !Sigma 0
		call hf1(17,tkin,weight)
		call hf1(117,pmom,weight)
		call hf2(217,pmom,theta_polar,weight)
	     else if(abs(lund).eq.45)then !Sigma -
		call hf1(18,tkin,weight)
		call hf1(118,pmom,weight)
		call hf2(218,pmom,theta_polar,weight)
c	else if(abs(lund).eq.51)then			!L3H
c		call hf1(21,tkin,weight)
c		call hf1(121,pmom,weight)
c		call hf2(221,pmom,theta_polar,weight)
	     else if(abs(lund).eq.52)then               !Deuteron
		call hf1(22,tkin,weight)
		call hf1(122,pmom,weight)
		call hf2(222,pmom,theta_polar,weight)
		pdeuteronlab = pmom
	end if
c
	call hf1(900,prec*1000.,weight)
	call hf1(901,qfmass*1000.,weight)
100	continue
ccc
ccc       Write special output file to compare mc_generator with GSIM
ccc
cc	if(iloop.eq.1)then
cc	   if(iwinner.eqv..true.)then
cccc	      call clas_out(ppbeam,pkgen,plgen)  !save events that pass cuts
cc	   else
cc	      call clas_out(ppbeam,pkgen,plgen)  !save events that fail cuts
cc	   endif
cc	endif
c
c       Once per event, sample random numbers
c
	call hf1(902,xflat,weight)
	call hf1(903,xnormal,weight)
	call hf1(904,xbw,weight)
c
c       Special code for pi+ photoproduction
c
	do 105 i=2,imark
	lundi = plist(6,i)
	if(lundi.eq.17)then	!pi+
c	   if(i.eq.imark)write(6,*)'Bloody murder!! a',lundi,imark
c
c          Velocity of lab frame relative to CM
c          Boost the pion into the total system c.m. frame
c          Compute angle between c.m. pion and beam axis
c       
c	   e0 = etot(plist(1,1))	!defines the beam plus target
	   e0 = plist(5,1)	!defines the beam plus target
	   if(e0.le.0)then
	      write(6,*)'How the heck did this happen?',
     1	                                           e0,plist
	      call exit()
c              e0 = 1.e32
	   end if
	   do j=1,3
	      beta(j) = plist(j,1)/e0
	   end do
	   call boost(plist(1,i),beta,ptemp)
	   costhpi = dot(ptemp,nbeam)/sqrt(psq(ptemp))
	   thmeson = acos(costhpi)*rade
c            write(6,*)costhpi,ptemp,nbeam
	   phimeson   = atan2(ptemp(2),ptemp(1)) * rade
	   if(phimeson.lt.0.0)phimeson = phimeson + 360. !azimuthal angle from 0 to 360.
	   phisector = mod(phimeson+30.,60.)-30. !phi within one sector +- 30 deg.
c       
c          Compute missing mass off pion, assuming we have a
c          two-body final state.  
c       
	   do k=1,4
	      ptemp(k) = plist(12+k,i)
	   enddo
	   call missingmass(plist(1,1),ptemp,rmiss,rmiss2)
	   pmom = sqrt(plist(17,i)**2 - plist(16,i)**2) * 1000.
c
c	   call hf1(28,ppbeam,weight)
	   call hf1(29,costhpi,weight)
	   call hf1(517,rmiss,weight)
	   call hf2(2226,ppbeam,costhpi,weight) 
	   call hf2(228 ,thmeson,phimeson,weight) !mixed c.m. and lab angles
	endif
 105	continue
c
c       Code for events with protons 
c       (eg pi0 and eta photoproduction)
c
	do 107 i=2,imark
	lundi = plist(6,i)
	if(lundi.eq.41)then	!proton
c	   if(i.eq.imark)write(6,*)'Bloody murder!! b',lundi,imark
c
c          Velocity of lab frame relative to CM
c          Boost the proton into the total system c.m. frame
c          Compute angle between c.m. proton and the beam axis
c       
c 	    e0 = etot(plist(1,1)	        !defines the beam plus target
	   e0 = plist(5,1)	        !defines the beam plus target
	   ww = plist(4,1)              !W of system 
	   if(e0.le.0)then
	      write(6,*)'How the heck did this happen 2 ?',
     1	                                           e0,plist
	      call exit()
c          e0 = 1.e32
	   end if
	   do j=1,3
	      beta(j) = plist(j,1)/e0
	   end do
	   call boost(plist(1,i),beta,ptemp)
	   costhpr    = dot(ptemp,nbeam)/sqrt(psq(ptemp))
	   costhmeson = -costhpr          !This is easy in the c.m. frame
	   thmeson    = acos(costhmeson)*rade
	   phimeson   = atan2(ptemp(2),ptemp(1)) * rade
	   if(phimeson.lt.0.0)phimeson = phimeson + 360. !azimuthal angle from 0 to 360.
	   phisector = mod(phimeson+30.,60.)-30. !phi within one sector +- 30 deg.
	   do j=1,3                       !direction of meson recoiling off proton in c.m. frame
	      dirmeson(j) = -ptemp(j)
	   enddo
	   dirmeson(4) = 0
c	   write(6,*)costhpr,phimeson,phisector,ptemp,nbeam
c       
c          Compute missing mass off proton, assuming we have a
c          two-body final state.  
c       
	   do k=1,4
	      ptemp(k) = plist(12+k,i)  !"decay time" (with energy loss) variables
c	      ptemp(k) = plist(k,i)
	   enddo
	   call missingmass(plist(1,1),ptemp,rmiss,rmiss2)
	   pmom = sqrt(plist(17,i)**2 - plist(16,i)**2) * 1000.
cc	   if(rmiss.lt.0.3)then
cc	      write(6,*)'>',rmiss,rmiss2
cc	      write(6,*)ptemp
cc	      write(6,*)(plist(k,i),k=1,4)
cc	      write(6,*)(plist(k,1),k=1,4)
cc	   endif
cc
	   call hf1(30,costhmeson,weight)
c
c          Histograms for pi_zero case
c
	   call hf1(518,rmiss,weight)
	   call hf1(519,rmiss2,weight)
	   call hf2(3226,ppbeam,costhmeson,weight) 
	   call hf2(228 ,thmeson,phimeson,weight) !mixed c.m. and lab angles
c
c          Histograms for eta case
c
	   call hf1(1518,rmiss,weight)
	   call hf1(1519,rmiss2,weight)
	   call hf2(1226,ppbeam,costhmeson,weight) 
c
c          Histograms for eta' and x(1280) cases
c       
	   call hf1(2518,rmiss,weight)
	   call hf1(2519,rmiss2,weight)
	   if(rmiss.gt.0.8 .and. rmiss.lt. 1.1)then !etaprime case
	      call hf2(4226,ww,costhmeson,weight) 
	      call hf1(621,ww,weight)
	      call hf1(623,ww,weight)
	      call hf1(627,costhmeson,weight)
c
	      if(costhmeson.ge.-1.0 .and. costhmeson.lt.-0.8)then
		 call hf2(5601,ww,pmompip,weight)
	      elseif(costhmeson.ge.-0.8 .and. costhmeson.lt.-0.6)then
		 call hf2(5602,ww,pmompip,weight)
	      elseif(costhmeson.ge.-0.6 .and. costhmeson.lt.-0.4)then
		 call hf2(5603,ww,pmompip,weight)
	      elseif(costhmeson.ge.-0.4 .and. costhmeson.lt.-0.2)then
		 call hf2(5604,ww,pmompip,weight)
	      elseif(costhmeson.ge.-0.2 .and. costhmeson.lt.-0.0)then
		 call hf2(5605,ww,pmompip,weight)
	      elseif(costhmeson.ge. 0.0 .and. costhmeson.lt. 0.2)then
		 call hf2(5606,ww,pmompip,weight)
	      elseif(costhmeson.ge. 0.2 .and. costhmeson.lt. 0.4)then
		 call hf2(5607,ww,pmompip,weight)
	      elseif(costhmeson.ge. 0.4 .and. costhmeson.lt. 0.6)then
		 call hf2(5608,ww,pmompip,weight)
	      elseif(costhmeson.ge. 0.6 .and. costhmeson.lt. 0.8)then
		 call hf2(5609,ww,pmompip,weight)
	      elseif(costhmeson.ge. 0.8 .and. costhmeson.lt. 1.0)then
		 call hf2(5610,ww,pmompip,weight)
	      endif
c
	   elseif(rmiss.gt.1.1.and.rmiss.lt.1.5)then !x1280 case
	      call hf2(4227,ww,costhmeson,weight) 
	      call hf1(622,ww,weight)
	      call hf1(624,ww,weight)
	      call hf1(628,costhmeson,weight)
	      call hf2(8120,costhmeson,ww,weight)
	   endif
c
c          TOF counter distribution of the protons
c       
	   call hf1(420+isectorp,float(numtofbarp),weight)
	   call hf2(430,float(numtofbarp),ppbeam,weight)
	   call hf2(431,costhmeson,float(numtofbarp),weight)
	   call hf2(432,pmom,    float(numtofbarp),weight)
c	   call hf2(5000+isectorp,float(numtofbarp),pmomp,weight)
	endif
 107	continue
c
c	For gamma + p --> K+ + Hyperon, pick out the K+ to define the
c	polarization direction.  Compute the center-of-mass scattering
c	angle of the kaon.  Compute boost parameters to get to Lambda
c       rest frame as well.
c
c
c	Invariant Mass of a detected proton and a detected pi-
c
	rmasspimp = 0
	do 140 i=2,imark
	   lundi = plist(6,i)
	   do 140 j=2,imark
	      lundj = plist(6,j)
	      if(lundi.eq.-17 .and. lundj.eq.41)then !got a pion and a proton
		 call epcm(plist(1,i),plist(1,j),ptemp)
		 rmasspimp = ptemp(4)
		 call hf1(514,rmasspimp,weight)
		 call hf2(515,phyperon(4),rmasspimp,weight)
	      endif
 140	continue
c
c	Invariant Mass of a detected proton and a detected pi+
c
        do 150 i=1,imark
	   lundi = plist(6,i)
	   do 150 j=1,imark
	      lundj = plist(6,j)
	      if(lundi.eq.17 .and. lundj.eq.41)then !got a pion and a proton
c          costh = dot(plist(1,i),plist(1,j))/
c     1  	sqrt(psq(plist(1,i))*psq(plist(1,j)))
c	   tht = acos(costh) * rade
c	   call hf1(??,tht,weight)
cc	   if(dpoverp_out.gt.0.0)then
cc	      call psmear(plist(1,i),ppione)
cc	      call psmear(plist(1,j),ppitwo)
cc	      call epcm(ppione,ppitwo,ptemp)
cc	   else
		 call epcm(plist(1,i),plist(1,j),ptemp)
cc	   end if
		 rmasspipp = ptemp(4)
		 call hf1(516,rmasspipp,weight)
	      endif
 150	continue
c
c
c	Opening angle between two Lambdas
c
c	do 200 i=1,imark-1
c	lundi = plist(6,i)
c	do 200 j=i+1,imark
c	lundj = plist(6,j)
c	if(lundi.eq.57 .and. lundj.eq.57)then	!got two lambdas
c		costh = dot(plist(1,i),plist(1,j))/
c     1  			sqrt(psq(plist(1,i))*psq(plist(1,j)))
c		call hf1(70,costh,weight)
c		goto 220
c	end if
c200	continue
c
c	Opening angle between two Protons
c
c220	do 240 i=1,imark-1
c	lundi = plist(6,i)
c	do 240 j=i+1,imark
c	lundj = plist(6,j)
c	if(lundi.eq.41 .and. lundj.eq.41)then	!got two protons
c		ipass = iaccept(plist(1,i),plist(1,j))
c		if(ipass.eq. .true.)then
c			costh = dot(plist(1,i),plist(1,j))/
c     1  			sqrt(psq(plist(1,i))*psq(plist(1,j)))
c			call hf1(71,costh,weight)
c		end if
c		goto 250
c	end if
c240	continue
c
c	Opening angle between two Pi-'s
c
c250	do 255 i=1,imark-1
c	lundi = plist(6,i)
c	do 255 j=i+1,imark
c	lundj = plist(6,j)
c	if(lundi.eq.-17 .and. lundj.eq.-17)then	!got two pi minuses
c		ekin1  = plist(5,i) - plist(4,i)
c		ekin2  = plist(5,j) - plist(4,j)
c		ipass1 = iaccept(plist(1,i),plist(1,j))
c		ipass2 = iseepi(ekin1)
c		ipass3 = iseepi(ekin2)
c		ipass  = ipass1.and.ipass2.and.ipass3
cc		write(6,*)ipass,ipass1,ipass2,ipass3,ekin1,ekin2
c		if(ipass.eqv..true.)then
c			costh = dot(plist(1,i),plist(1,j))/
c     1  			sqrt(psq(plist(1,i))*psq(plist(1,j)))
c			call hf1(72,costh,weight)
c		end if
c		goto 260
c	end if
c255	continue
cc
cc	Opening angle between two photons from one pi_zero
cc
c	ii = 0
c260	do 280 i=1,imark-1
c	ii = ii + 1
c	if(ii.ge.imark)goto 300
c	lundi = plist(6,ii)
c	j = ii + 1
c	lundj = plist(6,j)
c	if(lundi.eq.1 .and. lundj.eq.1)then	!got two photons
c		ipass = iaccept(plist(1,ii),plist(1,j))
c		if(ipass.eqv..true.)then
c			costh = dot(plist(1,ii),plist(1,j))/
c     1  			sqrt(psq(plist(1,ii))*psq(plist(1,j)))
c			call hf1(73,costh,weight)
c		ii = ii + 1
c		end if
c	end if
c280	continue
cc
cc	Opening angle between two photons from any pi_zeros
cc	This should include all the combinatorial pairs of photons
cc
c300	do 320 i=1,imark-1
c	lundi = plist(6,i)
c	do 320 j=i+1,imark
c	lundj = plist(6,j)
c	if(lundi.eq.1 .and. lundj.eq.1)then	!got two photons
c		ipass = iaccept(plist(1,i),plist(1,j))
c		if(ipass.eqv..true.)then
c			costh = dot(plist(1,i),plist(1,j))/
c     1  			sqrt(psq(plist(1,i))*psq(plist(1,j)))
c			call hf1(74,costh,weight)
c		end if
c	end if
c320	continue
c
c	Reconstruct Lambdas.
c	Compute proton angle in Lambda c.m. with respect to
c	polarization axis.
c       This section of code assumes that the K+ was already found.
c
c       Skip this code for the present...
c
	GOTO 371
 	do 370 i=2,imark
	   lundi = plist(6,i)
	   do 370 j=2,imark
	      lundj = plist(6,j)
	      if(lundi.eq.-17 .and. lundj.eq.41)then !got pi- and proton
		 call epcm(plist(1,i),plist(1,j),ptemp)	!4vec of ideal pair
c                4vec of straggled pi- and p
		 call epcms(plist(1,i),plist(1,j),plam2) 
		 if(kmark.gt.0)then
c                   4vec of straggled K+ and pi-
		    call epcms(plist(1,kmark),plist(1,i),pik) 
		 else
		    pik(4) = 0
		 endif
		 write(6,*)i,j
	 	 write(6,*)'Pion  :',(plist(k,i),k=1,5)
		 write(6,*)'Proton:',(plist(k,j),k=1,5)
		 write(6,*)'Lambda:',(ptemp(k),k=1,5)
c		 write(6,*)'Lam2  :',(plam2(k),k=1,5)
c		 write(6,*)'Lam2  :',(plam2(k),k=1,5)
c
		 rmasspip2 = plam2(4)**2.
		 rmasskp2  = pik(4)**2.
		 call hf2(84,rmasspip2,rmasskp2,weight)	!Dalitz plot
c
		 rmassinv = ptemp(4)
		 delmass = rmassinv - 1.115
c
c 		 Check to see that the Lambda kinematics is consistent.
c		 We may have gotten a bad pair, e.g. bad pion from kaon decay
c
		 if(abs(delmass).lt. 0.1)then !got a Lambda
		    call crossprod(nhatpol,ptemp,ptemp2,angle)
		    if(abs(angle-piover2).gt.0.01)then
c		       ibogus = ibogus + 1
c		       write(6,*)'Lambda not normal to pol. axis',ibogus
c		       write(6,*)'LAMBDA',ptemp
c		       write(6,*)'AXIS  ',nhatpol
c		       write(6,*)'ANGLE',angle
		       goto 360
c		     else
c		       write(6,*)'Good one:'
c		       write(6,*)'LAMBDA',ptemp
c		       write(6,*)'AXIS  ',nhatpol
c		       write(6,*)'ANGLE',angle
		    end if
c
c		    Define the z' and x' axes for Cz' and Cx'
c
		    zhatprimenorm = 1.0/absp(ptemp)
		    do k=1,3	
		       zhprime(k) = ptemp(k)*zhatprimenorm
		    end do
		    call crossprod(nhatpol,zhprime,xhprime,angle)
c
c                   Boost the proton to the Lambda center of mass frame
c                   Determine direction of proton in Lambda c.m. 
c		    with respect to the P polarization axis
c
c		    e0 = etot(ptemp)	
		    e0 = ptemp(5)
		    do k=1,3
		       beta(k) = ptemp(k)/e0
		    end do
		    call boost(plist(1,j),beta,ptemp2)
		    costhp = dot(ptemp2,nhatpol)/
     1  				sqrt(psq(ptemp2)*psq(nhatpol))
c		     call hf1(77,costhp,weight)
		    call hf1(78,costhk,weight)
		    call hf2(79,costhk,costhp,weight)
c
c                   Fill polarization acceptance bins as a 
c                   function of beam energy.  The bin was
c                   determined above, when the kaon was found.
c
		    call hf2(86,costhp,rkinbin,weight)
		    call hf2(860,costhk,costhp,weight) !same as #79
c
c                   Compute the proton direction w.r.t. the x' and
c                   z' axes, along which we expect no (single) polarization 
c                   but which can have Cx' and Cz' beam-recoil 
c                   polarizations.
c
		    costhpx = dot(ptemp2,xhprime)/
     1  		 sqrt(psq(ptemp2)*psq(xhprime))
c		    call hf1(80,costhpx,weight)
		    call hf2(82,costhk,costhpx,weight)
		    call hf2(992,costhpx,rkinbin,weight)
c
		    costhpz = dot(ptemp2,zhprime)/
     1  		 sqrt(psq(ptemp2)*psq(zhprime))
c		     call hf1(81,costhpz,weight)
		    call hf2(83,costhk,costhpz,weight)
		    call hf2(993,costhpz,rkinbin,weight)
c                    goto 400		!assume only one Lambda
c                     else
c		        write(6,*)'Failed Lambda mass test'
		 end if
 360		 continue
		 call hf1(503,rmassinv,weight) !as computed from ideal values
	      end if
 370	   continue
 371	   continue
c
c       Accumulate proton angular distributions as a 
c       function of an index which tells us which energy
c       and proton angle bin we are in.
c
c       In this section of code we assume that if we have a proton
c       then it came from Lambda decay, not some background reaction
c
c
c 400	continue
 	continue
	do 410 i=2,imark
	lundi = plist(6,i)
	if(lundi.eq.41)then	!got the proton (again)
c
c          Compute the decay proton's angle with respect to 
c          Lambda polarization direction.
c
	   call boost(plist(1,i),betahyp,pproton)
	   costhp = dot(pproton,nhatpol)/
     1  	sqrt(psq(pproton)*psq(nhatpol))
c	   write(6,*)'>loading h#85)',iloop,costhp,rkinbin
	   call hf2(85,costhp,rkinbin,weight)
  	   call hf2(850,costhk,costhp,weight)
cc	   if(ifiddle)then
c	      write(6,*)'RECONSTRUCTED: costhk,costhp',costhk,costhp
c	      write(6,*)'RECONSTRUCTED: nhatpol',nhatpol
c	      write(6,*)'RECONSTRUCTED: pproton',pproton
c	      write(6,*)'RECONSTRUCTED: betahyp',betahyp
cc	   endif
c
c   	   Compute the proton direction w.r.t. the x' and
c	   z' axes, along which we expect no (single) polarization 
c	   but which can have Cx' and Cz' beam-recoil 
c	   polarizations.
c
	   costhpx = dot(pproton,xhatprime)/
     1  	sqrt(psq(pproton)*psq(xhatprime))
	   call hf2(990,costhpx,rkinbin,weight)
  	   call hf2(9900,costhk,costhpx,weight)
c       
	   costhpz = dot(pproton,zhatprime)/
     1  	sqrt(psq(pproton)*psq(zhatprime))
	   call hf2(991,costhpz,rkinbin,weight)
c  	   call hf2(9910,costhk,costhpz,weight)  !careful, this histo is defined twice
  	   call hf2(9960,costhpx,costhpz,weight)
	   call hf1(77,costhp,weight)
	   call hf1(80,costhpx,weight)
	   call hf1(81,costhpz,weight)
c	   write(6,*)'Analyze:'
c	   write(6,*)'PLIST  :',i,plist(1,i),plist(2,i),plist(3,i),
c       1		                  plist(4,i)
c	   write(6,*)'NBEAM    :',nbeam
c	   write(6,*)'NHATPOL  :',nhatpol,rhatnorm
c	   write(6,*)'XHATPRIME:',xhatprime,xhatprimenorm
c	   write(6,*)'ZHATPRIME:',zhatprime,zhatprimenorm
c	   write(6,*)'PHYPERON :',phyperon
c	   write(6,*)'BETAHYPE :',betahyp
c	   write(6,*)'PPROTON  :',pproton
c	   write(6,*)'COSTHP   :',costhp
c	   write(6,*)'COSTHPX  :',costhpx
c	   write(6,*)'COSTHPZ  :',costhpz
c	   write(6,*)i,iebin,jkbin,rkinbin,ppbeam,costhk,costhp
	endif
 410	continue
c
c
c       Accumulate pi- angular distributions as a 
c       function of an index which tells us which energy
c       and proton angle bin we are in.
c
c       In this section of code we assume that if we have a proton
c       then it came from Lambda decay, not some background reaction
c
c
	continue
	do 420 i=2,imark
	lundi = plist(6,i)
	if(lundi.eq.-17)then	!got the pi- (again)
c
c          Compute the decay proton's angle with respect to 
c          Lambda polarization direction.  Reverse it relative to proton
c          case.
c
c          There is some contamination here from the decay of 
c          the K+ to pi+ pi+ pi-,  Try to reject this by looking
c          at the creation point of the pion: reject if more than 2 cm 
c          from the target (in X Y)
c
	   decaypoint = sqrt(plist(7,i)**2. + 
     1	          plist(8,i)**2.)     !   + plist(9,i)**2.)
	   if(decaypoint.gt.2.0)then
c	      write(6,*)'Got one'
	      goto 420
	   end if
	   call boost(plist(1,i),betahyp,ptemp2)
	   costhpi = -dot(ptemp2,nhatpol)/
     1  	sqrt(psq(ptemp2)*psq(nhatpol))
	   call hf2(87,costhpi,rkinbin,weight)
  	   call hf2(870,costhk,costhpi,weight)
c	   write(6,*)' '
c	   write(6,*)i,iebin,jkbin,rkinbin,ppbeam,costhk,costhp
c	   write(6,*)ptemp2,nhatpol,betahyp
	endif
 420	continue
c
c	Reconstruct Cascades via (K-,K+) missing mass
c	assuming a proton target, even if the real target was a nucleus.
c
c	do 480 i=2,imark
c	lundi = plist(6,i)
c	if(lundi.eq.+18)then	!got a K+
c	   beam(1) = 0
c	   beam(2) = 0
c	   beam(3) = ppbeam	!smeared value
c	   beam(4) = rmbeam
c	   targ(1) = 0
c	   targ(2) = 0
c	   targ(3) = 0
c	   targ(4) = 0.93828	!proton mass
c	   call epcm(beam,targ,psystem)
c	   call missingmass(psystem,plist(1,i),rmiss,rmiss2)
c
c	   call hf1(508,rmiss,weight)
c	   if(istop_cas)then
c	      call hf1(509,rmiss,weight)
c	   end if
c	end if
c 480	continue

c
c	Reconstruct pi- (or whatever) from Kaon and proton
c
 	do 490 i=2,imark
	   lundi = plist(6,i)
	   do 490 j=2,imark
	      lundj = plist(6,j)
	      if(lundi.eq.+18 .and. lundj.eq.41)then !got K+ and proton
		 call epcm(plist(1,i),plist(1,j),ptemp)
c                4vec of straggled K+ and p
		 call epcms(plist(1,i),plist(1,j),ptemp2) 
		 call missingmass(plist(1,1),ptemp2,rmiss,rmiss2)
c                write(6,*)i,j
c	 	 write(6,*)'Kaon  :',(plist(k,i),k=1,4)
c		 write(6,*)'Proton:',(plist(k,j),k=1,4)
c		 write(6,*)'Pi-   :',(ptemp(k),k=1,4),rmiss,rmiss2
		 call hf1(520,rmiss,weight)
		 call hf1(521,rmiss2,weight)
	      end if
 490    continue
c
c       Study semileptonic hyperon decay
c
c	Reconstruct neutrino from Kaon and proton and lepton
c       lepton:  -9 = mu-, +7 = electron
c
 	do i=2,imark
	   lundi = plist(6,i)  !Kaon
	   acceptweight_kaon = acceptlist(i)
	   do j=2,imark
	      lundj = plist(6,j)  !proton
	      acceptweight_proton = acceptlist(j)
	      do k=2,imark
		 lundk = plist(6,k)   !lepton
		 acceptweight_lepton = acceptlist(k)
		 if(lundi.eq.+18 .and. lundj.eq.41 .and.
c     1              lundk.eq.-9)then         !got all three of K+,proton,muon
     1              lundk.eq.7)then         !got all three of K+,proton,electron
		    if(iloop.eq.3 .and. itrigmask.eq.4)then
		       weight = acceptweight_proton
     1                        * acceptweight_kaon 
     1                        * acceptweight_lepton
c		       write(6,*)'>>>',weight,acceptweight_proton,
c     1                           acceptweight_kaon,acceptweight_lepton
		    endif
c
c                   Preliminary calculations: look at Delta TOF
c
		    path = 400.0 !cm (we lost the ability to get exact path length to TOF bars)
c		    call deltatof(plist(1,k),plist(13,k),path,rmmuon,dtofmuon)
		    pmuonlabs = absp(plist(13,k))
		    call hf2(6300,pmuonlabs,dtofmuon,weight)
c		    call deltatof(plist(1,k),plist(13,k),path,rmelectron,dtofelectron)
		    pelectronlabs = absp(plist(13,k))
		    call hf2(6301,pelectronlabs,dtofelectron,weight)
c
c
c                   Missing mass off Kaon alone
c
		    call missmass(plist(1,1),plist(1,i), plambda ,rmissl ,rmissl2)
		    call missmass(plist(1,1),plist(13,i),plambdas,rmissls,rmissls2) !smeared
c		    write(6,*)'Iloop         ',iloop
c		    write(6,*)'Lambda pure   ',plambda
c		    write(6,*)'Lambda smeared',plambdas
c
c                   Missing mass off Kaon, proton, and lepton  (i.e. should be the neutrino mass)
c
		    call epcm(plist(1,i),plist(1,j),ptemp) !combine Kaon and proton
		    call epcm(plist(1,k),ptemp,ptemp2)     !combined Kaon, proton, and lepton
		    call missmass(plist(1,1),ptemp2,pneutrino,rmissn,rmissn2)
		    pneutmom = absp(pneutrino) !horrible naming convention...scalar vs. 4-vector
c
		    call epcms(plist(1,i),plist(1,j),ptemp) !combine Kaon and proton !smeared
		    call epcm(plist(13,k),ptemp,ptemp2)     !combined Kaon, proton, and lepton (can't use epcms here: crappy code)
		    call missmass(plist(1,1),ptemp2,pneutrinos,rmissns,rmissns2)
		    pneutmoms = absp(pneutrinos)
c
		    call missmass(plist(1,1),ptemp,ptemp2,rmasstemp,rmasstemp2)  !should be the pion mass for lambda-->pi- p 
c
		    call hf1(532,rmissl,weight)
           	    call hf1(533,rmissls,weight)
		    call hf1(632,rmissn2,weight)
           	    call hf1(633,rmissns2,weight)
		    call hf1(634,pneutmom, weight)
		    call hf1(635,pneutmoms,weight)
		    call hf1(640,rmasstemp2,weight)   
		    call hf2(642,rmissns2,rmasstemp2,weight)
c
c                   Histograms for the lepton (electron or muon) in the lab frame
c
		    call getangles(plist(1,k),plepton,thetapart,cospart,phipart)
		    call hf1(103,plepton,weight)
		    call hf1(104,plepton,weight)
		    call hf2(204,plepton,thetapart,weight)
c
		    call hf2(636,rmissns2,pneutmoms,weight)
c
c                   Invariant mass of proton and candidate pions
c
		    call epcms(plist(1,j),plist(1,k),ptemp)  !combine proton and lepton
		    rimplepton = ptemp(4)  !mass value
		    call hf1(6302,rimplepton,weight)
c
c                   Temporarily give the lepton a pion mass
c
		    temp = plist(4,k)
		    plist(4,k) = rmpion
		    call epcms(plist(1,j),plist(1,k),ptemp)  !combine proton and lepton
		    rimppion = ptemp(4)  !mass value
		    plist(4,k) = temp  !undo the transgression
		    call hf1(6303,rimppion,weight)

c
c                   Boost particles to the overall CM frame (gamma + proton)
c
c                   If we boost the unsmeared variables, use plist(1,n), but if
c                   we boost the smeared variables, use plist(13,n).
c
c                   If we boost the smeared variables to the overall cm, then the later
c                   second boost to the hyperon rest frame will automatically also 
c                   be boosting smeared variables: no need to do things twice.
c
		    call getbeta(plist(1,1),beta)  !defines the beam plus target frame
c
		    call boost(plist(13,i),beta,pkaoncm)    !Kaon
		    call boost(plist(13,j),beta,pprotoncm)  !proton
		    call boost(plist(13,k),beta,pleptoncm)  !lepton
		    call boost(plambdas,   beta,plambdacm)  !Lambda
		    call boost(pneutrinos, beta,pneutrinocm)!missing neutrino
c		    call boost(plist(1,i),beta,pkaoncm)	!Kaon
c		    call boost(plist(1,j),beta,pprotoncm)  !proton
c		    call boost(plist(1,k),beta,pleptoncm)  !lepton
c		    call boost(plambda,   beta,plambdacm)  !Lambda
c		    call boost(pneutrino, beta,pneutrinocm)!missing neutrino
c
		    if(iloop.eq.0)then
		    pxcm = plambdacm(1)-pprotoncm(1)-pleptoncm(1)-pneutrinocm(1)
		    pycm = plambdacm(2)-pprotoncm(2)-pleptoncm(2)-pneutrinocm(2)
		    pzcm = plambdacm(3)-pprotoncm(3)-pleptoncm(3)-pneutrinocm(3)
		    pxlab = plambda(1)-plist(1,j)-plist(1,k)-pneutrino(1)
		    pylab = plambda(2)-plist(2,j)-plist(2,k)-pneutrino(2)
		    pzlab = plambda(3)-plist(3,j)-plist(3,k)-pneutrino(3)
		    write(6,*)'------ Next Semileptonic Event ------------------'
		    write(6,*)i,j,k
c		    write(6,*)'Overall CM frame'
		    write(6,*)'Overall CM frame (smeared variables)'
		    write(6,491)'Kaon     :',(plist(m,i),m=1,5),pkaoncm
		    write(6,491)'Proton   :',(plist(m,j),m=1,5),pprotoncm
		    write(6,491)'lepton   :',(plist(m,k),m=1,5),pleptoncm
		    write(6,491)'Lambda   :',plambda,plambdacm
		    write(6,491)'neutrino :',pneutrino,pneutrinocm
		    write(6,491)'Lab Lambda decay:',pxlab,pylab,pzlab
		    write(6,491)'CM  Lambda decay:',pxcm,pycm,pzcm
 491		    format(1x,a,5f12.7,'; ',5f12.7)
		    endif
c
c                   Boost proton, lepton, and neutrino to the hyperon rest frame
c                   from the overall cm frame
c
c                   If plambdacm is already smeared, then this next boost adds a 
c                   kind of smearing to all variables as well.
c
		    call getbeta(plambdacm,betahyp)   !defines the hyperon rest frame
		    call boost(plambdacm,  betahyp,plambdahyp)  !Lambda (sanity check)
		    call boost(pprotoncm,  betahyp,pprotonhyp)  !proton
		    call boost(pleptoncm,  betahyp,pleptonhyp)  !lepton
		    call boost(pneutrinocm,betahyp,pneutrinohyp)!missing neutrino

		    if(iloop.eq.0)then
		    pxhyp = pprotonhyp(1)+pleptonhyp(1)+pneutrinohyp(1)
		    pyhyp = pprotonhyp(2)+pleptonhyp(2)+pneutrinohyp(2)
		    pzhyp = pprotonhyp(3)+pleptonhyp(3)+pneutrinohyp(3)
		    write(6,*)' '
		    write(6,*)'Hyperon Rest frame reached from overall CM frame'
		    write(6,491)'Lambda   :',plambdahyp
		    write(6,491)'Proton   :',pprotonhyp
		    write(6,491)'lepton   :',pleptonhyp
		    write(6,491)'neutrino :',pneutrinohyp
		    write(6,491)'HYP Lambda decay:',pxhyp,pyhyp,pzhyp
		    endif
c
c                   Boost proton, lepton, and neutrino to the hyperon rest frame directly
c                   from the lab frame;  this is a second way to get to the same
c                   rest frame.  Do the Lorentz transformations show non-commutativity?
c                   Yes! Energy and momentum are conserved both ways, but there is 
c                   an overall rotation that changes the components (but not the magnitude)
c                   of the particles' momentum.
c
c                   If the boost parameter is computed using plambda, the boost is unsmeared,
c                   but if betahyp is computed from plambdas, the boost includes the 
c                   smearing distortion to all subsequently boosted particles.
c
c                   Distinguish this boost from the previous one with the suffix "2" on the
c                   variable names.
c
		    call getbeta(plambdas, betahyp2)  !defines the hyperon rest frame
		    call boost(plambdas,   betahyp2,plambdahyp2)  !Lambda (sanity check)
		    call boost(plist(13,j),betahyp2,pprotonhyp2)  !proton
		    call boost(plist(13,k),betahyp2,pleptonhyp2)  !lepton
		    call boost(pneutrinos ,betahyp2,pneutrinohyp2)!missing neutrino
c		    call getbeta(plambda, betahyp2)  !defines the hyperon rest frame
c		    call boost(plambda,   betahyp2,plambdahyp2)  !Lambda (sanity check)
c		    call boost(plist(1,j),betahyp2,pprotonhyp2)  !proton
c		    call boost(plist(1,k),betahyp2,pleptonhyp2)  !lepton
c		    call boost(pneutrino ,betahyp2,pneutrinohyp2)!missing neutrino
c
		    if(iloop.eq.0)then
		    pxhyp2 = pprotonhyp2(1)+pleptonhyp2(1)+pneutrinohyp2(1)
		    pyhyp2 = pprotonhyp2(2)+pleptonhyp2(2)+pneutrinohyp2(2)
		    pzhyp2 = pprotonhyp2(3)+pleptonhyp2(3)+pneutrinohyp2(3)
		    cosproton = dot(pprotonhyp,pprotonhyp2)/
     1                          (absp(pprotonhyp)*absp(pprotonhyp2))
		    thetaproton = acos(cosproton)*rade

		    write(6,*)' '
		    write(6,*)'Hyperon Rest frame reached directly from lab frame'
		    write(6,491)'Lambda   :',plambdahyp2
		    write(6,491)'Proton   :',pprotonhyp2
		    write(6,491)'lepton   :',pleptonhyp2
		    write(6,491)'neutrino :',pneutrinohyp2
		    write(6,491)'HYP Lambda decay:',pxhyp2,pyhyp2,pzhyp2
		    write(6,491)'Lorentz Rotation angle:',cosproton,thetaproton
		    endif
c
c                   Now we decide which boost to use for plotting variables in the hyperon 
c                   rest frame.  
c                   Choose to use the boost **from the overall cm frame** (no suffix "2")
c                   Choose to use the boost **directly from the lab frame** ( suffix "2")
c
		    call getangles(pprotonhyp,pprotonhypmag,dum1,dum2,dum3)
		    call getangles(pneutrinohyp,pneutmomhypmag,dum1,dum2,dum3)
		    call getangles(pleptonhyp,pleptonhypmag,thetapart,cospart,phipart)
		    costhprotonlepton = dot(pleptonhyp,pprotonhyp)
     1                          /sqrt(psq(pleptonhyp)*psq(pprotonhyp))

		    call crossprod(pkaoncm  ,pprotoncm,ptemp5,angle)
		    call crossprod(pleptoncm,pprotoncm,ptemp6,angle)
		    coplanarity = dot(ptemp5,ptemp6)/(absp(ptemp5)*absp(ptemp6))

		    call hf1(1039,pprotonhypmag,weight)
		    call hf1(1040,pleptonhypmag,weight)
		    call hf1(6340,pneutmomhypmag, weight)
		    call hf1(6342,costhprotonlepton,weight)
		    call hf2(2040,pleptonhypmag,costhprotonlepton,weight)

		    call hf2(2042,pleptonhypmag,coplanarity,weight)

		    costhleptonneutrino = dot(pleptonhyp,pneutrinohyp)
     1                          /sqrt(psq(pleptonhyp)*psq(pneutrinohyp))
		    elepton = pleptonhyp(5)
		    call hf2(2050,costhleptonneutrino,pleptonhypmag,weight)
		    call hf1(1042,costhleptonneutrino,weight)
		 end if
	      enddo
	   enddo
	enddo
 495	continue
c
c	write(6,*)'Got to here 1'
	GOTO 899  !jump closer to what we are presently interested in 
c       
c       The standard Lambda decay to pi- p
c	Reconstruct "zero" or "nothing"  from Kaon and proton and pion events
c
 	do i=2,imark
	   lundi = plist(6,i)  !Kaon
	   acceptweight_kaon = acceptlist(i)
	   do j=2,imark
	      lundj = plist(6,j)  !proton
	      acceptweight_proton = acceptlist(j)
	      do k=2,imark
		 lundk = plist(6,k)   !pi-
		 acceptweight_lepton = acceptlist(k)
		 if(lundi.eq.+18 .and. lundj.eq.41 .and.
     1              lundk.eq.-17)then         !got all of K+,proton,pion
		    if(iloop.eq.-17 .and. itrigmask.eq.4)then
		       weight = acceptweight_proton
     1                        * acceptweight_kaon 
     1                        * acceptweight_piminus
c		       write(6,*)'>>>',weight,acceptweight_proton,
c     1                           acceptweight_kaon,acceptweight_pion
		    endif
c
c                   Preliminary calculations: look at Delta TOF
c
		    path = 400.0 !cm (we lost the ability to get exact path length to TOF bars)
c		    call deltatof(plist(1,k),plist(13,k),path,rmmuon,dtofpion)
		    ppionlabs = absp(plist(13,k))
		    call hf2(6300,ppionlabs,dtofpion,weight)
c
c                   Missing mass off Kaon alone
c
		    call missmass(plist(1,1),plist(1,i),plambda,rmissl,rmissl2)
		    call missmass(plist(1,1),plist(13,i),plambdas,rmissls,rmissls2) !smeared
c
c                   Missing mass off Kaon, proton, and pion
c
		    call epcm(plist(1,i),plist(1,j),ptemp) !combine Kaon and proton
		    call epcm(plist(1,k),ptemp,ptemp2)     !combined Kaon, proton, and pion
		    call missmass(plist(1,1),ptemp2,pnothing,rmissnothing,rmissnothing2)
		    pnothingmag = absp(pnothing)
c
		    call epcms(plist(1,i),plist(1,j),ptemp) !combine Kaon and proton !smeared
		    call epcm(plist(13,k),ptemp,ptemp2)     !combined Kaon, proton, and pion (can't use epcms here: crappy code)
		    call missmass(plist(1,1),ptemp2,pnothings,
     1                                   rmissnothings,rmissnothings2)
		    pnothingmags = absp(pnothings)
c
		    call missmass(plist(1,1),ptemp,ptemp2,rmasstemp,rmasstemp2)  !should be the pion mass for lambda-->pi- p 
c		    
c                   write(6,*)i,j,k
c		    write(6,*)'Kaon     :',(plist(m,i),m=1,4)
c		    write(6,*)'Proton   :',(plist(m,j),m=1,4)
c		    write(6,*)'pi-      :',(plist(m,k),m=1,4)
c		    write(6,*)'Lambda   :',rmissl,rmissl2
c		    write(6,*)'nothing  :',rmissnothing,rmissnothing2

		    call hf1(637,pnothingmag, weight)
		    call hf1(638,pnothingmags,weight)
		    call hf1(641,rmasstemp2,weight)
		    call hf1(6320,rmissnothings,weight)
		    call hf1(6330,rmissnothings2,weight)
c
c                   Histograms for the nothing (recoiling against proton pion)
c
		    call getangles(plist(1,k),ppion,thetapart,cospart,phipart)
		    call hf1(103,ppion,weight)
		    call hf1(104,ppion,weight)
		    call hf2(204,ppion,theta_polar,weight)
c
		    call hf2(639,rmissnothings2,pnothingmags,weight)
		    call hf2(643,rmissnothings2,rmasstemp2,weight)
c
		    pion_momentum = absp(plist(1,k))
c		    write(6,*)pmomgev_pion,pion_momentum
c		    call hf1(114,pion_momentum,weight)
c
c                   Invariant mass of proton and candidate pions
c
		    call epcms(plist(1,j),plist(1,k),ptemp)  !combine proton and lepton
		    rimppion = ptemp(4)  !mass value
		    call hf1(6304,rimppion,weight)
c
c                   Temporarily give the pion a lepton mass
c
		    temp = plist(4,k)
		    plist(4,k) = rmmuon
		    call epcms(plist(1,j),plist(1,k),ptemp)  !combine proton and lepton
		    rimplepton = ptemp(4)  !mass value
		    plist(4,k) = temp      !undo the transgression
		    call hf1(6305,rimplepton,weight)
c       
c                   Boost particles to the overall CM frame (gamma + proton)
c
		    call getbeta(plist(1,1),beta)  !defines the beam plus target frame
c
		    call boost(plist(13,i),beta,pkaoncm)    !Kaon
		    call boost(plist(13,j),beta,pprotoncm)  !proton
		    call boost(plist(13,k),beta,ppiminuscm) !pi-
		    call boost(plambdas,   beta,plambdacm)  !Lambda
		    call boost(pnothings,  beta,pnothingcm) !Nothing
c		    call boost(plist(1,i),beta,pkaoncm)    !Kaon
c		    call boost(plist(1,j),beta,pprotoncm)  !proton
c		    call boost(plist(1,k),beta,ppiminuscm) !pi-
c		    call boost(plambda,   beta,plambdacm)  !Lambda
c        	    call boost(pnothing,  beta,pnothingcm) !Nothing
c
		    if(iloop.eq.0)then
		    pxcm = plambdacm(1)-pprotoncm(1)-ppiminuscm(1)
		    pycm = plambdacm(2)-pprotoncm(2)-ppiminuscm(2)
		    pzcm = plambdacm(3)-pprotoncm(3)-ppiminuscm(3)
		    pxlab = plambda(1)-plist(1,j)-plist(1,k)
		    pylab = plambda(2)-plist(2,j)-plist(2,k)
		    pzlab = plambda(3)-plist(3,j)-plist(3,k)
		    write(6,*)'------ Next Hadronic Decay Event ------------------'
		    write(6,*)i,j,k
		    write(6,*)'Overall CM frame'
		    write(6,491)'Kaon     :',(plist(m,i),m=1,5),pkaoncm
		    write(6,491)'Proton   :',(plist(m,j),m=1,5),pprotoncm
		    write(6,491)'pi-      :',(plist(m,k),m=1,5),ppiminuscm
		    write(6,491)'Lambda   :',plambda,plambdacm
		    write(6,491)'Lab Lambda decay:',pxlab,pylab,pzlab
		    write(6,491)'CM  Lambda decay:',pxcm,pycm,pzcm
		    endif
c
c                   Boost proton, and pion to the hyperon rest frame (HYP)
c                   from the overall cm frame
c
		    call getbeta(plambdacm,betahyp)  !defines the hyperon rest frame
		    call boost(plambdacm,  betahyp,plambdahyp)  !Lambda (sanity check)
		    call boost(pprotoncm,  betahyp,pprotonhyp)  !proton
		    call boost(ppiminuscm, betahyp,ppiminushyp) !pion
		    call boost(pnothingcm, betahyp,pnothinghyp) !nothing

		    if(iloop.eq.0)then
		    pxhyp = pprotonhyp(1)+ppiminushyp(1)
		    pyhyp = pprotonhyp(2)+ppiminushyp(2)
		    pzhyp = pprotonhyp(3)+ppiminushyp(3)
		    write(6,*)' '
		    write(6,*)'Hyperon Rest frame reached from overall CM frame'
		    write(6,491)'Lambda   :',plambdahyp
		    write(6,491)'Proton   :',pprotonhyp
		    write(6,491)'pi-      :',ppiminushyp
		    write(6,491)'HYP Lambda decay:',pxhyp,pyhyp,pzhyp
		    endif
c
c                   do the same boosts again for the smeared tracks
c
		    call getbeta(plambdacms,betahyps)  !defines the hyperon rest frame
		    call boost(plambdacm,  betahyps,plambdahyps)  !Lambda (sanity check)
		    call boost(pprotoncm,  betahyps,pprotonhyps)  !proton
		    call boost(ppiminuscm, betahyps,ppiminushyps) !pion
		    call boost(pnothingcm, betahyps,pnothinghyps) !nothing
c

c                   Boost proton and pion to the hyperon rest frame
c                   from the lab frame;  this is a second way to get to the same
c                   rest frame.  Do the Lorentz transformations show non-commutativity?
c                   Yes! Energy and momentum are conserved both ways, but there is 
c                   an overall rotation that changes the components (but not the magnitude)
c                   of the particles' momentum.
c
		    call getbeta(plambdas, betahyp)  !defines the hyperon rest frame
		    call boost(plambdas,   betahyp,plambdahyp2)  !Lambda (sanity check)
		    call boost(plist(13,j),betahyp,pprotonhyp2)  !proton
		    call boost(plist(13,k),betahyp,ppiminushyp2) !pion
		    call boost(pnothings,  betahyp,pnothinghyp2) !nothing
c		    call boost(plist(1,j),betahyp,pprotonhyp2)  !proton
c		    call boost(plist(1,k),betahyp,ppiminushyp2) !pion
c		    call boost(pnothing,  betahyp,pnothinghyp2) !nothing
c
		    if(iloop.eq.0)then
		    pxhyp2 = pprotonhyp2(1)+ppiminushyp2(1)
		    pyhyp2 = pprotonhyp2(2)+ppiminushyp2(2)
		    pzhyp2 = pprotonhyp2(3)+ppiminushyp2(3)
		    cosproton = dot(pprotonhyp,pprotonhyp2)/
     1                          (absp(pprotonhyp)*absp(pprotonhyp2))
		    thetaproton = acos(cosproton)*rade
		    write(6,*)' '
		    write(6,*)'Hyperon Rest frame reached directly from lab frame'
		    write(6,491)'Lambda   :',plambdahyp2
		    write(6,491)'Proton   :',pprotonhyp2
		    write(6,491)'pi-      :',ppiminushyp2
		    write(6,491)'HYP Lambda decay:',pxhyp2,pyhyp2,pzhyp2
		    write(6,491)'Lorentz Rotation angle:',cosproton,thetaproton
		    endif
c
c                   Now we decide which boost to use for plotting variables in the hyperon 
c                   rest frame.  
c                   Choose to use the boost **from the overall cm frame** (no suffix)
c                   Choose to use the boost **directly from the lab frame** (suffix "2")
c
c

		    call getangles(pprotonhyp,pprotonhypmag,dum1,dum2,dum3)
		    call getangles(pnothinghyp,pnothinghypmag,dum1,dum2,dum3)
		    call getangles(ppiminushyp,ppiminushypmag,thetapart,cospart,dum)
		    costhprotonpion = dot(ppiminushyp,pprotonhyp)
     1                          /sqrt(psq(ppiminushyp)*psq(pprotonhyp))


		    call crossprod(pkaoncm   ,pprotoncm,ptemp5,angle)
		    call crossprod(ppiminuscm,pprotoncm,ptemp6,angle)
		    coplanarity = dot(ptemp5,ptemp6)/(absp(ptemp5)*absp(ptemp6))

		    call hf1(6343,costhprotonpion,weight)
		    call hf1(1038,pprotonhypmag,weight)
		    call hf1(1041,ppiminushypmag,weight)
		    call hf1(6341,pnothinghypmag,weight)
		    call hf2(2041,ppiminushypmag,costhprotonpion,weight)
		    call hf2(2043,ppiminushypmag,coplanarity,weight)

		 endif
	      enddo
	   enddo
	enddo
 496	continue
c
c	Reconstruct K+ p pi- events ("Case 2")
c
 	do 500 i=2,imark
	   lundi = plist(6,i)         !Kaon
	   acceptweight_kaon = acceptlist(i)
	   do 500 j=2,imark
	      lundj = plist(6,j)      !proton
	      acceptweight_proton = acceptlist(j)
	      do 500 k=2,imark
		 lundk = plist(6,k)   !pi-
c		 write(6,*)'>>>',imark,i,j,k,lundi,lundj,lundk
		 acceptweight_piminus = acceptlist(k)
		 if(lundi.eq.+18 .and. lundj.eq.41 
     1           .and. lundk.eq.-17)then !got K+ and proton and pi-
c                   write(6,*)i,j,k
c                   write(6,*)'Kaon  :',(plist(m,i),m=1,4)
c                   write(6,*)'Proton:',(plist(m,j),m=1,4)
c                   write(6,*)'Pi-   :',(plist(m,k),m=1,4),rmiss,rmiss2
		    if(iloop.eq.3 .and. itrigmask.eq.9)then
		       weight = acceptweight_proton
     1                        * acceptweight_kaon 
     1                        * acceptweight_piminus
c		       write(6,*)'>>>',weight,acceptweight_proton,
c     1                           acceptweight_kaon,acceptweight_piminus
		    endif

		    call missingmass(plist(1,1),plist(1,i),rmissy,rmiss2)!Y
c		    call hf1(532,rmissy,weight) !MM off K+ alone
		    call missingmass(plist(1,1),plist(13,i),rmissysmear,rmiss2)!Y
		    call hf1(533,rmissysmear,weight) !MM off K+ alone but smeared
c
		    call epcm(plist(1,i),plist(1,j),ptemp)
		    call epcm(ptemp,     plist(1,k),ptemp2)
		    call missingmass(plist(1,1),ptemp2,rmiss,rmiss2)  !pi0
		    call hf1(530,rmiss,weight) !MM off K+ , p and pi- --> pi0
		    call epcm(plist(1,i),plist(13,j),ptemp)
		    call epcm(ptemp,     plist(13,k),ptemp2)
		    call missingmass(plist(1,1),ptemp2,rmiss,rmiss2)  !pi0 smeared
		    call hf1(531,rmiss,weight) !MM off K+ , p and pi- --> pi0 or K0
		    if(0.460.le.rmiss .and. 0.560.ge.rmiss)then  !select K0
		       wvalue = plist(4,1)
		       call hf1(575,costhmeson,weight)
		       call hf2(592,wvalue,costhmeson,weight) !acceptance...
		    endif
c
		    call epcm(plist(1,i),plist(1,j),ptemp) !K+ and proton
		    call missingmass(plist(1,1),ptemp,rmisspim,rmiss2)
		    call hf1(536,rmisspim,weight) !MM off K+ and proton
		    call epcm(plist(13,i),plist(13,j),ptemp) !K+ and proton
		    call missingmass(plist(1,1),ptemp,rmisspim,rmiss2)
		    call hf1(537,rmisspim,weight) !MM off K+ and proton
c
		    call epcm(plist(1,k),plist(1,j),plam)  !pi- and proton 
		    rmlam = plam(4)
		    call hf1(540,rmlam,weight) ! Lambda
		    call epcm(plist(13,k),plist(13,j),plam)  !pi- and proton smeared 
		    rmlamsmear = plam(4)
		    call hf1(541,rmlamsmear,weight) !Lambda smeared
		    call hf2(550,rmlamsmear,rmiss,weight)
c       
c                   Select events with pi0 and a gamma to pick up Case 2D
c
		    wvalue = plist(4,1)
		    if(rmiss.gt.0.180                         .and. 
     1                (rmlam .ge. 1.105 .and. rmlam.le.1.125))then !  .and.
c     1                (rmissysmear.ge. 1.250 .and. rmissysmear.le.1.60))then
		       call hf2(840,costhk,wvalue,weight)
		       call hf2(841,rmissysmear,wvalue,weight)
		    endif
c
		    call epcm(plist(1,i),plist(1,k),ptemp)  !K+ and pi-
		    call missingmass(plist(1,1),ptemp,rmsigplus,rmiss2)  !Sigma+
		    call hf1(534,rmsigplus,weight)
		    rinvmasskpi = ptemp(4)
		    call hf1(538,rinvmasskpi,weight)
		    call epcm(plist(13,i),plist(13,k),ptemp)  !K+ and pi- smeared
		    call missingmass(plist(1,1),ptemp,rmsigplus,rmiss2)  !Sigma+
		    call hf1(535,rmsigplus,weight)
		    rmkpismear = ptemp(4)
		    call hf1(539,rmkpismear,weight)
		    call hf2(839,rmsigplus,rmlamsmear,weight)
c
c                   Plots of IM(K+ pi-) lineshape for bins of photon energy
c
		    if(ppbeam.gt.1.56 .and. ppbeam.le.1.77)then
		       call hf1(5381,rinvmasskpi,weight) 
		    elseif(ppbeam.gt.1.77 .and. ppbeam.le.1.99)then
		       call hf1(5382,rinvmasskpi,weight) 
		    elseif(ppbeam.gt.1.99 .and. ppbeam.le.2.23)then
		       call hf1(5383,rinvmasskpi,weight) 
		    elseif(ppbeam.gt.2.23 .and. ppbeam.le.2.47)then
		       call hf1(5384,rinvmasskpi,weight) 
		    elseif(ppbeam.gt.2.47 .and. ppbeam.le.2.73)then
		       call hf1(5385,rinvmasskpi,weight) 
		    elseif(ppbeam.gt.2.73 .and. ppbeam.le.3.00)then
		       call hf1(5386,rinvmasskpi,weight) 
		    elseif(ppbeam.gt.3.00 .and. ppbeam.le.3.27)then
		       call hf1(5387,rinvmasskpi,weight) 
		    elseif(ppbeam.gt.3.27 .and. ppbeam.le.3.56)then
		       call hf1(5388,rinvmasskpi,weight) 
		    elseif(ppbeam.gt.3.56 .and. ppbeam.le.3.83)then
		       call hf1(5389,rinvmasskpi,weight) 
		    endif

c
c                   Select Sigma+ events by selecting its mass, rejecting events
c                   containing Lambda candidates, and making a loose cut on 
c                   (gamma,K+) missing mass.  Then fill acceptance table
c                   Case 2B
c
c		    write(6,*)'>> ',iloop,weight,rmsigplus,rmlamsmear
		    if((rmsigplus.gt.1.180 .and. rmsigplus.lt.1.220) .and.
     1                   .not.
     1  		 (rmlamsmear.ge.1.105 .and.rmlamsmear.le.1.125))then !  .and.
c     1                (rmissysmear.ge. 1.250 .and. rmissysmear.le.1.60))then

c		       write(6,*)'>',iloop,weight,costhk,wvalue
		       call hf2(830,costhk,wvalue,weight)
		       call hf2(831,rmissysmear,wvalue,weight)
c		       write(6,*)'   ',iloop,weight
		    endif
		 
		    call hf2(700,rmissysmear,rmkpismear,weight) 

c		    call hf1(601,plist(7,1),weight) ! target vertex distribution
c		    call hf1(602,plist(8,1),weight) ! target vertex distribution
c		    call hf1(603,plist(9,1),weight) ! target vertex distribution
c
c                   Plots of (gamma,K) lineshape for bins of photon energy
c
		    if(ppbeam.gt.1.56 .and. ppbeam.le.1.77)then
		       call hf1(5321,rmissy,weight) !MM off K+ alone
		    elseif(ppbeam.gt.1.77 .and. ppbeam.le.1.99)then
		       call hf1(5322,rmissy,weight) !MM off K+ alone
		    elseif(ppbeam.gt.1.99 .and. ppbeam.le.2.23)then
		       call hf1(5323,rmissy,weight) !MM off K+ alone
		    elseif(ppbeam.gt.2.23 .and. ppbeam.le.2.47)then
		       call hf1(5324,rmissy,weight) !MM off K+ alone
		    elseif(ppbeam.gt.2.47 .and. ppbeam.le.2.73)then
		       call hf1(5325,rmissy,weight) !MM off K+ alone
		    elseif(ppbeam.gt.2.73 .and. ppbeam.le.3.00)then
		       call hf1(5326,rmissy,weight) !MM off K+ alone
		    elseif(ppbeam.gt.3.00 .and. ppbeam.le.3.27)then
		       call hf1(5327,rmissy,weight) !MM off K+ alone
		    elseif(ppbeam.gt.3.27 .and. ppbeam.le.3.56)then
		       call hf1(5328,rmissy,weight) !MM off K+ alone
		    elseif(ppbeam.gt.3.56 .and. ppbeam.le.3.83)then
		       call hf1(5329,rmissy,weight) !MM off K+ alone
		    endif
		 end if
 500	continue
c
c	Reconstruct p pi+ pi- events ("Case 1")
c
c       Items for x(1280) study
c
 	do 600 i=2,imark
	   lundi = plist(6,i)         !pi+
	   acceptweight_piplus = acceptlist(i)
	   do 600 j=2,imark
	      lundj = plist(6,j)      !proton
	      acceptweight_proton = acceptlist(j)
	      do 600 k=2,imark
		 lundk = plist(6,k)   !pi-
		 acceptweight_piminus = acceptlist(k)
		 if(lundi.eq.+17 .and. lundj.eq.41 
     1           .and. lundk.eq.-17)then !got pi+ and proton annd pi-
c
c                   Include a hack to adjust acceptance weight for the event because
c                   above we gave it only the weighting due to the proton
c                   but not yet the pi+ and pi-
c                   ILOOP=1  first time through: all events plotted
c                   ILOOP=2  weighting is just the fiducial weighting
c                   ILOOP=3  weighting includes the decay probability of the particles
c
c		    write(6,*)'>Case 1',iloop,itrigmask
		    if(iloop.eq.2 .and. itrigmask.eq.3)then !assume we are using trigger mask 3 (proton only)
		       weight = acceptweight_proton              
		       if(acceptweight_piplus .eq.0)weight = 0.  
		       if(acceptweight_piminus.eq.0)weight = 0.
		    endif
		    if(iloop.eq.3 .and. itrigmask.eq.3)then !assume we are using trigger mask 3 (proton only)
		       weight = acceptweight_proton
     1                        * acceptweight_piplus 
     1                        * acceptweight_piminus
c		       write(6,510)'>weights',i,j,k,weight,acceptweight_proton,
c     1                        acceptweight_piplus,acceptweight_piminus 
c 510		       format(1x,a,3i3,4f6.3)
		    endif
c
c		    write(6,*)'>> loop= ',iloop,weight,
c     1              acceptweight_proton,acceptweight_piplus,acceptweight_piminus
c
		    if(weight.eq.0.0)goto 600  !if this combo did not survive, don't bother plotting
c
c		    call missingmass(plist(1,1),plist(1,j),rmissp,rmiss2) !mesons off proton
		    call missmass(plist(1,1),plist(1,j),pxmeson,rmissp,rmiss2) !mesons off proton
c		    write(6,*)'MC_analyze plist(3,1)',plist(3,1),rmissp 
		    call hf1(570,rmissp,weight) !MM off proton alone
		    call missingmass(plist(1,1),plist(13,j),rmisspsmear,rmiss2) !smeared
		    call hf1(571,rmisspsmear,weight) !MM off proton alone but smeared
c		    call hf1(28,ppbeam,weight)
c
 		    call epcm(plist(1,i),plist(1,j),ptemp) 
 		    call epcm(ptemp,     plist(1,k),ptemp2)
c 		    call missingmass(plist(1,1),ptemp2,rmiss,rmiss2)  
 		    call missmass(plist(1,1),ptemp2,precval,rmiss,rmiss2)  
 		    call hf1(572,rmiss2,weight)        !MM off p pi+ pi-
		    call hf2(580,rmissp,rmiss2,weight) !MM gamma p pi+ pi-  vs MM gamma p
c
 		    call epcm(plist(13,i),plist(13,j),ptemp)
 		    call epcm(ptemp,     plist(13,k),ptemp2)
 		    call missmass(plist(1,1),ptemp2,precoil,rmiss,rmiss2)  
 		    call hf1(573,rmiss2,weight) !MM off p pi+ pi- smeared
		    call hf2(581,rmisspsmear,rmiss2,weight)     !smeared variable
		    call hf2(583,rmisspsmear,rmiss ,weight)     !smeared variable
		    pperp = sqrt(precoil(1)**2.+precoil(2)**2) !momentum perpendicular to beam
c
 		    call epcm(plist(1,i),plist(1,k),ptemp)  !pi+ pi- invariant mass 
		    rhomass = ptemp(4)
		    call hf1(595,rhomass,weight)            !look for rho's  
		    if(rmiss2.lt.0.1)then
		       call hf1(596,rhomass,weight)         !look for rho's  
		    endif
c
		    call epcm(plist(1,k),plist(1,j),ptemp) !invariant mass of pi- and proton
		    call missingmass(plist(1,1),ptemp,rmpipg,rmpipg2) !MM off p pi- or IM of pi+ gamma 
		    call epcm(plist(1,i),plist(1,j),ptemp) !invariant mass of pi+ and proton
		    call missingmass(plist(1,1),ptemp,rmpimg,rmpimg2) !MM off p pi+ or IM of pi- gamma 
c
		    px = plist(13,j)  !GeV/c
		    py = plist(14,j) 
		    pz = plist(15,j) 
		    pmompro = sqrt(px*px +py*py +pz*pz) !GeV/c
		    px = plist(13,i)  !GeV/c
		    py = plist(14,i) 
		    pz = plist(15,i) 
		    pmompip = sqrt(px*px +py*py +pz*pz) !GeV/c
		    px = plist(13,k) 
		    py = plist(14,k) 
		    pz = plist(15,k) 
		    pmompim = sqrt(px*px +py*py +pz*pz) !GeV/c 
		    call hf2(452,rmiss2,pmompip,weight)  
		    call hf2(453,rmiss2,pmompim,weight)  
c
c                   Rough division into eta and gamma recoils.
c                   This lets us compute acceptance for eta pi pi vs rho gamma
c
		    if(rmisspsmear.gt.0.92.and.rmisspsmear.lt.1.00)then  !select etaprimes
		       call hf1(593,rmiss2,weight)    !MM off p pi+ pi- smeared
		       if(0.620.lt.rhomass.and.0.920.gt.rhomass)then
c			  call hf1(597,rmiss2,weight) !MM2 off p pi+ pi- smeared
		       endif
c
c                      Rough reaction ID cuts to select eta' --> eta pi pi OR rho gamma
c
		       if(-0.04.le.rmiss2 .and. 0.04.ge.rmiss2)then !select gamma
			  call hf1(599,rhomass,weight) 
 			  call hf2(454,rmpipg2,rmpimg2,weight)
			  if(0.620.lt.rhomass.and.0.920.gt.rhomass)then !select rho
			     call hf1(716,pperp,weight)  !photon perp momentum
			     if(pperp.gt.0.100)then                     !reject low gamma pperp
				if(iloop.eq.3)then
c				   write(6,*)'Got a rho gamma event...',weight
				   if(weight.gt.0)then
c				      iwriteit = .true. !write this event to output file
c				      call prevent(plist,imark+1,ikeep)
				   endif
				endif 
c
				call hf2(455,rmpipg2,rmpimg2,weight)
				call hf1(626,ppbeam,weight)
				call hf2(4236,ww,costhmeson,weight)
				call hf1(705,pmompip*1000,weight)
				call hf1(709,pmompro*1000,weight)
				call hf1(714,pmompim*1000,weight)
				if(costhmeson.ge.-1.0 .and. costhmeson.lt.-0.8)then
				   call hf2(4601,ww,pmompip,weight)
				   call hf2(4611,ww,rmiss2,weight)
				elseif(costhmeson.ge.-0.8 .and. costhmeson.lt.-0.6)then
				   call hf2(4602,ww,pmompip,weight)
				   call hf2(4612,ww,rmiss2,weight)
				elseif(costhmeson.ge.-0.6 .and. costhmeson.lt.-0.4)then
				   call hf2(4603,ww,pmompip,weight)
				   call hf2(4613,ww,rmiss2,weight)
				elseif(costhmeson.ge.-0.4 .and. costhmeson.lt.-0.2)then
				   call hf2(4604,ww,pmompip,weight)
				   call hf2(4614,ww,rmiss2,weight)
				elseif(costhmeson.ge.-0.2 .and. costhmeson.lt.-0.0)then
				   call hf2(4605,ww,pmompip,weight)
				   call hf2(4615,ww,rmiss2,weight)
				elseif(costhmeson.ge. 0.0 .and. costhmeson.lt. 0.2)then
				   call hf2(4606,ww,pmompip,weight)
				   call hf2(4616,ww,rmiss2,weight)
				elseif(costhmeson.ge. 0.2 .and. costhmeson.lt. 0.4)then
				   call hf2(4607,ww,pmompip,weight)
				   call hf2(4617,ww,rmiss2,weight)
				elseif(costhmeson.ge. 0.4 .and. costhmeson.lt. 0.6)then
				   call hf2(4608,ww,pmompip,weight)
				   call hf2(4618,ww,rmiss2,weight)
				elseif(costhmeson.ge. 0.6 .and. costhmeson.lt. 0.8)then
				   call hf2(4609,ww,pmompip,weight)
				   call hf2(4619,ww,rmiss2,weight)
				elseif(costhmeson.ge. 0.8 .and. costhmeson.lt. 1.0)then
				   call hf2(4610,ww,pmompip,weight)
				   call hf2(4620,ww,rmiss2,weight)
				endif
			     endif
			  endif
			  
		       elseif(0.26.le.rmiss2.and.0.34.ge.rmiss2)then !select eta
c			  if(iloop.eq.3)then
c			     write(6,*)'Got an eta pi pi event...',weight
c			     if(weight.gt.0)then
c				iwriteit = .true. !write this event to output file
c				call prevent(plist,imark+1,ikeep)
c			     endif
c			  endif
			  call hf1(625,ppbeam,weight)
			  call hf2(4246,ww,costhmeson,weight) 
			  call hf1(715,pmompip*1000,weight)
			  call hf1(719,pmompro*1000,weight)
			  call hf1(724,pmompim*1000,weight)
			  if(costhmeson.ge.-1.0 .and. costhmeson.lt.-0.8)then
			     call hf2(4701,ww,pmompip,weight)
			     call hf2(4711,ww,rmiss2,weight)
			  elseif(costhmeson.ge.-0.8 .and. costhmeson.lt.-0.6)then
			     call hf2(4702,ww,pmompip,weight)
			     call hf2(4712,ww,rmiss2,weight)
			  elseif(costhmeson.ge.-0.6 .and. costhmeson.lt.-0.4)then
			     call hf2(4703,ww,pmompip,weight)
			     call hf2(4713,ww,rmiss2,weight)
			  elseif(costhmeson.ge.-0.4 .and. costhmeson.lt.-0.2)then
			     call hf2(4704,ww,pmompip,weight)
			     call hf2(4714,ww,rmiss2,weight)
			  elseif(costhmeson.ge.-0.2 .and. costhmeson.lt.-0.0)then
			     call hf2(4705,ww,pmompip,weight)
			     call hf2(4715,ww,rmiss2,weight)
			  elseif(costhmeson.ge. 0.0 .and. costhmeson.lt. 0.2)then
			     call hf2(4706,ww,pmompip,weight)
			     call hf2(4716,ww,rmiss2,weight)
			  elseif(costhmeson.ge. 0.2 .and. costhmeson.lt. 0.4)then
			     call hf2(4707,ww,pmompip,weight)
			     call hf2(4717,ww,rmiss2,weight)
			  elseif(costhmeson.ge. 0.4 .and. costhmeson.lt. 0.6)then
			     call hf2(4708,ww,pmompip,weight)
			     call hf2(4718,ww,rmiss2,weight)
			  elseif(costhmeson.ge. 0.6 .and. costhmeson.lt. 0.8)then
			     call hf2(4709,ww,pmompip,weight)
			     call hf2(4719,ww,rmiss2,weight)
			  elseif(costhmeson.ge. 0.8 .and. costhmeson.lt. 1.0)then
			     call hf2(4710,ww,pmompip,weight)
			     call hf2(4720,ww,rmiss2,weight)
			  endif
		       endif

		    elseif(rmisspsmear.gt.1.1.and.rmisspsmear.lt.1.4)then !the x(1280)'s
		       call hf1(594,rmiss2,weight) !MM off p pi+ pi- smeared
		       if(0.620.lt.rhomass.and.0.920.gt.rhomass)then
c			  call hf1(598,rmiss2,weight) !MM off p pi+ pi- smeared
		       endif
		    endif
c
c                   Do another round of selections, this time starting with
c                   selection of the recoiling eta or recoiling photon.
c
c                   Rough reaction ID cuts to select x(1285)--> eta pi pi
c
		    if(1.225.le.rmisspsmear .and. 1.335.ge.rmisspsmear)then   !selects f1
c		       write(6,*)'f1       ',ww,costhmeson
		       call hf1(597,rmiss2,weight) ! MM2 off p pi+ pi- smeared
		       call hf1(5970,rmiss,weight) ! MM  off p pi+ pi- smeared
		       if(costhmeson.le.+0.9)then
			  if    (ww.gt.2.00 .and. ww.le.2.10)then
			     call hf2(10001,rmiss,costhmeson,weight)
			  elseif(ww.gt.2.10 .and. ww.le.2.20)then
			     call hf2(10002,rmiss,costhmeson,weight)
			  elseif(ww.gt.2.20 .and. ww.le.2.30)then
			     call hf2(10003,rmiss,costhmeson,weight)
			  elseif(ww.gt.2.30 .and. ww.le.2.40)then
			     call hf2(10004,rmiss,costhmeson,weight)
			  elseif(ww.gt.2.40 .and. ww.le.2.50)then
			     call hf2(10005,rmiss,costhmeson,weight)
			  elseif(ww.gt.2.50 .and. ww.le.2.60)then
			     call hf2(10006,rmiss,costhmeson,weight)
			  elseif(ww.gt.2.60 .and. ww.le.2.70)then
			     call hf2(10007,rmiss,costhmeson,weight)
			  elseif(ww.gt.2.70 .and. ww.le.2.80)then
			     call hf2(10008,rmiss,costhmeson,weight)
			  endif
		       endif 
		    endif
		    if(0.940.le.rmisspsmear .and. 0.975.ge.rmisspsmear)then   !selects eta'
c		       write(6,*)'etaprime ',ww,costhmeson
		       call hf1(598,rmiss2,weight) ! MM2 off p pi+ pi- smeared
		       call hf1(5980,rmiss,weight) ! MM  off p pi+ pi- smeared
		       if(costhmeson.le.+0.9)then
			  if    (ww.gt.2.00 .and. ww.le.2.10)then
			     call hf2(20001,rmiss,costhmeson,weight)
			  elseif(ww.gt.2.10 .and. ww.le.2.20)then
			     call hf2(20002,rmiss,costhmeson,weight)
			  elseif(ww.gt.2.20 .and. ww.le.2.30)then
			     call hf2(20003,rmiss,costhmeson,weight)
			  elseif(ww.gt.2.30 .and. ww.le.2.40)then
			     call hf2(20004,rmiss,costhmeson,weight)
			  elseif(ww.gt.2.40 .and. ww.le.2.50)then
			     call hf2(20005,rmiss,costhmeson,weight)
			  elseif(ww.gt.2.50 .and. ww.le.2.60)then
			     call hf2(20006,rmiss,costhmeson,weight)
			  elseif(ww.gt.2.60 .and. ww.le.2.70)then
			     call hf2(20007,rmiss,costhmeson,weight)
			  elseif(ww.gt.2.70 .and. ww.le.2.80)then
			     call hf2(20008,rmiss,costhmeson,weight)
			  endif
		       endif 
		    endif
		    if(0.28.le.rmiss2 .and. 0.32.ge.rmiss2)then  !select eta range first
		       call hf1(5760,rmissp,weight)     ! MM off proton alone
		       call hf1(576,rmisspsmear,weight) ! MM off proton alone but smeared
		       call hf2(780,rmisspsmear,costhmeson,weight)  !all W
		       if(ww.ge.1.8 .and. ww.lt.1.9)
     1                    call hf2(770,rmisspsmear,costhmeson,weight) 
		       if(ww.ge.1.9 .and. ww.lt.2.0)
     1                    call hf2(771,rmisspsmear,costhmeson,weight) 
		       if(ww.ge.2.0 .and. ww.lt.2.1)
     1                    call hf2(772,rmisspsmear,costhmeson,weight) 
		       if(ww.ge.2.1 .and. ww.lt.2.2)
     1                    call hf2(773,rmisspsmear,costhmeson,weight) 
		       if(ww.ge.2.2 .and. ww.lt.2.3)
     1                    call hf2(774,rmisspsmear,costhmeson,weight) 
		       if(ww.ge.2.3 .and. ww.lt.2.4)
     1                    call hf2(775,rmisspsmear,costhmeson,weight) 
		       if(ww.ge.2.4 .and. ww.lt.2.5)
     1                    call hf2(776,rmisspsmear,costhmeson,weight) 
		       if(ww.ge.2.5 .and. ww.lt.2.6)
     1                    call hf2(777,rmisspsmear,costhmeson,weight) 
		       if(ww.ge.2.6 .and. ww.lt.2.7)
     1                    call hf2(778,rmisspsmear,costhmeson,weight) 
		       if(ww.ge.2.7 .and. ww.lt.2.8)
     1                    call hf2(779,rmisspsmear,costhmeson,weight) 
c
		       if(1.15.le.rmissp .and. 1.35.ge.rmiss)then !select x(1285)
			  call hf1(574,costhmeson,weight)
			  call hf1(5740,costhmeson,weight)
			  call hf2(591,ww,costhmeson,weight) !acceptance...
c
			  call epcm(plist(1,i),plist(1,j),ptemp) !pi+ proton
			  call missingmass(plist(1,1),ptemp,rmissa0m,rmissa0m2) !a0- recoil  
			  call epcm(plist(1,k),plist(1,j),ptemp) !pi- proton
			  call missingmass(plist(1,1),ptemp,rmissa0p,rmissa0p2) !a0+ recoil  
			  call hf2(450,rmissa0p2,rmissa0m2,weight)
			  call epcm(plist(13,i),plist(13,j),ptemp) !pi+ proton
			  call missingmass(plist(1,1),ptemp,rmissa0m,rmissa0m2) !a0- recoil  
			  call epcm(plist (13,k),plist(13,j),ptemp) !pi- proton
			  call missingmass(plist(1,1),ptemp,rmissa0p,rmissa0p2) !a0+ recoil  
			  call hf2(451,rmissa0p2,rmissa0m2,weight)
c    
c                         Now that we have x(1280)-->eta pi pi, do some kinematic tests
c                         Boost directly from the lab frame to the x rest frame
c
			  exmeson = etot(pxmeson) !energy of x in the lab frame
			  do m=1,3
			     beta(m) = pxmeson(m)/exmeson
			  end do
c
c                         Boost unsmeared variables
c
			  call boost(plist(1,i),beta,pxcmpip)  !4-vector of pi+ in x meson rest frame
			  call boost(plist(1,k),beta,pxcmpim)  !4-vector of pi- in x meson rest frame
			  call boost(precval,   beta,pxcmeta)  !4-vector of eta in x meson rest frame
			  call boost(pbeamlab,  beta,pxcmbeam) !4-vector of incoming photon in x meson rest frame
c			  if(iloop.eq.1)write(6,*)pbeamlab,beta,pxcmbeam
c
c                         Boost smeared variables
c
c			  call boost(plist(13,i),beta,pxcmpip)  !4-vector of pi+ in x meson rest frame
c			  call boost(plist(13,k),beta,pxcmpim)  !4-vector of pi- in x meson rest frame
c			  call boost(precoil,    beta,pxcmeta)  !4-vector of eta in x meson rest frame
c			  call boost(pbeamlab,   beta,pxcmbeam) !4-vector of incoming photon in x meson rest frame
c
c                         Check that momentum is conserved among pi+ pi- and eta in the x(1280) frame
c
c			  pxsum = pxcmpip(1) + pxcmpim(1) + pxcmeta(1)
c			  pysum = pxcmpip(2) + pxcmpim(2) + pxcmeta(2)
c			  pzsum = pxcmpip(3) + pxcmpim(3) + pxcmeta(3)
c			  psummag=sqrt(pxsum**2.+pysum**2.+pzsum**2)
c			  write(6,*)'>  ',psummag,pxsum,pysum,pzsum
c
c                         The angle between the pions in the x(1280) rest frame
c
			  costhpipi = dot(pxcmpip,pxcmpim)/sqrt(psq(pxcmpip)*psq(pxcmpim))
c
c                   Relative momentum of the pi+ pi- in the x rest frame
c
			  pxdif = pxcmpip(1) - pxcmpim(1)
			  pydif = pxcmpip(2) - pxcmpim(2)
			  pzdif = pxcmpip(3) - pxcmpim(3)
			  pdifmag = sqrt(pxdif**2.+pydif**2.+pzdif**2.)
c			  write(6,*)costhpipi,pdifmag
			  call hf2(9000,costhpipi,pdifmag,weight) ! in the x rest frame
c
c                         Compute the relative momentum of the pions in their pair rest frame
c
			  call epcm(pxcmpip,pxcmpim,ppipi)
			  epipi = etot(ppipi) 
			  do m=1,3
			     beta(m) = ppipi(m)/epipi
			  end do
			  call boost(pxcmpip,beta,ppipipip)  !4-vector of pi+ in pipi rest frame
			  call boost(pxcmpim,beta,ppipipim)  !4-vector of pi- in pipi rest frame
c			  write(6,*)'>>+ ',pxcmpip
c			  write(6,*)'>>- ',pxcmpim
c			  write(6,*)'>>T ',ppipi,epipi
c
c                         Check that momentum totals to 0 for pi+ pi- in the pi pi rest frame
c
			  pxsum = ppipipip(1) + ppipipim(1)
			  pysum = ppipipip(2) + ppipipim(2)
			  pzsum = ppipipip(3) + ppipipim(3)
			  psummag=sqrt(pxsum**2.+pysum**2.+pzsum**2)
c			  write(6,*)'>  ',psummag,pxsum,pysum,pzsum
c
c                         Relative momentum of the pi+ pi- in the pipi rest frame
c
			  pxdif = ppipipip(1) - ppipipim(1)
			  pxdif = ppipipip(2) - ppipipim(2)
			  pxdif = ppipipip(3) - ppipipim(3)
			  pdifmag = sqrt(pxdif**2.+pydif**2.+pzdif**2.)
			  call hf2(9001,costhpipi,pdifmag,weight) ! in the pipi rest frame
c
c                         Rebuild the eta pi pi Dalitz plot in the x rest frame
c
			  call epcm(pxcmpip,pxcmeta,pa0p) !pi+ eta
			  call epcm(pxcmpim,pxcmeta,pa0m) !pi- eta
			  rmassa0p2 = pa0p(4)**2.
			  rmassa0m2 = pa0m(4)**2.
			  call hf2(461,rmassa0p2,rmassa0m2,weight)
c
c			  Angle between the pion and the x 'direction' in the x rest frame
c
 			  costhpipx = dot(pxcmpip,dirmeson)/sqrt(psq(pxcmpip)*psq(dirmeson))
 			  costhpimx = dot(pxcmpim,dirmeson)/sqrt(psq(pxcmpim)*psq(dirmeson))
c
c                         Angle between the pion and the normal to the reaction plane
c                         (Leads to same result as using the other axis defined above
c
c			  call crossprod(nbeam,plist(1,j),nhatpol,angle)
c			  costhpipx = dot(pxcmpip,nhatpol )/sqrt(psq(pxcmpip)*psq(nhatpol ))
c			  costhpimx = dot(pxcmpim,nhatpol )/sqrt(psq(pxcmpim)*psq(nhatpol ))
c			  
c			  write(6,*)'dirmeson',dirmeson
c			  write(6,*)'pxcmpip ',pxcmpip
c			  write(6,*)'pxcmpim ',pxcmpim
c			  write(6,*)'costhpipx,costhpimx',costhpipx,costhpimx
c
			  cosweightpip = costhpipx   !weight from P-wave and S-wave pieces
			  cosweightpim = costhpimx   !
			  sinweightpip = sqrt(1.- costhpipx**2.) 
			  sinweightpim = sqrt(1.- costhpimx**2.) 
c			  write(6,*)sinweightpip,sinweightpim,cosweightpip,cosweightpim
c
c                         Compute the Breit-Wigner amplitudes for the process
c                         Do the real and imaginary parts for each eta pi combination and add...
c
			  denominp = (rmanot**2.-rmassa0p2)**2. + (rmanot*gammaa0)**2.
			  amprealp = rmanot*gammaa0*(rmanot**2.-rmassa0p2)/denominp 
			  ampimagp = (rmanot*gammaa0)**2./denominp
			  denominm = (rmanot**2.-rmassa0m2)**2. + (rmanot*gammaa0)**2.
			  amprealm = rmanot*gammaa0*(rmanot**2.-rmassa0m2)/denominm
			  ampimagm = (rmanot*gammaa0)**2./denominm
c
c                         Fold in the angular weighting in the decay.
c                         The a0+ Breit-Wigner amplitude gets weighted by the angle of the pi-
c                         with respect to the given direction axis, and vice-versa.
c                         That is, the eta pi+ amplitude is weighted by sin or cos of the pi-
c                         decay direction with respect to the eta pi+ (a.k.a a0+)
c
			  amprealp = amprealp   * cosweightpim  
			  ampimagp = ampimagp   * cosweightpim
			  amprealm = amprealm   * cosweightpip
			  ampimagm = ampimagm   * cosweightpip

			  amptotmagp = amprealp**2 + ampimagp**2.
			  amptotmagm = amprealm**2 + ampimagm**2.
c			  amptotmag= amptotmagp + amptotmagm           !incoherent version
			  amptotmag= (amprealp+amprealm)**2 + (ampimagp+ampimagm)**2. !coherent sum
			  decayweight = amptotmag
c
c                         Back to basics:  just let the decay weight be unity
c
			  decayweight = 1.0
c
			  call hf2(462,rmassa0p2,rmassa0m2,weight*decayweight)
c
c                         Create the 'folded over' version of the Dalitz plot
c                         Take a horizontal slice of the distribution to project
c                         onto the X axis.  The slice boundaries are at 
c                         0.86 GeV^2
c                         1.17 GeV^2
c                         which are adjusted to match the slice taken of the real data
c                         in routine ~/g11/f1eta/etapipi_dalitz_process.f
c
c                         At present the sliced version is histogram 459
c                         while the unsliced projection is histogram 458.
c
			  if(rmassa0m2.gt.rmassa0p2)then 
			     call hf2(456,rmassa0p2,rmissa0m2,weight*decayweight)
			     call hf1(457,rmassa0m2,weight*decayweight)
			     call hf1(458,rmassa0p2,weight*decayweight)
			     if(rmassa0m2.gt.0.86 .and. rmassa0m2.lt.1.17)then
				call hf1(459,rmassa0p2,weight*decayweight)
			     endif
			  else
			     call hf2(456,rmassa0m2,rmsssa0p2,weight*decayweight)
			     call hf1(457,rmassa0p2,weight*decayweight)
			     call hf1(458,rmassa0m2,weight*decayweight)
			     if(rmissa0p2.gt.0.86 .and. rmissa0p2.lt.1.17)then
				call hf1(459,rmassa0m2,weight*decayweight)
			     endif
			  endif

			  if(.FALSE.)then !SKIP the step of writing events out one at a time
c
c                         Proceed to write out these x(1280)events to a file that we can
c                         subsequently use to weight the Monte Carlo phase space Dalitz distribution
c
c                         indexes of which bin this event will fall in the 50x50 Dalitz plot in eta pi+ and eta pi-
c                         mass of x-state that is decaying to eta pi pi
c                         mass of eta pi+ ^2
c                         mass of eta pi- ^2
c                         cosine of angle between a0+ and boost direction
c                         cosine of angle between a0- and boost direction
c                         eta 3-momentum
c                         pi+ 3-momentum
c                         pi- 3-momentum
c       
c                         Trap a handful of bad events for which the masses are way out of range
c
			  if(iloop.eq.1)then
			     indexa0p = int((rmassa0p2 - 0.40)/0.02 + 1) !between 1 and 50
			     indexa0m = int((rmassa0m2 - 0.40)/0.02 + 1) !between 1 and 50
c			     indexa0p = int((rmassa0p2 +.001 - 0.41)/0.02) + 1 !between 1 and 50
c			     indexa0m = int((rmassa0m2 +.001 - 0.41)/0.02) + 1
c			     write(6,*)rmassa0p2,indexa0p,indexa0px
c
c
c                            Angle of the a0 with respect to the x boost direction ("s-channel helicity" system)
c
			     cosxa0p  = -dot(beta,pxcmpim)/sqrt(psq(beta)*psq(pxcmpim)) !a0 with respect to boost direction
			     cosxa0m  = -dot(beta,pxcmpip)/sqrt(psq(beta)*psq(pxcmpip))
c
c                            Angle of the a0 with respect to the photon beam in the x rest frame  ("Gottfried Jackson" system)
c
c			     cosxa0p  = -dot(pxcmbeam,pxcmpim)/sqrt(psq(pxcmbeam)*psq(pxcmpim)) !a0 with respect to photon
c			     cosxa0m  = -dot(pxcmbeam,pxcmpip)/sqrt(psq(pxcmbeam)*psq(pxcmpip))
c
c                            Angle with respect to reaction plane normal
c
c			     call crossprod(nbeam,plist(1,j),nhatpol,angle)
c			     cosxa0p  = -dot(nhatpol,pxcmpim)/sqrt(psq(nhatpol)*psq(pxcmpim)) !a0 with respect to plane normal
c			     cosxa0m  = -dot(nhatpol,pxcmpip)/sqrt(psq(nhatpol)*psq(pxcmpip))
c
c
c                            Angle with respect to the beam in the c.m., ie. the beam axis z-hat ("Adair" system)
c
c			     cosxa0p  = -dot(nbeam,pxcmpim)/sqrt(psq(nbeam)*psq(pxcmpim)) !a0 with respect to plane normal
c			     cosxa0m  = -dot(nbeam,pxcmpip)/sqrt(psq(nbeam)*psq(pxcmpip))
c
c                            Here is where we produce the event by event output file for the times we want the Dalitz plot events one by one
c
			     if(indexa0p.lt.0 .or. indexa0p.gt.50 .or.
     1                          indexa0m.lt.0 .or. indexa0m.gt.50)then
				write(6,*)'How did we get here?',indexa0p,indexa0m
			     else
				write(4,550)float(indexa0p),float(indexa0m),pxmeson(4),
     1			            rmassa0p2,rmassa0m2,
     1			            cosxa0p,cosxa0m,
     1                              pxcmeta(1),pxcmeta(2),pxcmeta(3),
     1                              pxcmpip(1),pxcmpip(2),pxcmpip(3),
     1                              pxcmpim(1),pxcmpim(2),pxcmpim(3)
 550				format(1x,2f4.0,14f10.5)
			     endif
 			     ENDIF  !We have skipped the event-by-event output writing

			  endif
		       endif
		    endif
c
c                   Rough reaction ID cuts to select x(1285)--> rho gamma
c
		    if((0.1.gt.rmiss)  .and.                       !select gamma
     1		       (1.15.le.rmissp .and. 1.45.ge.rmissp))then  !select x(1285)
c
		       call epcm(plist(1,i),plist(1,k),ptemp) !pi+ pi-
		       rhomass = ptemp(4)
		       call hf1(582,rhomass,weight)
		       call hf1(5760,rmissp,weight) !MM off proton alone
		    endif
		 end if
 600	continue
c
c       BAIL OUT OF FURTHER ANALAYSIS OF THIS EVENT
c
c	GOTO 1000

c
c	Reconstruct K+ pi+ pi- events ("Case 3")
c
 	do 700 i=2,imark
	   lundi = plist(6,i)         !Kaon
	   acceptweight_kaon = acceptlist(i)
	   do 700 j=2,imark
	      lundj = plist(6,j)      !pi+
	      acceptweight_piplus = acceptlist(j)
	      do 700 k=2,imark
		 lundk = plist(6,k)   !pi-
		 acceptweight_piminus = acceptlist(k)
		 if(lundi.eq.+18 .and. lundj.eq.17 
     1           .and. lundk.eq.-17)then !got K+ and pi+ and pi-
c
c                   Include a hack to adjust acceptance weight for the event because
c                   above we gave it only the weighting due to the K+ and the pi-,
c                   but not yet the pi+
c
c		    write(6,*)'>>>',i,j,k,' itrigmask  iloop',itrigmask,iloop
c		    write(6,*)'Kaon  :',(plist(m,i),m=1,4),acceptweight_kaon
c		    write(6,*)'Pi+   :',(plist(m,j),m=1,4),acceptweight_piplus
c		    write(6,*)'Pi-   :',(plist(m,k),m=1,4),acceptweight_piminus
		    if(iloop.eq.3 .and. itrigmask.eq.9)then !assume we are using trigger mask 9
		       weight = acceptweight_piminus * 
     1		                acceptweight_piplus  * acceptweight_kaon
c		       write(6,*)'Weight:',weight
		    endif
c
c		    call hf1(28,ppbeam,weight)  !record beam energy
c
		    call missingmass(plist(1,1),plist(1,i),rmissy,rmiss2)!Y
c	  	     call hf1(532,rmissy,weight)      !MM off K+ alone NOT smeared
		    call missingmass(plist(1,1),plist(13,i),rmissysmear,rmiss2)!Y
c
c                   Hack to select "background" region
c
c		    if(rmissysmear.lt.1.6)goto 700
c
		    call hf1(533,rmissysmear,weight) !MM off K+ alone but smeared
c
		    call epcm(plist(1,i),plist(1,j),ptemp) !reconstruct recoil neutron or Lambda
c		    write(6,*)'1a) kaon',(plist(nn,i),nn=1,4)
c		    write(6,*)'1b) pi+ ',(plist(nn,j),nn=1,4)
c		    write(6,*)'1c) ptemp',ptemp
		    call epcm(ptemp,     plist(1,k),ptemp2)
c		    write(6,*)'2a) pi- ',(plist(nn,k),nn=1,4)
c		    write(6,*)'2b)ptemp2',ptemp2
		    call missmass(plist(1,1),ptemp2,pneut,rmiss,rmiss2) !MM off K+ pi- pi+  
		    call hf1(800,rmiss,weight) 
c		    write(6,*)'3a) beam ',(plist(nn,1),nn=1,4)
c		    write(6,*)'3b) rmiss',rmiss,rmiss2
c		    write(6,*)'3c) pneut',pneut
c
		    call epcms(plist(1,i),plist(1,j),ptemp)  !smeared version epcmS
c		    write(6,*)'3d) pneut',pneut
		    call epcm(ptemp,      plist(1,k),ptemp2)
c		    write(6,*)'3e) pneut',pneut
		    call missmass(plist(1,1),ptemp2,pneuts,rmiss,rmiss2) !MM off K+ pi- pi+  
c		    write(6,*)'3f) pneut',pneut
		    call hf1(801,rmiss,weight)
c		    write(6,*)'3g) pneut',pneut
c
c                   Select events with a missing neutron
c                   Then look for Sigma+- events and make Dalitz plots
c
		    if(rmiss.gt.0.9 .and. rmiss.lt.1.0)then
		       pneut(4) = rmneut !artificially correct the neutron mass
c		       write(6,*)'3h) pneut',pneut
		       call epcm(pneut,plist(1,j),psigp) ! Sigma+
		       rmsigplus = psigp(4)
c		       write(6,*)' '
c		       write(6,*)'pneut',pneut
c		       write(6,*)'pi+  ',(plist(nn,j),nn=1,4)
c		       write(6,*)'pi-  ',(plist(nn,k),nn=1,4)
c		       write(6,*)'Sig+ ',psigp,rmsigplus

		       call epcm(pneut,plist(1,k),psigm) ! Sigma-
		       rmsigminus = psigm(4)

		       call hf1(802,rmsigplus,weight)
		       call hf1(803,rmsigminus,weight)
		       call hf2(804,rmsigplus,rmsigminus,weight)
		       call hf2(807,rmsigplus**2.,rmsigminus**2.,weight)

		       call epcm(psigp,plist(1,k),ptemp)  !Sigma+ and pi-
		       rmsigppim = ptemp(4)
		       call epcm(psigm,plist(1,j),ptemp)  !Sigma- and pi+
		       rmsigmpip = ptemp(4)
		       call epcm(plist(1,i),plist(1,k),ptemp) !Kaon and pi-
		       rmkppim = ptemp(4)
		       call epcm(plist(1,i),plist(1,j),ptemp) !Kaon and pi+
		       rmkppip = ptemp(4)
		       call epcm(psigp,plist(1,i),ptemp)  !Sigma+ and K+
		       rmsigpkp = ptemp(4)
		       call epcm(psigp,plist(1,i),ptemp)  !Sigma+ and K+
		       rmsigpkp = ptemp(4)
		       call epcm(psigm,plist(1,i),ptemp)  !Sigma- and K+
		       rmsigmkp = ptemp(4)
		       call hf2(808,rmkppip,rmsigppim,weight)
		       call hf2(809,rmkppim,rmsigppip,weight)
		       call hf2(810,rmkppip**2,rmsigppim**2,weight)
		       call hf2(811,rmkppim**2,rmsigppip**2,weight)
		       call hf2(81200,rmsigpkp**2,rmsigppim**2,weight)
		       call hf2(81300,rmsigmkp**2,rmsigmpip**2,weight)
		    endif
c
c                   Select events with a missing Lambda;  assume the parent is the Lambda(1520)
c
		    if(rmiss.gt.1.08 .and. rmiss.lt.1.15)then
c		       pneuts(4) = rmasslambda  !artificially correct the Lambda mass (puts structure in plots)
		       call epcm(pneuts,plist(1,j),ptemp) ! Sigma(1385)+
c		       call epcm(pneut,plist(1,j),ptemp) ! Sigma(1385)+
		       rmsigplusstar = ptemp(4)
		       call epcm(pneuts,plist(1,k),ptemp) ! Sigma(1385)-   (smeared version)
c		       call epcm(pneut,plist(1,k),ptemp) ! Sigma(1385)-   (unsmeared version)
		       rmsigminusstar = ptemp(4)
		       call hf1(8020,rmsigplusstar,weight)
		       call hf1(8030,rmsigminusstar,weight)
c		       call hf2(8040,rmsigplusstar,rmsigminusstar,weight)
		       call hf2(8040,rmsigplusstar**2,rmsigminusstar**2,weight)
		       if(rmissysmear.gt.1.50 .and. rmissysmear.lt.1.54)then
c hack		       if(rmissysmear.gt.1.60 )then
			  call hf2(8070,rmsigplusstar**2.,rmsigminusstar**2.,weight)
		       endif
		       if(rmissysmear.gt.1.50 .and. rmissysmear.lt.1.520)then
			  call hf2(8071,rmsigplusstar**2.,rmsigminusstar**2.,weight)
		       elseif(rmissysmear.gt.1.52 .and. rmissysmear.lt.1.540)then
			  call hf2(8072,rmsigplusstar**2.,rmsigminusstar**2.,weight)
		       endif
		       call epcm(plist(1,i),plist(1,k),ptemp) ! K* from K+ pi-
		       rmkppim = ptemp(4)
		       call hf2(8073,rmissysmear,rmkppim,weight)
		       call hf2(8074,rmsigplusstar ,rmkppim,weight)
		       call hf2(8075,rmsigminusstar,rmkppim,weight)
		    endif
c       
		    call epcm(plist(1,i),plist(1,k),ptemp) !K+ and pi-
		    call missingmass(plist(1,1),ptemp,rmiss,rmiss2) !MM off K+ pi-  
		    call hf1(805,rmiss,weight) 
		    call epcms(plist(1,i),plist(1,k),ptemp) !K+ and pi- smeared
		    call missingmass(plist(1,1),ptemp,rmiss,rmiss2) !MM off K+ pi-  
		    call hf1(806,rmiss,weight) 
c
		    call missingmass(plist(1,1),plist(1,i),rmissy,rmiss2)!Y
		    call missingmass(plist(1,1),plist(13,i),rmissysmear,rmiss2)!Y
c
c                   Fill a bunch of "basic" momentum and angle histograms for comparison to real data
c                   We want the overall c.m quantities
c		    write(6,*)'> Got to here: plot basic c.m. quantities'
c
c		    call hf1(28,ppbeam,weight)
		    e0=etot(plist(1,1))	!defines the beam plus target
		    if(e0.le.0)then
		       write(6,*)'How the heck did this happen?',
     1	                                           e0,plist
		       call exit()
		    end if
		    do n=1,3
		       beta(n) = plist(n,1)/e0
		    end do
		    call boost(plist(1,i),beta,ptemp)           !kaon
		    pmomkaon = absp(ptemp)
		    costhk   = dot(ptemp,nbeam)/sqrt(psq(ptemp))
		    thkaon   = acos(costhk)*rade
		    call hf1(7601,pmomkaon,weight)
c		    call hf1(7611,thkaon,weight)
		    call hf1(7611,costhk,weight)
		    call boost(plist(1,j),beta,ptemp)           !pi+
		    pmompiplus    = absp(ptemp)
		    costhpiplus   = dot(ptemp,nbeam)/sqrt(psq(ptemp))
		    thpiplus      = acos(costhpiplus)*rade
		    call hf1(7602,pmompiplus,weight)
c		    call hf1(7612,thpiplus,weight)
		    call hf1(7612,costhpiplus,weight)
c
c		    if(iloop.eq.3) write(6,*)'>> ',iloop,thpiplus,weight
		    call boost(plist(1,k),beta,ptemp)           !pi-
		    pmompiminus    = absp(ptemp)
		    costhpiminus   = dot(ptemp,nbeam)/sqrt(psq(ptemp))
		    thpiminus      = acos(costhpiminus)*rade
		    call hf1(7603,pmompiminus,weight)
c		    call hf1(7613,thpiminus,weight)
		    call hf1(7613,costhpiminus,weight)

c	            write(6,*)costhpi,ptemp,nbeam

c
c                   Plots of (gamma,K) lineshape for bins of photon energy
c                   (caution: these same histograms are filled for Case 2 as well
c
c		    if(iloop.eq.1)write(6,*) 'Filling 532'
		    call hf1(532,rmissy,weight) !MM off K+ alone
		    if(ppbeam.gt.1.56 .and. ppbeam.le.1.77)then
		       call hf1(5321,rmissy,weight) !MM off K+ alone
		    elseif(ppbeam.gt.1.77 .and. ppbeam.le.1.99)then
		       call hf1(5322,rmissy,weight) !MM off K+ alone
		    elseif(ppbeam.gt.1.99 .and. ppbeam.le.2.23)then
		       call hf1(5323,rmissy,weight) !MM off K+ alone
		    elseif(ppbeam.gt.2.23 .and. ppbeam.le.2.47)then
		       call hf1(5324,rmissy,weight) !MM off K+ alone
		    elseif(ppbeam.gt.2.47 .and. ppbeam.le.2.73)then
		       call hf1(5325,rmissy,weight) !MM off K+ alone
		    elseif(ppbeam.gt.2.73 .and. ppbeam.le.3.00)then
		       call hf1(5326,rmissy,weight) !MM off K+ alone
		    elseif(ppbeam.gt.3.00 .and. ppbeam.le.3.27)then
		       call hf1(5327,rmissy,weight) !MM off K+ alone
		    elseif(ppbeam.gt.3.27 .and. ppbeam.le.3.56)then
		       call hf1(5328,rmissy,weight) !MM off K+ alone
		    elseif(ppbeam.gt.3.56 .and. ppbeam.le.3.83)then
		       call hf1(5329,rmissy,weight) !MM off K+ alone
		    endif
c
c                   Select Sigma+ events and fill acceptance table
c
		    wvalue = plist(4,1)
		    if(rmsigplus.gt.1.180 .and. rmsigplus.lt.1.960)then
c		       call hf2(810,costhk,wvalue,weight)
c		       call hf2(811,rmissy,wvalue,weight)
		    endif
c
c                   Select Sigma- events and fill acceptance table
c
		    if(rmsigminus.gt.1.190 .and. rmsigplus.lt.1.205)then
		       call hf2(820,costhk,wvalue,weight)
		       call hf2(821,rmissy,wvalue,weight)
		    endif

		 endif
 700	continue
c
c	Reconstruct K+ p K- events (for Lambda (1520))
c
 	do 800 i=2,imark
	   lundi = plist(6,i)         !Kaon
	   acceptweight_kaon = acceptlist(i)
	   do 800 j=2,imark
	      lundj = plist(6,j) !proton
	      acceptweight_proton = acceptlist(j)
	      do 800 k=2,imark
		 lundk = plist(6,k) !K-
		 acceptweight_kminus = acceptlist(k)
		 if(lundi.eq.+18 .and. lundj.eq.41 
     1           .and. lundk.eq.-18)then !got K+ and proton and K-
c                   write(6,*)i,j,k
c                   write(6,*)'K+    :',(plist(m,i),m=1,4)
c                   write(6,*)'Proton:',(plist(m,j),m=1,4)
c                   write(6,*)'K -   :',(plist(m,k),m=1,4)
c		    if(iloop.eq.3 .and. itrigmask.eq.9)then
c		       weight = acceptweight_proton
c     1                        * acceptweight_kaon 
c     1                        * acceptweight_kminus
c		       write(6,*)'>>>',weight,acceptweight_proton,
c     1                           acceptweight_kaon,acceptweight_kminus
c		    endif
c
		    call missingmass(plist(1,1),plist(1,i),rmissy,rmissy2) !Y
		    call hf1(1532,rmissy,weight) !MM off K+ alone
		    call missingmass(plist(1,1),plist(13,i),rmissysmear,rmissysmear2)	!Y
		    call hf1(1533,rmissysmear,weight) !MM off K+ alone but smeared
		    call hf1(1534,ww,weight) !W distribution for these events
c        
		    call epcm(plist(1,i),plist(1,j),ptemp)
		    call epcm(ptemp,     plist(1,k),ptemp2)
		    call missingmass(plist(1,1),ptemp2,rmiss,rmiss2) !"nothing"
		    call hf1(1530,rmiss2,weight) !MM off K+ , p and pi- --> "nothing"
		    call epcm(plist(1,i),plist(13,j),ptemp)
		    call epcm(ptemp,     plist(13,k),ptemp2)
		    call missingmass(plist(1,1),ptemp2,rmiss,rmiss2) !"nothing" smeared
		    call hf1(1531,rmiss2,weight) !MM off K+ , p and pi- --> "nothing"
c
		 end if
800	continue

c

c	
c	Reconstruct p pi- pi+ pi- events (for K0 Lambda studies on deuterium in g14 CLAS data)
c
 	do i=2,imark
	 lundi = plist(6,i)	!proton
	 if(lundi.eq.41)then
c	  write(6,*)'Proton    :',(plist(ix,i),  ix=1,4)
	  acceptweight_proton = acceptlist(i)
	  do j=2,imark
             lundj = plist(6,j)	!pi+
	     if(lundj.eq.+17)then
c		write(6,*)'pi+       :',(plist(ix,j),  ix=1,4)
		acceptweight_piplus = acceptlist(j)
		do k=2,imark
		   lundk = plist(6,k) !pi-
		   if(lundk.eq.-17)then
c		      write(6,*)'pi-       :',(plist(ix,k),  ix=1,4)
		      acceptweight_piminus = acceptlist(k)
		      do kk=k+1,imark
c		      do kk=2,imark
			 lundkk = plist(6,kk) !pi- number 2
			 if(lundkk.eq.-17)then
c			    write(6,*)'pi- 2     :',(plist(ix,kk), ix=1,4)
			    acceptweight_piminus2 = acceptlist(kk)
c
c			    write(6,*)i,j,k,kk
c			    write(6,*)'Proton    :',(plist(ix,i),  ix=1,5)
c			    write(6,*)'pi+       :',(plist(ix,j),  ix=1,5)
c			    write(6,*)'pi-       :',(plist(ix,k),  ix=1,5)
c			    write(6,*)'pi- 2     :',(plist(ix,kk), ix=1,5)
			    if(iloop.eq.3 .and. itrigmask.eq.15)then
			       weight = acceptweight_proton
     1		               * acceptweight_piplus
     2		               * acceptweight_piminus*acceptweight_piminus2
c			       write(6,*)'>>>',weight,acceptweight_proton,
c     1                          acceptweight_piplus,acceptweight_piminus,
c     1                          acceptweight_piminus2
			    else
			       weight = 1.0
			    endif
c			   
c                           Some raw momentum distributions
c
			    pinitial(1) = 0.0
			    pinitial(2) = 0.0
			    pinitial(3) = ppbeam 
			    pinitial(4) = sqrt(2.*ppbeam*rmneut + rmneut**2.)
			    pinitial(5) = ppbeam + rmneut

			    pmomgev = absp(plist(1,i))   !lab proton
			    pmompip = absp(plist(1,j))   !pi+
			    pmompim1 = absp(plist(1,k))  !pi- first, from K0
			    pmompim2 = absp(plist(1,kk)) !pi- second, from Lambda
			    call hf1(109,pmomgev,weight) !lab proton
c			    call hf1(105,pmompip,weight) !lab piplus
c			    call hf1(114,pmompim1,weight)!lab piminus 1
c			    call hf1(114,pmompim2,weight)!lab piminus 2
			    call hf1(6260,ppbeam,weight)  !lab photon
c
c                           Invariant mass combinations
c
			    call invarmass(plist(1,i),plist(1,j) ,ptemp5, rmtemp)  !proton pi+
  			    call invarmass(ptemp5,    plist(1,k) ,ptemp6, rmtemp)  !proton pi+ pi-
			    call invarmass(ptemp6,    plist(1,kk),ptemp7, rmtemp)  !proton pi+ pi- pi-
c			    call missmass(plist(1,1),ptemp7,pspect,rmissspec,rmissspec2) !Spectator proton
			    call missmass(ptotrest,ptemp7,pspect,rmissspec,rmissspec2) !Spectator proton

			    recoilp = absp(pspect)

			    call hf1(2551,recoilp,weight)   !MM off pi+ pi- pi- p  --> spectator momentum

			    if(recoilp .le. 0.200) then       	!GeV/c recoil momentum cut

			       call hf1(2550,rmissspec,weight) !MM off pi+ pi- pi- p  --> spectator
			       
			       call invarmass(plist(1,j),plist(1,k), pK01 , rmpippim1) !IM of pi+ and pi-1  (K0)
			       call invarmass(plist(1,j),plist(1,kk),pK02 , rmpippim2) !IM of pi+ and pi-2
			       call hf1(2531,rmpippim1,weight) !look for the K0
			       call hf1(2532,rmpippim2,weight) !look for the K0
			       call hf1(2533,rmpippim1,weight) !look for the K0 combinatorics
			       call hf1(2533,rmpippim2,weight) !look for the K0
			       call invarmass(plist(1,i),plist(1,k) ,plambda1, rmppim1) !IM of proton and pi-1
			       call invarmass(plist(1,i),plist(1,kk),plambda2, rmppim2) !IM of proton and pi-2
			       call hf1(2541,rmppim1,weight) !look for the Lambda
			       call hf1(2542,rmppim2,weight) !look for the Lambda
			       call hf1(2543,rmppim1,weight) !look for the Lambda combinatorics
			       call hf1(2543,rmppim2,weight) !look for the Lambda
c
			       pmomK0lab1 = absp(pK01)
			       pmomK0lab2 = absp(pK02)
			       pmomlambdalab1 = absp(plambda1)
			       pmomlambdalab2 = absp(plambda2)
			       call hf1(2560,pmomK0lab1,weight)
			       call hf1(2561,pmomK0lab2,weight)
			       call hf1(2562,pmomlambdalab1,weight)
			       call hf1(2563,pmomlambdalab2,weight)
c
c                           Boost K0 and Lambda to the CM defined by beam and "stationary" neutron
c
			       e0gn=pinitial(5)	!defines the beam plus target
c                            e0gn = 1.e32
			       do nn=1,3
				  beta(nn) = pinitial(nn)/e0gn !beta factor for transformation to cm
			       end do
			       call boost(pK01,    beta,pK01cm)
			       call boost(pK02,    beta,pK02cm)
			       call boost(plambda1,beta,plambda1cm)
			       call boost(plambda2,beta,plambda2cm)
c
c                           Check whether momentum is conserved in cm
c
			       pcmx = pK01cm(1) + plambda2cm(1)
			       pcmy = pK01cm(2) + plambda2cm(2)
			       pcmz = pK01cm(3) + plambda2cm(3)
			       pcmtot = sqrt(pcmx**2. + pcmy**2. + pcmz**2.)
c			    write(6,*)pcmtot,pcmx,pcmy,pcmz
			       call hf1(2568,pcmx,weight)
			       call hf1(2569,pcmx,weight)
			       call hf1(2570,pcmx,weight)
			       call hf1(2571,pcmtot,weight)
c
			       pmomK0cm1 = absp(pK01cm)
			       pmomK0cm2 = absp(pK02cm)
			       pmomlambdacm1 = absp(plambda1cm)
			       pmomlambdacm2 = absp(plambda2cm)

			       costhk0lab = pK01(3)/pmomK0lab1
			       costhk0cm  = pK01cm(3)/pmomK0cm1
			       call hf2(8121,costhk0lab,ppbeam,weight)
			       call hf2(8123,costhk0cm ,ppbeam,weight)
c
			       call hf1(4531,rmpippim1,weight) !look for the K0
c                              Poor man's Monte Carlo:
			       pK01(4) = rmK0
			       call missmass(pinitial,pK01,ptemp5,rmisshyper,rmisshyper2) !Hyperon recoil off nucleon
			       call hf1(2544,rmisshyper,weight) !look for the Hyperon in missing mass
			       call hf2(2545,rmisshyper,rmissspec,weight)

			       call hf1(2564,pmomK0cm1,weight)
			       call hf1(2565,pmomK0cm2,weight)
			       call hf1(2566,pmomlambdacm1,weight)
			       call hf1(2567,pmomlambdacm2,weight)
			    endif
c
c
c                           Repeat most of the analysis using smeared variables
c                           ***************************************************
c
c                           Momentum of the recoiling proton reconstructed from the rest...
			    call invarmass(plist(13,i),plist(13,j) ,ptemp5, rmtemp)
			    call invarmass(ptemp5,    plist(13,k)  ,ptemp6, rmtemp)
			    call invarmass(ptemp6,    plist(13,kk) ,ptemp7, rmtemp)
c			    call missmass(plist(1,1),ptemp7,pspect,rmissspec,rmissspec2) !Spectator proton
			    call missmass(ptotrest,ptemp7,pspect,rmissspec,rmissspec2) !Spectator proton

			    recoilp = absp(pspect)
			    call hf1(3551,recoilp,weight)   !MM off pi+ pi- pi- p  --> spectator momentum
c
c			    if(recoilp .le. 0.200) then       	!GeV/c recoil momentum cut
			    if(recoilp .le. 1.200) then       	!GeV/c recoil momentum cut

			       call hf1(3550,rmissspec,weight) !MM off pi+ pi- pi- p  --> spectator

			       call invarmass(plist(13,j),plist(13,k), pK01, rmpippim1)	!IM of pi+ and pi-1
			       call invarmass(plist(13,j),plist(13,kk),pK02, rmpippim2)	!IM of pi+ and pi-2
			       call hf1(3531,rmpippim1,weight) !look for the K0
			       call hf1(3532,rmpippim2,weight) !look for the K0
			       call hf1(3533,rmpippim1,weight) !look for the K0 combinatorics
			       call hf1(3533,rmpippim2,weight) !look for the K0
			       call invarmass(plist(13,i),plist(13,k),ptemp5, rmppim1) !IM of proton and pi-1
			       call invarmass(plist(13,i),plist(13,kk),ptemp5, rmppim2) !IM of proton and pi-2
			       call hf1(3541,rmppim1,weight) !look for the Lambda
			       call hf1(3542,rmppim2,weight) !look for the Lambda
			       call hf1(3543,rmppim1,weight) !look for the Lambda combinatorics
			       call hf1(3543,rmppim2,weight) !look for the Lambda
c
			       call hf1(5531,rmpippim1,weight) !look for the K0
			       call hf1(5532,rmpippim2,weight) !look for the K0
			       call hf1(5533,rmpippim1,weight) !look for the K0 combinatorics
			       call hf1(5533,rmpippim2,weight) !look for the K0 combinatorics
			       call hf1(5541,rmppim1,weight) !look for the Lambda
			       call hf1(5542,rmppim2,weight) !look for the Lambda
			       call hf1(5543,rmppim1,weight) !look for the Lambda combinatorics
			       call hf1(5543,rmppim2,weight) !look for the Lambda
c
			       pinitial(1) = 0.0
			       pinitial(2) = 0.0
			       pinitial(3) = ppbeam 
			       pinitial(4) = sqrt(2.*ppbeam*rmneut + rmneut**2.)
c                              Poor man's Monte Carlo:
			       pK01(4) = rmK0
			       call missmass(pinitial,pK01,ptemp5,rmisshyper,rmisshyper2) !Hyperon recoil off nucleon
			       call hf1(3544,rmisshyper,weight) !look for the Hyperon in missing mass
			       call hf2(3545,rmisshyper,rmissspec,weight)
			    endif
c
c                           MORE analysis histgrams
c                           Compare momenta of the two pi- tracks, both pure and smeared.
			    recoilpim1 = absp(plist(1,k))
			    recoilpim2 = absp(plist(1,kk))
			    call hf2(3552,recoilpim1,recoilpim2,weight)
			    recoilpim1 = absp(plist(13,k))
			    recoilpim2 = absp(plist(13,kk))
			    call hf2(3553,recoilpim1,recoilpim2,weight)
c
			 endif
		      enddo
		   endif
		enddo
	     endif
	  enddo
	endif
	enddo

c
c	Reconstruct (p) pi- pi+ n events (for pi Delta studies on deuterium in CLAS g14 data)
c
 	do i=2,imark
	 lundi = plist(6,i)	              !proton
	 if(lundi.eq.41)then
c	  write(6,*)'Proton    :',(plist(ix,i),  ix=1,4)
	  acceptweight_proton = acceptlist(i)
	  do j=2,imark
             lundj = plist(6,j)	              !pi+
	     if(lundj.eq.+17)then
c		write(6,*)'pi+       :',(plist(ix,j),  ix=1,4)
		acceptweight_piplus = acceptlist(j)
		do k=2,imark
		   lundk = plist(6,k)         !pi-
		   if(lundk.eq.-17)then
c		      write(6,*)'pi-       :',(plist(ix,k),  ix=1,4)
		      acceptweight_piminus = acceptlist(k)
		      do kk=2,imark
			 lundkk = plist(6,kk) !neutral (neutron or photon)
			 if(lundkk.eq.42 .or. lundkk.eq.1)then
c			    write(6,*)'neutral   :',(plist(ix,kk), ix=1,4)
			    acceptweight_neutral = acceptlist(kk)
c
c		 write(6,*)i,j,k,kk
c		 write(6,*)'Proton    :',(plist(ix,i),  ix=1,4),acceptweight_proton
c		 write(6,*)'pi+       :',(plist(ix,j),  ix=1,4),acceptweight_piplus
c		 write(6,*)'pi-       :',(plist(ix,k),  ix=1,4),acceptweight_piminus
c		 write(6,*)'neutral   :',(plist(ix,kk), ix=1,4),acceptweight_neutral,
c     1                                    lundkk
			   
			    if(iloop.eq.3 .and. itrigmask.eq.16)then
			       weight = acceptweight_piplus
     1		               * acceptweight_piminus*acceptweight_neutral
c			       write(6,*)'>>>',weight,
c     1                          acceptweight_piplus,acceptweight_piminus,
c     2                          acceptweight_neutral
			    else
			       weight = 1.0
			    endif
			    call invarmass(plist(1,j),plist(1,k),ptemp5, rmpippim)   !IM of pi+ and pi- (rho)
			    call invarmass(plist(1,j),plist(1,kk),ptemp5, rmpipn)  !IM of pi+ and neutron (Delta)
			    call invarmass(plist(1,k),plist(1,kk),ptemp6, rmpimn)  !IM of pi- and neutron (Delta)
			    call invarmass(ptemp5,    plist(1,k) ,ptemp6, rmtemp)
			    call missmass(plist(1,1),ptemp6,pspect,rmissspec,rmissspec2) !Spectator proton
			    call hf1(2550,rmissspec,weight) !MM off pi+ pi- n  --> spectator
			    recoilp = absp(pspect)
c			    call hf1(2551,recoilp,weight)   !MM off pi+ pi- n  --> spectator momentum
			    call hf1(6531,rmpipn,weight)    !look for the Delta+
			    call hf1(6532,rmpimn,weight)    !look for the Delta-
			    call hf1(6535,rmpippim,weight)  !look for the rho
			    rmpipn2 = rmpipn**2.
			    rmpimn2 = rmpimn**2.
			    call hf2(6631,rmpipn2,rmpimn2,weight)  !Dalitz plot of the two Deltas
			    recoilneut = absp(plist(1,kk))   !momentum of the decay neutron
			    betaneut = recoilneut/plist(5,kk)  !momentum over energy
c			    betaneut = recoilneut/sqrt(recoilneut**2 + rmneut**2)
			    call hf1(6533,recoilneut,weight)  !neutron momentum
			    call hf1(6534,betaneut,weight) !neutron speed (beta)
			    
c                           Repeat using smeared variables
c       
c
c                           Momentum of the recoiling proton reconstructed from the rest...
			    call invarmass(plist(13,kk),plist(13,j) ,ptemp5, rmtemp)  
			    call invarmass(ptemp5,      plist(13,k) ,ptemp6, rmtemp)
			    call missmass(plist(1,1),ptemp6,pspect,rmissspec,rmissspec2) !Spectator proton
			    call hf1(3550,rmissspec,weight) !MM off pi+ pi- n  --> spectator
			    recoilp = absp(pspect)
			    call hf1(3551,recoilp,weight)   !MM off pi+ pi- n  --> spectator momentum
c
c                           Use neutron to plot lab angular distribution
c
			    costhneut = plist(3,kk)/
     1          (sqrt(plist(1,kk)**2 + plist(2,kk)**2 + plist(3,kk)**2))
			    call hf2(8122,costhneut,ppbeam,weight)
			 endif
		      enddo
		   endif
		enddo
	     endif
	  enddo
	endif
	enddo
c
c	Reconstruct  d pi+ pi-  events (for dibaryon studies in g13 and g14)
c
 	do i=2,imark
	lundi = plist(6,i)	!deuteron
	if(lundi.eq.52)then
c	  write(6,*)'Deuteron    :',(plist(ix,i),  ix=1,4)
	   acceptweight_deuteron = acceptlist(i)
	   do j=2,imark
	      lundj = plist(6,j) !pi+
	      if(lundj.eq.+17)then
c		write(6,*)'pi+       :',(plist(ix,j),  ix=1,4)
		 acceptweight_piplus = acceptlist(j)
		 do k=2,imark
		    lundk = plist(6,k) !pi-
		    if(lundk.eq.-17)then
c                   write(6,*)'pi-       :',(plist(ix,k),  ix=1,4)
		       acceptweight_piminus = acceptlist(k)
c
c	       write(6,*)i,j,k
c	       write(6,*)'Deuteron  :',(plist(ix,i),  ix=1,4),acceptweight_deuteron
c	       write(6,*)'pi+       :',(plist(ix,j),  ix=1,4),acceptweight_piplus
c	       write(6,*)'pi-       :',(plist(ix,k),  ix=1,4),acceptweight_piminus
c
c		       if(iloop.eq.3 .and. itrigmask.eq.17)then
c			  weight = acceptweight_piplus
c     1		                * acceptweight_piminus*acceptweight_deuteron
c			       write(6,*)'>>>',weight,
c     1                          acceptweight_piplus,acceptweight_piminus,
c     2                          acceptweight_deuteron
c		       else
c			  weight = 1.0
c		       endif
		       call invarmass(plist(1,i),plist(1,k),ptemp, rmdpim)   !IM of deuteron and pi-
		       call invarmass(plist(1,i),plist(1,j),ptemp, rmdpip)   !IM of deuteron and pi+
		       call invarmass(plist(1,j),plist(1,k),ptemp, rmpippim) !IM of pi+ and pi-
		       call hf1(8801,rmdpip,weight) !look for the d pi+
		       call hf1(8802,rmdpim,weight) !look for the d pi-
		       call hf2(8803,rmdpim,rmdpip,weight) !correlation
		       call hf1(8804,rmpippim,weight) !look for the "rho" in pi pi
c		       call hf1(28,ppbeam,weight)
		       call hf2(8617,rmpippim,ppbeam,weight)
c
c                      Fill a bunch of "basic" momentum and angle histograms for comparison to real data
c                      We want the overall c.m quantities
c
		       e0=etot(plist(1,1)) !defines the beam plus target
		       do n=1,3
			  beta(n) = plist(n,1)/e0
c			  write(6,*)'beta',beta
		       end do
		       call boost(plist(1,i),beta,ptemp) !deuteron
		       pmomdeut = absp(ptemp)
		       costhd   = dot(ptemp,nbeam)/sqrt(psq(ptemp))
		       thdeut   = acos(costhd)*rade
		       call hf1(8601,pmomdeut,weight)
		       call hf1(7611,thdeut,weight)
		       call hf1(8611,costhd,weight)
		       call hf2(8618,pmomdeut,pdeuteronlab,weight)

		       call boost(plist(1,j),beta,ptemp) !pi+
		       pmompiplus    = absp(ptemp)
		       costhpiplus   = dot(ptemp,nbeam)/sqrt(psq(ptemp))
		       thpiplus      = acos(costhpiplus)*rade
		       call hf1(7602,pmompiplus,weight)
		       call hf1(7612,thpiplus,weight)
		       call hf1(8612,costhpiplus,weight)

		       call boost(plist(1,k),beta,ptemp) !pi-
		       pmompiminus    = absp(ptemp)
		       costhpiminus   = dot(ptemp,nbeam)/sqrt(psq(ptemp))
		       thpiminus      = acos(costhpiminus)*rade
		       call hf1(7603,pmompiminus,weight)
		       call hf1(7613,thpiminus,weight)
		       call hf1(8613,costhpiminus,weight)

		       call hf2(8614,costhpiplus,costhpiminus,1.0)
		       call hf2(8615,pmompiplus,pmompiminus,1.0)
		       call hf2(8616,costhd,costhpiminus,1.0)
		       call hf2(8616,costhd,costhpiplus,1.0)

		    endif
		 enddo
	      endif
	   enddo
	endif
	enddo
c
c	Reconstruct  p pi-  events (for calibration reaction p pi- in g14)
c       (Fresh code after loss of some files, starting September 2014)
c
 	do i=2,imark
	lundi = plist(6,i)	!proton
	if(lundi.eq.41)then
c	  write(6,*)'proton    :',(plist(ix,i),  ix=1,4)
	   acceptweight_proton = acceptlist(i)
c	   write(6,*)'>Acceptlist',
c     1           acceptlist(1),acceptlist(2),acceptlist(3),acceptlist(4)
	   do j=2,imark
	      lundj = plist(6,j) !pi-
	      if(lundj.eq.-17)then
c		write(6,*)'pi-       :',(plist(ix,j),  ix=1,4)
		 acceptweight_piminus = acceptlist(j)
c
c		 write(6,*)i,j
c		 write(6,*)'proton    :',(plist(ix,i),  ix=1,4),acceptweight_proton
c		 write(6,*)'pi-       :',(plist(ix,j),  ix=1,4),acceptweight_piminus
c
		 if(iloop.eq.3 .and. itrigmask.eq.18)then
		    weight = acceptweight_proton * acceptweight_piminus
c		    write(6,*)'>>>In the final loop',weight,
c     1                          acceptweight_proton,acceptweight_piminus
		 else
		    weight = 1.0
		 endif
c		 call hf1(601,vertex(1),weight) ! target vertex distribution
c		 call hf1(602,vertex(2),weight) ! target vertex distribution
c		 call hf1(603,vertex(3),weight) ! target vertex distribution

		 call epcm(plist(1,i),plist(1,j),ptemp)
		 call missmass(plist(1,1),ptemp,pspect,rmissspec,rmissspec2) !Spectator proton
		 call hf1(2650,rmissspec,weight) !MM off p pi-  --> spectator
		 recoilp = absp(pspect)
c		 write(6,*)'Recoil momentum',recoilp
		 call hf1(2651,recoilp,weight) !MM off pi+ pi- n  --> spectator momentum

		 call epcm(plist(1,i),plist(13,j),ptemp)
		 call missmass(plist(1,1),ptemp,pspect,rmissspec,rmissspec2) !Spectator proton smeared
		 call hf1(3650,rmissspec,weight) !MM off p pi-  --> spectator
		 recoilp = absp(pspect)
		 call hf1(3651,recoilp,weight) !MM off pi+ pi- n  --> spectator momentum smeared

		 call hf1(2626,ppbeam,weight)  !distribution of beam energies
		 costhpi = plist(3,j)/
     1             (sqrt(plist(1,j)**2 + plist(2,j)**2 + plist(3,j)**2))   !pion lab angle
		 costhp  = plist(3,i)/
     1             (sqrt(plist(1,i)**2 + plist(2,i)**2 + plist(3,i)**2))   !proton lab angle
		 phip = atan2(plist(2,i),plist(1,i)) * rade   !recompute this due to 2-proton confusion
		 if(phip.lt.0.0)phip = phip + 360. !azimuthal angle from 0 to 360.
		 phipim = atan2(plist(2,j),plist(1,j)) * rade !recompute this for guaranteed consistency
		 if(phipim.lt.0.0)phipim = phipim + 360. !azimuthal angle from 0 to 360.

		 call hf2(3122,costhpi,ppbeam,weight)
c
c                Check of lab angle distribution calculations for combinations 
c                consistent with low recoil momentum (i.e. not the combinatoric case
c                of selecting the recoil proton as the active player in the kinematics)
c
		 if(recoilp.le.0.200)then
		    theta_polarpim = acos(costhpi)*rade
		    theta_polarp   = acos(costhp)*rade
		    call hf2(1606,theta_polarpim,phipim,weight)
		 endif

c		 call invarmass(plist(1,i),plist(1,j),ptemp, rmppim) !IM of proton and pi-
c		 call invarmass(plist(1,k),plist(1,kk),ptemp6, rmpimn) !IM of pi- and neutron (Delta)
c		 call invarmass(ptemp5,    plist(1,k) ,ptemp6, rmtemp)
c		 call missmass(plist(1,1),ptemp6,pspect,rmissspec,rmissspec2) !Spectator proton
c		 call hf1(6531,rmpipn,weight) !look for the Delta+
c		 call hf1(6532,rmpimn,weight) !look for the Delta-
c		 call hf1(6535,rmpippim,weight)	!look for the rho
c		 call invarmass(plist(1,i),plist(1,j),ptemp, rmdpip) !IM of deuteron and pi+
c		 call invarmass(plist(1,j),plist(1,k),ptemp, rmpippim) !IM of pi+ and pi-
c		 call hf1(8801,rmdpip,weight) !look for the d pi+
c		 call hf1(8802,rmdpim,weight) !look for the d pi-
c		 call hf2(8803,rmdpim,rmdpip,weight) !correlation
c		 call hf1(8804,rmpippim,weight) !look for the "rho" in pi pi
c		 call hf1(28,ppbeam,weight)
c		 call hf2(8617,rmpippim,ppbeam,weight)
cc
cc                      Fill a bunch of "basic" momentum and angle histograms for comparison to real data
cc                      We want the overall c.m quantities
cc
c		 e0=etot(plist(1,1)) !defines the beam plus target
c		 do n=1,3
c		    beta(n) = plist(n,1)/e0
cc                    write(6,*)'beta',beta
c		 end do
c		 call boost(plist(1,i),beta,ptemp) !deuteron
c		 pmomdeut = absp(ptemp)
c		 costhd   = dot(ptemp,nbeam)/sqrt(psq(ptemp))
c		 thdeut   = acos(costhd)*rade
c		 call hf1(8601,pmomdeut,weight)
c		 call hf1(7611,thdeut,weight)
c		 call hf1(8611,costhd,weight)
c		 call hf2(8618,pmomdeut,pdeuteronlab,weight)

c		 call boost(plist(1,j),beta,ptemp) !pi+
c		 pmompiplus    = absp(ptemp)
c		 costhpiplus   = dot(ptemp,nbeam)/sqrt(psq(ptemp))
c		 thpiplus      = acos(costhpiplus)*rade
c		 call hf1(7602,pmompiplus,weight)
c		 call hf1(7612,thpiplus,weight)
c		 call hf1(8612,costhpiplus,weight)

c		 call boost(plist(1,k),beta,ptemp) !pi-
c		 pmompiminus    = absp(ptemp)
c		 costhpiminus   = dot(ptemp,nbeam)/sqrt(psq(ptemp))
c		 thpiminus      = acos(costhpiminus)*rade
c		 call hf1(7603,pmompiminus,weight)
c		 call hf1(7613,thpiminus,weight)
c		 call hf1(8613,costhpiminus,weight)

c		 call hf2(8614,costhpiplus,costhpiminus,1.0)
c		 call hf2(8615,pmompiplus,pmompiminus,1.0)
c		 call hf2(8616,costhd,costhpiminus,1.0)
c		 call hf2(8616,costhd,costhpiplus,1.0)
c
	      endif
	   enddo
	endif
	enddo

 899	CONTINUE
c	write(6,*)'Got to here 2'

c	GOTO 9999
	
c
	
c	Reconstruct  p p pbar events (for summer 2016+ baryonium exploration for GlueX)
c 
c 	Apply a hack up front to analyze this reaction only if exactly 3 tracks 
c 	are present in the event (+1 for initial state +1 for intermediate state)
c
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

c		       if(iloop.eq.3 .and. itrigmask.eq.17)then
c		 	  weight = acceptweight_proton1
c      1		                * acceptweight_proton2*acceptweight_antiproton
c	               else
c                         weight = 1.0
c		       endif
c       
c                  We want the overall c.m quantities
c               
		 e0 = etot(plist(1,1)) !defines the beam plus target
		 ww = plist(4,1) !W of system 
		 do n=1,3
		    beta(n) = plist(n,1)/e0
c                 write(6,*)'beta',beta
		 end do
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
C		 GOTO 9999
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
C 9999		 continue
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
		 e0 = etot(pimp2apcm) !use the created-baryon pair
		 do n=1,3
		    beta(n) = pimp2apcm(n)/e0
		 end do
c                 write(6,*)'beta',beta
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
c	Reconstruct  Lambda anti-Lambda proton events (for summer 2017+)
c       The same code is supposed to handle Sigma0 anti-Sigma0 events.
c
c 	Skip events that don't have the right number of particles in the tree
c       This is really a hack to speed up the code; could be skipped
c	6    - pure 5 body phase space
c	9    - intermediate resonance decaying to L Lbar
c	13   - intermediate resonance decaying to S0 S0 bar
c
910	continue
c	write(6,*)'Got to baryon-anti-baryon code:', imark
	if(.not.(imark.eq.9 .or. imark.eq.13 .or. imark.eq.6 .or. imark.eq.8))then
c		write(6,*)'Skipping event with wrong number of tracks'
		goto 920
	endif
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
c       Hacky cheat to see whether we are simulating Mechanisms 1, 2, 3 or 4
c
	if(iloop.eq.1)then
c          ipart1 = plist(6,1)  !the LundID of the first  particle in the decay table for this event
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
				do m=2,imarkweight
				   lundm = plist(6,m) !piminus
c 				   write(6,*)'m',m,lundm
				   if(lundm.eq.-17)then
c                                  write(6,*)'pi-minus   :',(plist(ix,m),  ix=1,5)
				      acceptweight_piminus = acceptlist(m)
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
		e0 = etot(plist(1,1)) !defines the beam plus target
		ww = plist(4,1) !W of system 
		do n=1,3
		   beta(n) = plist(n,1)/e0
c		    write(6,*)'beta',beta
		end do
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
c               this next section of code is quite a hack:  to work correctly the input definition file must create particle in the
c               correct order
c
		call invarmass(apcm, pipcm, alamcm,  rmalam) !IM of anti-proton and piplus --> anti-Lambda
		call invarmass(aplab,piplab,alamlab, rmalamlab) !IM of anti-proton and piplus --> anti-Lambda
		if(imechanism.eq.2)then  ! Mechanism 2 (p-anti-Lam) created at upper vertex
		   rvertex(1) = plist(7,j)  !recoil proton (primary) vertex
		   rvertex(2) = plist(8,j)
		   rvertex(3) = plist(9,j)
		   dvertex(1) = plist(7,i)  !decay proton production vertex
		   dvertex(2) = plist(8,i)
		   dvertex(3) = plist(9,i)
		   call invarmass(plist(1,i),plist(1,m),plamlab,rmplamlab)
		   call invarmass(p1cm,pimcm,plamcm, rmplam) !IM of proton1 and piminus    --> Lambda 
		   call invarmass(p2cm,plamcm,pp1plam, rmp1plam) !IM of proton2 (recoil) and Lambda 
		   call invarmass(p2cm,alamcm,pp1alam, rmp1alam) !IM of proton2 (recoil) and anti-Lambda
c		   write(6,*)'Two'
		elseif(imechanism.eq.1)then !Mechanism 1 (Lam-anti-Lam) created at upper vertex
		   rvertex(1) = plist(7,i)  !recoil proton vertex  
		   rvertex(2) = plist(8,i)
		   rvertex(3) = plist(9,i)
		   dvertex(1) = plist(7,j)  !decay proton production vertex
		   dvertex(2) = plist(8,j)
		   dvertex(3) = plist(9,j)
c		   write(6,*)'Primary: rvertex',rvertex
c		   write(6,*)'Decay  : dvertex',dvertex
		   call invarmass(plist(1,j),plist(1,m),plamlab,rmplamlab)
		   call invarmass(p2cm,pimcm,plamcm, rmplam) !IM of proton2 and piminus    --> Lambda 
		   call invarmass(p1cm,plamcm,pp1plam, rmp1plam) !IM of proton1 (recoil) and Lambda 
		   call invarmass(p1cm,alamcm,pp1alam, rmp1alam) !IM of proton1 (recoil) and anti-Lambda 
c		   write(6,*)'One'
		elseif(imechanism.eq.3)then !Mechanism 3 (anti-Lam p) created
		   rvertex(1) = plist(7,i)  !recoil proton vertex   - check for correctness!
		   rvertex(2) = plist(8,i)
		   rvertex(3) = plist(9,i)
		   dvertex(1) = plist(7,j)  !decay proton production vertex
		   dvertex(2) = plist(8,j)
		   dvertex(3) = plist(9,j)
c		   write(6,*)'Primary: rvertex',rvertex
c		   write(6,*)'Decay  : dvertex',dvertex
		   call invarmass(plist(1,j),plist(1,m),plamlab,rmplamlab)
		   call invarmass(p2cm,pimcm,plamcm, rmplam) !IM of proton2 and piminus    --> Lambda 
		   call invarmass(p1cm,plamcm,pp1plam, rmp1plam) !IM of proton1 (recoil) and Lambda 
		   call invarmass(p1cm,alamcm,pp1alam, rmp1alam) !IM of proton1 (recoil) and anti-Lambda 
c		   write(6,*)'Three'
		elseif(imechanism.eq.4)then !Mechanism 4 (anti-Lam Lambda) created at lower vertex
		   rvertex(1) = plist(7,i)  !recoil proton vertex   - check for correctness!
		   rvertex(2) = plist(8,i)
		   rvertex(3) = plist(9,i)
		   dvertex(1) = plist(7,j)  !decay proton production vertex
		   dvertex(2) = plist(8,j)
		   dvertex(3) = plist(9,j)
c		   write(6,*)'Primary: rvertex',rvertex
c		   write(6,*)'Decay  : dvertex',dvertex
		   call invarmass(plist(1,j),plist(1,m),plamlab,rmplamlab)
		   call invarmass(p2cm,pimcm,plamcm, rmplam) !IM of proton2 and piminus    --> Lambda 
		   call invarmass(p1cm,plamcm,pp1plam, rmp1plam) !IM of proton1 (recoil) and Lambda 
		   call invarmass(p1cm,alamcm,pp1alam, rmp1alam) !IM of proton1 (recoil) and anti-Lambda 
c		   write(6,*)'Four'
		else
		   write(6,*)'This hack is messed up at top...',i
		   call exit()
		endif
		call invarmass(plamcm,alamcm,plamalam, rmplamalam)     !IM of Lambda and anti-Lambda 
		call gett(pbeamcm,plamalam,ttplamalam,ttprimeplamalam) !Mandelstam t for beam and Lambda anti-Lambda pair
		call gett(pbeamcm,pp1alam,ttp1alam,ttprimep1alam)      !Mandelstam t for beam and proton anti-Lambda pair
		if(imechanism.eq.1)then  !dreadful hack to make all reaction sequences work
		   call invarmass(plamalam,p1cm,ptemp, rmtemp) !IM of Lambda, anti-Lambda, and remaining proton 
		elseif(imechanism.eq.2)then
		   call invarmass(plamalam,p1cm,ptemp, rmtemp) !IM of Lambda, anti-Lambda, and remaining proton 
		elseif(imechanism.eq.3)then 
		   call invarmass(plamalam,p1cm,ptemp, rmtemp) !IM of Lambda, anti-Lambda, and remaining proton 
		elseif(imechanism.eq.4)then 
		   call invarmass(plamalam,p1cm,ptemp, rmtemp) !IM of Lambda, anti-Lambda, and remaining proton 
		else
		   write(6,*)'This hack is messed up at bottom...',i
		   call exit()
		endif
            	call missingmass(psystem,ptemp,rmissplamalam,rmissplamalam2)
c                write(6,*)'Overall missing mass^2 :',rmissplamalam2
		call hf1(9887,rmissplamalam2,weight) !look at the overall missing mass distribution
		rmissenergy = psystem(5) - ptemp(5)  
c		 write(6,*)'Missing energy :',rmissenergy
c		 write(6,*)psystem(5)
c                write(6,*)ptemp(5)
c 		 write(6,*)'Got to here', rmissplamalam2,rmissenergy
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
c		if(iloop.eq.1)then
c		   if(abs(rmplam-1.115683).gt. 0.000001)then
c		      write(6,*)'Bloody murder',rmplam
c		      write(6,*)i,j,k,l,m
c		      write(6,*)'Proton 1    :',(plist(ix,i),ix=1,5),acceptweight_proton1
c		      write(6,*)'Proton 2    :',(plist(ix,j),ix=1,5),acceptweight_proton2
c		      write(6,*)'Anti-proton :',(plist(ix,k),ix=1,5),acceptweight_antipro
c		      write(6,*)'Pion plus   :',(plist(ix,l),ix=1,5),acceptweight_piplus
c		      write(6,*)'Pion minus  :',(plist(ix,m),ix=1,5),acceptweight_piminus
c		      call exit()
c		   else
cc		      write(6,*)'All OK', icount
cc		      write(6,*)i,j,k,l,m
cc		      write(6,*)'Proton 1    :',(plist(ix,i),ix=1,5),acceptweight_proton1
cc		      write(6,*)'Proton 2    :',(plist(ix,j),ix=1,5),acceptweight_proton2
cc		      write(6,*)'Anti-proton :',(plist(ix,k),ix=1,5),acceptweight_antipro
cc		      write(6,*)'Pion plus   :',(plist(ix,l),ix=1,5),acceptweight_piplus
cc		      write(6,*)'Pion minus  :',(plist(ix,m),ix=1,5),acceptweight_piminus
c		   endif
c		endif
		call hf1(9801,rmplam,weight) !look at the Lambda  mass distribution 
		call hf1(9802,rmalam,weight) !look at the anti-Lambda  mass distribution 
c		if(costhplamcm.gt.0.0)then
		   call hf1(9851,rmplamalam,weight) !look at the anti-Lambda Lambda mass distribution
		   call hf1(9856,-ttprimeplamalam,weight) !look at the t distribution for the produced pair 
c		else
		   call hf1(9852,rmp1alam,weight)   !look at the anti-Lambda proton mass distribution
		   call hf1(9857,-ttprimep1alam,weight)   !look at the t distribution for the produced pair 
c		endif
		
		call hf2(9850,rmp1plam,rmp1alam,weight)   !correlation (not Dalitz plot)
		call hf2(9855,rmplamalam,rmp1alam,weight) !correlation (not Dalitz plot)
		rmp1plam2 = rmp1plam**2.
		rmp1alam2 = rmp1alam**2.
		call hf2(9854,rmp1plam2,rmp1alam2,weight) ! Dalitz plot
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
		   call makebeta(plamalam,betapair)
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
		   if(costhplamcm.gt.0.0)then !see what Lambda distribution looks like in pair frame for forward Lambdas in CM frame.
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
c		   endif
c       
c               Construct the proton and anti-proton 4-vectors in the rest frame of each respective
c               hyperon.   This requires boosting them first from c.m. to the mutual pair rest frame and then to the
c               rest frame of each hyperon separately.
c
		   call boost(p1cm, betapair, p1pair) !boost both protons to hyperon pair rest frame
		   call boost(p2cm, betapair, p2pair) 
		   call boost(apcm, betapair, appair) !the anti-proton that is part of the anti-lambda
		   call makebeta(plampair,betatophyp)
		   if(imechanism.eq.2)then  !hack to resolve combinatoric issue
		      call boost(p1pair, betatophyp, pphyp) !the Lambda's proton in the Lambda rest frame
		   elseif(imechanism.eq.1 .or. imechanism.eq.3 .or. imechanism.eq.4)then
		      call boost(p2pair, betatophyp, pphyp) !the Lambda's proton in the Lambda rest frame
		   else
		      write(6,*)'very confused...'
		      call exit()
		   endif
		   call makebeta(alampair,betatoahyp)
		   call boost(appair, betatoahyp, aphyp) !the anti_Lambda's anti-proton in anti-Lambda rest frame
c		write(6,*)'Proton in Hyperon rest frame          :',pphyp
c		write(6,*)'Anti-proton in Anti-Hyperon rest frame:',aphyp
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
c                write(6,*)'Got to here...'
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
	
c       
c       End of event-examination loop
c       
 1000   continue
c
 2000	return
	end

c
c	Shut down analysis 
c
	subroutine analysis_finish
c
	parameter (nwpawc = 20000000)
	character*256 prg_name, histnm, filenm, prefix
	real corrsum(3,5),correrr(3,5)
	data alpha,alphabar/0.642,-0.642/  !weak decay asymmetry
c	data alpha,alphabar/1.0 ,-1.0/  !HACK FOR TESTS 
	common /pawc/hmemor(nwpawc)
	common /stuff/poldir(4),xx,idocuts,icount,ikeep,ihackflag,
     1                iattempt
	common /beam/pbeam,ppbeam,rmbeam,rmtarg,jcode,xsigma,ysigma,ncode
	common /run_control/n_arg,prg_name,filenm,histnm,prefix,
     1                      irestart,lfl,lprefix,maxbosoutput
	common /correlations/corrsum,correrr
c     1     polsumplamX,polsumplamY,polsumplamZ,
c     1     polsumalamX,polsumalamY,polsumalamZ,
c     1	   corrsumXX,corrsumYY,corrsumZZ,corrsumZY,corrsumYZ,
c     1	   corrsumXY,corrsumYX,corrsumXZ,corrsumZX
c
c	write(6,*)'For beam momentum',pbeam,' MeV/c,'
	percent  = float(ikeep)/float(icount) * 100.
	fraction = float(icount)/float(iattempt) * 100.
	write(6,10)icount,ikeep,percent,fraction
 10	format(1x,'Analyzed',i8,' events',
     1         1x,';  Passed',i8,' events',
     1         1x,';  Events accepted =',f7.2,' %',/,
     1         1x,
     1  'Analyzed events percent of all generation attempts',f7.2,' %')
        open(unit=8,file='mc_generator.out',status='unknown',
     1                                       access='append')
	write(8,20)pbeam,icount,ikeep,percent
 20	format(1x,f7.2,2i10,f7.2)
	close(8)
c
c       Special output for the case of Lambda-anti-Lambda process
c
	write(6,*)' '
	write(6,*)'For Lambda anti-Lambda analysis: polarization results (%):'

	do ii=1,3
	   corrsum(ii,4) = corrsum(ii,4) / float(icount)
	   corrsum(ii,5) = corrsum(ii,5) / float(icount)
	   correrr(ii,4) = sqrt(correrr(ii,4) / float(icount) - corrsum(ii,4)**2.)
	   correrr(ii,5) = sqrt(correrr(ii,5) / float(icount) - corrsum(ii,5)**2.)
	   corrsum(ii,4) = (3./alpha) * corrsum(ii,4) * 100.
	   corrsum(ii,5) = (3./alpha) * corrsum(ii,5) * 100.
	   correrr(ii,4) = (3./alpha) * correrr(ii,4) / sqrt(float(icount)) * 100.  !error on the mean
	   correrr(ii,5) = (3./alpha) * correrr(ii,5) / sqrt(float(icount)) * 100.
	enddo
	
	write(6,40)'     Lambda Polarization (X,Y,Z): ',
     1      (corrsum(ii,4),ii=1,3),(correrr(ii,4),ii=1,3)
	write(6,40)'anti-Lambda Polarization (X,Y,Z): ',
     1      (corrsum(ii,5),ii=1,3),(correrr(ii,5),ii=1,3)
 40	format(1x,a,3f8.2,5x,3f8.2)

	write(6,*)' '
	write(6,*)'Correlation Matrix:       '
	rnorm = 9./(alpha*alphabar)
	do ii=1,3
	   do jj=1,3
	      corrsum(ii,jj) = corrsum(ii,jj) / float(icount)
	      correrr(ii,jj) =
     1	            sqrt(correrr(ii,jj) / float(icount) - corrsum(ii,jj)**2.)
	      corrsum(ii,jj) = rnorm * corrsum(ii,jj) * 100.
	      correrr(ii,jj) = rnorm * correrr(ii,jj) / sqrt(float(icount)) * 100. !error on the mean
	      correrr(ii,jj) = abs(correrr(ii,jj))
	   enddo
	enddo
	do ii=1,3
	   write(6,50)(corrsum(ii,jj),jj=1,3),(correrr(ii,jj),jj=1,3)
	enddo
 50	format(35x,3f8.2,5x,3f8.2)
c
c       Compute the singlet fraction
c
	sf = 25.0*(1 + (corrsum(1,1) - corrsum(2,2) + corrsum(3,3))/100.)
	write(6,60)'Singlet Fraction (%)              ',sf
 60	format(1x,/,1x,a,f8.2)
	write(6,*)' '
c
c       End of Lambda anti-Lambda section
c
	if(idocuts.eq.1)then
	   write(6,*)'Writing histogram file: ',histnm(1:lfl)
	   call hropen(1,'MCDIR',histnm(1:lfl),
     1  	'N',1024,istat)
	   call hcdir('//MCDIR',' ') !Memory
	   call hmdir('RAWEVENT',' ')	
	   call hmdir('ACCEPTED',' ')
	   call hmdir('WEIGHTED',' ')
	   call hcdir('//PAWC/RAWEVENT',' ')
	   call hcdir('//MCDIR/RAWEVENT',' ')
	   call hrout(0,icycle,0)
	   call hcdir('//PAWC/ACCEPTED',' ')
	   call hcdir('//MCDIR/ACCEPTED',' ')
	   call hrout(0,icycle,0)
	   call hcdir('//PAWC/WEIGHTED',' ')
	   call hcdir('//MCDIR/WEIGHTED',' ')
	   call hrout(0,icycle,0)
	   call hrend('MCDIR')
	endif
c       
	open(unit=20,file='randomseed.num',status='unknown')
	iseq = 1
	call grndmq(iseed1,iseed2,iseq,'G')
	write(20,*)iseq,iseed1,iseed2
	write(6,*)'Saving random number seeds to "randomseed.num":'
	write(6,*)'Sequence:',iseq
	write(6,*)'Seed1   :',iseed1
	write(6,*)'Seed2   :',iseed2
	close(20)
c
c	close(4)
c	write(6,*)'Closing file for saving Dalitz plot data for x(1280)'
c
	write(6,*)'Output files written.'
	return
	end 

c
c	We use the usual coordinates for the experiment
c	1 = X = to the left as viewed in beam direction
c	2 = Y = up
c	3 = Z = beam direction
c
	function iaccept(p1,p2)
	logical iaccept
	real p1(5),p2(5)
	data tantest/0.9/
c
c	iaccept = .true.
c	return
c
	r1     = sqrt(p1(2)**2 + p1(3)**2)
	tanth1 = r1/p1(1)
	r2     = sqrt(p2(2)**2 + p2(3)**2)
	tanth2 = r2/p2(1)
	if(abs(tanth1).lt.tantest .and. abs(tanth2).lt.tantest)then
		iaccept = .true.
	else
		iaccept = .false.
	end if
	return
	end

c
c	We use the usual coordinates for the experiment
c	1 = X = to the left as viewed in beam direction
c	2 = Y = up
c	3 = Z = beam direction
c	angles are in degrees
c
	function iaccept_l3h(theta_polar,theta_xz,theta_yz,pmom,lund)
	logical iaccept_l3h
	data pmin,pmax/0,450./	!momentum cuts
	iaccept_l3h = .false.
	if(pmom.gt.pmin .and. pmom.lt.pmax)then
		iaccept_l3h = .true.
	end if
	return
	x = lund  !fool compiler
	x = theta_polar
	x = theta_xz
	x = theta_yz
	end

c
c	We use the usual coordinates for the experiment
c	1 = X = to the left as viewed in beam direction
c	2 = Y = up
c	3 = Z = beam direction
c	angles are in degrees
c
	function iaccept_k(theta_polar,theta_xz,theta_yz,pmom,lund)
	logical iaccept_k
	if(abs(theta_polar).lt.10.0)then
		iaccept_k = .true.
	else
		iaccept_k = .false.
	end if
	return
	x = lund  !fool compiler
	x = theta_polar
	x = theta_xz
	x = theta_yz
	x = pmom
	end
c
c	We use the usual coordinates for the experiment
c	1 = X = to the left as viewed in beam direction
c	2 = Y = up
c	3 = Z = beam direction
c	angles are in degrees
c
	function iaccept_p(theta_polar,theta_xz,theta_yz,pmom,lund)
	logical iaccept_p
	iaccept_p = .false.
	if(abs(theta_polar).gt.90.0)then
		if(pmom.gt.250.0)then
			iaccept_p = .true.
		end if
	end if
	return
	x = lund  !fool compiler
	x = theta_polar
	x = theta_xz
	x = theta_yz
	x = pmom
	end

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
	acceptweight  = 1.0
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
c	We use the usual coordinates for the experiment
c	1 = X = to the left as viewed in beam direction
c	2 = Y = up
c	3 = Z = beam direction
c
	function iseepi(ekin)
	logical iseepi
 	data ekmin /0.030/	!GeV
c
c	iseepi = .true.
c	return
c
	if(ekin.gt.ekmin)then
		iseepi = .true.
	else
		iseepi = .false.
	end if
	return
	end

c
c	Compute invariant rest mass
c	1 = px; 2 = py; 3=pz; 4=rmass; 5=total energy
c	Input is in GeV, output is in GeV/c2
c
	subroutine invarmass(p1,p2,pt,rmassinv)
	real p1(5),p2(5),pt(5)
	do i=1,3
	   pt(i) = p1(i) + p2(i)
	end do
	pt(5) = p1(5) + p2(5)
	rmassinv =pt(5)*pt(5)-pt(1)*pt(1)-pt(2)*pt(2)-pt(3)*pt(3)
	if(rmassinv.gt.0)then
		rmassinv=sqrt(rmassinv)	!GeV
		pt(4) = rmassinv
	else
		write(6,*)'INVARMASS Ack! '
		write(6,*)p1
		write(6,*)p2
		write(6,*)pt
		call exit()
	end if
	return
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

c
c	Routine to smear the particle's momentum according to a parameter.
c	Smears each component the same way, thereby smearing both
c	energy and angle.  The "y * smearscale" term is an ad hoc smearing that is
c       independent of momentum, and is supposed to represent some kind
c       of "absolute" smearing not dependent upon momentum.  This would be 
c       related to energy-loss straggling, for example.
c
c       CLAS version
c
c       5-15-2015  - turn off the ad hoc smearing if we are not doing any other
c                    smearing as well.  All or nothing...
c                  - also, add the energy parameter to this calculation
c       11-11-2015 - tests with various forms
c
	subroutine psmear(ppure,pspect)
	real ppure(5),pspect(5),ptotrest(5)
        common /analysis/dpoverp_in,dpoverp_out,ptotrest
	common /monte/ iseed,xflat,xnormal,xbw
c
	smearscale = dpoverp_out   ! 0.005  !0.005
	pspect(4) = ppure(4)		!mass
	if(dpoverp_out.gt.0.0)then
c	   call rannor(x,dummy)
c	   call rannor(y,dummy)
c	   pmag      = absp(ppure)
c	   pmagsmear = pmag*(1.0 + dpoverp_out * x) + smearscale * y
c	   pmagsmear = pmag +   smearscale * y
c	   pscale    = pmagsmear/pmag
c	   pspect(1) = ppure(1)*pscale
c	   pspect(2) = ppure(2)*pscale
c	   pspect(3) = ppure(3)*pscale
	   do i=1,3
	      call rannor(x,dummy)
	      call rannor(y,dummy)
	      pspect(i) = ppure(i)*(1.0 + dpoverp_out * x) + y * smearscale
	   enddo
c
c          Rescale to keep the original magnitude, but now the angle is smeared
c
	   pmagpure   = absp(ppure)
	   pmagsmear  = absp(pspect)
	   ratio      = pmagpure/pmagsmear
	   do i = 1,3
	      pspect(i) = pspect(i)*ratio
	   enddo
	else
	   do i=1,3
	      pspect(i) = ppure(i)
	   enddo
	endif
	pspect(5) = etot(pspect)
c
	xnormal = x !for diagnostic plotting
c	write(6,*)dpoverp_out,x,y
c	write(6,*)ppure
c	write(6,*)pspect
	return
	end

c
c	Routine to smear the particle's momentum according to a parameter.
c	Smears each component the same way, thereby smearing both
c	energy and angle.  
c
c       GlueX version
c
	subroutine gsmear(ppure,pspect,lund)
	real ppure(5),pspect(5),ptotrest(5)
	data rade /57.29577951/

	common /analysis/dpoverp_in,dpoverp_out,ptotrest
	common /monte/ iseed,xflat,xnormal,xbw
c
c	write(6,*)'Using GlueX smearing...' 
	rmas   = ppure(4)	!GeV
	px     = ppure(1)	!GeV/c 
	py     = ppure(2)	!GeV/c
	pz     = ppure(3)	!GeV/c
	pmom   = sqrt(px**2.+py**2.+pz**2)  !GeV/c
	evalue = sqrt(pmom**2. + rmas**2.)  !GeV
	if(pmom.ne.0)then
	   arg=pz/pmom
	   if(arg.gt. 1.0)arg= 1.0
	   if(arg.lt.-1.0)arg=-1.0
	   theta_polar = acos(arg) * rade !lab angles, in degrees
	   theta_xz    = atan2(px,pz) * rade
	   theta_yz    = atan2(py,pz) * rade
	   phi         = atan2(py,px) * rade
	   if(phi.lt.0.0)phi = phi + 360. !azimuthal angle from 0 to 360.
	   costhetapolar = arg
	else
	   theta_polar = 1000.
	   theta_xz    = 1000.
	   theta_yz    = 1000.
	   phi         = 1000.
	   costhetapolar = -999.
	end if
	pspect(4) = ppure(4)    !mass
	if(lund.ne.1)then       !charged particle: momentum smearing including angles
	   if(dpoverp_out.ne.0)then
	      dpop = sqrt(dpoverp_out**2. + (0.005*px)**2.)
	      call rannor(v,dummy)
	      pspect(1) = ppure(1)*(1.0 + dpop * v)
	      dpop = sqrt(dpoverp_out**2. + (0.005*py)**2.)
	      call rannor(v,dummy)
	      pspect(2) = ppure(2)*(1.0 + dpop * v)
	      dpop = sqrt(dpoverp_out**2. + (0.005*pz)**2.)
	      call rannor(v,dummy)
	      pspect(3) = ppure(3)*(1.0 + dpop * v)
	   else
	      pspect(1) = ppure(1)
	      pspect(2) = ppure(2)
	      pspect(3) = ppure(3)
	   endif
	else                   !photon:  energy smearing only, no angle smearing
	   if(theta_polar.lt.11.)then        !FCAL region
	      deoe = 0.05/sqrt(evalue)
	   else                              !BCAL region
	      deoe = 0.135/sqrt(evalue) 
	   endif
	   call rannor(v,dummy)
	   ratio = 1 + deoe*v   !fractional change in photon energy
	   pspect(1) = ppure(1) * ratio
	   pspect(2) = ppure(2) * ratio
	   pspect(3) = ppure(3) * ratio
	endif
	pspect(5) = pspect(1)**2.+pspect(2)**2+pspect(3)**2
	pspect(5) = sqrt(pspect(5) + rmas*rmas)
c	write(6,10)'>>',ppure,rmas,lund,pmom,evalue
c	write(6,10)'> ',pspect,v
c 10	format(1x,a,6f8.4,i4,2f8.4)
	return
	end


c
c       Compute the Time-of-Flight paddle hit by a track, given its
c       momentum and angle at the target.  Based on a translation table
c       built by Rob Feuerbach on 5-18-99.
c
c       2-28-01 Fixed a few small bugs for the negative tracks
c       9-8-01  Change the matching algorithm to be more precise.  This was
c               needed to help deal with TOF bar knockouts in the forward
c               directions.
c       9-8-01  This version branched from an older one: the table it uses
c               to match TOF counters to tracks has twice the resolution both
c               in momentum and in angle.
c
	subroutine tofbar(pmomgev,theta,icharge,numtofbar,pathtobar)
c
c       Table encoding is: (momentum index, angle index, 1=momentum value
c       2=angle value, 3=tof number hit, 4=pathlength in cm 
	real ptable(221,285,4) 
	real ntable(221,285,4) 
	real pvalues(221)
	real avalues(285)
c
	common /tofbarlocal1/ptable    !separate commons to make g77 happy
	common /tofbarlocal2/ntable
	common /tofbarlocal3/first,ifail,pvalues,avalues
	common /monte/ iseed,xflat,xnormal,xbw
	data ifirst,ifail/1,0/
c
	numtofbar = -1
	pathtobar = 0.0
	if(ifirst.eq.1)then
	   ifirst = 0
	   write(6,*)'Reading ptable of (p,angle)--->time-of-flight paddle'
c	   open(unit=1,file='scid_p_theta.txt',status='old')
c	   open(unit=1,file='scid_p2.5_theta.txt',status='old')
c	   open(unit=2,file='scid_p_theta_neg.txt',status='old')
c	   open(unit=1,file='scid_p_142theta_pos.txt',status='old')
c	   open(unit=2,file='scid_p_142theta_neg.txt',status='old')
c	   open(unit=1,file='scid_wires_p_theta_run11821.txt',status='old')
c	   open(unit=2,file='scid_wires_n_theta_run11821.txt',status='old')
	   open(unit=1,file='scid_wires_p_run11821_s1.txt',status='old')
	   open(unit=2,file='scid_wires_n_run11821_s1.txt',status='old')
	   do i=1,221      !momenta
	      do j=1,285   !angles
		 do k=1,4
		    ptable(i,j,k) = -1
		    ntable(i,j,k) = -1
		 enddo
	      enddo
	   enddo
	   write(6,*)'    Positive particles...'
	   do i=1,221
	      do j=1,285
		 read(1,*,end=20)pval,aval,itof,plen
		 ptable(i,j,1) = pval
		 ptable(i,j,2) = aval
		 ptable(i,j,3) = itof
		 ptable(i,j,4) = plen
c		 write(6,*)i,j,pval,aval,itof,plen,pkeep
		 if(j.eq.1)pkeep = pval
		 if(pval.ne.pkeep)then
		    write(6,*)'Deep trouble:',i,j,pval,pkeep,aval,itof,plen
		    call exit()
		 end if
	      enddo
	   enddo
 20	   close(1)
	   do i=1,221
	      do j=1,285
		 do k=1,4
		    if(ptable(i,j,k) .eq. -1.)then
		       write(6,*)'Empty table entry: ',i,j,k
		       call exit()
		    end if
		 enddo
	      enddo
	   enddo
	   write(6,*)'    Negative particles...'
	   do i=1,221
	      do j=1,285
		 read(2,*,end=30)pval,aval,itof,plen
		 ntable(i,j,1) = pval
		 ntable(i,j,2) = aval
		 ntable(i,j,3) = itof
		 ntable(i,j,4) = plen
c		 write(6,*)i,j,pval,aval,itof,plen,pkeep
		 if(j.eq.1)pkeep = pval
		 if(pval.ne.pkeep)then
		    write(6,*)'Deep trouble:',i,j,pval,pkeep,aval,itof,plen
		    call exit()
		 end if
	      enddo
	   enddo
 30	   close(2)
	   do i=1,221
	      do j=1,285
		 do k=1,4
		    if(ntable(i,j,k) .eq. -1.)then
		       write(6,*)'Empty table entry: ',i,j,k
		       call exit()
		    end if
		 enddo
	      enddo
	   enddo
c       
c          Generate a table of momentum values that gives the bin edges
c          at the low momentum end of each bin in momentum
c       
	   do i=2,221
	      pvalues(i) = 0.5*(ptable(i,1,1) + ptable(i-1,1,1))
c	      write(6,*)i,pvalues(i),ptable(i,1,1)
	   enddo
	   pvalues(1) = ptable(1,1,1) - 0.01
c       
c          Generate a table of lab angle values (in degrees) that gives
c          the bin edges at the high angle side of each bin.  Since we
c          are using 1 degree binning this is just a 0.5 degree offset.
c       
	   do j = 1,285
	      avalues(j) = ptable(1,j,2) + 0.5
c	      write(6,*)j,avalues(j),ptable(1,j,2)
	   enddo
	endif
c
c       Stupid search algorithm to see which TOF was hit
c
	x = ran(iseed)
	ivalu = 4*x
	if(ivalu.eq.0)then    !50% chance of no shift, 25% shift to either side
	   idelt = -1
	elseif (ivalu.eq.3)then
	   idelt = +1
	else
	   idelt = 0
	end if
	fuzz  = 6.*x - 3.  !guesstimated width in cm if TOF hit is uncertain by 1
c
	pvalue = pmomgev
	if(theta.gt.142.0 .or. theta.lt.8.0)then
	   numtofbar = -1
	   pathtobar = 0.0
	   return
	endif
	if(pvalue.lt.ptable(1,1,1))then
	   numtofbar=-1
	   pathtobar = 0.0
	   return
	endif
	if(icharge.gt.0)then
	   do i=221,1,-1
	      if(pvalue.gt.pvalues(i))then
		 do j=1,285
		    if(theta.lt.avalues(j))then
cc		       numtofbar = ptable(i,j,3) + idelt
		       numtofbar = ptable(i,j,3)
		       pathtobar = ptable(i,j,4) + fuzz
		       goto 40
		    endif
		 enddo
	      endif
	   enddo
	else
	   do i=221,1,-1
	      if(pvalue.gt.pvalues(i))then
		 do j=1,285
		    if(theta.lt.avalues(j))then
cc		       numtofbar = ntable(i,j,3) + idelt
		       numtofbar = ntable(i,j,3)
		       pathtobar = ntable(i,j,4) + fuzz
		       goto 40
		    endif
		 enddo
	      endif
	   enddo
	endif
	ifail = ifail + 1
       	write(6,*)ifail, ' Failed to find a TOF match:' ,pmomgev,theta
	numtofbar = -1
	pathtobar = 0.0
	return
 40	continue
	if(numtofbar.lt.1)numtofbar  = 0   !undo idelt for low  angles.
	if(numtofbar.gt.48)numtofbar = 48  !undo idelt for high angles.
c	if(numtofbar.eq.0)
c	1    write(6,*)icharge,pmomgev,theta,numtofbar,pathtobar
	return
	end

c
c       Compute the Time-of-Flight paddle hit by a track, given its
c       momentum and angle at the target.  Based on a translation table
c       built by Rob Feuerbach on 5-18-99.
c
c       2-28-01 Fixed a few small bugs for the negative tracks
c       9-8-01  Change the matching algorithm to be more precise.  This was
c               needed to help deal with TOF bar knockouts in the forward
c               directions.
c
	subroutine tofbarb(pmomgev,theta,icharge,numtofbar,pathtobar)
c
c       THIS IS THE 'OLDER' VERSION WITH LESS FINE ANGLE AND MOMENTUM BINS
c
c       Table encoding is: (momentum index, angle index, 1=momentum value
c       2=angle value, 3=tof number hit, 4=pathlength in cm 
	real ptable(111,143,4) 
	real ntable(111,143,4) 
	real pvalues(111)
	real avalues(143)
c
	common /tofbarlocal1b/ptable 
	common /tofbarlocal2b/ntable
	common /tofbarlocal3b/first,ifail,pvalues,avalues
	common /monte/ iseed,xflat,xnormal,xbw
	data ifirst,ifail/1,0/
c
	numtofbar = -1
	pathtobar = 0.0
	if(ifirst.eq.1)then
	   ifirst = 0
	   write(6,*)'Reading ptable of (p,angle)--->time-of-flight paddle'
c	   open(unit=1,file='scid_p_theta.txt',status='old')
c	   open(unit=1,file='scid_p2.5_theta.txt',status='old')
c	   open(unit=2,file='scid_p_theta_neg.txt',status='old')
	   open(unit=1,file='scid_p_142theta_pos.txt',status='old')
	   open(unit=2,file='scid_p_142theta_neg.txt',status='old')
	   do i=1,111      !momenta
	      do j=1,143   !angles
		 do k=1,4
		    ptable(i,j,k) = -1
		    ntable(i,j,k) = -1
		 enddo
	      enddo
	   enddo
	   write(6,*)'    Positive particles...'
	   do i=1,111
	      do j=1,143
		 read(1,*,end=20)pval,aval,itof,plen
		 ptable(i,j,1) = pval
		 ptable(i,j,2) = aval
		 ptable(i,j,3) = itof
		 ptable(i,j,4) = plen
c		 write(6,*)i,j,pval,aval,itof,plen,pkeep
		 if(j.eq.1)pkeep = pval
		 if(pval.ne.pkeep)then
		    write(6,*)'Deep trouble:',i,j,pval,pkeep,aval,itof,plen
		    call exit()
		 end if
	      enddo
	   enddo
 20	   close(1)
	   do i=1,111
	      do j=1,143
		 do k=1,4
		    if(ptable(i,j,k) .eq. -1.)then
		       write(6,*)'Empty table entry: ',i,j,k
		       call exit()
		    end if
		 enddo
	      enddo
	   enddo
	   write(6,*)'    Negative particles...'
	   do i=1,111
	      do j=1,143
		 read(2,*,end=30)pval,aval,itof,plen
		 ntable(i,j,1) = pval
		 ntable(i,j,2) = aval
		 ntable(i,j,3) = itof
		 ntable(i,j,4) = plen
c		 write(6,*)i,j,pval,aval,itof,plen,pkeep
		 if(j.eq.1)pkeep = pval
		 if(pval.ne.pkeep)then
		    write(6,*)'Deep trouble:',i,j,pval,pkeep,aval,itof,plen
		    call exit()
		 end if
	      enddo
	   enddo
 30	   close(2)
	   do i=1,111
	      do j=1,143
		 do k=1,4
		    if(ntable(i,j,k) .eq. -1.)then
		       write(6,*)'Empty table entry: ',i,j,k
		       call exit()
		    end if
		 enddo
	      enddo
	   enddo
c       
c          Generate a table of momentum values that gives the bin edges
c          at the low momentum end of each bin in momentum
c       
	   do i=2,111
	      pvalues(i) = 0.5*(ptable(i,1,1) + ptable(i-1,1,1))
c	      write(6,*)i,pvalues(i),ptable(i,1,1)
	   enddo
	   pvalues(1) = ptable(1,1,1) - 0.01
c       
c          Generate a table of lab angle values (in degrees) that gives
c          the bin edges at the high angle side of each bin.  Since we
c          are using 1 degree binning this is just a 0.5 degree offset.
c       
	   do j = 1,143
	      avalues(j) = ptable(1,j,2) + 0.5
c	      write(6,*)j,avalues(j),ptable(1,j,2)
	   enddo
	endif
c
c       Stupid search algorithm to see which TOF was hit
c
	x = ran(iseed)
	ivalu = 4*x
	if(ivalu.eq.0)then    !50% chance of no shift, 25% shift to either side
	   idelt = -1
	elseif (ivalu.eq.3)then
	   idelt = +1
	else
	   idelt = 0
	end if
	fuzz  = 6.*x - 3.  !guesstimated width in cm if TOF hit is uncertain by 1
c
	pvalue = pmomgev
	if(theta.gt.142.0 .or. theta.lt.8.0)then
	   numtofbar = -1
	   pathtobar = 0.0
	   return
	endif
	if(pvalue.lt.ptable(1,1,1))then
c	   write(6,*)'TOFBAR: momentum too low...',pvalue
	   numtofbar=-1
	   pathtobar = 0.0
	   return
	endif
	if(icharge.gt.0)then
	   do i=111,1,-1
	      if(pvalue.gt.pvalues(i))then
		 do j=1,143
		    if(theta.lt.avalues(j))then
cc		       numtofbar = ptable(i,j,3) + idelt
		       numtofbar = ptable(i,j,3)
		       pathtobar = ptable(i,j,4) + fuzz
		       goto 40
		    endif
		 enddo
	      endif
	   enddo
	else
	   do i=111,1,-1
	      if(pvalue.gt.pvalues(i))then
		 do j=1,143
		    if(theta.lt.avalues(j))then
cc		       numtofbar = ntable(i,j,3) + idelt
		       numtofbar = ntable(i,j,3)
		       pathtobar = ntable(i,j,4) + fuzz
		       goto 40
		    endif
		 enddo
	      endif
	   enddo
	endif
	ifail = ifail + 1
       	write(6,*)ifail, ' Failed to find a TOF match:' ,pmomgev,theta
	numtofbar = -1
	pathtobar = 0.0
	return
 40	continue
	if(numtofbar.lt.1)numtofbar  = 0   !undo idelt for low  angles.
	if(numtofbar.gt.48)numtofbar = 48  !undo idelt for high angles.
c	if(numtofbar.eq.0)
c     1      write(6,*)icharge,pmomgev,theta,numtofbar,pathtobar
	return
	end

c
c       Compute the drift chamber wires hit by this track, given its
c       momentum and angle at the target.  Based on a translation table
c       built by Rob Feuerbach on 6-29-01.
c
	subroutine dcwire(pmomgev,theta,icharge,idcwirelist)
c
c       Table encoding is: (momentum index, angle index, 
c       1...34 = number of the 1...34 wires)
c 
	integer ptable(111,143,34) 
	integer ntable(111,143,34) 
	real pvalues(111)
	real avalues(143)
c
	common /dclocal1/ptable    !separate commons to make g77 happy
	common /dclocal2/ntable
	common /dclocal3/first,ifail,pvalues,avalues
	common /monte/ iseed,xflat,xnormal,xbw
	integer idcwirelist(34),numberlist(34)
	data ifirst,ifail/1,0/
c
	if(ifirst.eq.1)then
	   ifirst = 0
	   write(6,*)'Reading ptable of (p,angle)--->wire lists'
	   open(unit=1,file='wires_p_theta_run11821.txt',status='old')
	   open(unit=2,file='wires_n_theta_run11821.txt',status='old')
	   do i=1,111      !momenta
	      do j=1,143   !angles
		 do k=1,34
		    ptable(i,j,k) = -1
		    ntable(i,j,k) = -1
		 enddo
	      enddo
	   enddo
	   write(6,*)'    Positive particles...'
	   do i=1,111
	      do j=1,143
		 read(1,*,end=20)pval,aval,numberlist
		 pvalues(i) = pval
		 avalues(j) = aval
		 do k = 1,34
		    ptable(i,j,k) = numberlist(k)
		    if(j.eq.1)pkeep = pval
		    if(pval.ne.pkeep)then
		       write(6,*)'Deep trouble:',i,j,k,pval,pkeep,aval
		       call exit()
		    end if
		 enddo
	      enddo
	   enddo
 20	   close(1)
	   do i=1,111
	      do j=1,143
		 do k=1,34
		    if(ptable(i,j,k) .eq. -1.)then
		       write(6,*)'Empty table entry: ',i,j,k
		       call exit()
		    end if
		 enddo
	      enddo
	   enddo
	   write(6,*)'    Negative particles...'
	   do i=1,111
	      do j=1,143
		 read(2,*,end=30)pval,aval,numberlist
		 do k = 1,34
		    ntable(i,j,k) = numberlist(k)
		    if(j.eq.1)pkeep = pval
		    if(pval.ne.pkeep)then
		       write(6,*)'Deep trouble:',i,j,k,pval,pkeep,aval
		       call exit()
		    end if
		 enddo
	      enddo
	   enddo
 30	   close(1)
	   do i=1,111
	      do j=1,143
		 do k=1,34
		    if(ntable(i,j,k) .eq. -1.)then
		       write(6,*)'Empty table entry: ',i,j,k
		       call exit()
		    end if
		 enddo
	      enddo
	   enddo
	endif
c
	do i=1,34
	   idcwirelist(i) = -1
	enddo
c
c       Stupid search algorithm to see which wires were hit
c
	pvalue = pmomgev
	if(theta.gt.142.0 .or. theta.lt.8.0)then
	   return   !track is not within with polar acceptance
	endif
	if(pvalue.lt.pvalues(1))then
	   return
	endif
	if(icharge.gt.0)then
	   do i=111,1,-1
	      if(pvalue.gt.pvalues(i))then
		 do j=1,143
		    if(theta.lt.avalues(j))then
		       do k=1,34
			  idcwirelist(k) = ptable(i,j,k)
		       enddo
		       goto 40
		    endif
		 enddo
	      endif
	   enddo
	else
	   do i=111,1,-1
	      if(pvalue.gt.pvalues(i))then
		 do j=1,143
		    if(theta.lt.avalues(j))then
		       do k=1,34
			  idcwirelist(k) = ntable(i,j,k)
		       enddo
		       goto 40
		    endif
		 enddo
	      endif
	   enddo
	endif
	ifail = ifail + 1
       	write(6,*)ifail, ' Failed to find a DC match:' ,pmomgev,theta
	return
 40	continue
	return
	end


c
c       This code is used to model the CLAS scintillator (TOF) 
c       trigger inefficiency for the g1a data set 
c       that was noted and studied by Burin A.
c
	function trigger_CLAS(isector,numtofbar)
	logical trigger_CLAS
	integer isec_count(6)
	real effys(6),toftrig(48,6)  !48 TOF counters by 6 sectors
	common /monte/ iseed,xflat,xnormal,xbw
	common /trigger_CLASlocal/toftrig,isec_count
	data isec_count/0,0,0,0,0,0/
	data ifirst/1/
c
	trigger_CLAS = .true.
	return                      !short circuit this whole code 9-4-14
c
	if(ifirst.eq.1)then
	   ifirst = 0
	   write(6,*)'Load table of TOF paddle trigger effics'
	   open(unit=1,file='tof_paddle_trigger.txt',status='old')
	   do j=1,48
	      do i=1,6
		 toftrig(j,i) = -1
	      enddo
	   enddo
 10	   read(1,*,end=50)idnum,effys
	   if(idnum.ge.1.and.idnum.le.48)then
	      do i=1,6
		 if(effys(i).ge.0 .and. effys(i).le.100)then
		    toftrig(idnum,i)=effys(i)/100.
		 else
		    write(6,*)'Table error:',idnum,effys
		    call exit()
		 endif
	      enddo
	   else
	      write(6,*)'Table error:',idnum,effys
	      call exit()
	   endif
	   goto 10
 50	   close(1)
	endif
c
	if(isector.lt.1 .or. isector.gt.6)then
	   write(6,*)'TOF Sector screw up',isector,numtofbar
	   trigger_CLAS = .false.
	   return
c	   call exit()
	endif
	isec_count(isector) = isec_count(isector) + 1
	if(numtofbar.le.0)then
	   trigger_CLAS = .false.
	   return
	endif
c
c       We have lost the trigger efficiency file, so...
c
c	triggereff = toftrig(numtofbar,isector)
	triggereff = 1.0            
	xrandom = 0
	if(triggereff.gt.0)then
	   xrandom = ran(iseed)
	   if(xrandom.le.triggereff)then
	      trigger_CLAS = .true.
	   else
	      trigger_CLAS = .false.
	   endif
	endif
c	write(6,*)'trigger_CLAS:',isector,numtofbar,
c     1      triggereff,trigger_CLAS,xrandom !,isec_count
	return
	end


	function knockout_CLAS_tof(isector,numtofbar)
	logical knockout_CLAS_tof
	integer effys(6),toflist(48,6)  !48 TOF counters by 6 sectors
	common /knockout_CLAS_toflocal/toflist
	data ifirst/1/
c
	knockout_CLAS_tof = .true.
c
	if(ifirst.eq.1)then
	   ifirst = 0
	   write(6,*)'Load table of TOF paddle knockouts'
c
cc	   open(unit=1,file='tof_paddle_knockout_g1a.txt',status='old')
cc	   open(unit=1,file='tof_paddle_knockout_g1c.txt',status='old')
cc	   open(unit=1,file='tof_paddle_knockout_g11.txt',status='old')
	   open(unit=1,file='tof_paddle_knockout_g14.txt',status='old')
c
	   do j=1,48
	      do i=1,6
		 toflist(j,i) = -1
	      enddo
	   enddo
 10	   read(1,*,end=50)idnum,effys
	   if(idnum.ge.1.and.idnum.le.48)then
	      do i=1,6
		 if(effys(i).ne.0 .or. effys(i).ne.1)then
		    toflist(idnum,i)= effys(i)
		 else
		    write(6,*)'Table error:',idnum,effys
		    call exit()
		 endif
	      enddo
	   else
	      write(6,*)'Table error:',idnum,effys
	      call exit()
	   endif
	   goto 10
 50	   close(1)
	endif
c
	if(isector.lt.1 .or. isector.gt.6)then
	   write(6,*)'Sector screw up',isector,numtofbar
	   knockout_CLAS_tof = .false.
c	   call exit()
	endif
	if(numtofbar.eq.0)then         !changed 9-4-14;  until 2/7/01, this line said "eq"
	   knockout_CLAS_tof = .true. 
c	   write(6,*)'knockout_CLAS:',isector,numtofbar
	   return
	endif
c
c       We have lost the trigger efficincy file, so...
c
c	knockout = toflist(numtofbar,isector)
	knockout = 1
	if(knockout.eq.1)then
	   knockout_CLAS_tof = .true.
	else
	   knockout_CLAS_tof = .false.
	endif
	write(6,*)'knockout_CLAS:',isector,numtofbar,knockout
	return
	end

c
c       Decide if this track in the drift chambers will be kept or not
c
c       The big table that holds the info about which wires are alive is:
c       LIVEWIRELIST(1...34 layers, 1..192 wires, 1..6 sectors)
c
	function knockout_CLAS_dc(isector,idcwirelist)
	parameter (nwpawc = 20000000)
	character*7 nsname(6)
	logical knockout_CLAS_dc,livewirelist(34,192,6)
	integer idcwirelist(34),iwireshift(34)
	common /pawc/hmemor(nwpawc)
	common /knockout_CLAS_dclocal/livewirelist
	common /monte/ iseed,xflat,xnormal,xbw
	data ifirst/1/
	data nsname/'ns1_g1c','ns2_g1c','ns3_g1c',
     1              'ns4_g1c','ns5_g1c','ns6_g1c'/
c
	write(6,*)'knockout_CLAS_dc entered...'
	knockout_CLAS_dc = .true.
c
	if(ifirst.eq.1)then
	   ifirst = 0
	   write(6,*)'Load table of Drift Chamber dead/live wires...'
c
c          Read the list of live wires for g1c prepared by 
c          Luminita Todor on 6-29-01.
c
	   do isec = 1,6
	      open(unit=1,file=nsname(isec),status='old')
	      jlayer = 0
	      do ilayer = 1,36
		 if(.not.(ilayer.eq.5 .or. ilayer.eq.6))
     1                jlayer = jlayer + 1
		 do iwire = 1,192
		    read(1,*)istatus
		    if(ilayer.eq.5 .or. ilayer.eq.6)then
		       if(istatus.eq.0)then
			  write(6,*)'Deep shit:',ilayer,iwire
			  call exit()
		       end if
		    else
		       if(istatus.ne.0)then
			  livewirelist(jlayer,iwire,isec) = .false.  !bad wire
		       else
			  livewirelist(jlayer,iwire,isec) = .true.   !good wire
		       end if
		    end if
		 enddo
	      enddo
	      close(1)
	   enddo
	endif
c
	if(isector.lt.1 .or. isector.gt.6)then
	   write(6,*)'Sector screw up',isector
	   call exit()
	endif
c
	isum = 0
	isl1 = 0
	isl2 = 0
  	isl3 = 0
 	isl4 = 0
	isl5 = 0
	isl6 = 0
c
c       Set all axial and stereo wire shifts to zero for midplane tracks
c
	do ilayer = 1,34
	   iwireshift(ilayer) = 0
	enddo
	tanval = tan(6.*3.141592/180.)   !tangent of 6 deg stereo angle
	x = 2.*ran(iseed)-1
c
	do ilayer = 1,34
	   idcwiretest = idcwirelist(ilayer) + iwireshift(ilayer)
	   if(idcwiretest.le.0)  idcwiretest = 1
	   if(idcwiretest.ge.143)idcwiretest = 142
c	   write(6,*)ilayer,idcwirelist(ilayer),
c     1                           iwireshift(ilayer),idcwiretest
	   if(idcwiretest .gt. 0)then
	      if(livewirelist(ilayer,idcwiretest,isector))then
		 isum = isum + 1
		 call hf2(6000+isector,float(idcwiretest),
     1                    float(ilayer),1.0)
		 if(ilayer.ge.1  .and. ilayer.le.4)then
		    isl1 = isl1 + 1
		 else if(ilayer.ge.5  .and. ilayer.le.10)then
		    isl2 = isl2 + 1
		 else if(ilayer.ge.11 .and. ilayer.le.16)then
		    isl3 = isl3 + 1
		 else if(ilayer.ge.17 .and. ilayer.le.22)then
		    isl4 = isl4 + 1
		 else if(ilayer.ge.23 .and. ilayer.le.28)then
		    isl5 = isl5 + 1
		 else if(ilayer.ge.29 .and. ilayer.le.34)then
		    isl6 = isl6 + 1
		 else
		    write(6,*)'Could not determine super-layer number!!?'
		    call exit()
		 endif
	      endif
	   endif
	enddo
c
c       Three hits per superlayer are required, according to 
c       Jim Mueller, for our data.
c
	inumsl = 0
	if(isl1.ge.3)inumsl = inumsl+ 1
	if(isl2.ge.3)inumsl = inumsl+ 1
	if(isl3.ge.3)inumsl = inumsl+ 1
	if(isl4.ge.3)inumsl = inumsl+ 1
	if(isl5.ge.3)inumsl = inumsl+ 1
	if(isl6.ge.3)inumsl = inumsl+ 1
c
	if(inumsl .lt. 5)then
c	if(isum.lt.29)then
	   knockout_CLAS_dc = .false.    !track is rejected
c	else
c	   call hf1(6007,float(isum),1.)      !DC track hit multiplicity
	endif
c	call hf1(6008,float(isum),1.)         !DC track hit multiplicity
c
	write(6,*)'knockout_CLAS_dc:',isector,inumsl,isum,knockout_CLAS_dc
	write(6,*)'                :',isl1,isl2,isl3,isl4,isl5,isl6
	return
	end


	subroutine epcms(p1,p2,pcm)
c       calculates total momentum and cm energy of two particles
	REAL P1(20),P2(20),PCM(5)
	AP1S=P1(13)**2+P1(14)**2+P1(15)**2
	AP2S=P2(13)**2+P2(14)**2+P2(15)**2
	APCMS=0.
	DO 100 I=1,3
	   PCM(I)=P1(I+12)+P2(I+12)
	   APCMS=APCMS+PCM(I)**2
 100	CONTINUE
	E      = SQRT(AP1S+P1(4)**2)+SQRT(AP2S+P2(4)**2)
	PCM(4) = SQRT(E*E-APCMS)
	pcm(5) = e
	RETURN
	END

c
c	Output to file in format suitable for txt2part
c
c	Coding for the list of particles generated in the decay chain:
c	All variables in this list are measured in the lab.
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
c	16	TOTAL ENERGY AT DECAY POINT
c	17	X OF DECAY
c	18	Y OF DECAY
c	19	Z OF DECAY
c	20	(CODE)
c
 	subroutine clas_out(egamma,pkgen,plgen)
	logical ifirst
	real pkgen(5),plgen(5)
 	data ifirst/.true./
	common /com_clas_out/ icount
	data icount /0/
c
c	On first call just open the data file
c
	if(ifirst)then
		ifirst = .false.
		write(6,*)'Opening TXT2PART format output file'
		open(unit=2,file='mc_clas.dat',status='UNKNOWN')
	end if
	icount = icount + 1
c
c	decode PLIST 4-vectors into CLAS format
c       (number of particles)
c       (ebeam energy) (photon energy) (time offset=0)
c       (PID 11=kaon 18=lambda) (E) (px) (py) (pz)
c       (vertex location x y z = 0 0 0)
c       (PID 11=kaon 18=lambda) (E) (px) (py) (pz)
c       (vertex location x y z = 0 0 0)
	numpart = 2
	ebeam   = 2.5
	toff    = 0.0
	nkaon   = 11
	nlam    = 18
	vx  = 0.0
	vy  = 0.0
	vz  = 0.0
c
	write(2,10)numpart
 10	format(i1)
	write(2,20)ebeam,egamma,toff
 20	format(f8.6,2f9.6)
	write(2,30)nkaon,pkgen
 30	format(i2,4f10.6)
	write(2,40)vx,vy,vz
 40	format(f8.6,2f9.6)
	write(2,30)nlam,plgen
	write(2,40)vx,vy,vz
	return 
	end 

c
c	Input from to file in format of txt2part
c
c	Coding for the list of particles generated in the decay chain:
c	All variables in this list are measured in the lab.
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
c	16	TOTAL ENERGY AT DECAY POINT
c	17	X OF DECAY
c	18	Y OF DECAY
c	19	Z OF DECAY
c	20	(CODE)
c
	subroutine clas_in(ppbeam,pkgsim,plgsim,r0)
	logical ifirst
	real pkgen(5),plgen(5)
	real pkgsim(5),plgsim(5),r0(3)
 	data ifirst/.true./
	common /com_clas_in/ jcount
	data jcount /0/
c
c	On first call just open the data file
c
	if(ifirst)then
		ifirst = .false.
		write(6,*)'Opening TXT2PART format input file'
		open(unit=2,
c     1       file='/dsk55/mcnabb/acceptance/kaon/kl.1001-1010.txt',   !with decays
c     1  file='/dsk55/mcnabb/acceptance/kaon/kl.3001-3010.no3005.txt', !no decays
c     1  file='/dsk55/mcnabb/acceptance/kaon/kl4.txt', !with decays and long target
     1  file=
     1  '/raid13/kmoriya/1520ascii/mc_generator_CLAS_1520_001.ascii', !with decays and long target
     1	status='OLD')
	end if
	jcount = jcount + 1
c
c	decode PLIST 4-vectors into CLAS format
c       (number of particles)
c       (ebeam energy) (photon energy) (time offset=0)
c       (PID 11=kaon 18=lambda) (E) (px) (py) (pz)
c       (vertex location x y z = 0 0 0)
c       (PID 11=kaon 18=lambda) (E) (px) (py) (pz)
c       (vertex location x y z = 0 0 0)
c
	read(2,*,err=10)numpart  !should be 2
	read(2,*,err=10)ebeam,egamma,toff
	read(2,*,err=10)nkaon,pkgen
	read(2,*,err=10)vxk,vyk,vzk
	read(2,*,err=10)nlam,plgen
	read(2,*,err=10)vxl,vyl,vzl
	write(6,*)'>>',jcount
	write(6,*)'  ',numpart
	write(6,*)'  ',ebeam,egamma,toff
	write(6,*)'  ',nkaon,pkgen
	write(6,*)'  ',vxk,vyk,vzk
	write(6,*)'  ',nlam,plgen
	write(6,*)'  ',vxl,vyl,vzl
c
	pkgsim(1) = pkgen(2)
	pkgsim(2) = pkgen(3)
	pkgsim(3) = pkgen(4)
	plgsim(1) = plgen(2)
	plgsim(2) = plgen(3)
	plgsim(3) = plgen(4)
	r0(1)     = vxk
	r0(2)     = vyk
	r0(3)     = vzk
	ppbeam    = egamma
	return 
c
 10	pkgsim(1) = 0
	pkgsim(2) = 0
	pkgsim(3) = 0
	plgsim(1) = 0
	plgsim(2) = 0
	plgsim(3) = 0
	r0(1)     = 0
	r0(2)     = 0
	r0(3)     = 0
	ppbeam    = 1.0
	return 
	end 

	subroutine pseudo_bellis(icharge,pmomgev,theta,phisector,isector,acc)

c-      Inputs:
c-          icharge - INTEGER 0=electron, 1=positive hadron, -1=negative hadron,
c-          pmomgev - REAL particle momentum GeV/c
c-          theta   - REAL particle polar angle in degrees (0 to 180.)
c-          phi     - REAL particle azimuthal angle in degrees in the sector -30. to 30.
c-          isector - INTEGER sector number.    
c-
c-      Output:
c-          acc     - REAL acceptance 0. to 1.
c
c       For g11 data set.
c
	data dera/0.017453293/
	acc = 0.
	costhetarad = cos(theta*dera)
	phirad   = phisector*dera
	c3 = 1.87603
	c2 = -4.54334
	c1 = 0.022356
	c0 = 0.993861
	costhetacut = c3*phirad**6 + c2*phirad**4 + c1*phirad**2 + c0
	phicut = 0.4
	forwardcosthetacut = 0.985
	if((costhetarad .lt. costhetacut)        .and.
     1     (abs(phirad) .lt. phicut)             .and.
     1     (costhetarad .lt. forwardcosthetacut) .and.
     1     (costhetarad .gt. -0.5)               .and.
     1     .not.(isector.eq.1 .and. costhetarad.lt.-0.25) .and.
     1     .not.(isector.eq.3 .and. costhetarad.lt. 0.15) .and.
     1     .not.(isector.eq.5 .and. costhetarad.lt. 0.00))then
	   acc = 1.0
	endif
c	write(6,*)icharge,pmomgev,theta,phisector,isector,acc
	return
	x = pmomgev !fool the gfortran compiler
	x = icharge
	end


c****************************************************************************
c       Convert a particle 4-vector to angles in standard coordinates
c       Z - beam axis, when given a vector in lab coordinates
c       Y - up in lab frame
c       X - to left at seen along beam
c
c       R. A. Sch. CMU
c       5-1-2015 - first version.  I have no idea why I didn't write this
c                  code a very long time ago.
c
c       Input:
c       particle - input giving px,py,pz,mass,energy
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
c       Compute the boost beta that gets is to the rest frame of whatever
c       particle 4-vector we are given
c
c       R.A. Sch,  CMU,  5-2015
c
c       PVECTOR (input) is defined as PX, PY, PZ, MASS, ENERGY
c       BETA (output) is the boost 3-vector
c 
        subroutine getbeta(pvector,beta)
c
	real pvector(5),beta(3)
c
	e0=etot(pvector)
	do j=1,3
	   beta(j) = pvector(j)/e0
	end do
	return
	end


c*******************************************************************************
c       Compute the boost beta that gets is to the rest frame of whatever
c       particle 4-vector we are given
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
c       P3 is the "produced" particle
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
	
