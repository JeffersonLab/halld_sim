      PROGRAM MC_GEN
c*********************************************************************************
c	MC_GEN
c	Monte Carlo generation of events with sequential N-body decays.
c
c	The topology of the events is this:
c	1 + 2(target) ---> 3 + 4
c			   3 ---> 5 + 6 + ... + n
c			       4 --->  n+1 + ...polv
c				etc...
c	The masses, etc, for particles are read in from a carefully-structured
c       file whose generic name is MC_GEN.DEF.  It is the "definition" file for the
c       process we will be simulating, a.k.a. a "configuration" file.
c
c	The top-level program MC_GEN is responsible for the generation of
c       events under the control of the input-file.  There are several subtroutines,
c       some of which are attached directly to this file and others that are 
c       separate.  A (very large) subroutine called MC_ANALYZE is (optionally)
c       called for each generated event to allow for fast on-the-spot analysis and
c       histogramming of the generated events.  
c
c       This program writes raw event output to external files in ASCII format  
c       for input to the JLab/GlueX standard GEANT package.
c       The optional analysis subroutine MC_ANALYZE creates and stores histograms 
c       in the HBOOK format for displaying with PAW or ROOT.  
c
c	R. A. Schumacher, CMU, October 1991.
c	Original:
c		10-28-91 - Started from program DATGEN to develop this program
c	Modified:
c		11-18-91 - for at-rest decay the spin vector can point
c			   in any direction
c		11-25-92 - add three-body decay modes
c	  	 12-4-92 - add resonance production (incoherent only)
c		 7-19-95 - big upgrades in three-body decay modes
c			 - develop weak decay code to see angular correlations
c		 8-10-95 - reorganized the data input at the beginning 
c			   of the program; all data in the input file.
c		  3-1-96 - added an energy loss calculation in the target 
c		 3-13-96 - add distributed primary vertex
c	 	 8-13-97 - add input parameter for CLAS spectrometer B field
c                8-18-97 - UNIX version
c                9-29-97 - Include the Geant random number generator rather
c                          than "ran".
c                        - Pick up the input file name as the "first"
c                          command line argument.
c                 4-1-99 - various fixes related to mixed-mode arithmetic
c                          such that the program runs on Linux systems 
c                   3-99 - improve hyperon photoproduction model
c                   6-00 - add stuff for polarization of Hyperons
c                 1-5-01 - apply the instrumental smearing of each track
c                          here in mc_generator, at the time the event is
c                          generated.
c                 2-9-01 - expanded PLIST such that smeared tracks have their
c                          mass stored in slot 20, for smoother analysis flow
c                6-28-01 - expanded to allow pion cross sections to be used
c                          as input, if needed.
c                 7-3-01 - recover via backups and hardcopies from accidental
c                          deletion of this file.
c                 4-2-02 - fix polarization info calculations; now get Cx and Cz
c                8-20-02 - introduce the ability to read TXT2BOS files as input
c                 6-3-05 - upgrade the cross section modeling to use fit results to 
c                          the measured values
c               11-20-05 - add code to allow text output for TXT2PART input.
c               12-01-05 - increase maxdecays from 10 to 12 (for the Lambda(1520))
c                1-18-06 - correct Z vertex time offset calculation to match GSIM.
c                8-10-06 - comment out all need to use 'mc_default.f';
c                          add bremsstrahlung shape for beam momentum generation
c                8-14-06 - increased maximum particles in table from 40 to 55
c                8-21-06 - make TXT2BOS output record formation point of tracks 
c                        - add a line to input definition file to give file name prefix
c                          for TXT2BOS output files
c                8-30-06 - Add a rule for generating events that follows "t-channel"
c                          dominance
c                2-02-07 - added logic to prevent a technical snafu:
c                          allow decay via phase space if and only if the production 
c                          angle (of the meson) was not pre-selected
c                3-30-07 - fix bug in vertex time computation in txt2bos_output
c                1-28-08 - upgrade code to compile under gfortran; many modules affected...
c                5-05-08 - allow for 4- and 5-body phase space by incorporating genbod.f
c                9- 3-08 - let the TXT2BOS file have a 4-digit extension
c                11-5-08 - let the lineshape of a produced particle be a 
c                          relativistic Breit-Wigner a la J. D. Jackson.
c               11-13-08 - found and fixed an old bug whereby txt2bos output sometimes
c                          ended after writing fewer than the specified number of events
c               11-18-08 - introduce a relativistic Breit-Wigner lineshape option
c                          this involves some programming chicanery with the sign of the 
c                          "width" of the particle to pick either a pure Lorentzian
c                          or a full relativistic BW
c                1-23-09 - Add a Sigma(1385) L=1 selection.  Some day these L -value
c                          selections should be coded into the .def files....
c                2-13-09 - Add "Particle 70" which is generated as a flat mass distribution
c                          between limits set by the nominal mass (low edge) and width 
c                          (high edge) parameters
c                3-18-09 - Added Lambda(1405) lineshape selection: when ctau=-1 read
c                          and use a table of values indexed into photon energy bins and
c                          decay cases Sigma+pi-, Sigma0pi0, Sigma-pi+.  The logic of the
c                          code in built in near the place where we also can select the 
c                          relativistic BW.  The new subroutine is 'lineshapepick'
c               10-30-09 - Added a 2-body phase-space weighting of the decay mass
c                          distribution for the case of eta_prime --> rho0 gamma.  The weighting
c                          is in proportion to the available decay momentum.  It is implemented
c                          using ndynamics=9
c               11-02-09 - Expanded the calculation of 2-body phase space momentum weighting
c                          to the cases of production of K+ Y* and K* Y.  This involves 
c                          weighting by the ratio of output c.m. momentum to input c.m. momentum.
c                          It is turned on using ndynamics = 10 or 11. 
c                3-26-10 - lengthen name of output files.  Now it can be up to 150 characters (!)
c                3-18-11 - insert an extra call to RAN to tweak the sequence of random numbers
c                        - fix bug in Lambda(1405) special case: ctau = -1 flag that was used for
c                          making lineshape selection correct also caused 1405 to have ctau = 1cm!                    
c                4-15-11 - added an IPRULE=11 for the x(1280) production case
c                8-09-13 - added 4 and 5 -body phase space decays for which P1 is forced to 
c                          have <PFERMI recoil momentum
c               11-05-13 - tweak the algorithm for having one "broad" state decay into another
c                          "broad" state.  Relates to the study of L(1520) --> S(1385) pi
c               12-03-13 - for IPRULE=7, change the recoil particle from being always a ground 
c                          state Lambda to being whatever it really is in the decay list.
c                4-07-14 - add routine REGGEPICK_RHO to mimic rho photoproduction a la 
c                          Donnachie and Landshof;  goes with IPRULE=13
c                9-04-15 - start the recovery process from a disk crash and loss of changes since
c                          the end of May of this year.
c                9-05-14 - introduce Hulthen momentum distribution for quasi-free breakup of deuteron
c               12-18-14 - allow the initial state to include a neutron in Fermi motion.  This 
c                          is the first time the code allows a non-stationary initial target.
c                          The subroutine getdeuteronmomentum(p) was written
c                5-15-14 - Make several 4-vectors of length 5 to put total energy in element 5
c                9-30-15 - modify write-out list to write out photons and pi0 4-vectors (for BoGa runs)
c                        - some code tweaks to satisfy the gfortran compiler on the Mac 
c         	 6-26-17 - tweaks to allow Lambda anti-Lambda events to be generated
c                7-27-17 - big move to compute "Regge" decay in the particle-decay part of the code
c                          rather than in the particle-generation part of the code
c               11-28-17 - Read in the "shape" of the raw beam photon spectrum from a file
c                2-08-18 - Add bits to allow Lambda p --> Lambda p scattering
c                5-25-18 - Updates to Lambda-anti-Lambda code for event generation
c                3-06-19 - some adjustments to generate Delta-anti-Delta, increase MAXNUMLIS to 60
c                6-16-19 - Code to read in a GlueX beam flux file from outside
c       Fission of MC_GENERATOR into two programs;  this one is MC_GEN which generates events but
c                          only optionally analyzes events
c                6-26-19 - Develop MC_GEN, with the goal of making a stand-alone generator
c                          with no event analysis included.        
c                7-12-19 - reexport the program to the CMU central computers
c                          to run larger-scale simulations on computer clusters
c                7-19-19 - Introduce a decay mode with 3-body decay followd by Regge slope selection
c                8-26-19 - Small bug fix in writing out ASCII data for Mechanism G2 - order of particles
c               10-29-19 - Tiny bug fix in Reggepick to trap rare machine precision issue    
c               11-19-19 - Code to allow simulation of the "double Regge" process vi a "Mechanism 5"
c               02-26-20 - Small code upgrades to handle "double-Regge" better
c               03-23-20 - create a 3-body mechanism with and overall b1 t-slope to the produced pair(GeV^-2) and
c                          and an invariant P23 mass cluster pairing slope b2 (GeV^-1)
c                          We are calling this "Mechanism 6".     
c     
        parameter (version=11.60) !Version of the program... keep this updated
c      
	external ran  !forces the use of "my" version of random number generator (see RNDNOR.F)
	parameter (maxnumpar=65) !Maximum number of particles in the table
	parameter (maxnumlis=60) !Maximum number of particles in the decay list
	parameter (maxdecays=12) !Maximum number of decay modes per particle
	character*8 names(maxnumpar)
	character   dummy*50,letter*1
	character*80  title 
        character*200 fluxfilename
	integer lundid(maxnumpar),numdecays(maxnumpar)
	real rm(maxnumpar),width(maxnumpar),ct(maxnumpar)
	real zpart(maxnumpar)
	integer idecay(maxnumpar,maxdecays,7)
	real branch(maxnumpar,maxdecays,3)
	real plist(24,maxnumlis)
	real pp(5,maxnumpar)
	real ptot(5),pv(5),pvsmear(5),p1(5),p2(5),p3(5),p4(5),p5(5) !4-vectors in KINPACK
	real pc(5),ptotrest(5),p23(5),pvtemp(5),ptarget(5)
	real r0(3),r1(3)
	real polv(3),pol1(3),pol2(3),pol3(3),pol4(3),pol5(3)
	logical idolund,idogluex,idocuts,ihackflag,iwidthflag,imess
	data rmpipm    /.13957/     !pion mass
	data rmkaon    /.493677/    !K+ mass
        data rmrho     /.770/       !mean rho mass 
	data rmkstar   /.892/       !mean K* mass
	data rmetaprime/.95778/     !etaprime mass
	data rmproton  /.938272081/ !proton mass 
	data rmlambda  /1.115683/   !lambda mass 
	data rmsigma   /1.198/      !highest sigma mass 
	data rmdeuteron/1.875/      !deuteron mass 
	data rmx1280   /1.282/      !x(1280) mass 
        data idcountmax/0/
c        data iseed /123456/
	common /decaytable/idecay,lundid,rm
	common /monte/ iseed,xflat,xnormal,xbw
	common /stuff/poldir(4),xx,icount,ikeep,ihackflag,iattempt
        common /fermi/prec,qfmass,pfermi,ndynamics
	common /beam/pbeam,ppbeam,rmbeam,rmtarg,jcode,xsigma,ysigma,ncode
	common /analysis/dpoverp_in,dpoverp_out,ptotrest
        common /bremselectname/fluxfilename,lfluxfilename
	common /target/icode,size1,size2,size3,size4,maxtstep
	common /misc/cthcm_gen,itoken          !for tests and hacks
	common /temp/eloss_clas,eloss_prop,pcmag_clas,pvmag_prop
c
c       Iargc: Unix function to return the number of command line args.
c       Getarg: Unix function to get a specific argument.
c       Lenocc: cernlib length of a string.
c
	character*256  prg_name, histnm, filenm, prefix
	common /run_control/n_arg,prg_name,filenm,histnm,prefix,
     1                      irestart,lfl,lprefix,maxbosoutput
c
	write(6,*)'****************************************************'
	write(6,*)'Starting Carnegie Mellon Monte Carlo Event Generator'
        write(6,*)'  Program MC_GEN, developed from MC_GENERATOR'
        write(6,*)'  R. A. Sch.,  schumacher@cmu.edu'
        write(6,1)'  Version ',version
	write(6,*)'****************************************************'
 1      format(1x,a,f6.2)
c     
c       Start out with some unix hackery to get command line
c       Argument "0" is the program name
c
	N_arg = iargc( )
        call getarg(0,Prg_name)
	LFL = LENOCC(Prg_name)
	if ( N_arg .lt. 1 ) then
           write(6,2)Prg_name(1:LFL)
 2         format(' Usage: ',A,' input_file1 ')
	   Stop
	endif
	call getarg(1,filenm)
	lfl = lenocc(filenm)
c
	do i=1,maxnumpar
	   do j=1,maxdecays
	      do k=1,7
		 idecay(i,j,k)=-1
	      enddo
	   enddo
	enddo
	imess = .false.
c     
c	Get the configuartion information for what we are going to be calculating.
c
	open(unit=1,file=filenm(1:LFL),status='old') !,readonly)
c
c       First skip all the comments;  they are delimited by a string of "=" signs
c
 3      read(1,'(a)')dummy
	if(dummy(1:5).ne.'=====')then
c		write(6,'(1x,a)')dummy
		goto 3
	end if
c        
	read(1,'(a)')title
 4      continue
	read(1,*,err=4) nevent
c	write(6,*)nevent
 5      continue
	read(1,*,err=5) ngevent
c	write(6,*)ngevent
 6      continue
       	read(1,*,err=6) pbeammin
c	write(6,*)pbeammin
 7      continue
	read(1,*,err=7) pbeammax
c	write(6,*)pbeammax
 8      continue
	read(1,*,err=8) iprule  !rule for generating photons
c	write(6,*)iprule
 9      continue
	read(1,'(a)',err=9)fluxfilename   !used if we are reading an external file for photon beam shape 
	lfluxfilename = iscan(fluxfilename,' ')-1  !KERNLIB routine
 10     continue
	read(1,*,err=10) nprint
c	write(6,*)nprint
 11     continue
	read(1,*,err=11)item
c	write(6,*)item
	if(item.eq.0)then
	   ihackflag = .false.
	elseif(item.eq.1)then
	   ihackflag = .true.
	end if    
 12     continue
	read(1,*,err=12)item
c	write(6,*)item
	if(item.eq.0)then
	   idogluex  = .false.
	endif
	if(item.eq.1)then
	   idogluex  = .true.
	endif
 13     continue
	read(1,*,err=13)item
c	write(6,*)item
	if(item.eq.0)then
	   idocuts  = .false.
	endif
	if(item.eq.1)then
	   idocuts  = .true.
	endif
 14     continue
	read(1,'(a)',err=14)prefix !prefix for output files
	lprefix = iscan(prefix,' ')-1  !KERNLIB routine
 15     continue
	read(1,*,err=15)maxbosoutput !max number of events in output file
c       write(6,*)maxbosoutput
 16     continue
	read(1,*,err=16) dpoverp_out
c	write(6,*)dpoverp_out
 17     continue
	read(1,*,err=17) pfermi
c       write(6,*)pfermi
 18     continue
	read(1,*,err=18) icode
c	write(6,*)icode
 19     continue
	read(1,*,err=19) size1,size2,size3,size4
c	write(6,*)size1,size2,size3,size4
 20     continue
	read(1,*,err=20) jcode,xsigma,ysigma
c	write(6,*)jcode,xsigma,ysigma
 21     continue
        read(1,*,err=21) irestart
c        write(6,*)irestart
 22     continue
	read(1,*,err=22) reserved		
c	write(6,*)reserved
c        write(6,*)'Read in the initial parameters...'
        
c	Read over a fixed number of dummy lines that are comments in the input file
c
	do i=1,13
           read(1,'(a)')dummy
c          write(6,*)dummy
	end do
c
c	Loop over the particle definition table 
c
c	The syntax includes that the first line is the beam particle,
c	the second line is the target particle, and the third line is
c	a "pseudo" particle that defines how the beam-target system will
c	decay.   The rest of the table defines all other particles needed
c	to define the decay chain for the reaction.
c       
c	IDECAY(A,B,C)
c		A - position index in list of particles
c		B - decay mode index
c		C = 1  number of bodies in decay mode
c		    2, 3, 4, 5, 6  Lund ID''s of up to five decay products
c		    7  code to specify how to treat decay dynamics: 
c						0= pure phase space
c						1= hyperon decay 
c						2= nuclear Fermi gas
c						3= Lambda photoprod.
c                                               4= Sigma0 photoprod.
c                                               et cetera...
c     BRANCH(A,B,C)
c		A - position index in list of particles
c		B - decay mode index
c		C = 1  decay branching fraction for this decay mode
c		    2  first  parameter specifying the decay dynamics (e.g. b-slope     in Regge model)
c		    3  second parameter specifying the decay dynamics (e.g. W_threshold in Regge model)
c
	i = 0		!Counts particles in the table, including Beam and Targ.
 23     continue
	i = i + 1
 25	continue
	read(1,*,err=25)idummy,
     1    names(i),rm(i),width(i),zpart(i),ct(i),rlundid,numdecays(i)
	lundid(i) = rlundid    !a hack to avoid changing all the input files for Linux
c	write(6,*)idummy,i,' ',names(i),int(rlundid),numdecays(i)
	if(idummy.eq.-1)goto 40
c       width(i) = width(i) / 2.354 /1000.	!Width sigma in GeV
	width(i) = width(i) / 2.0   /1000.      !Gamma/2     in GeV
	if(idummy.ne.i)write(6,*)'Definition table has messed up syntax',idummy,i
	if(i.le.2)goto 23			!beam and target can't decay
	do 35 j=1,abs(numdecays(i)) 	!decay modes of particle (may be zero).
 26        continue
           read(1,*,err=26)nbody, ndynamics, dyncode1, dyncode2
           idecay(i,j,1) = nbody !flags N=2 3, 4 or 5 body decay
           idecay(i,j,7) = ndynamics !specify decay type: phase space, etc.
           branch(i,j,2) = dyncode1
           branch(i,j,3) = dyncode2
           if(nbody.eq.2)then
 27           continue
              read(1,*,err=27)id1,id2,decayfrac
              idecay(i,j,2) = id1 !Lund ID of first decay particle
              idecay(i,j,3) = id2 !Lund ID of 2nd   decay particle
              branch(i,j,1) = decayfrac !sum over J's should be 1.0
c	      write(6,*)'>2-body decay mode',id1,id2,decayfrac
           else if(nbody.eq.3)then
 28           continue
              read(1,*,err=28)id1,id2,id3,decayfrac
              idecay(i,j,2) = id1 !Lund ID of first decay particle
              idecay(i,j,3) = id2 !Lund ID of 2nd   decay particle
              idecay(i,j,4) = id3 !Lund ID of 3rd   decay particle
              branch(i,j,1) = decayfrac !sum over J's should be 1.0
c             write(6,*)'>3-body decay mode',id1,id2,id3,decayfrac
           else if(nbody.eq.4)then
 29           continue
              read(1,*,err=29)id1,id2,id3,id4,decayfrac
              idecay(i,j,2) = id1 !Lund ID of first decay particle
              idecay(i,j,3) = id2 !Lund ID of 2nd   decay particle
              idecay(i,j,4) = id3 !Lund ID of 3rd   decay particle
              idecay(i,j,5) = id4 !Lund ID of 3rd   decay particle
              branch(i,j,1) = decayfrac !sum over J's should be 1.0
c	      write(6,*)' 4-body decay mode',id1,id2,id3,id4,decayfrac
           else if(nbody.eq.5)then
 30           continue
              read(1,*,err=30)id1,id2,id3,id4,id5,decayfrac
              idecay(i,j,2) = id1 !Lund ID of first decay particle
              idecay(i,j,3) = id2 !Lund ID of 2nd   decay particle
              idecay(i,j,4) = id3 !Lund ID of 3rd   decay particle
              idecay(i,j,5) = id4 !Lund ID of 4th   decay particle
              idecay(i,j,6) = id5 !Lund ID of 5th   decay particle
              branch(i,j,1) = decayfrac !sum over J's should be 1.0
c              write(6,*)' 5-body decay mode',id1,id2,id3,id4,id5,decayfrac
           else
              write(6,*)'N>5 body decay not implemented',nbody
              call exit()
           end if
 35	continue
	goto 23
c
c       All the configuration info in now in hand.  Start feeding it out again...
c
 40	close(1)

	write(6,'(1x,/,1x,79(''*''))')
	write(6,*)title
	write(6,*)'Configuration file :              ',filenm(1:LFL)
	write(6,*)'Generating',nevent,' events'
	write(6,*)'or        ',ngevent,' events that pass all cuts'
        write(6,52)'Minimum beam momentum              = ',pbeammin,' GeV/c'
        write(6,52)'Maximum beam momentum              = ',pbeammax,' GeV/c'
 52     format(1x,a,f8.2,a)
	write(6,*)'Beam generated using code',iprule
        write(6,*)'Beam flux distribution :          ',
     1      fluxfilename(1:lfluxfilename)
	write(6,*)'The incident        particle is a ',names(1)
	write(6,*)'The target          particle is a ',names(2)
	write(6,*)'The beam/target system has',numdecays(3),
     1                                              ' decay modes.'
	do 60 idecaynum = 1,numdecays(3)
	   write(6,*)'Decay mode',idecaynum,':'
	   write(6,55)
     1  	idecay(3,idecaynum,1),idecay(3,idecaynum,2),
     1  	idecay(3,idecaynum,3),idecay(3,idecaynum,4),
     1  	idecay(3,idecaynum,5),idecay(3,idecaynum,6),
     1          idecay(3,idecaynum,7),
     1  	branch(3,idecaynum,1)
 55	   format(1x,'Number of bodies, Lund IDs, Branching_Fraction',/,
     1            1x, 7i5, f6.3)
 60	continue
        if(irestart.eq.0)then
           write(6,*)
     1     'The random number generator is restarting from its original seed'
        elseif(irestart.eq.1)then
           write(6,*)
     1     'The random number generator is starting from the previous sequence'
        endif
c
c	Coding for the list of particles generated in the decay chain:
c	All variables in this list are given in the lab frame.
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
c	16	MASS (same as element (4))
c	17	TOTAL ENERGY AT DECAY POINT
c	18	X OF DECAY
c	19	Y OF DECAY
c	20	Z OF DECAY
c	21	CTAU for THIS PARTICLE
c       22      Charge of the particle
c       23      spare
c       24      spare 
c
c	KINPACK parameters used in kinematics calculations
c	PP(1,i) = px, PP(2,i) = py, PP(3,i) = pz, PP(4,i) = mass, PP(5,i) = energy
c
c	LUND parameters (probably not used any more in this code)
c	p(i,1) = px, p(i,2) = py, p(i,3) = pz, p(i,4) = total energy, 
c	p(i,5) = mass, k(i,1) = ??, k(i,2) = ID code
c
c	The beam axis is hereby declared to be the Z axis.
c
c	Get the cm energy and total momentum of the system
c	assuming the average momentum.
c
	pbeam   = 0.5*(pbeammin+pbeammax)
	pp(1,1) = 0	!Beam particle
	pp(2,1) = 0
	pp(3,1) = pbeam
	pp(4,1) = rm(1)
	pp(1,2) = 0	!Target
	pp(2,2) = 0
	pp(3,2) = 0
	pp(4,2) = rm(2)
	call epcm(pp(1,1),pp(1,2),ptot)
	wvalue = ptot(4)
	etotal = sqrt(ptot(4)*ptot(4) + pbeam*pbeam)
	tbeam  = etotal - rm(1) - rm(2)  !lab kinetic energy
	pcmvaluein = (wvalue**2.-(rm(1)**2.+rm(2)**2.))**2.-
     1	                         (2.*rm(1)*rm(2))**2.
	if(pcmvaluein.gt.0.0)then
	   pcmvaluein = sqrt(pcmvaluein)/(2.*wvalue) !c.m. initial momentum
	else
	   write(6,*)'Yikes!',wvalue,rm(1),rm(2)
	   call exit
	endif
	rmbeam = rm(1)		!load common block for later
	rmtarg = rm(2)		!load common block for later
        if(rmtarg.ge.0.939 .and. rmtarg.le.0.940)then
           ncode = 1   !signals a neutron target in common block
        else
           ncode = 0   !signals NOT a neutron target
        endif
	write(6,75)wvalue,pcmvaluein,etotal,tbeam,dpoverp_out*100,
     1   pfermi,ihackflag,idogluex,idocuts
75	format(1x,/,1x,'For the average beam momentum we have:',/,
     1        1x,'Invariant rest mass of system      = ',f8.3,' GeV/c2',
     1         /,
     1        1x,'C.M. momentum in initial state     = ',f8.3,' GeV/c',/,
     1        1x,'Total energy in lab frame          = ',f8.3,' GeV',/,
     1        1x,'Typical kinetic energy of beam     = ',f8.3,' GeV',/,
     1        1x,'Momentum spread of final particles = ',f8.3,' %',/,
     1        1x,'Fermi momentum for nuclear recoil  = ',f8.3,' GeV/c',/,
     1        1x,'Enable various coding hacks ON/OFF(1/0) = ',7x,l1,/,
     1        1x,'Writing event file to disk  ON/OFF(1/0) = ',7x,l1,/,
     1        1x,'Calling the event analyzer  ON/OFF(1/0) = ',7x,l1)

	if(icode.eq.1)then
           write(6,77)size1,size2,size3,size4
 77        format(3x,'Rectangular geometry with sizes',/,
     1   3x,'X = ',f8.3,' cm',/,
     1   3x,'Y = ',f8.3,' cm',/,
     1   3x,'Zmin = ',f8.3,' cm',/,
     1   3x,'Zmax = ',f8.3,' cm')
        else if(icode.eq.2)then
           write(6,78)size1,size2,size3,size4
 78        format(1x,'Cylindrical geometry with full widths',/,
     1   3x,'D (diameter 1) = ',f8.3,' cm',/,
     1   3x,'D (diameter 2) = ',f8.3,' cm',/,
     1   3x,'Zmin = ',f8.3,' cm',/,
     1   3x,'Zmax = ',f8.3,' cm')
        else
           write(6,*)'No target specified'
        end if
	if(jcode.eq.1)then
           write(6,79)xsigma,ysigma
 79        format(1x,'Primary vertex uniform in Z',/,
     1               3x,'X_sigma         = ',f8.3,' cm',/,
     1               3x,'Y_sigma         = ',f8.3,' cm')
        elseif(jcode.eq.0)then
           write(6,*)'Point-like primary vertex at (0,0,0) used'
	end if
        if(ncode.eq.0)then
c           write(6,*)'Non-neutron target particle'
        elseif(ncode.eq.1)then
           write(6,*)'Neutron target particle'
        else
           write(6,*)'TARGET?'
           call exit()
        endif
	write(6,'(1x,79(''*''))')
	write(6,*)' '
c     
c       Generate a filename for the hbook output file:
c
	lfl = lenocc(filenm)
	do j = 1,lfl
	   if(filenm(j:j).eq.'.') goto 80
	end do
 80     do k = j,1,-1
	   if(filenm(k:k).eq.'/') goto 81
        enddo
 81     k = k + 1
        n = j - k + 1
        histnm(1:n) = filenm(k:j)
	histnm = histnm(1:n)//'hbk'
        lfl = n + 3
        write(6,*)'Histogram file name will be:',histnm(1:lfl),j,k,lfl
c
c       Initialize the random number generator:
c       As of Jan 2010, reuse the external functions 
c       provided by RNDNOR.F, which include a mimic of the
c       Fortran instrinsic function "ran".  
c
	if(irestart.eq.0)then
	   write(6,*)'Repeating standard random number sequence:'
	   call grndmq(iseed1,iseed2,iseq,' ')
	else
	   open(unit=20,file='randomseed.num',status='unknown')
	   read(20,*,err=82)iseq,iseed1,iseed2
	   write(6,*)'Continuing previous random number sequence:'
	   write(6,*)'Sequence:',iseq
	   write(6,*)'Seed1   :',iseed1
	   write(6,*)'Seed2   :',iseed2
	   call grndmq(iseed1,iseed2,iseq,'S')
	   goto 84
 82        call grndmq(iseed1,iseed2,iseq,' ')
 84	   close(20)
	endif
c
c       ...and the other one...
c       As of June 2005, the intrinsic RAN function does not care 
c       about the value of ISEED.  If we want a different sequence 
c       we can only shift it by hand by some number of steps. Use 123456
c
c       As of January 2008, the RAN function doesn''t work with gfortran 
c       unless the value of iseed is zero.
c
c       As of January 2010, we use the external function RAN that is part
c       of the RNDNOR.F package that we got from... somewhere.  This allows
c       use to once again restart the random number sequence for a run where
c       we left off at the of the previous run.
c
	iseed = 0
	do i=1,4
	   x=ran(iseed)
	   write(6,*)'  Random Number Test: ',x
	enddo
c        
        numwrite = 0
	ievent   = 0
        iattempt = 0
c
c
c	TOP OF LOOP for letting the system decay
c
100	continue
c	write(6,*)' '
c	write(6,*)'> ************************************************'
c	write(6,*)'> New decay chain'
c	do iii=1,maxnumpar
c	   do jjj = 1,maxdecays
c	      write(6,999)iii,jjj,(idecay(iii,jjj,kkk),kkk=1,7)
c 999	      format(1x,9i9)
c	   enddo
c	enddo
c
c       Run control
c
c	miterate = 0
	numwrite   = numwrite + 1
	ievent     = ievent   + 1
	iattempt   = iattempt + 1
        imechanism = 999         ! used later for (optionally) writing out GlueX-style data (obsolete)
c	write(6,*)'> Generating event',ievent,ikeep,iattempt
 	test = mod(float(ievent),5000.)*1000.
	if(test.eq.0.0)then
	   write(6,102)ievent,ikeep,nevent,ngevent,iattempt
 102	   format(1x,'Number of events generated so far:',2i10,3i12)
	end if
 105	if(ievent.gt.nevent)goto 500
	if(ikeep.gt.ngevent)goto 500
c
c	Initialize the particle list for this event
c
	do 110 mark = 2,maxnumlis
        do 110 i = 1,24
	   plist(i,mark) = 0.0
 110	continue
c
c       If the incident particles have momentum spread, pick the momentum now.
c       IPRULE = 0 means pick beam momentum at random and 'flat'
c       IPRULE = 1 means generate beam momentum a la bremsstrahlung
c       IPRULE = 2 means generate beam momentum a la bremsstrahlung with ad hoc sculpting
c       IPRULE = 3 means select beam energy according to an external beam-energy flux file
c
 	if(iprule.eq.0)then
	   pdelta = pbeammax-pbeammin
	   if(pdelta.gt.0.0)then
	      ppbeam = pbeammin + ran(iseed)*pdelta
	   else
	      ppbeam = pbeammin
	   endif
	   iflat = 0
        else if(iprule.eq.1)then    !Bremsstrahlung distribution
	   call bremsstrahlung(pbeammin,pbeammax,ppbeam)
	   iflat = 1
c           write(6,*)'IPRULE',iprule,ppbeam, pbeammin,pbeammax
        else if(iprule.eq.2)then !Bremsstrahlung distribution with ad hoc hacking
           ibremcount = 0
 144       call bremsstrahlung(pbeammin,pbeammax,ppbeam)
	   call bremsculpt(ppbeam,ibremkeep)
           if(ibremkeep.eq.0)then
              ibremcount = ibremcount + 1
c              write(6,*)'Reject photon',ibremcount,ppbeam
              goto 144
           endif
	   iflat = 1
        else if(iprule.eq.3)then !random selection using external file of actual beam energy profile
c     
c          We still use a locally-determined minimum and maximum beam energy to maintain flexibility
c
           pdelta = pbeammax-pbeammin
 149       if(pdelta.gt.0.0)then
              ppbeam = pbeammin + ran(iseed)*pdelta
           else
              ppbeam = pbeammin
           endif
c           fluxfilename = "photon_flux_sp17.ascii"  !hardcode hack for preliminary version of the program
c           lfluxfilename = 23
	   call bremselect(ppbeam,ibremkeep) !we use a given distribution in this case
           if(ibremkeep.eq.0)then
              ibremcount = ibremcount + 1
c              write(6,*)'Reject photon',ibremcount,ppbeam
              goto 149
           endif
c          write(6,*)'>>> Photon energy generated : ',ppbeam
c           ppbeam = 9.0  !! for tests...
	   iflat = 1
c
c        Some legacy code that may still be useful
c           
c        else if(iprule.eq.3)then !get input from a txt2part GSIM file
c	   call txt2bos_input(plist,imark)
c	   iflat    = 0
c	   goto 150  !skip the whole event generation code and go to analyzer
c
c           Some comments from older code with some possibly-useful references...
c
c           apomeron = 6.0       ! See Donnachie and Landshoff, hep-ph/9912312 (12 Dec 1999)
c           aregge   = 15.9      ! these values are for protons, not deuterons.
c           alpha0pomeron = 1.08
c           alpha0regge   = 0.55
c           betapomeron   = 0.25
c           betaregge     = 0.93
c	    rmrecoil = rmproton
c	    rmtarget = rmproton
c	    bslope   = 6.6  !from Ballam et al.  Phys Rev D 5, 545 (1972).
c 	    call exponentialpick(ppbeam,pbeammax,bslope,sminf,rmtarget,rmmeson,
c     1                     rmrecoil,cosinethetacm,iregge)
c
        else
           write(6,*)
     1           'Unknown rule requested for selecting beam momentum'
           call exit()
        endif
c
c       Done selecting the beam energy.
c       Next we actually generate the particles in the event.
c
        ww   = sqrt(2.*ppbeam*rmproton + rmproton**2.)
        cosinethetacm = 2.*ran(iseed)-1 !may overwrite this later...
	pp(1,1) = 0		!Beam particle
	pp(2,1) = 0
	pp(3,1) = ppbeam
	pp(4,1) = rm(1)
        pp(5,1) = sqrt(ppbeam**2. + rm(1)**2)
c
c       Generally the target particle is at rest, unless we have a deuteron 
c       target, in which case we want a deuteron-momentum distribution of
c       initial state momenta.  This case is signaled by selecting (in the
c       .def file) the target to be a neutron
c
        pp(4,2) = rm(2)         !Target particle 
        if(ncode.eq.0)then      !NOT a neutron/deuteron target (normal case)
           pp(1,2) = 0
           pp(2,2) = 0
           pp(3,2) = 0
           pp(5,2) = pp(4,2)    !Particle is at rest in lab; energy is just mass
           call epcm(pp(1,1),pp(1,2),ptot)
           call epcm(pp(1,1),pp(1,2),ptotrest) !redundant statement for non-neutron target
           wvalue = ptot(4)
        else                    !neutron/deuteron target
           call getdeuteronmomentum(pp(1,2))
           call epcm(pp(1,1),pp(1,2),ptot)
           wvalue = ptot(4)
           pp(1,2) = 0          !recompute initial state with neutron at rest in a deuteron
           pp(2,2) = 0
           pp(3,2) = 0
           pp(4,2) = rmdeuteron
           call epcm(pp(1,1),pp(1,2),ptotrest)
           wvalue = ptotrest(4)
        endif

	etotal = sqrt(ptot(4)*ptot(4) + ppbeam*ppbeam) !in lab
	tbeam  = etotal - rm(1) - rm(2)  !lab kinetic energy
	pcmvaluein = (wvalue**2.-(rm(1)**2.+rm(2)**2.))**2.-
     1	                         (2.*rm(1)*rm(2))**2.
	if(pcmvaluein.gt.0.0)then
	   pcmvaluein = sqrt(pcmvaluein)/(2.*wvalue) !c.m. initial momentum
	else
	   write(6,*)'Yikes!',wvalue,rm(1),rm(2)
	   call exit
	endif
c	write(6,*)'Incident c.m. momentum',pcmvaluein,wvalue,rm(1),rm(2)
c
c	Pick the primary vertex within the target.
c
	call getvertex(r0)
c
c	Load the "pseudo" particle consisting of the incident and target 
c	particles into the output list.  We then treat this object as
c	just another decaying state.
c
c        write(6,*)'Pseudo particle:',ptot
	plist(1,1) = ptot(1)	!x momentum component of lab frame
	plist(2,1) = ptot(2)	!y momentum component
	plist(3,1) = ptot(3)	!z momentum component
	plist(4,1) = ptot(4)	!total invariant rest of initial system (W)
	plist(5,1) = ptot(5)    !total energy "pseudo" particle in lab. 
	plist(6,1) = 999        !LUND particle ID,flag for parent pseudo-particle
	plist(7,1) = r0(1)      !X location of creation
	plist(8,1) = r0(2)      !Y location of creation
	plist(9,1) = r0(3)      !Z location of creation
	plist(10,1) = 0.0       !X polarization
	plist(11,1) = 0.0       !Y polarization
	plist(12,1) = 0.0       !Z polarization
	plist(13,1) = ptot(1)	!x momentum component of lab frame at decay
	plist(14,1) = ptot(2)	!y momentum component
	plist(15,1) = ptot(3)	!z momentum component
	plist(16,1) = ptot(4)	!total invariant rest mass of system
	plist(17,1) = ptot(5)	!total energy "pseudo" particle in lab.
	plist(18,1) = r0(1)	!X location of decay; no propagation
	plist(19,1) = r0(2)	!Y location of decay
	plist(20,1) = r0(3)	!Z location of decay
	plist(21,1) = 0		!ctau for this particle type
	plist(22,1) = 999       !charge of the particle
c
c	Top of loop of decay chain
c
c	write(6,*)'> Decay the event...'
	inext = 2	!marks the next free entry in the decay list
	imark = 0	!marks the next particle to be processed
150	continue
	imark = imark + 1
	if(plist(6,imark).ne.0)then !any undecayed particles left?
	   lund = plist(6,imark)
c           write(6,*)'> Look for the next particle',imark,lund
	   if(lund.eq.999)then	!flags the initial-state "particle"
	      index = 3
	   else
	      do index=4,maxnumpar !scan list of known particles
		 if(abs(lund).eq.lundid(index))goto 160
	      end do
	      write(6,*)'Could not find particle in table',lund
	      goto 500
	   end if
c       
c          Define the next decay particle and initialize the arrays
c          for its daughter particles.
c     
 160       do m  = 1,5 
	      pv(m) = plist(m,imark)
	      pc(m) = pv(m)
	      p1(m) = 0. 
	      p2(m) = 0.
	      p3(m) = 0.
	      p4(m) = 0.
	      p5(m) = 0.
	   end do
c           write(6,*)'Loading PV',imark,pv
c
c          Propagate the particle prior to decay
c          Also apply the instrumental smearing to the momentum
c          The (crude) smearing calculation is applied to all particles
c          regardless of charge or range.
c
	   r0(1) = plist(7,imark) !creation point
	   r0(2) = plist(8,imark)
	   r0(3) = plist(9,imark)
	   if(ct(index).gt.0.0)then !ctau<0 is a special flag for the Lambda(1405) case
              ctau  = ct(index)
           else
              ctau = 0.0        !in the special case of the Lambda(1405): force ctau to be zero 
           endif
	   zpar  = zpart(index) !charge of particle 
	   if(lund.lt.0)zpar = -zpar !this is how we deal with antiparticles
	   if(lund.ne.999)then  !special treatment for initial state
              call propagate(pv,r0,r1,ctau,zpar,lund)
c	      write(6,*)'>Instrumental smearing'
c	      write(6,*)'Using GlueX smearing'
              call gsmear(pv,pvsmear,lund) !instrumental smearing a la GlueX
c             If commenting out the smearing subroutine, do:
c	      do i=1,3
c		 pvsmear(i) = pv(i)
c	      enddo
c	      pvsmear(4) = pv(4)
	   else
	      do i=1,3
		 r1(i)      = r0(i)
		 pvsmear(i) = pv(i)
	      enddo
	      pvsmear(4) = pv(4)
	      pvsmear(5) = pv(5)
	   endif
	   plist(13,imark) = pvsmear(1)	!momentum after dE/dx, at decay
	   plist(14,imark) = pvsmear(2)
	   plist(15,imark) = pvsmear(3)
	   plist(16,imark) = plist(4,imark) !copy of mass; NO smearing done
	   plist(17,imark) = pvsmear(5)  !etot(pvsmear) !energy at decay point
	   plist(18,imark) = r1(1) !decay point for straight-line track
	   plist(19,imark) = r1(2)
	   plist(20,imark) = r1(3)
	   plist(21,imark) = ctau !c_tau of this particle type
	   plist(22,imark) = zpar !charge of this particle type
c
c	   Pick a decay mode.
c	   A negative value for NDEC signifies that we should take 
c	   the first available decay mode in the list.  Later, if the
c	   decay masses don''t work out we can move to the next mode.
c	   This trick is used to make the f0 decay to KK if the
c	   f0 mass is above the KK threshold and to pi_pi if the
c	   mass is below  KK threshold.
c
	   ndec = numdecays(index) !number of decay modes for this particle
	   if(ndec.eq.0)then	!this particle is stable
	      goto 150
	   end if
	   if (ndec.eq.1)then
	      mode = 1
	   else if(ndec.gt.0)then	!pick decay mode a la branching frac.
	      which= ran(iseed)
	      rsum = 0
	      do mode = 1,abs(ndec)
		 bran = branch(index,mode,1)
		 if(which.ge.rsum .and. which.le.
     1  	      (rsum+bran))goto 170
		 rsum = rsum + bran
	      end do
	      write(6,*)'Failed to select a decay mode',index,rsum
	      write(6,*)'Is the sum of decay branches equal to 1.0?'
	      write(6,*)which,ndec,bran,lund
	      call prevent(plist,inext,ievent)
	      goto 100		!try another event
c              goto 500  !shut down analysis
	   else
	      mode = 1		!force the first mode in the list for ndec < 0
	   end if
 170	   nbody     = idecay(index,mode,1)
	   ndynamics = idecay(index,mode,7)
           dyncode1  = branch(index,mode,2)  !usually the b-slope for Regee model
           dyncode2  = branch(index,mode,3)  !usually the thresold W for the Regge model
c	   write(6,*)'> Particle decay dynamics code:',ndynamics
c	   write(6,*)'> Particle decay dynamics code:',ndynamics
c
c          Prevent unwanted phase space decay when we have pre-selected
c          a decay angle.  
c
c	   if(index.eq.3)then	!parent particle treated in special way
c	      if(ndynamics.eq.0 .and. iflat.ne.0)then
c		 write(6,*)'Phase space decay of parent not allowed'
c		 write(6,*)'when decay angle was pre-selected'
c		 call exit()
c	      endif
cc	      if(ndynamics.ne.0 .and. iflat.eq.0)then
cc		 write(6,*)'Incompatible combination of input parameters:'
cc		 write(6,*)'decay angle pre-selection (IPRULE) ',iprule,iflat
cc		 write(6,*)'specified decay dynamics (NDYNAMICS) ',ndynamics
cc		 call exit()
cc	      endif
c	   endif
c
	   id1 = idecay(index,mode,2) !Lund ID of decay particle
	   id2 = idecay(index,mode,3)
	   id3 = idecay(index,mode,4)
	   id4 = idecay(index,mode,5)
	   id5 = idecay(index,mode,6)
	   if(lund.lt.0)then	!for anti-particles, change the
	      id1 = -id1	!sign of the decay products, (decay
	      id2 = -id2	!table is for particles only)
	      id3 = -id3
	      id4 = -id4
	      id5 = -id5
	   end if
           ii = 1
           mm = 1
           nn = 1
	   do kk=4,maxnumpar	!scan list of known particles
	      if(abs(id1).eq.lundid(kk))goto 180
	   end do
	   write(6,*)'Could not find decay par1 in table',id1
	   call exit()
c	   goto 500
 180	   do jj=4,maxnumpar	!scan list of known particles
	      if(abs(id2).eq.lundid(jj))goto 190
	   end do
	   write(6,*)'Could not find decay par2 in table',id2
	   call exit()
c	   goto 500
 190	   if(nbody.ge.3)then !we have a 3- or 4- or 5-body decay
	      do ii=4,maxnumpar !scan list of known particles
		 if(abs(id3).eq.lundid(ii))goto 195
	      end do
	      write(6,*)'Could not find decay par3 in table',id3
              call exit()
c	      goto 500
	   end if
 195	   if(nbody.ge.4)then	!we have a 4- or 5-body decay
	      do nn=4,maxnumpar !scan list of known particles
		 if(abs(id4).eq.lundid(nn))goto 196
	      end do
	      write(6,*)'Could not find decay par4 in table',id4
              call exit()
c	      goto 500
	   end if
 196	   if(nbody.eq.5)then	!we have a 5-body decay
	      do mm=4,maxnumpar !scan list of known particles
		 if(abs(id5).eq.lundid(mm))goto 200
	      end do
	      write(6,*)'Could not find decay par5 in table',id5
              call exit()
c	      goto 500
	   end if
 200	   continue
c           write(6,*)'> Decay IDs   of particles:',id1,id2,id3,id4,id5
c           write(6,*)'> Indices     of particles:',kk,jj,ii,nn,mm
c           write(6,*)'> Widths:',
c     1               width(kk),width(jj),width(ii),width(nn),width(mm)
c           write(6,*)'> Masses:',
c     1               rm(kk),rm(jj),rm(ii),rm(nn),rm(mm)
c
c	   If we are decaying to a particle which has a significant width,
c	   now is the time to pick its rest mass.  This approach
c	   neglects many complications due to phase space distortions
c	   and the presence of decay thresholds for some cases.
c	   The branching ratios are not adjusted according to the masses
c	   that are picked nor their relation to decay thresholds.
c
c	   For Gaussian, WIDTH should be the sigma value.
c	   For Breit-Wigner, WIDTH should be Gamma/2.
c
c	   If switching from Gaussian to Breit-Wigner, be sure to
c	   change the definition of WIDTH, above.
c
 	   iwidthflag = .false.	!tells us whether ANY decay particle has a width
c
c          First treat the special case that the particle is ID number "100",
c          in which case we create its width as a "box" of constant
c          density and defined boundaries.  Use the nominal mass and 
c          width parameters to define the boundaries of the box.
c
c          We have no good way to control the L of the decay lineshape
c          This is a piece of pure hackery:  if the decaying particle is
c          the Lambda(1520) (indexed as the 59th particle in the decay table), 
c          set L=2, else make it L=0. 
c
c          Also treat the special case of the Lambda(1405) decaying to
c          Sigma pi with a *different lineshape* in each decay mode
c
c          On top of everything else, we have to do this for each of the two particles
c          that are decaying at this stage, the KKth and the JJth particles.  Use IX as the index.
c
 204       continue
           do itwo = 1,2
              if(itwo.eq.1)then
                 ix = kk        !treat the KK-th particle
                 iy = jj
              else
                 ix = jj        !treat the JJ-th particle
                 iy = kk
              endif
c              write(6,*)'Mass for decay particle',itwo,ix,
c     1                   rm(ix),lundid(ix),ct(ix),width(ix)
              if(lundid(ix).eq.100)then !special case of constant or "rectangular" mass distribution
                 rmlo  = rm(ix)
                 rmhi  = rm(ix) + 2.*width(ix) !the inelegant way we encoded the high mass limit
                 range = rmhi-rmlo
                 x     = ran(iseed)
                 if(itwo.eq.1)then
                    p1(4) = rmlo + x*range !uniform mass distribution within range
c                    write(6,*)'"Box" mass 1', rmlo,rmhi,p1(4),width(ix) 
                 else
                    p2(4) = rmlo + x*range !uniform mass distribution within range
c                    write(6,*)'"Box" mass 2', rmlo,rmhi,p2(4),width(ix) 
                 endif
              elseif(width(ix).eq.0)then
                 if(itwo.eq.1)then
                    p1(4) = rm(ix)
                    if(lundid(ix).eq.60)p1(4) = p1(4) + rmproton + 0.001 ! magic Lambda case... a hack
                 else
                    p2(4) = rm(ix)
                    if(lundid(ix).eq.60)p2(4) = p2(4) + rmproton + 0.001! magic Lambda case
                 endif
c                 write(6,*)'Zero-width particle',ix,rm(ix)
              else 
                 if(ct(ix).eq.-1)then !special case of Lambda(1405) lineshapes
cc                    call getmasses(ix,rm1,rm2,1)
cc                    call lineshapepick(rm(ix),rm1,rm2,ppbeam,rmassvalue,iflag)
cc                    if(iflag.eq.1)goto 105 !abort this event
cc                    p2(4) = rmassvalue
cc                    write(6,*)'Lineshape mass value:',rmassvalue
                    write(6,*)
     1              'Selected lineshape option is not implemented'
                    call exit()
                 elseif(width(ix).gt.0)then !generate a pure Lorentzian decay lineshape
 205                continue
                    call bwran(x) !simple Breit-Wigner random number
                    if(itwo.eq.1)then
                       p1(4) = rm(ix) + width(ix)*x
c                       write(6,*)'Simple BW 1>',rm(ix),width(ix),p1(4),x
                    else
                       p2(4) = rm(ix) + width(ix)*x
c                       write(6,*)'Simple BW 2>',rm(ix),width(ix),p2(4),x
                    endif
                    if(p1(4).lt.0)then
c                       write(6,*)'Bad mass 1: ',rm(ix),width(ix),p1(4),x
                       goto 205
                    endif
                    if(p2(4).lt.0)then
c                       write(6,*)'Bad mass 2: ',rm(ix),width(ix),p2(4),x
                       goto 205
                    endif
                 else           !negative "width" parameter signals that we want the relativistic BW
                    gamval = -2.*width(ix) !earlier we divided FWHM by 2; undo this here
                    if(ct(ix).eq.-10)then
                       lvalue = 0
                    elseif(ct(ix).eq.-11)then
                       lvalue = 1
                    elseif(ct(ix).eq.-12)then
                       lvalue = 2
                    else
                       lvalue = 0
                    endif
                     call getmasses(ix,rm1,rm2,1)
c                    rm3 = plist(iy,4)   !need code that figures this out automatically
                    rm3 = rmproton !hack
                    if(lundid(ix).eq.96)then !special case of proton anti-Lambda
                       rm3 = rmlambda
                    endif
c
c                    write(6,*)'> Calling rel. BW lineshape calculation',
c     1                   ix,lundid(ix),rm(ix),rm1,rm2,rm3,lundid(iy),ww
                    call relativisticlineshape(rm(ix),gamval,
     1                            rm1,rm2,lvalue,rm3,ww,bwvalue,ifail)
c                    write(6,*)'Orbital L> ',lvalue,itwo,bwvalue
                    if(itwo.eq.1)then
                       p1(4) = bwvalue  !we are treating the first of the two decay particles...
                    else
                       p2(4) = bwvalue  !...now the second one.
                    endif
c                    if(iloop.eq.1)write(6,*)'BW value >>>',itwo,bwvalue
                 endif
                 iwidthflag = .true.
              end if
c              write(6,*)'Bottom of two-body mass selection...',ix,itwo
              if(itwo.gt.2)then
                 write(6,*)'Yikes!'
                 call exit()
              endif
           enddo
c
c           write(6,*)'Many-body decay:',nbody
	   if(nbody.eq.3.or.nbody.eq.4 .or. nbody.eq.5)then
	      if(width(ii).eq.0)then
		 p3(4) = rm(ii)
	      else
c                 call rannor(x,dummy)
		 call bwran(x)	!Breit-Wigner random number
		 p3(4) = rm(ii) + width(ii)*x
		 iwidthflag = .true.
	      end if
c
	      if(nbody.ge.4)then
		 if(width(nn).eq.0)then
		    p4(4) = rm(nn)
		 else
c                 call rannor(x,dummy)
		    call bwran(x) !Breit-Wigner random number
		    p4(4) = rm(nn) + width(nn)*x
		    iwidthflag = .true.
		 end if
c
		 if(nbody.eq.5)then
		    if(width(mm).eq.0)then
		       p5(4) = rm(mm)
		    else
c                       call rannor(x,dummy)
		       call bwran(x) !Breit-Wigner random number
		       p5(4) = rm(mm) + width(mm)*x
		       iwidthflag = .true.
		    end if
		 endif
	      end if
	   else
	      p3(4) = 0
	      p4(4) = 0
	      p5(4) = 0
	   end if
c
c	   Check that the decay products are less massive than
c	   the particle that is decaying.  If not, either the
c	   decay table is messed up, or we are dealing with a
c	   "threshold" case like the f0 meson or the Lambda(1405).
c
c	   write(6,*)'>Mass check',pv(4),p1(4),p2(4),p3(4),p4(4),p5(4)
	   if(pv(4).le.(p1(4)+p2(4)+p3(4))+p4(4)+p5(4))then
	      miterate = miterate + 1
              iattempt = iattempt + 1
 	      if(miterate.eq.1)then
                 write(6,*)'Decay masses too large; iterating masses'
                 write(6,210)lund,id1,id2,id3,id4,id5,
     1                    pv(4),p1(4),p2(4),p3(4),p4(4),p5(4)
              endif
	      if(miterate.ge.10)then   !too many iterations... give up
                 write(6,*)'Give up on this event', miterate
c                 write(6,210)lund,id1,id2,id3,id4,id5,
c     1                pv(4),p1(4),p2(4),p3(4),p4(4),p5(4)
 210             format(1x,6i7,/,1x,6f8.3)
                 miterate = 0
                 goto 105
              endif
c
	      if(miterate.gt.20)then !flush the decay mass but not the parent mass
c		  if(.not.imess)then
                 write(6,*)'>>> Mass selection inefficient!', miterate
c                 write(6,210)lund,id1,id2,id3,id4,id5,
c     1                pv(4),p1(4),p2(4),p3(4),p4(4),p5(4)
c 210             format(1x,'Sub-threshold? Must pick new masses or decay mode',
c     1                  6i4,6f8.3)
c                 end if
c                 if(width(kk).ne.0.0 .or. width(jj).ne.0.0)then
                    goto 204    !save this event with new decay masses selection
c                 else
c                    write(6,*)'Give up on decay  with delta-fn widths'
c                    miterate = 0
c                    goto 105
c                 endif
	      endif
c	      if(.not.imess)then
c		 write(6,210)lund,id1,id2,id3,id4,id5,
c     1                       pv(4),p1(4),p2(4),p3(4),p4(4),p5(4)
c		 imess = .true.
c	      end if
	      if(iwidthflag)then 
		 goto 105       !abort if any particles have finite width
	      else		!try other decay mode
		 mode = mode + 1
		 if(mode.le.abs(ndec))then !next mode
c		    write(6,*)'Try next decay mode'
		    goto 170	
		 else		!no more modes to try
c		    write(6,220)ndec,lund,id1,id2,id3,id4,
c     1  		      kk,width(kk)
c 220		    format(1x,'Bloody Murder',7i7,f10.3)
c		    call prevent(plist,inext,ievent)
c                    call exit()
c		     stop
		    goto 105	!abort this event and try again
		 end if
	      end if
	   end if
           
c
c   	   Now let particle called PV decay
c
c	   write(6,*)'>Top of decay selection loop',nbody,ndynamics
	   if(nbody.eq.2)then
c
c	      Load up the direction in which particle PV may 
c	      have some polarization.  This direction was defined
c	      at the time this particle was (strongly) created.
c
 	      polv(1) = plist(10,imark) !X
	      polv(2) = plist(11,imark) !Y
	      polv(3) = plist(12,imark) !Z
	      polnorm = psq(polv)	
c
c	      If the decay is at rest we have no axis for weak
c	      decay asymmetry.  If we were asking for a weak
c	      decay but we are at rest, change to pure phase
c	      space decay.
c
	      if(ndynamics.eq.1)then
		 pmag2 = psq(pv)
		 if(pmag2.eq.0)then
		    ndynamics = 0
                    write(6,*)'Zero momentum PV decay'
		 end if
		 if(polnorm.eq.0)then
		    ndynamics = 0
c                   write(6,*)'Zero pol. vector PV decay'
		 end if
	      end if
c     
c             Hack to save kaon momentum in c.m. frame
c             Used later for 3-body phase space calculation involving Sigma(1385)
c
              if (id1.eq.+18)then
                 pcminter = pcm
              else
                 pcminter = 1.0
              endif
c
c             Begin here a long series of possible decay modes, as specified for each
c             particle in the definition file. NDYNAMICS is the variable that specifies
c             the manner of the decay.
c
c             For 2-body decays we allow:
c              NDYNAMICS = 0 means isoptropic 2-body decay + a polarization assigned to products
c              NDYNAMICS = 1 means weak decay, as per Lambda -->  pi N
c              NDYNAMICS = 2 means nuclear decay with small momentum (check if this still works)
c              NDYNAMICS = 3 decay at a pre-selected polar angle (random by default)
c              NDYNAMICS = 4 means isotropic 2-body decay with polarization transfer hacked in
c              NDYNAMICS = 5 means isotropic 2-body decay & accepted with 2-body phase space weighting
c              NDYNAMICS = 6 means at a pre-selected polar angle & accepted with 2-body phase space weighting
c              NDYNAMICS = 7 means isotropic 2-body decay & accepted with 3-body phase space weighting
c              NDYNAMICS = 8 means Regge-type of decays with t-slope and W-dependence parameters given
c             
c             For 3-body decays we allow:
c              NDYNAMICS = 0 means pure three-body phase space
c              NDYNAMICS = 1 means Fermi gas nuclear model
c              NDYNAMICS = 2 means nuclear Hulthen wavefunction model
c              NDYNAMICS = 3 means 3-body decay with subsequent Regge-type direction selection
c
	      if(ndynamics.eq.0)then !decay isotropically & polarize products
		 call kin1r(pv,p1,p2,cthcm,pcm,ierr)
c                 write(6,*)pv         
c                 write(6,*)p1
c                 write(6,*)p2
c                 write(6,*)pcm
c                 write(6,*)cthcm
		 if(ierr.ne.0)then
		    write(6,*)'KIN1R is messed (nd=0)',pv,p1,p2,cthcm
		    call prevent(plist,inext,ievent)
		 endif
		 call polarize(pv,p1,p2,id1,id2,pol1,pol2,cthcm)
c                 write(6,*)' '
c                 write(6,*)
c     1             '> Returning from two-body phase space decay:'
c                 write(6,*)'PV',pv
c                 write(6,*)'P1',p1, id1
c                 write(6,*)'P2',p2, id2
c                 write(6,*)'POL1',pol1
c                 write(6,*)'POL2',pol2
c                 write(6,*)'PCM',pcm
c                 write(6,*)'CTHCM',cthcm 
c     
	      else if(ndynamics.eq.1)then !e.g. Lambda decay case
c                 write(6,*)'Calling weakdecay'
		 call weakdecay(pv,p1,p2,polv,pol1,pol2,cthcm,lund,ierr)
c                 write(6,*)'Retruning from weakdecay'
		 if(ierr.ne.0)write(6,*)'WEAKDECAY is messed...'

 	      else if(ndynamics.eq.2)then !nuclear recoil case
		 isum = 0
 230		 call kin1r(pv,p1,p2,cthcm,pcm,ierr)
		 if(ierr.ne.0)write(6,*)'KIN1R is messed (nd=2)',pv,p1,p2,cthcm
c
c	 	 Assume p2 is the recoiling nucleus
c		 which started at rest in the lab.
c		 For now, do a most naive box cut,
c		 at 250 MeV/c.
c
		 precoil = absp(p2)
		 if(precoil .gt. 0.500)then  !GeV/c
		    isum = isum + 1
		    if(isum.gt.500)then
		       write(6,*)'Ack! No luck generating small recoil events'
		       isum = 0
		    end if
		    goto 230
		 else
		    continue
		 end if
		 call polarize(pv,p1,p2,id1,id2,pol1,pol2,cthcm)
                 
	      else if(ndynamics.eq.3)then !photoproduction at a pre-selected decay polar angle
		 call product_stuffit(pv,p1,p2,pol1,pol2,
     1  	      id1,id2,cosinethetacm,lund,ierr)
                 if(ierr.ne.0)then
c                    write(6,*)'Returning from product_stuffit'
c                    write(6,*)lund,id1,id2,cosinethetacm
c                    write(6,*)pv
c                    write(6,*)p1
c                    write(6,*)p2
                 endif	
                 if(p1(1)==p1(1))then !test for NaN status
                    continue
                 else
                    write(6,*)'Bloody Hell',p1,p2
                 endif
                 if(pol1(1)==pol1(1))then !test for NaN status
                    continue
                 else
                    write(6,*)'Bloody Hell: POL1 is bad',pol1,pol2
                 endif
		 if(ierr.ne.0)write(6,*)'product_stuffit is messed...'
                 
	      else if(ndynamics.eq.4)then !isotropic angle & polarization transfer
		 call kin1r(pv,p1,p2,cthcm,pcm,ierr)
		 if(ierr.ne.0)write(6,*)'KIN1R is messed (nd=6)',pv,p1,p2,cthcm
		 do i = 1,3
		    pol1(i) = polv(i)
		    pol2(i) = polv(i)
		 enddo
                 
	      else if(ndynamics.eq.5)then !isotropic angle & phase space weighting
c                Keep the event with a weight proportional to the
c                available 2-body phase space, i.e. according to the available momentum.
c
c                We should not use this 'dynamics' if there is more than one path
c                in the decay chain of the simulation that has the need for 
c                weighting of this sort.  Routine PHASESPACEWEIGHT dynamically seeks
c                the highest c.m. momentum and adjusts its weighting accordingly.  
c                Warning: If two or more decay paths use the routine, this algorithm 
c                will get confused. 
c
		 call kin1r(pv,p1,p2,cthcm,pcm,ierr)            !pure phase space
		 pcmvaluein = 1.0  !dummy input momentum for decay-weighting alone
		 call phasespaceweight(pv,p1,p2,pcmvaluein,ipskeep)    !apply momentum weighting
		 if(ipskeep.eq.0)then !the event fails the phase-space weighting: abort
c		    write(6,*)'Event fails phase space weighting...abort'
		    goto 105
		 endif

	      else if(ndynamics.eq.6)then !generic particle photoproduction + phase space weighting
c		 write(6,*)'>>Calling product_stuffit'
		 call product_stuffit(pv,p1,p2,pol1,pol2,
     1  	      id1,id2,cosinethetacm,lund,ierr)
                 if(ierr.ne.0)write(6,*)'product_stuffit is messed...2'
		 call phasespaceweight(pv,p1,p2,pcmvaluein,ipskeep)    !apply momentum weighting
		 if(ipskeep.eq.0)then !the event fails the phase-space weighting: abort
c		    write(6,*)'Event fails phase space weighting...abort'
		    goto 105
		 endif

	      else if(ndynamics.eq.7)then !isotropic angle & 3-body phase-space weighting
c                Keep the event with a weight proportional to the
c                available 3-body phase space.
c
c                We should not use this 'dynamics' if there is more than one path
c                in the decay chain of the simulation that has the need for 
c                weighting of this sort.  Routine PHASESPACEWEIGHT3 dynamically seeks
c                the highest c.m. momentum and adjusts its weighting accordingly.  
c                Warning: If two or more decay paths use the routine, this algorithm 
c                will get confused. 
c
c                There is a necessary hack here:  the 3-body phase space needs to know
c                the momentum of the first of the three bodies that was generated already 
c                previously.  We have to reach back and get that momentum for the intermediate state,
c                called PCMINTER below.  This code was written for the case of a K+ Lambda(1520)
c                state (with momentum PCMINTER) decaying into a Sigma(1385) and a pion.
c
c                 pcminter = 1.0  !hack test
		 call kin1r(pv,p1,p2,cthcm,pcm,ierr)             !isotropic 2-body phase space
		 call phasespaceweight3(pv,p1,p2,pcminter,pcmvaluein,wvalue,
     1                                             threeweight,ipskeep) !apply 3-body event weighting
		 if(ipskeep.eq.0)then !the event fails the phase-space weighting: abort
c		    write(6,*)'Event fails 3-body phase space weighting...abort'
		    goto 105
		 endif

	      else if(ndynamics.eq.8)then !decay of parent particle according to Regge-like angle weighting
c     
c                The *.def file must have the recoiling object FIRST in the list of decay particles,
c                and the object to which t is transferred is the SECOND of the two objects
c                given in the 2-body decay list
c
c                indexparent = 3 !the third particle in the list is the c.m. particle
c                call getmasses(indexparent,rm1,rm2,1)  !works only for the first mode in the decay list
                 rm1 = p1(4)  !we already know the specific (possibly "wide") decay masses of the products
                 rm2 = p2(4)
c                write(6,*)'Decay masses:', rm1,rm2
                 bslope1    = dyncode1 ! was set in definition file
                 rmproduce  = rm2 ! the "produced" thing is the forward-going object at the photon vertex
                 rmrecoil   = rm1 ! the recoiling thing is listed first in the definition file
                 rmtarget   = rmproton
                 wthreshold = dyncode2 !specified in the definition file
                 sminf      = 0.001 !parameter controling the "flatness" at back angles
                 call reggepick(ppbeam,pbeammax,wthreshold,bslope1,sminf,
     1                rmtarget,rmproduce,rmrecoil,cosinethetacm,iregge)

c                 cosinethetacm = .5      !Test hack
                 
c                 write(6,*)'RP> ',ppbeam,pbeammax,cosinethetacm,iregge,bslope1
c                 write(6,*)'    ',rmtarget,rmrecoil,rmproduce,wthreshold
                 if(iregge.lt.0)then !abandon this event
                    write(6,*)'Reggepick error'
                    write(6,*)ndynamics,rmproduce,rmrecoil,rmtarget,rm1,rm2
                    write(6,*)p1
                    write(6,*)p2
                    goto 105
                 endif
                 iflat = 1
		 call product_stuffit(pv,p1,p2,pol1,pol2,
     1  	      id1,id2,cosinethetacm,lund,ierr)
		 call polarize(pv,p1,p2,id1,id2,pol1,pol2,cosinethetacm)
                 if(ierr.ne.0)then
                    write(6,*)'Returning from product_stuffit'
                    write(6,*)lund,id1,id2,cosinethetacm
                    write(6,*)'PV  ', pv
                    write(6,*)'P1  ', p1
                    write(6,*)'P2  ', p2
                    write(6,*)'POL1', pol1
                    write(6,*)'POL2', pol2
                 endif            

              else
		 write(6,*)'Unspecified decay dynamics ',ndynamics
		 call exit()
	      end if
c     
c             End of options for letting two-body final states decay
c ********************************************************************
              
	   else if(nbody.eq.3)then ! three-body decay modes
	      if(ndynamics .eq. 0)then !decay via pure phase space
c                 write(6,*)' '
c	 	  write(6,*)'>Calling  3-body decay code',pv(4),p1(4),p2(4),p3(4)
		 call kin3body(pv,p1,p2,p3,ierr)
c		 write(6,*)'>Return from 3-body decay code'
c                 write(6,*)'parent',pv
c                 write(6,*)'P1    ',p1
c                 write(6,*)'P2    ',p2
c                 write(6,*)'P3    ',p3
	      else if(ndynamics.eq.1)then !nuclear Fermi gas decay
c		 write(6,*)'>Calling 3-body decay code',pv(4),p1(4),p2(4),p3(4)
		 call kin3bodyqf(pv,p1,p2,p3,ierr)
                 if(ierr.gt.0)then
                    write(6,*)'QF generation failed: discard event'
                    goto 105  !abort event
                 endif
	      else if(ndynamics.eq.2)then !Hulthen momentum distribution (for deuteron)
		 write(6,*)'>Calling Hulthen decay code',pv(4),p1(4),p2(4),p3(4)
		 call kin3bodyhulthen(pv,p1,p2,p3,ierr)
                 if(ierr.gt.0)then
                    write(6,*)'QF generation failed: discard event'
                    goto 105  !abort event
                 endif
c
c             First generate pure 3-body phase space.   Then create a pseudo particle out of P2 and P3.
c             Arrange that this particle P23 satisfies the specified t-slope requirement.
c             Thus, we do not have to pre-select a mass of the intermediate state particle.
c             NDYNAMICS = 3  - rotate the 3-body phase-space state so the P23 combo has given t-slope
c             NDYNAMICS = 4  - rotate the 3-body phase-space state so the P23 combo has given t-slope and
c                              also enforce the t-slope at the bottom vertex:  "double - Regge" mechanism
c             NDYNAMICS = 5  - rotate the 3-body phase-space state so the P23 combo has given t-slope and
c                              keep event where the invariant mass of P23 follows a desired slope
c
              else if(ndynamics.ge.3 .and. ndynamics.le.5)then !combine 3-body production & Regge picture
                 idcount   = 0
 240             pvtemp(1) = 0
                 pvtemp(2) = 0
                 pvtemp(3) = 0
                 pvtemp(4) = pv(4)
                 pvtemp(5) = pv(4)
		 call kin3body(pvtemp,p1,p2,p3,ierr) !create the three final state particles in CM frame
c		 write(6,*)'>Return from 3-body decay code for CM'
c                 write(6,*)'parent',pvtemp
c                 write(6,*)'P1    ',p1
c                 write(6,*)'P2    ',p2
c                 write(6,*)'P3    ',p3
c                Group particles 2 and 3 together (P23), recoiling against particle 1 P(1)
c                Let the *.def file have the recoiling object P(1) FIRST in the list of decay particles,
c                Particles 2 and 3 taken together form the SECOND particle to which 4-mometum t transfers
c
c                indexparent = 3 !the third particle in the list is the c.m. particle
c                call getmasses(indexparent,rm1,rm2,1)  !works only for the first mode in the decay list
                 rm1 = p1(4)  !we already know the decay masses of the products
                 rm2 = p2(4)
                 rm3 = p3(4)
                 call invarmass(p2,p3,p23,rmassinv)
c                 write(6,*)'3 body decay in the CM frame >>'
c                 write(6,*)'PV >>',pv
c                 write(6,*)'P1 >>',p1
c                 write(6,*)'P2 >>',p2
c                 write(6,*)'P3 >>',p3
c                 write(6,*)'P23>>',p23
                 rm23 = p23(4)
                 bslope1     = dyncode1 ! was set in definition file
                 bslope2     = dyncode2 ! was set in definition file
                 if(ndynamics.eq.3)then ! Mechanism 4, single Regge with momentum transfer to produced pair
                    rmproduce  = rm23 ! the "produced" thing is the forward-going object at the photon vertex
                    rmrecoil   = rm1  ! the recoiling thing is listed first in the definition file (e.g. target proton)
                 elseif(ndymanics.eq.4)then ! Mechanism 5, double Regge 
                    rmproduce  = rm1  ! the "produced" thing is the forward-going object at the photon vertex
                    rmrecoil   = rm23 ! the recoiling thing here the pair of other particles taken together
                 elseif(ndymanics.eq.5)then  ! Mechanism 6,  single Regge with P23 invariant mass tuning
                    rmproduce  = rm23  ! the "produced" thing is the forward-going object at the photon vertex
                    rmrecoil   = rm1   ! the recoiling thing is listed first in the definition file (e.g. target proton)
                 endif
c                
c                Now test if the P23 pair has an invriant mass in line with the selected
c                slope.  This is our way of modeling possible clustering of the particle pair.
c                "Mechanism 6" 
c                 
                 if(ndynamics.eq.5)then
c                    write(6,*)'Trimming invariant mass'
c                    write(6,*)p2
c                    write(6,*)p3
c                    write(6,*)p23
c                    write(6,*)rm23
                    rmassmin = p2(4) + p3(4) !least value the invariant mass can have
                    call sculptmass(rm23,rmassmin,rmassmax,bslope2,iyesno)
                    if(iyesno.eq.0)then
c                       write(6,*)'Mass rejected',rm23
                       goto 240 !re-throw the 3-body phase space particles
                    else
c                       write(6,*)'Mass accepted',rm23
                    endif
                 endif
c
                 rmtarget   = rmproton
                 wdummy     = xdummy !(actually not used at this time)
                 sminf      = 0.001 !parameter controling the "flatness" at back angles
                 call reggepick(ppbeam,pbeammax,wdummy,bslope1,sminf,
     1                rmtarget,rmproduce,rmrecoil,cosinethetacm,iregge)
c                 write(6,*)'RP> ',ppbeam,pbeammax,cosinethetacm,iregge,bslope1
c                 write(6,*)'    ',rmtarget,rmrecoil,rmproduce,wdummy
                 if(iregge.lt.0)then !abandon this event
                    write(6,*)'Reggepick error'
                    write(6,*)ndynamics,rmproduce,rmrecoil,rmtarget,rm1,rm2
                    write(6,*)p1,p2,p3,p23
                    goto 105
                 endif
                 iflat = 1
                 if(ndynamics.eq.3)then
c                    write(6,*)'Rotating P23 to the desired t'
                 elseif(ndynamics.eq.4)then
c                    write(6,*)'Rotating P1 to the desired t'
                    cosinethetacm = -cosinethetacm   !kinematic trick to save coding steps
                 elseif(ndynamics.eq.5)then
c                    write(6,*)'Rotating P23 to the desired t'
                 else
                    write(6,*)'Inconsistent dynamics specified'
                    call exit()
                 endif
c                 write(6,*)'>> CM angle',cosinethetacm
		 call product_rotate(pv,p1,p2,p3,p23,cosinethetacm,ierr)
c
c                Return the tracks to the Lab frame.
c                 
                 call boost_to_lab(pv,p1,p2,p3,p23,ierr)
                 if(ierr.ne.0)then !abandon this event
c                    write(6,*)' > Returning from product_rotate'
c                    write(6,*)cosinethetacm
c                    write(6,*)'PV  ', pv
c                    write(6,*)'P1  ', p1
c                    write(6,*)'P2  ', p2
c                    write(6,*)'P3  ', p3
c                    write(6,*)'P23 ', p23
                    goto 105
                 endif            
c
c                "Mechanism 5"
c                Now test whether the momentum transfer at the bottom vertex is small
c                if we are doing the Double Regge mechanism (ndynammics = 4)
c                Remember we are in the lab now.  The momentum transfer is between the 
c                target proton (must construct it first) and the recoiling target particle.
c
c                Note: in the lab frame, tprime is always 0, by construction, because tmin and t 
c                are the same.
c
c                Warning: we must hard-code the target/recoil particle to be number 2 in the decay list. 
c                Be sure the input definition file has the proton in the second slot.
c                 
                 if(ndynamics.eq.4)then
                    ptarget(1) = 0.0
                    ptarget(2) = 0.0
                    ptarget(3) = 0.0
                    ptarget(4) = rmproton
                    ptarget(5) = rmproton
                    call gett(ptarget,p2,tt,ttmin,ttprime)
c                 write(6,*)'t at bottom vertex',tt,ttmin,ttprime
                    call doublereggepick(ppbeam,rmtarget,tt,bslope2,iyesno,ierr)
                    if(iyesno.eq.0)then
                       idcount = idcount + 1
c                       write(6,*)'Lower vertex t is too shallow; try again'
                       if(idcount.ge.idcountmax)then
                          idcountmax = idcount
                          if(mod(idcountmax,1000).eq.0)then
                             write(6,*)'Slow bottom vertex selection',idcount
                          endif
                       endif
                       goto 240 ! re-throw the 3-body phase space particles
                    endif
                 endif
	      else
		 write(6,*)'Unspecified 3-body dynamics ',ndynamics
		 call exit()
	      end if
c              call polarize(pv,p1,pv2,id1,id2,pol1,pol2,cthcm)
c              call polarize(pv,p1,pv3,id1,id3,pol1,pol3,cthcm) !this code not checked
c              
	   else if (nbody.eq.4 .or. nbody.eq.5)then ! four/five-body decay modes
	      if(ndynamics .eq. 0)then !decay via pure phase space
c                 write(6,*)'>Calling 4- or 5-body decay code'
		 call kin45body(nbody,pv,p1,p2,p3,p4,p5,ierr)
c                 call exit()
c                 write(6,*)'>Return from 4- 5-body decay code'
		 if(ierr.eq.1)then
		    goto 105    !abort event
		 endif
              elseif(ndynamics.eq.2)then !demand small recoil momentum for particle p1
c		 write(6,*)'>Calling 4- or 5-body decay code'
                 isum = 0
 250             call kin45body(nbody,pv,p1,p2,p3,p4,p5,ierr)
		 if(ierr.eq.1)then
		    goto 105    !abort event
		 endif
c
c	 	 Assume p1 is the spectator. Do a most naive box cut,
c		 at PFERMI GeV/c.
c     
		 precoil = absp(p1)
		 if(precoil .gt. pfermi)then !GeV/c
		    isum = isum + 1
		    if(isum.gt.500)then
		       write(6,*)'Ack! No luck generating small recoil events'
		       isum = 0
		    end if
		    goto 250
                 end if
c                 write(6,*)isum, precoil
              else
		 write(6,*)'Unspecified N-body dynamics ',ndynamics
		 call exit()
	      end if
c              call polarize(pv,p1,pv2,id1,id2,pol1,pol2,cthcm)
	   else
	      write(6,*)'Too many final-state particles requested',nbody
	      call exit()
	   end if
c
c	   write(6,*)'> Load particle into event table PLIST',inext
c           write(6,*)id1,p1
c           write(6,*)id2,p2
c           write(6,*)id3,p3
	   plist(1,inext)  = p1(1)    !X momentum
	   plist(2,inext)  = p1(2)    !Y momentum
	   plist(3,inext)  = p1(3)    !Z momentum
	   plist(4,inext)  = p1(4)    !mass
	   plist(5,inext)  = p1(5)    !total energy
	   plist(6,inext)  = id1
	   plist(7,inext)  = plist(18,imark) !X creation point
	   plist(8,inext)  = plist(19,imark) !Y creation point
	   plist(9,inext)  = plist(20,imark) !Z creation point
	   plist(10,inext) = pol1(1) !X polarization
	   plist(11,inext) = pol1(2) !Y polarization
	   plist(12,inext) = pol1(3) !Z polarization
	   inext = inext + 1
	   plist(1,inext)  = p2(1) !X momentum
	   plist(2,inext)  = p2(2)
	   plist(3,inext)  = p2(3)
	   plist(4,inext)  = p2(4)    !mass
	   plist(5,inext)  = p2(5)    !total energy
	   plist(6,inext)  = id2
	   plist(7,inext)  = plist(18,imark) !X creation point
	   plist(8,inext)  = plist(19,imark) !Y creation point
	   plist(9,inext)  = plist(20,imark) !Z creation point
	   plist(10,inext) = pol2(1) !X polarization
	   plist(11,inext) = pol2(2) !Y polarization
	   plist(12,inext) = pol2(3) !Z polarization
	   inext = inext + 1
	   if(nbody.eq.2)goto 150  !Quit here for 2-body modes
	   plist(1,inext)  = p3(1) !X momentum
	   plist(2,inext)  = p3(2)
	   plist(3,inext)  = p3(3)
	   plist(4,inext)  = p3(4) !mass
	   plist(5,inext)  = p3(5) !total energy
	   plist(6,inext)  = id3
	   plist(7,inext)  = plist(18,imark) !X creation point
	   plist(8,inext)  = plist(19,imark) !Y creation point
	   plist(9,inext)  = plist(20,imark) !Z creation point
	   plist(10,inext) = pol3(1) !X polarization
	   plist(11,inext) = pol3(2) !Y polarization
	   plist(12,inext) = pol3(3) !Z polarization
	   inext = inext + 1
	   if(nbody.eq.3)goto 150  !Quit here for 3-body modes
	   plist(1,inext)  = p4(1) !X momentum
	   plist(2,inext)  = p4(2)
	   plist(3,inext)  = p4(3)
	   plist(4,inext)  = p4(4) !mass
	   plist(5,inext)  = p4(5) !total energy
	   plist(6,inext)  = id4
	   plist(7,inext)  = plist(18,imark) !X creation point
	   plist(8,inext)  = plist(19,imark) !Y creation point
	   plist(9,inext)  = plist(20,imark) !Z creation point
	   plist(10,inext) = pol4(1) !X polarization
	   plist(11,inext) = pol4(2) !Y polarization
	   plist(12,inext) = pol4(3) !Z polarization
	   inext = inext + 1
	   if(nbody.eq.4)goto 150  !Quit here for 4-body modes
	   plist(1,inext)  = p5(1) !X momentum
	   plist(2,inext)  = p5(2)
	   plist(3,inext)  = p5(3)
	   plist(4,inext)  = p5(4) !mass
	   plist(5,inext)  = p5(5) !total energy
	   plist(6,inext)  = id5
	   plist(7,inext)  = plist(18,imark) !X creation point
	   plist(8,inext)  = plist(19,imark) !Y creation point
	   plist(9,inext)  = plist(20,imark) !Z creation point
	   plist(10,inext) = pol5(1) !X polarization
	   plist(11,inext) = pol5(2) !Y polarization
	   plist(12,inext) = pol5(3) !Z polarization
	   inext = inext + 1
	   goto 150
	else
c
c          Here is a hack to select a second photon after the event
c          was already generated by the first "real" photon.  We studied this in
c          11-2009 to simulate events with "wrong" photons in the tagger.
c          100% of the events generated with the following lines turned on will 
c          be messed up by having the wrong photon energy.  The kinematics of 
c          all the decay products will correspond to the original "correct" photon,
c          but the initial state will be wrongly computed using the false photon.
c
cc   	   ppbeamold = ppbeam
cc   	   call bremsstrahlung(pbeammin,pbeammax,ppbeam)
ccc	   write(6,*)'Changing the photon energy ',ppbeamold,ppbeam
cc	   pp(1,1) = 0		!Beam particle
cc	   pp(2,1) = 0
cc	   pp(3,1) = ppbeam
cc	   pp(4,1) = rm(1)
cc	   pp(1,2) = 0		!Target
cc	   pp(2,2) = 0
cc	   pp(3,2) = 0
cc	   pp(4,2) = rm(2)
cc	   call epcm(pp(1,1),pp(1,2),ptot)
cc	   wvalue = ptot(4)
cc	   plist(1,1) = ptot(1)	!x momentum component of lab frame
cc	   plist(2,1) = ptot(2)	!y momentum component
cc	   plist(3,1) = ptot(3)	!z momentum component
cc	   plist(4,1) = ptot(4)	!total invariant rest mass of system (W)
cc	   plist(5,1) = wvalue	!total energy "pseudo" particle in lab.
c
c          end of hack
c     
	   if(ievent.le.nprint)then
	      if(ievent.eq.1)	write(6,*)'The first event(s):'
	      call prevent(plist,inext,ievent)
	   endif
c
c          Calling the event analyzer must be optional, since we are planning to
c          split off the analyzer as a stand-alone program that reads in the
c          data from the disk file.
c          IDOCUTS turns analysis on or off.   But to remove the analyzer completely
c          it is necessary to comment it out here AND below in the closeout routines.
c
c           if(idocuts)write(6,*)'Event analyzer is not linked into the code.' 
           if(idocuts) call mc_analyze(plist,imark-1)
	   if(idogluex)call write_output(plist,imark-1,imechanism)
	   if(numwrite.eq.ngevent)then
              write(6,*)'Hit the limit of number of events to write:',
     1             numwrite,ngevent
              goto 500
           endif
	   goto 100
	end if
c
c       All done except for the final accounting.
c
 500	write(6,*)'Total number of events generated =',ievent-1
        if(idocuts)then
c           write(6,*)'No histograms created'
           call histogram_finish
           call analysis_finish
        else
           write(6,*)'No event analysis or histogramming was done.'
        endif
c
c       Save the random number seed(s) in case we later want to pick up
c       the sequence from here.
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
	write(6,*)'Successful end of event generation.'
        call exit()
	end
c******************** End of main program mc_gen ***********************************

      
c
c	Shut down the event generation
c
	subroutine analysis_finish
c
	character*256 prg_name, histnm, filenm, prefix
	real corrsum(3,5),correrr(3,5)
	data alpha,alphabar/0.642,-0.642/  !weak decay asymmetry
c	data alpha,alphabar/1.0 ,-1.0/  !HACK FOR TESTS 
	common /stuff/poldir(4),xx,icount,ikeep,ihackflag,iattempt
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
c        open(unit=8,file='mc_generator.out',status='unknown',
c     1                                       access='append')
c	write(8,20)pbeam,icount,ikeep,percent
c 20	format(1x,f7.2,2i10,f7.2)
c	close(8)
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
	return
	end 

c********************************************************************************
c
c	Close the histogramming function 
c
	subroutine histogram_finish
c
	parameter (nwpawc = 20000000)
	character*256 prg_name, histnm, filenm, prefix
	common /pawc/hmemor(nwpawc)
	common /run_control/n_arg,prg_name,filenm,histnm,prefix,
     1                      irestart,lfl,lprefix,maxbosoutput
c
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
c       
	write(6,*)'Histrogam file closed.'
	return
	end 
      
c********************************************************************************
c
c       Output to file in text format suitable for conversion to full-scale Monte Carlo modeling
c        -MC_GENERATOR_2_HDDM (for GlueX) 
c
 	subroutine write_output(plist,imark,imechanism)
	parameter (maxnumlis=60) !Maximum number of particles in the decay list
	character nam*4,fname*200
	real plist(24,maxnumlis),pgsim(4)
	integer igsimid(maxnumlis)
	logical ifirst
 	data ifirst/.true./
	common /beam/pbeam,ppbeam,rmbeam,rmtarg,jcode,xsigma,ysigma,ncode
	common /eventiolocal/ifirst,ifile,ieventout
	common /monte/ iseed,xflat,xnormal,xbw
c
	character*256  prg_name, histnm, filenm, prefix
	common /run_control/n_arg,prg_name,filenm,histnm,prefix,
     1                      irestart,lfl,lprefix,maxbosoutput
c
c	On first call just open the data file
c
	if(ifirst)then
		ifirst = .false.
		ifile  = 1                !specifies starting number of output files
		ieventout = 0
		write(nam,'(i4)'),ifile   !fancy internal write to convert # to ascii
		if(ifile.lt.10)then
		   nam(3:3)=char(48)
		endif
		if(ifile.lt.100)then
		   nam(2:2)=char(48)
		endif
		if(ifile.lt.1000)then
		   nam(1:1)=char(48)
		endif
c		fname = '/raid5/schumach/'//histnm(1:lfl-4)//'_'//nam//'.ascii'
		fname = prefix(1:lprefix)//histnm(1:lfl-4)//'_'//nam//'.ascii'
c                write(6,*)'New output file: ',prefix(1:lprefix),lprefix
c                write(6,*)'New output file: ',histnm(1:lfl-4),fl-4
c                write(6,*)'New output file: ',nam
		write(6,*)'New Monte Carlo output ',fname
		open(unit=2,file=fname,status='UNKNOWN')
c                call exit()
             end if
	ieventout = ieventout+1
	if(ieventout.gt.maxbosoutput)then  !open new file after <maxbosoutput> events
	   close(2)
	   ifile = ifile+1
	   ieventout = 1
	   write(nam,'(i4)'),ifile
	   if(ifile.lt.10)then
	      nam(3:3)=char(48)
	   endif
	   if(ifile.lt.100)then
	      nam(2:2)=char(48)
	   endif
	   if(ifile.lt.1000)then
	      nam(1:1)=char(48)
	   endif
	   fname = prefix(1:lprefix)//histnm(1:lfl-4)//'_'//nam//'.ascii'
	   write(6,*)'Opening raw data output file ',fname
	   open(unit=2,file=fname,status='UNKNOWN')
	endif
c
c       Repackage tracks in a format suitable for input to TXT2PART or the GlueX HDDM converter
c
	vz  = plist(9,2)    !take vertex from first track (the same for all tracks)
c	toff    = vz/30.0   !photon time offset in nanoseconds (old, incorrect)
	toff    = (-10.0 - vz)/30.0   !photon time offset in nanoseconds, for g11
	ebeam   = 4.023     !in GeV, a hack... for g11 data set
c
	icount = 0    !counts number of GSIM "known" particles in event
        iproton = 0   !counts the number of protons in the final state
	do i = 2,imark
c
c          Map Lund ID's to GEANT ID's
c          Note that GEANT doesn't do resonances, so some entries are "zero"
c          Also, put "0" for the ID of any particle we don't want GSIM to
c          decay for us, since we are producing the decay products ourselves.
c
	   igsimid(i) = -1
	   if(plist(6,i).eq.  1)igsimid(i) = 1 !0 !1   !photon
	   if(plist(6,i).eq. -1)igsimid(i) = 1 !0 !1   !photon
	   if(plist(6,i).eq. -7)igsimid(i) = 2   !positron
	   if(plist(6,i).eq.  7)igsimid(i) = 3   !electron
	   if(plist(6,i).eq.  2)igsimid(i) = 4   !neutrino
	   if(plist(6,i).eq. -9)igsimid(i) = 5   !mu+
	   if(plist(6,i).eq.  9)igsimid(i) = 6   !mu-
	   if(plist(6,i).eq. 23)igsimid(i) = 7 !0 !7   !pi0
	   if(plist(6,i).eq. 17)igsimid(i) = 8   !pi+
	   if(plist(6,i).eq.-17)igsimid(i) = 9   !pi-
	   if(plist(6,i).eq. 18)igsimid(i) = 11  !K+
	   if(plist(6,i).eq.-18)igsimid(i) = 12  !K-
	   if(plist(6,i).eq. 42)igsimid(i) = 13  !neutron
	   if(plist(6,i).eq. 41)igsimid(i) = 14  !proton
	   if(plist(6,i).eq.-41)igsimid(i) = 15  !anti-proton
	   if(plist(6,i).eq. 19)igsimid(i) = 0 !16  !K short (only! short, not long) 
	   if(plist(6,i).eq. 24)igsimid(i) = 0 !17  !eta
	   if(plist(6,i).eq. 57)igsimid(i) = 0 !18  !Lambda
	   if(plist(6,i).eq.-57)igsimid(i) = 0 !18  !Lambda
	   if(plist(6,i).eq. 43)igsimid(i) = 0 !19  !Sigma+
	   if(plist(6,i).eq. 44)igsimid(i) = 0 !20  !Sigma0
	   if(plist(6,i).eq.-44)igsimid(i) = 0 !20  !Sigma0
	   if(plist(6,i).eq. 45)igsimid(i) = 0 !21  !Sigma-
	   if(plist(6,i).eq. 46)igsimid(i) = 0 !??  !Xi- (Cascade)
	   if(plist(6,i).eq.-46)igsimid(i) = 0 !??  !Xi- (Cascade)
      	   if(plist(6,i).eq. 33)igsimid(i) = 0 !57  !rho0
	   if(plist(6,i).eq. 27)igsimid(i) = 0 !58  !rho+
	   if(plist(6,i).eq.-27)igsimid(i) = 0 !59  !rho-
	   if(plist(6,i).eq. 34)igsimid(i) = 0 !60  !omega meson
	   if(plist(6,i).eq. 36)igsimid(i) = 0 !61  !eta prime
	   if(plist(6,i).eq. 35)igsimid(i) = 0 !62  !phi

	   if(plist(6,i).eq. 58)igsimid(i) =  0  !Lambda(1405)
	   if(plist(6,i).eq. 59)igsimid(i) =  0  !Lambda(1520)
	   if(plist(6,i).eq. 56)igsimid(i) =  0  !Sigma+(1385)
	   if(plist(6,i).eq. 48)igsimid(i) =  0  !Sigma0(1385)
	   if(plist(6,i).eq. 60)igsimid(i) =  0  !Sigma-(1385)
	   if(plist(6,i).eq. 20)igsimid(i) =  0  !K*0
	   if(plist(6,i).eq. 21)igsimid(i) =  0  !K*+
	   if(plist(6,i).eq. 61)igsimid(i) =  0  !Delta++
	   if(plist(6,i).eq.-61)igsimid(i) =  0  !Delta++
	   if(plist(6,i).eq. 63)igsimid(i) =  0  !Delta+
	   if(plist(6,i).eq. 62)igsimid(i) =  0  !Delta0
	   if(plist(6,i).eq.-62)igsimid(i) =  0  !anti-Delta0
	   if(plist(6,i).eq. 70)igsimid(i) =  0  !special "box" particle
	   if(plist(6,i).eq. 74)igsimid(i) =  0  !f1(1285)
	   if(plist(6,i).eq. 78)igsimid(i) =  0  !f0(1380)
	   if(plist(6,i).eq. 91)igsimid(i) =  0  !Baryonium
	   if(plist(6,i).eq. 92)igsimid(i) =  0  !di-proton
	   if(plist(6,i).eq. 93)igsimid(i) =  0  !Lambda-anti-Lambda
	   if(plist(6,i).eq. 94)igsimid(i) =  0  !
	   if(plist(6,i).eq. 95)igsimid(i) =  0  !anti-Lambda - proton
	   if(plist(6,i).eq. 96)igsimid(i) =  0  !anti-Sigma0 - Sigma0
	   if(plist(6,i).eq. 97)igsimid(i) =  0  !anti-Delta0 - Delta0
	   if(plist(6,i).eq.101)igsimid(i) =  0  !anti-Sigma0 - Lambda
	   if(plist(6,i).eq.-101)igsimid(i) =  0 !Sigma0 - anti-Lambda
	   if(plist(6,i).eq. 99)igsimid(i) =  0  !Non-Such
	   if(igsimid(i).eq.-1)then
	      write(6,*)'Unknown Lund-->GSIM mapping for Lund ID:',plist(6,i)
c	      call exit()
	   endif
	   if(igsimid(i).gt.0)then
	      icount = icount + 1
	   endif
	   if(igsimid(i).eq.14)then !count protons
	      iproton = iproton + 1
              if(iproton.eq.1)ipmark1 = i
              if(iproton.eq.2)ipmark2 = i
	   endif
	enddo
c        write(6,*)'imark',imark,'icount',icount
c	do i = 2,imark
c           write(6,*)plist(6,i),igsimid(i)
c       enddo
c        
c       Special hack the the p-anti-p case which has 3 particles in the final state.
c       If there are two protons in the event, swap them in half of all events
c       in order to randomize their occurrance in the file that gets passes to HDGEANT
c        
c       This was done for the p-anti-p case with 3 particles in the final state.
c
c        write(6,*)'There are ',iproton,' protons in this event.',
c     1       ipmark1,ipmark2
        if(iproton.eq.2.and.icount.eq.3)then
           x = ran(iseed)
           if(x.lt.0.5)then
c              write(6,*)'Protons will be swapped'
              do j = 1,24
                 plisttemp = plist(j,ipmark1)
                 plist(j,ipmark1) = plist(j,ipmark2)
                 plist(j,ipmark2) = plisttemp
              enddo
           endif
        endif
c
c       Another special hack to Lambda-anti-Lambda events which have 5  particles in the final state.
c       Because HDGEANT's HDDM input file wants to "know" the full decay structure of the events,
c       we must put in a hack that always puts the proton that did not come from Lambda decay first.
c
        if(iproton.eq.2.and.icount.eq.5)then
c           write(6,*)'Lambda - anti - Lambda'
c           write(6,*)(plist(6,j), j=2,imark)
           if(plist(6,3).eq.93)then
c              write(6,*)'Found proton in imark 2 case' !nothing for more be done
           elseif(plist(6,3).eq.96)then   !(this is for "Mechanism 2" Kaon exchange)
c              write(6,*)'Found proton in imark 6 case'
c              write(6,*)'PLIST rearranged 6-->2, 2-->3 ... 5-->6'
              do j = 1,24
                 plisttemp  = plist(j,6)
                 plist(j,6) = plist(j,5)
                 plist(j,5) = plist(j,4)
                 plist(j,4) = plist(j,3)
                 plist(j,3) = plist(j,2)
                 plist(j,2) = plisttemp
              enddo
              igsimidtmp = igsimid(6)
              igsimid(6) = igsimid(5)
              igsimid(5) = igsimid(4)
              igsimid(4) = igsimid(3)
              igsimid(3) = igsimid(2)
              igsimid(2) = igsimidtmp
           else
c              write(6,*)
c     1        'This event is of another ilk.... '
c              call exit()
           endif
        endif
c     
c
c       Test for "NaN"
c       This is a hack that works in gfortran
c       (Back in the old days, when the program messed up here, 
c       it would simply crash...)
c
	do i = 2,imark
	   if(igsimid(i).gt.0)then
	      pgsim(1) = plist(5,i)
	      do k=1,3
		 pgsim(k+1) = plist(k,i)
	      enddo
	      vx  = plist(7,i)	!creation location for this track
	      vy  = plist(8,i)
	      vz  = plist(9,i)
c     
	      icrash = 0
              ichoke = 0
	      if(vx==vx)then
		 continue
	      else
		 icrash = 1
                 ichoke = ichoke +10
	      endif
	      if(vy==vy)then
		 continue
	      else
		 icrash = 1
                 ichoke = ichoke +100
	      endif
	      if(vz==vz)then
		 continue
	      else
		 icrash = 1
                 ichoke = ichoke +1000
	      endif
	      do j=1,4
		 if(pgsim(j)==pgsim(j))then
		    continue
		 else
		    icrash = 1
                    ichoke = ichoke +1000*j
		 endif
	      enddo
c
c             Also test for ridiculously large values
c
              toobig = 1.E4  !centimeters
              if(abs(vx) .gt. toobig)then
                 icrash = 1
                 ichoke = ichoke +100000
              endif
              if(abs(vy) .gt. toobig)then
                 icrash = 1
                 ichoke = ichoke +1000000
              endif
              if(abs(vy) .gt. toobig)then
                 icrash = 1
                 ichoke = ichoke +10000000
              endif
c
	      if(icrash.eq.1)then
		 write(6,*)'Aborting event output due to NaN status',ieventout, ichoke
                 write(6,*)vx,vy,vz,pgsim,plist(4,i),plist(6,i)
                 call prevent(plist,8,100)
		 ieventout = ieventout-1
c                 call exit()
		 return
	      endif
	   endif
        enddo
c
c       This event is OK, so write it out...
c       GEANT version
c        
	write(2,10)icount
 10	format(i2)
c        write(6,*)'Output stream photon energy:'
c	write(6,20)ppbeam,imechanism
	write(2,20)ppbeam,imechanism
 20	format(f10.6,i4)
c
	do i = 2,imark
	   if(igsimid(i).gt.0)then
	      pgsim(1) = plist(5,i)
	      do k=1,3
		 pgsim(k+1) = plist(k,i)
	      enddo
	      vx  = plist(7,i)	!creation location for this track
	      vy  = plist(8,i)
	      vz  = plist(9,i)
c
	      write(2,30)igsimid(i),pgsim
 30	      format(i2,4f10.6)
	      write(2,40)vx,vy,vz
 40	      format(3f11.6)
c
c	      write(6,*)plist(6,i),igsimid(i)
	   endif
	enddo
        return 
	end 

c********************************************************************************
c
c	Input from file in TXT2BOS format 
c
c       plist holds the tracks for this event
c       imark indexes the last track in the event
c
 	subroutine txt2bos_input(plist,imark)
	parameter (maxnumlis=60) !Maximum number of particles in the decay list
	character nam*3,fname*150
	real plist(24,maxnumlis),pgsim(4)
	integer igsimid(maxnumlis)
	logical ifirst
c 	data ifirst/.true./
	common /beam/pbeam,ppbeam,rmbeam,rmtarg,jcode,xsigma,ysigma,ncode
	common /eventiolocal/ifirst,ifile,ieventin
c
	character*256  prg_name, histnm, filenm, prefix
	common /run_control/n_arg,prg_name,filenm,histnm,prefix,
     1                      irestart,lfl,lprefix,maxbosoutput
	data rmprot/0.93827/  !proton rest mass
c       
c	On first call just open the data file
c
	if(ifirst)then
	   ifirst = .false.
	   ifile  = 1		!specifies starting number of input files
	   ieventin = 0
	   write(nam,'(i3)'),ifile !fancy internal write to convert # to ascii
	   if(ifile.lt.10)then
	      nam(2:2)=char(48)
	   endif
	   if(ifile.lt.100)then
	      nam(1:1)=char(48)
	   endif
c           fname = '/raid5/schumach/'//histnm(1:lfl-4)//'_'//nam//'.ascii'
	   fname = prefix(1:lprefix)//histnm(1:lfl-4)//'_'//nam//'.ascii'
	   write(6,*)'New TXT2BOS input ',fname
	   open(unit=2,file=fname,status='OLD')
	end if
	ieventin = ieventin+1
	if(ieventin.gt.maxbosoutput)then  !open new file after <maxbosoutput> events
	   close(2)
	   ifile = ifile+1
	   ieventout = 1
	   write(nam,'(i3)'),ifile
	   if(ifile.lt.10)then
	      nam(2:2)=char(48)
	   endif
	   if(ifile.lt.100)then
	      nam(1:1)=char(48)
	   endif
c	   fname = '/raid5/schumach/'//histnm(1:lfl-4)//'_'//nam//'.ascii'
	   fname = prefix(1:lprefix)//histnm(1:lfl-4)//'_'//nam//'.ascii'
	   write(6,*)'Opening TXT2BOS input ',fname
	   open(unit=2,file=fname,status='OLD')
	endif
c
	read(2,10)imark
 10	format(i2)
	read(2,20)ebeam,ppbeam,toff
 20	format(f8.6,2f10.6)
c	write(6,*)'> ',imark,ebeam,ppbeam,toff
	plist(1,1) = 0.0                !create pseudo-particle in initial state
	plist(2,1) = 0.0
	plist(3,1) = ppbeam             !total lab momentum
	plist(4,1) = sqrt(2.*ppbeam*rmprot + rmprot**2.)
	plist(5,1) = ppbeam+rmprot      !total lab energy, for photoproduction
	plist(6,1) = 999
	plist(22,1)= +1                 !charge
	pbeam      = ppbeam
	rmtarg     = rmprot
c
	do i = 2,imark+1
	   read(2,30)igsimid(i),pgsim
 30	   format(i2,4f10.6)
	   read(2,40)vx,vy,vz
 40	   format(3f11.6)
	   plist(5,i) = pgsim(1)       !total energy
	   do k=1,3
	      plist(k,i)    = pgsim(k+1)  !momentum
	      plist(k+12,i) = pgsim(k+1)  !momentum at decay point
	   enddo
	   rmass = sqrt(pgsim(1)**2.-(pgsim(2)**2.+pgsim(3)**2.+pgsim(4)**2.))
	   plist(4,i)  = rmass
	   plist(7,i)  = vx	!creation location for this track
	   plist(8,i)  = vy
	   plist(9,i)  = vz
	   plist(16,i) = rmass
	   plist(18,i) = vx	!decay location for this track
	   plist(19,i) = vy
	   plist(20,i) = vz
c
c          Map GSIM ID's to LUND ID's
c          ...yes, yes, this code is inefficient...
c
	   if(igsimid(i).lt.1 .or. igsimid(i).gt.14)then
	      write(6,*)'Unknown GSIM-->Lund mapping for particle ID:',igsimid(i)
	      call exit()
	   endif
	   if(igsimid(i).eq.  1)plist(6,i) = 1     !photon
	   if(igsimid(i).eq.  2)plist(6,i) = -7    !positron
	   if(igsimid(i).eq.  3)plist(6,i) = 7     !electron
	   if(igsimid(i).eq.  4)plist(6,i) = 2     !neutrino
	   if(igsimid(i).eq.  5)plist(6,i) = -9    !mu+
	   if(igsimid(i).eq.  6)plist(6,i) = 9     !mu-
	   if(igsimid(i).eq.  7)plist(6,i) = 23    !pi0
	   if(igsimid(i).eq.  8)plist(6,i) = 17    !pi+
	   if(igsimid(i).eq.  9)plist(6,i) = -17   !pi-
	   if(igsimid(i).eq. 11)plist(6,i) = 18    !K+
	   if(igsimid(i).eq. 12)plist(6,i) = -18   !K-
	   if(igsimid(i).eq. 13)plist(6,i) = 42    !neutron
	   if(igsimid(i).eq. 14)plist(6,i) = 41    !proton
	   plist(22,i) = 0                          !charge
	   if(plist(6,i).gt.1)then
	      plist(22,1) = +1	!charge
	   else
	      plist(22,1) = -1	!charge
	   endif
c	   write(6,*)'>>',i,pgsim,igsimid(i),int(plist(6,i))
	enddo
c
	return 
	end 

c*****************************************************************************
c	POLARIZE
c
c	Assign polarization directions to the strongly produced 
c	particles.  If there is no incident vector direction (decay
c	at rest), then pick a zero length polarization axis.
c
c	Introduce some rudimentary dynamics by making the size of the
c	polarization vector proportional to a simple function of theta,
c	where theta is the production angle of the particle.
c
c       PV :    input - 4-vector of the parent particle that decayed (lab frame)
c       PDECAY: input - 4-vector of the decay particle that is to be polarized (lab frame)
c	CTHCM:  input - cosine of center of mass angle at which particle
c			pdecay was generated.  Used to assign the 
c			"strength" of the polarization (toy model).
c	POLDECAY output - the direction along which the particle PDECAY
c			may have some polarization in strong decays
c
c     7-26-2018 - Rework this code so both decay particles are computed together
c                 Depending on some ad hoc rules, the polarization axis is selected
c     
	subroutine polarize(pv,pdecay1,pdecay2,lundid1,lundid2,pol1,pol2,cthcm)
	real pv(5),pdecay1(5),pdecay2(5),pol1(3),pol2(3),khat(5)
	logical ihackflag,ihavespinhalf
	common /stuff/poldir(4),xx,icount,ikeep,ihackflag,iattempt
c
c        Decide which polarization axis to use depending on an ad hoc rule:
c        For two Lambda's, use their mutual normal direction.
c        For one Lambda, use the normal to the parent and the Lambda.
c
	if(ihavespinhalf(lundid1).and.ihavespinhalf(lundid2))then
c           if(abs(lundid1).eq.57.and.abs(lundid2).eq.57)then
              call crossprod(pdecay1,pdecay2,khat,angle)
              iflag = 2   ! means polarize both particles along their normal
c           else
c              write(6,*)'Unknown ID combination',lundid1,lundid2
c              call exit()
c           endif
        else
           if(ihavespinhalf(lundid1))then
              call crossprod(pv,pdecay1,khat,angle) !Note: its wants the first one
              iflag = 11  !means polarize only particle 1 using parent and daughter
           elseif(ihavespinhalf(lundid2))then
              call crossprod(pv,pdecay2,khat,angle) !Note: its wants the first one
              iflag = 12
           else
              do i=1,3
                 pol1(i) =  0
                 pol2(i) =  0
                 poldir(i)   =  0
              enddo
              return
           endif
        endif
        rkhatnorm = absp(khat)
        if(rkhatnorm.eq.0)then
           rkhatnorm = 1.e36
        end if
c     
c          Trap bad cases...
c
        if(cthcm.gt.1.0)then
           write(6,*)'Polarize angle bad',cthcm
           cthcm = 1.0
        endif
c
        fact       = sqrt(1.-cthcm**2.)*(+.3 + 1.0*cthcm)  !mimic Kaon photoproduction (toy model)
        if(iflag.eq.2)then
           rkhatnormf = 1.0/rkhatnorm 
           do i=1,3
              pol1(i) =  khat(i)*rkhatnormf
              pol2(i) =  -pol1(i)
              poldir(i)    =  pol1(i) ! load common block
           enddo
        elseif(iflag.eq.11)then
           rkhatnormf = 1.0/rkhatnorm * fact
           do i=1,3
              pol1(i) =  khat(i)*rkhatnormf
              pol2(i) =  0.0
              poldir(i)    =  pol1(i) ! load common block
           enddo
        elseif(iflag.eq.12)then
           rkhatnormf = 1.0/rkhatnorm * fact
           do i=1,3
              pol2(i) =  khat(i)*rkhatnormf
              pol1(i) =  0.0
              poldir(i)    =  pol2(i) ! load common block
           enddo
        else
           write(6,*)'POLARIZE is very confused'
           call exit()
        endif
c
c           write(6,*)'POLARIZE:'
c           write(6,*)pol1
c           write(6,*)pol2
	return
	end

c
c	Gotta figure out which LUND id''s correspond to spin 1/2 particles
c
c       7-22-18 - Add the special cases of Lambda-anti-Lambda, and proton-anti-Lambda
c
      function ihavespinhalf(id)
	parameter (imax=14)
	integer ilist(imax)
	logical ihavespinhalf
	data ilist/999,57,41,42,43,44,45,46,47,49,94,7,2,9/
	it = abs(id)
	do i=1,imax
	   if(it.eq.ilist(i))then
	      ihavespinhalf = .true.
	      return
	   end if
	end do
	ihavespinhalf = .false.
	return	
	end

c****************************************************************************
c	WEAKDECAY
c	Compute decay which has a weak decay asymmetry, like for the 
c	Lambda, such that the direction of the decaying particles
c	in the rest frame of the decaying particle is correlated 
c	according to S-P wave interference via an asymmetry parameter
c       alpha.
c
c       Input: p0(5)   - the decaying particle (Hyperon)
c              polv(3) - the degree of polarization of the parent particle
c              lund    - ID of the parent particle
c
c       Output: p1(5), p2(5)     - products of the decay
c               pol1(3), pol2(3) - polarization vectors of the two decay products
c               cthcm            - cosine of the decay product with respect to the initial
c               polarization direction
c
c	R.A.Schumacher, CMU, 1991
c	7-31-95		make it more general than for just Lambda
c       7-04-17         sort out why anti-Lambda are not treated equally to the Lambdas
c       7-26-18         lots of testing for the Lambda anti-Lambda case
c      
	subroutine weakdecay(p0,p1,p2,polv,pol1,pol2,cthcm,lund,ierr)
	real p0(5),p1(5),p2(5),polv(3),rot(3),pol1(3),pol2(3)
	real p0tmp(5),p1tmp(5),p2tmp(5),dtmp(5),beta(3),polvnorm(3)
	logical ihackflag
	data alphalam/.642/	!asymmetry parameter for Lambdas
	data twopi /6.283185308/
	data piover2/1.5707963/
c
	common /monte/ iseed,xflat,xnormal,xbw
	common /stuff/poldir(4),xx,icount,ikeep,ihackflag,iattempt
c
	if(lund.eq.57)then	!Lambda case
	   alpha = alphalam 
	elseif(lund.eq.-57)then
           alpha = alphalam !antiparticle has the opposite asymmetry, BUT polarization POLV is flipped already
        else
	   write(6,*)'Decay asymmetry is not specified',lund
	   alpha = 1.0
	end if
c        write(6,*)'Weakdecay 0:',lund,alpha,alphalam
c        write(6,*)'            ',polv
c     
c	pick decay angle of Lambda relative to polarization direction  POLV
c
	polnorm = absp(polv)
        if(polnorm.gt.1.0)then
c           write(6,*)'Weakdecay polarization factor is too big',polnorm
           polnorm = 1.0
        endif
c	a = 1.0 		!as a test, assume 100% polarization of Lambda with 100% asymmetry
c	a = alpha * 1.0		!as a test, assume 100% polarization of Lambda
	a = alpha * polnorm	!defines degree of polarization along POLV
c
	y = ran(iseed)
c	x = 2.*ran(iseed) - 1.  !flat distribution in cosine theta
	x = (sqrt((a-1)*(a-1)+4.0*a*y) - 1.0)/a !linearly weighted cosine of polar angle wrt polv

c        x = .0        !!!!!!!!!!!!!   @@@@@@@@@@@@@@@@
        
c	write(6,*)x,y,a,polnorm
        if(x.gt. 1.0)x= 1.0
	if(x.lt.-1.0)x=-1.0
        theta = acos(x)         !polar decay angle wrt POLV in radians
	phik  = twopi*ran(iseed)		!phi angle wrt POLV

c        PHIK = 0.0  !3./2. * 3.1415 ! twopi/4.   !!!!!!!!!!! @@@@@@@@@@@@
c
c        write(6,*)'Weakdecay 1  lund p0:',lund,p0
c        write(6,*)'Weakdecay 2  alpha,polv,polnorm,theta,phik;',
c     1         alpha,polv,polnorm,theta,phik
c
c	Generate the direction vector for the decay nucleon
c	Start out in POLV direction, rotate by polar angle THETA in the
c	POLV/P0 plane, then rotate about POLV by PHIK radians
c	DTMP will be the decay direction of the product in the c.m. frame
c 	of the decaying particle P0.
c
        if(polnorm.eq.0)polnorm = 1.0
	do i=1,3
	   dtmp(i)     = polv(i)
	   polvnorm(i) = polv(i)/polnorm
	end do
	call crossprod(polv,p0,rot,angle)
	if(abs(angle-piover2).gt.0.001)then
           temp = angle-piover2
c           write(6,*)'POLV is not normal to P0',temp
	end if
c	write(6,*)'POLV ',polv
c	write(6,*)'P0   ',p0
c	write(6,*)'ANGLE',angle
c	write(6,*)'ROT  ',rot
c	write(6,*)'THETA,PHIK,X',theta,phik,x

	call carrot(dtmp,rot,theta)
	call carrot(dtmp,polv,phik)
c
c	Get the polar and azimuthal angles of the decay with respect to P0
c
	cthcm = dot(dtmp,p0)/(absp(dtmp)*absp(p0))
c
c	Velocity of lab frame relative to CM
c
	e0=p0(5)
        do i=1,3
	   p0tmp(i)=0.
	   beta(i)=-p0(i)/e0
	end do
	p0tmp(4)=p0(4)
	p1tmp(4)=p1(4)
	p2tmp(4)=p2(4)
c
c	Let the decaying particle decay in its rest frame with particle 2 
c	in heading in direction DTMP
c
	call kin12(p0tmp,p1tmp,p2tmp,dtmp,ierr)
c	write(6,*)'WEAK GENERATED    : polvnorm',polvnorm
c	write(6,*)'WEAK GENERATED    : p0(hyp) ',p0
c	write(6,*)'WEAK GENERATED    : pproton ',p1tmp
c	write(6,*)'WEAK GENERATED    : betahyp ',beta
c	
c	Boost the decay products back to the original (often lab) frame
c
	call boost(p1tmp,beta,p1)
	call boost(p2tmp,beta,p2)
c
c	The decay products of this particle are not themselves polarized
c	Thus, we are ignoring the longitudinal polarization of the
c	nucleons from Lambda decay, for example.
c
	do i=1,3
	   pol1(i) = 0
	   pol2(i) = 0
	end do
c        write(6,*)'Weakdecay p0',p0
c        write(6,*)'Weakdecay p1',p1
c        write(6,*)'Weakdecay p2',p2
	return
	end

c****************************************************************************
c	PRODUCT_STUFFIT
c
c       Some previous routine set up the event parameters and selected the cm
c       energy and decay angle.  This routine selects the phi angle and manipulates
c       the 4-vectors of the decay particles.
c
c     Input:
c       P0 - parent particle 4-vector (px,py,pz,M,E) in the LAB frame
c       ID1 -Lund ID of the first  decay particle
c       ID2 -Lund ID of the second decay particle
c       CTHCM - cosine of decay angle of particle P1 in the rest frame of particle P0,
c               which is the overall CM frame
c
c       Output:
c       P1 - decay particle 4-vector (px,py,pz,M,E) in the LAB frame
c       P2 - decay particle 4-vector (px,py,pz,M,E) in the LAB frame
c       POL1 - polarization vector of P1
c       POL2 - polarization vector of P2
c     
c       R.A.Schumacher, CMU
c       11-2006
c       8-3-16 put in a bit to let a "dibaryon" decay (ID 91 or 92)
c       6-28-17 rename the subroutine from "meson_stuffit" to "product_stuffit"
c               and genealize some of the variable names      
c
	subroutine product_stuffit(p0,p1,p2,pol1,pol2,id1,id2,
     1             cthcm,lund,ierr)
	real p0(5),p1(5),p2(5),rot(5),pol1(3),pol2(3)
	real p0tmp(5),p1tmp(5),p2tmp(5),dtmp(5),beta(3),ex(5)
	data twopi,piover2 /6.283185308,1.570796327/
c
	common /monte/ iseed,xflat,xnormal,xbw
c
c       Hackery to check which particle is a LUND meson (or a decaying dibaryon)
c
	idummy = lund !fool the compiler...
        iproduct = 2  !Make sure .def file has recoil product listed first
c     
c     Test parameter passing
c
c        write(6,*)'Product_stuffit',id1,id2,cthcm
c        write(6,*)p0
c        write(6,*)p1
c        write(6,*)p2
c
	theta = acos(cthcm)	        !polar decay angle wrt beam
	phik  = twopi*ran(iseed)	!phi angle wrt x-axis

ccc        phik  = 0   ! Test Hack
c
c	Generate the direction vector for the first decay particle
c	Start out in P0 (beam) direction, rotate by polar angle THETA in the
c	EX/P0 plane, then rotate about P0 by PHIK radians
c	DTMP will be the decay direction of the parent particle in the rest frame
c 	of the decaying particle P0.
c
	ex(1) = 1
	ex(2) = 0
	ex(3) = 0
	ex(4) = 0
	ex(5) = 0
	do i=1,5
           dtmp(i) = p0(i)
	end do
	call crossprod(p0,ex,rot,angle)	!defines a rotation axis ROT
	call carrot(dtmp,rot,theta)	!polar angle rotation
	call carrot(dtmp,p0,phik)	!azimuthal angle rotation
c
c	Get the polar and azimuthal angles of the decay with respect to P0
c
	cthcm = dot(dtmp,p0)/(absp(dtmp)*absp(p0))
	if(abs(cthcm-cos(theta)).gt.0.001)then
		write(6,*)'Deep trouble:',cthcm,cos(theta)
	end if
c
c	Velocity of parent (lab) frame relative to CM
c
	do i=1,3
	   p0tmp(i) = 0.
	   beta(i)  = -p0(i)/p0(5)
	enddo
        p0tmp(5)=p0(4) !at rest, so energy = mass
	p0tmp(4)=p0(4)
	p1tmp(4)=p1(4)
	p2tmp(4)=p2(4)
c
c	Let the decaying object decay in its rest frame, with the first
c	decay particle heading in direction DTMP
c
	if(iproduct.eq.1)then
           call kin12(p0tmp,p1tmp,p2tmp,dtmp,ierr)
	else
c           write(6,*)'calling kin12',p0tmp
           call kin12(p0tmp,p2tmp,p1tmp,dtmp,ierr)
c           write(6,*)'return  kin12',p1tmp
c           write(6,*)'return  kin12',p2tmp
        end if
c	
c	Boost the decay products back to the lab frame
c
	call boost(p1tmp,beta,p1)
	call boost(p2tmp,beta,p2)
c
c       Don't give the decay products any polarization.
c       This will be done by a separate call from the parent routine.
c
	do i=1,3
	   pol1(i) = 0
	   pol2(i) = 0
	enddo
	return
	end

c****************************************************************************
c       PRODUCT_ROTATE
c           
c       We have Particle 0 decaying into Particle 1 and Particle 23.
c       We also have the center-of-mass decay polar angle of Particle 23 that
c       gives us the desired momentum transfer (computed by previous code).
c       This routine puts Particle 23 on a line at the desired polar angle.
c       Particle 23 is composed of Particle 2 and Particle 3, and we already
c       have all the 4-vectors in the CM frame.   As a last step, boost
c       particles 1, 2, 3, and 23 to the lab frame.
c
c     Input:
c       P0    - parent particle 4-vector (px,py,pz,M,E) in the LAB frame
c       COSTHETACM - cosine of decay angle of particle P1 in the rest frame of particle P0,
c               which is the overall CM frame
c       P1    - the recoiling particle decay 4-vector in the zero-momemtum frame (not lab frame)
c       P23   - the created particle 4-vector in the zero-momentum frame that we want to rotate to COSTHETACM 
c       P2    - the decay product of particle P23 in the zero-momentum frame (not lab frame)
c       P3    - the decay product of particle P23 in the zero-momentum frame (not lab frame)

c     Output: (modified 4-vectors of the produced particles, all in the LAB frame
c       P1  - decay particle 4-vector (px,py,pz,M,E) in the LAB frame ("recoil')
c       P23 - decay of the composite particle 4-vector (px,py,pz,M,E) in the LAB frame
c       P2  - decay particle 4-vector (px,py,pz,M,E) in the LAB frame (decay of Particle 23)
c       P3  - decay particle 4-vector (px,py,pz,M,E) in the LAB frame (decay of Particle 23)
c     
c       R.A.Schumacher, CMU
c       11-2006
c       6-28-17 rename the subroutine from "meson_stuffit" to "product_stuffit"
c               and generalize some of the variable names      
c       7-19-19 create this program for the case of quasi-three-body decays
c      
	subroutine product_rotate(p0,p1,p2,p3,p23,costhetacm,ierr)
	real p0(5),p1(5),p2(5),p3(5),p23(5),rot(5)
	real p0tmp(5),p1tmp(5),p2tmp(5),p3tmp(5),p23tmp(5),dtmp(5)
        real beta(3),ex(5),ze(5)
	data twopi,piover2 /6.283185308,1.570796327/
        data ex/1,0,0,0,0/      !the X axis
        data ze/0,0,1,0,0/      !the Z (beam) axis
	data rade /57.29577951/ !convert radians to degrees
c     
	common /monte/ iseed,xflat,xnormal,xbw
c
c       Hackery to check which particle is a LUND meson (or a decaying dibaryon)
c
	idummy = lund !fool the compiler...
        iproduct = 2  !Make sure .def file has recoil product listed first
        ierr = 0
c     
        theta = acos(costhetacm) !polar angle wrt beam that results in correct t distribution
	phik  = twopi*ran(iseed)	!phi angle wrt x-axis
c        phik  = 0   ! Test Hack
c
c     write(6,*)' '
c        write(6,*)'Product_rotate',costhetacm,theta*rade
c        write(6,*)'P0   lab',p0
c        write(6,*)'P1    cm',p1
c        write(6,*)'P2    cm',p2
c        write(6,*)'P3    cm',p3
c        write(6,*)'P23   cm',p23
c
c       Rotate particle P23 so that it points along the beam axis.  Rotate
c       all the other particles by the same amount, so the whole system is
c       oriented with the intermediate produced particle P23 along the beam.
c       This is still in the CM frame of the system.
c
        call crossprod(p23,ze,rot,angle)
        call carrot(p23,rot,angle)
        call carrot(p1,rot,angle)
        call carrot(p2,rot,angle)
        call carrot(p3,rot,angle)
c
c        write(6,*)'Rotate P23 to zero degrees',angle,angle*rade
c        write(6,*)'P0   lab',p0
c        write(6,*)'P1    cm',p1
c        write(6,*)'P2    cm',p2
c        write(6,*)'P3    cm',p3
c        write(6,*)'P23   cm',p23
c     
c       Rotate the system such that P23 angle results in the desired t value
c
	call carrot(p23,ex,theta)	!polar angle rotation
	call carrot(p1 ,ex,theta)	!polar angle rotation
	call carrot(p2 ,ex,theta)	!polar angle rotation
	call carrot(p3 ,ex,theta)	!polar angle rotation
c        write(6,*)'Rotate system for desired t value',theta,theta*rade
c        write(6,*)'P0   lab',p0
c        write(6,*)'P1    cm',p1
c        write(6,*)'P2    cm',p2
c        write(6,*)'P3    cm',p3
c        write(6,*)'P23   cm',p23
c
c       We rotated about the X axis to get the desired polar angle.
c       Now rotate by a random amount around the Z (beam) axis to scramble
c       the azimuthal information.
c       
	call carrot(p23,p0,phik) !azimuthal angle rotation
	call carrot(p1 ,p0,phik) !azimuthal angle rotation
	call carrot(p2 ,p0,phik) !azimuthal angle rotation
	call carrot(p3 ,p0,phik) !azimuthal angle rotation
c        write(6,*)'Rotate azimuthally',phik,phik*rade
c        write(6,*)'P0   lab',p0
c        write(6,*)'P1    cm',p1
c        write(6,*)'P2    cm',p2
c        write(6,*)'P3    cm',p3
c        write(6,*)'P23   cm',p23
c       Check that mometum is still conserved.
c        px = p1(1)+p2(1)+p3(1) 
c        py = p1(2)+p2(2)+p3(2) 
c        pz = p1(3)+p2(3)+p3(3)
c        pxx = p1(1) + p23(1)
c        pyy = p1(2) + p23(2)
c        pzz = p1(3) + p23(3)
c        write(6,*)'Momentum',px,py,pz
c        write(6,*)'Momentum',pxx,pyy,pzz
c     
c	Check the polar angle of P23 with respect to beam axis
c
	cthcm = dot(p23,ze)/(absp(p23)*absp(ze))
	if(abs(cthcm-cos(theta)).gt.0.001)then
c           write(6,*)'Trouble in rotating system:',
c     1         cthcm,cos(theta),costhetacm
           ierr = 1
           return
        end if
        return
        end

c*************************************************************************
c     BOOST_TO_LAB
c     Ad hoc porgram to boost the set of particles that have been
c     rotated in the CM frame back to the Lab frame.
c     Input:
c       P0    - parent particle 4-vector (px,py,pz,M,E) in the LAB frame
c       P1    - the recoiling particle decay 4-vector in the zero-momemtum frame (not lab frame)
c       P23   - the created particle 4-vector in the zero-momentum frame that we want to rotate to COSTHETACM 
c       P2    - the decay product of particle P23 in the zero-momentum frame (not lab frame)
c       P3    - the decay product of particle P23 in the zero-momentum frame (not lab frame)

c     Output: (modified 4-vectors of the produced particles, all in the LAB frame
c       P1  - decay particle 4-vector (px,py,pz,M,E) in the LAB frame ("recoil')
c       P23 - decay of the composite particle 4-vector (px,py,pz,M,E) in the LAB frame
c       P2  - decay particle 4-vector (px,py,pz,M,E) in the LAB frame (decay of Particle 23)
c       P3  - decay particle 4-vector (px,py,pz,M,E) in the LAB frame (decay of Particle 23)
c     
c       R.A.Schumacher, CMU
c       12-2019
c        
      subroutine boost_to_lab(p0,p1,p2,p3,p23,ierr)
      real p0(5),p1(5),p2(5),p3(5),p23(5)
      real p0tmp(5),p1tmp(5),p2tmp(5),p3tmp(5),p23tmp(5),dtmp(5)
      real beta(3),ex(5),ze(5)
c
c     Velocity of laboratory frame relative to CM
c     (Note that P0 has been in the lab frame throughout this routine.)
c        
      do i=1,3
         beta(i)  = -p0(i)/p0(5)
      enddo
c	
c     Boost the particles fro the CM to the lab frame.
c     Note that we are engaging in the dubious practice of
c     modifying the input variables to the program to make the
c     output variables.
c
      do i=1,5
         p1tmp(i) = p1(i)
         p2tmp(i) = p2(i)
         p3tmp(i) = p3(i)
         p23tmp(i) = p23(i)
      enddo
      call boost(p1tmp,beta,p1)
      call boost(p2tmp,beta,p2)
      call boost(p3tmp,beta,p3)
      call boost(p23tmp,beta,p23)
c        write(6,*)'Boost to the lab frame',beta
c        write(6,*)'P0   lab',p0
c        write(6,*)'P1   lab',p1
c        write(6,*)'P2   lab',p2
c        write(6,*)'P3   lab',p3
c        write(6,*)'P23  lab',p23
cc       Check that mometum is still conserved.
c        px = p1(1)+p2(1)+p3(1) 
c        py = p1(2)+p2(2)+p3(2) 
c        pz = p1(3)+p2(3)+p3(3)
c        pxx = p1(1) + p23(1)
c        pyy = p1(2) + p23(2)
c        pzz = p1(3) + p23(3)
c        write(6,*)'Momentum',px,py,pz
c        write(6,*)'Momentum',pxx,pyy,pzz
      return
      end

c****************************************************************************
c	Put the primary vertex somewhere in the target
c
c	Rectangular geometry.  SIZE1,2,3,4 are the X,Y,Zmin,Zmax full widths
c	of the target.
c
c       Modified:
c       12-6-2005 - displaced vertex added
c        7-5-2018 - fix up code so it works for GlueX       
c
	subroutine getvertex(r0)
	real r0(3),r1(3),direc(3)
	common /monte/ iseed,xflat,xnormal,xbw
	common /beam/pbeam,ppbeam,rmbeam,rmtarg,jcode,xsigma,ysigma,ncode
	common /target/icode,size1,size2,size3,size4,maxtstep
c
c        write(6,*)'Vertex:'
c        write(6,*)jcode,xsigma,ysigma,size1,size2,size3,size4
c
	if(jcode.ne.1)then
	   r0(1) = 0.0
	   r0(2) = 0.0
	   r0(3) = 0.0
	else
	   ibad = 0
 10	   ibad = ibad + 1
	   if(ibad.gt.10)then
	      write(6,*)'Beam spot is really bad'
	      call exit()
	   endif
c	   call rannor(x,y)
	   call rannor(x,q)
	   call rannor(y,q)
	   x = x*xsigma
	   y = y*ysigma
	   z = size3 + (size4-size3) * ran(iseed)
	   r0(1) = x
	   r0(2) = y
	   r0(3) = z
	   r1(1) = x
	   r1(2) = y
	   r1(3) = z
	   iflag = 0		!says vertex is inside the target...
	   rmag = sqrt(x*x+y*y+z*z)
	   if(rmag.gt.0)then
	      do i=1,3
		 direc(i) = -r0(i)/rmag
	      enddo
	   else
	      do i=1,3
		 direc(i) = -1.
	      enddo
	   endif
	   call checkout(r0,r1,direc,deltad,iflag) !...or is it outside?
	   if(iflag.ge.1)then
	      write(6,*)'Vertex outside target ',ibad,r0,r1,iflag
	      goto 10
	   end if
	end if
	return
	end

c****************************************************************************
c       DSDOMAKE
c
c       Compute the cross section from the Legendre coefficents
c
	function dsdomake(x,ipar,cval)
	real cval(5),p(5)
c
	if(ipar.eq.0)then
	   dsdomake=0.
	   return
	endif
	p(1) = 1.		! L = 0
	p(2) = x		! L = 1
	do i=3,ipar		! L = 2,3...
	   j = i-1
	   p(i) = (1./j)*((2.*j-1.)*x*p(i-1)-(j-1)*p(i-2))
	enddo
c       
	dsdo = cval(1)*p(1) 
	do i=2,ipar
	   dsdo = dsdo + cval(i)*cval(1)*p(i)
	enddo
	dsdomake = dsdo
	return
	end

c**************************************************************************
c       BREMSSHRAHLUNG
c
c       Select a photon energy according to a bremsstrahlung spectrum shape
c       (Box approx)
c       August 2006
c       Modified: 10-17-2006 - fix it!
c                  3-31-2014 - change (fix?) it according to 4-11-89 notes
c
	subroutine bremsstrahlung(pbeammin,pbeammax,ppbeam)
	common /bremlocal/ ifirst,pendpoint,pdelta,xmax,ratio,scale
	common /monte/ iseed,xflat,xnormal,xbw
	data ifirst/1/
c
	if(ifirst.eq.1)then
	   ifirst = 0
	   pendpoint = pbeammax/0.95 !endpoint specific for CLAS tagger
	   pdelta = pbeammax-pbeammin 
	   fracmin  = pbeammin/pendpoint
           ratio    = 2.*pendpoint/.000511
           xx       = pbeammin/pendpoint
           scale    = (1 - xx + 0.75*xx*xx)*(log(ratio*(1/xx - 1.))-0.5)
           scale    = scale/pbeammin
c	   xmax     = 1.0/fracmin  - 1.0
	   if(pdelta.le.0.0)then
	      write(6,*)'Big trouble...'
	      call exit()
	   endif
	endif
c 10	ppbeam   = pbeammin + ran(iseed)*pdelta
c	fraction = ppbeam/pendpoint
c	xcurve   = 1.0/fraction - 1.0
c	xtest    = xmax * ran(iseed)
c	if(xtest.le.xcurve)then
c 10	ppbeam   = pbeammin + ran(iseed)*pdelta
c	ycurve   = pbeammin/ppbeam   !the extreme approxiation of a 1/k spectrum
c	ytest    = ran(iseed)
 10	ppbeam   = pbeammin + ran(iseed)*pdelta
        xx       = ppbeam/pendpoint
        factor   = (1 - xx + 0.75*xx*xx)*(log(ratio*(1/xx - 1.))-0.5)
        factor   = factor/ppbeam
        ycurve   = factor/scale
c	ycurve   = pbeammin/ppbeam   !the extreme approxiation of a 1/k spectrum
	ytest    = ran(iseed)
        if(ytest.le.ycurve)then
c	   write(6,15)'Got one:',ppbeam,ratio,scale,xx,factor,ycurve,ytest,
c     1                           pbeammin,pbeammax
c 15	   format(1x,a,10f8.3)
	else
c	   write(6,15)'Failed :',ppbeam,ratio,scale,xx,factor,ycurve,ytest,
c     1                           pbeammin,pbeammax
	   goto 10
	endif
	return
	end

c**************************************************************************
c      Decide whether to keep an event with momentum transfer -t according
c     to the desired t-slope.
c
c     The algorithm assumes the parent t distribution is flat, which is not
c     strictly true, resulting in an effective t distribution skewed slightly
c     to larger values of t
c
c     Input:
c      PPBEAM   - beam momentum for computing W
c      RMTARGET - target rest mass for computing W
c      TT       - t of the vertex of interest
c      BSLOPE2  - exponential slope we require at this vertex
c      IYESNO   - flag to keep or reject this candidate
c
c      Output:
c      IYESNO - flag to keep or reject this candidate
c
c      R.A.Sch,  December 2019
c      
      subroutine doublereggepick(ppbeam,rmtarget,tt,bslope2,iyesno,ierr)
      data ifirst/1/
      data scale/1000./
      common /doublelocal/ifirst
      common /monte/ iseed,xflat,xnormal,xbw
c
      if(ifirst.eq.1)then
         ifirst=0
         ierr = 0
         write(6,*)'First call to DOUBLEREGGEPICK'
      endif
      ierr = 0
      ss = 2.0*ppbeam*rmtarget + rmtarget**2.
c      snot = 1./bslope2
c      snot        = 1./0.25     !For Pomeron trajectory 
      snot        = 1.0          !We have no reliable way to estimate this, so pick a default value
c      probdensity = scale * exp(2.*bslope2*tt) !key statement...
      probdensity = scale * exp(2.*bslope2*log(ss/snot)*tt) !key statement...
c     
c     Throw the dart
c
      ytest = scale*ran(iseed) 
      if(ytest.lt.probdensity)then
         iyesno = 1
c         write(6,*)'Double:',sqrt(ss),ss,snot,tt,bslope2,probdensity,ytest,iyesno
      else
         iyesno = 0
      endif
      return
      end
      
c***************************************************************************
c       Model the cross section as a rough Regge-like, i.e. t-channel-
c       dominant distribution.
c
c       Input: ppbeam        = previously selected beam momentum
c              ppbeammax     = maximum beam (photon) energy
c              wscale        = W scaling factor (not used in this version of the code)
c              bslope        = Regge slope parameter
c              sminf         = minimum fraction of largest cross section to pick
c              rmproduce     = mass of the object produced in the final state 
c              rmrecoil      = mass of system recoiling against the produced object, in GeV
c              rmtarget      = mass of target
c       Output:cosinethetacm = c.m. production angle selected with "t-channel" weighting
c              iregge        = error flag
c
c       The formula used is essentially Eq 4.86 from Perkins, 2nd Ed. with
c       invented coefficients.
c
c       This version allows two different values of the b-slope, which we may need when
c       studying the anti-Lambda Lambda channel.
c
c       R.A.Sch. June 2017
c       7-25-2017 - various reworkings to make it suitable for Lambda -anti-Lambda      
c
	subroutine reggepick(ppbeam,pbeammax,wscale,bslope,sminf,rmtarget,
     1                        rmproduce,rmrecoil,cosinethetacm,iregge)
	data rmg/0.0/		!photon mass  (use this term so the equations look symmetric)
	data ifirst/1/  
	common /monte/ iseed,xflat,xnormal,xbw
	common /reggelocal/ifirst,snot,dsigmadtnot,dsdtmax,dsdtmin,bslopelocal
c
c     The problem with s_not, the Regge scale parameter
c     See Donnachie and Landshoff, hep-ph/9912312 (12 Dec 1999) = Phys. Lett B, 478, 146, (2000)
c
c        wscale = 3.17           !proton plus 2 Lambdas
c        snot        = 1./bslope   !Their recipe in general, which leads to bad numerics for us
c        snot        = 1./0.25     !For Pomeron trajectory 
        snot        = 1.           !We have no way to estimate this reasonably, so pick a default value
c        snot        = 1./(wscale*wscale)
c        wmax        = sqrt(2.*pbeammax*rmtarget + rmtarget**2.)
c        smax        = wmax*wmax
c        wthreshold  = 2. * 0.93827 + rmrecoil !hack for tests
c
	if(ifirst.eq.1 .or. bslope.ne.bslopelocal)then
	   dsigmadtnot = 1000.0       !arbitrary scale
	   dsdtmax     = dsigmadtnot
	   dsdtmin     = dsigmadtnot*sminf !set "floor" for production cross section
           if(ifirst.eq.1)then
              ifirst = 0
              write(6,*)' ' 
              write(6,*)'Regge exchange parameters:' 
              write(6,5)bslope,snot,wscale,dsdtmax,dsdtmin
 5            format(1x,
     1         't-slope (b)                   = ',f8.3,' GeV^-2',/,1x,
     1         'Scale factor s_not            = ',f8.3,' GeV^2',/,1x,
     1         'W scale (if used)             = ',f8.3,' GeV',/,1x,
     1         'D(t=0)                        = ',f8.3,' microbarn/GeV^2',/,1x,
     1         't-independent minimum         = ',f8.3,' microbarn/GeV^2')
              write(6,*)' '
           endif
           bslopelocal = bslope
c          call exit()
	end if
c     
	iregge    = 0
	ww        = sqrt(2.*ppbeam*rmtarget + rmtarget**2.)
c 
c	if(bslope.ne.bslopelocal)then
c	   write(6,*)'b-slope changed',bslope,bslopelocal
c          write(6,*)
c     1         'Cannot allow multiple decays with Regge-type behavior'
c           write(6,*)ppbeam,rmtarget
c	   iregge = -1
c	   call exit()
c	endif
	eg           = ppbeam
	egcm         = (ww*ww - (rmtarget**2. - rmg**2.    ))/(2.*ww)
	eproducecm   = (ww*ww - (rmrecoil**2. - rmproduce**2.))/(2.*ww)
	erecoilcm    = (ww*ww + (rmrecoil**2. - rmproduce**2.))/(2.*ww)
	pgcm         = egcm
	pproducecm   = eproducecm*eproducecm - rmproduce*rmproduce
        if(pproducecm.lt.0.0 .and. pproducecm.gt.-2.0e-5)then
c           write(6,*)'Reggepick round-off error fix:',pproducecm
           pproducecm = 1.0e-6     ! put it right at threshold
        endif
        if(pproducecm.ge.1.0e-8)then
           pproducecm = sqrt(pproducecm)
        else
           write(6,*)'CM momentum incorrect:',eg,ww,rmtarget,pproducecm,
     1                              egcm,eproducecm,erecoilcm,
     2                              rmproduce,rmrecoil         
           iregge = -1
           return
        endif
	tmin = (rmproduce*rmproduce + rmg*rmg) - 
     1         2.*(egcm*eproducecm - pgcm*pproducecm)
	tmax = (rmproduce*rmproduce + rmg*rmg) - 
     1         2.*(egcm*eproducecm + pgcm*pproducecm)
	ss   = ww**2.
c        write(6,*)'>>>',ww,eg,egcm,eproducecm,erecoilcm,pproducecm,tmin,tmax
c        write(6,*)'   ',rmproduce,rmrecoil
c
c       Top of loop to find a t value, the associated 'cross section',
c       and then throw darts to find whether we should accept this 
c       Monte Carlo value of cosine(theta) of the production angle
c
	iflag  = 0
c
 10	iregge = iregge + 1
	if(iregge.ge.2000)then
	   iflag = 1
	endif
	tt = tmin + ran(iseed)*(tmax-tmin)
c	dsdt = dsdtmax * exp(2.*bslope*log(ss/snot)*tt)  !key statement...
	dsdt = dsigmadtnot * exp(2.*bslope*log(ss/snot)*tt)  !key statement...
c        write(6,*)'>>>>',dsdt,dsdtmax,bslope,tt,ss,snot
	dsdt = max(dsdt,dsdtmin)
	if(dsdt.gt.dsdtmax)then
	   write(6,*)'Raising dsdtmax:',dsdt,dsdtmax,bslope,tt,ss,snot
           dsdtmax = dsdt
	endif
	xtest = ran(iseed)*dsdtmax
c	write(6,*)'>',iregge,tt,dsdt,xtest
	if(xtest.gt.dsdt)goto 10
	cosinethetacm = (tt-(rmproduce*rmproduce + rmg*rmg)+2.*(egcm*eproducecm))/
     1                                              (2.*pgcm*pproducecm)
	if(iflag.eq.1)then
	   write(6,*)'Inefficient t-channel generation:', 
     1                                            iregge,tmin,tmax
	endif

c        cosinethetacm = 0.7071  !test value 
        if(cosinethetacm==cosinethetacm)then   !test for NaN
           continue
        else
           write(6,*)'How on Earth...?'
           write(6,*)ww,wscale,eg,egcm,eproducecm,erecoilcm,pgcm
           write(6,*)tt,rmproduce,rmg,pproducecm
           write(6,*)tmin,tmax
        endif
        if(cosinethetacm.gt.1.0)then   !yes, this can happen
           write(6,*)'REGGEPICK selected cosine theta =',cosinethetacm
           cosinethetacm = 1.0
        endif
        if(cosinethetacm.lt.-1.0)then
           write(6,*)'REGGEPICK selected cosine theta =',cosinethetacm
           cosinethetacm = -1.0
        endif
c
c       write(6,*)'Returning'
c	write(6,*)'>',iregge,cosinethetacm,tt,dsdt,xtest,dsdtmin
c	write(6,*)' ',tmin,tmax
c
	return
	end

c***************************************************************************
c       Model the cross section for as a simple exponential slope.
c       It is independent of W, though we compute it for completeness 
c
c       Used, for example, for rho photoproduction a 
c       la Ballam et al., PRD 5, 545 (1972).
c
c       Input: ppbeam        = previously selected beam momentum
c              ppbeammax     = maximum beam (photon) energy
c              bslope        = slope parameter
c              sminf         = minimum fraction of largest cross section to pick
c              rmmeson       = mass of final state meson in GeV
c              rmrecoil      = mass of system recoiling against meson in GeV
c       Output:cosinethetacm = c.m. production angle selected with "t-channel" weighting
c              iregge        = error flag 
c
c       R.A.Sch. April 2014
c
	subroutine exponentialpick(ppbeam,pbeammax,bslope,sminf,rmtarget,
     1                        rmmeson,rmrecoil,cosinethetacm,iregge)
	data rmg/0.0/		!photon mass  (use it so equations look symmetric)
	common /monte/ iseed,xflat,xnormal,xbw
	common /explocal/ifirst,snot,dsigmadtnot,dsdtmax,dsdtmin,wthreshold
	data ifirst/1/  
c       
	if(ifirst.eq.1)then
	   ifirst = 0
	   dsigmadtnot = 1.0       !arbitrary scale
	   wmax        = sqrt(2.*pbeammax*rmtarget + rmtarget**2.)
	   smax        = wmax*wmax
	   dsdtmax     = dsigmadtnot
	   dsdtmin     = dsigmadtnot*sminf !set "floor" for production cross section
	   wthreshold  = rmmeson + rmrecoil
c
	   write(6,*)' ' 
	   write(6,*)'Angular distribution parameters' 
	   write(6,*)'D(t=0) (scale)       = ',dsdtmax,' microbarn/GeV^2'
	   write(6,*)'t-independent minimum= ',dsdtmin,' microbarn/GeV^2'
	   write(6,*)'t-slope (b)          = ',bslope,' GeV^-2'
	   write(6,*)'W threshold          = ',wthreshold,' GeV'
	   write(6,*)' ' 
c	   call exit()
	end if
c
	iregge    = 0
	ww        = sqrt(2.*ppbeam*rmtarget + rmtarget**2.)
	if(ww.lt.wthreshold)then
c	   write(6,*)'Below reaction threshold',ww,wthreshold
	   iregge = -1
	   return 
	endif
	eg        = ppbeam
	egcm      = (ww*ww - (rmtarget**2. - rmg**2.    ))/(2.*ww)
	emesoncm  = (ww*ww - (rmrecoil**2. - rmmeson**2.))/(2.*ww)
	erecoilcm = (ww*ww + (rmrecoil**2. - rmmeson**2.))/(2.*ww)
	pgcm      = egcm
	pmesoncm  = emesoncm*emesoncm - rmmeson*rmmeson
        if(pmesoncm.gt.1.0e-6)then
           pmesoncm = sqrt(pmesoncm)
        else
           write(6,*)'Threshold event:',eg,ww,pmesoncm,egcm,
     1                                   emesoncm,erecoilcm,
     2                                   rmmeson,rmrecoil         
           iregge = -1
           return
        endif
	tmin = (rmmeson*rmmeson + rmg*rmg) - 
     1         2.*(egcm*emesoncm - pgcm*pmesoncm)
	tmax = (rmmeson*rmmeson + rmg*rmg) - 
     1         2.*(egcm*emesoncm + pgcm*pmesoncm)
	ss   = ww**2.
c	write(6,*)'>>>',ww,eg,egcm,emesoncm,erecoilcm,pmesoncm,tmin,tmax
c
c       Top of loop to find a t value, the associated 'cross section',
c       and then throw darts to find whether we should accept this 
c       Monte Carlo value of cosine(theta) of the production angle
c
	iflag  = 0
c
 10	iregge = iregge + 1
	if(iregge.ge.500)then
	   iflag = 1
	endif
	tt = tmin + ran(iseed)*(tmax-tmin)
	dsdt = max(dsigmadtnot * exp(2.*bslope*tt),dsdtmin)
	if(dsdt.gt.dsdtmax)then
	   write(6,*)'Exp.pick: Yikes!',dsdt,dsdtmax,dsigmadtnot,bslope,tt,ss,snot
	endif
	xtest = ran(iseed)*dsdtmax
c	write(6,*)'>',iregge,tt,dsdt,xtest
	if(xtest.gt.dsdt)goto 10
	cosinethetacm = (tt-(rmmeson*rmmeson + rmg*rmg) + 2.*(egcm*emesoncm))/
     1                                              (2.*pgcm*pmesoncm)
	if(iflag.eq.1)then
	   write(6,*)'Inefficient t-channel generation:', 
     1                                            iregge,tmin,tmax
	endif
        if(cosinethetacm==cosinethetacm)then   !test for NaN
           continue
        else
           write(6,*)'How on Earth...'
           write(6,*)ww,wthreshold,eg,egcm,emesoncm,erecoilcm,pgcm
           write(6,*)tt,rmmeson,rmg,pmesoncm
           write(6,*)tmin,tmax
        endif
        if(cosinethetacm.gt.1.0)then   !yes, this can happen
           cosinethetacm = 1.0
           write(6,*)'EXPPICK selected cosine theta =',cosinethetacm
        endif
        if(cosinethetacm.lt.-1.0)then
           cosinethetacm = -1.0
           write(6,*)'EXPPICK selected cosine theta =',cosinethetacm
        endif
c
c	write(6,*)'>',iregge,cosinethetacm,tt,dsdt,xtest,dsdtmin
c	write(6,*)' ',tmin,tmax
c
	return
	end

c************************************************************************
c       GETMASSES
c
c       This routine gets the masses of the first two-body decay mode of 
c       the particle indexed by 'index' in the particle table.  This involves
c       indexing carefully into the array that holds all the particle decay 
c       information.
c       
c	Bad programming practice: the PARAMETER variables are not in an 
c       include file, so there is a danger of screwing up the indexing
c       if the values of maxnumpar or maxdecays changes in the main program.
c
c       Inputs:
c          Index - the ordinal number specifying which particle in the input
c                  list of all particles is wanted.
c          iflag - 1 finds the particle's decay masses 
c                  0 on first call, finds the particle's decay masses
c                    while on later calls the routine merely returns what
c                    was previously found.
c       Outputs:
c          rm1   - mass of the FIRST two-body decay mode particle for the
c                  particle indexed by 'index'
c          rm2   - mass of the other particle in the first 2-body decay mode
c
c       R.A.Sch. 11-19-08
c
	subroutine getmasses(index,rm1,rm2,iflag)
	parameter (maxnumpar=65) !Maximum number of particles in the table
	parameter (maxdecays=12) !Maximum number of particles in the decay list
	integer idecay(maxnumpar,maxdecays,7)
	integer lundid(maxnumpar)
	real rm(maxnumpar)
	logical ifirst /.true./
	common /getmasseslocal/ifirst,rm1save,rm2save
	common /decaytable/idecay,lundid,rm
c
	if (ifirst .or. iflag.eq.1) then
	   if(iflag.ne.1)ifirst=.false.
	   nbody     = idecay(index,1,1)
	   ndynamics = idecay(index,1,7)
	   id1   = idecay(index,1,2)
	   id2   = idecay(index,1,3)
	   if(nbody.ne.2 .and. ndynamics.ne.0)then
	      write(6,*)'Decay to BW lineshape not defined for this kinematics'
	      write(6,*)nbody,ndynamics
	      call exit()
	   endif
	   id1 = abs(id1)
	   id2 = abs(id2)
	   rm1save = -1
	   rm2save = -1
	   do i=1,maxnumpar
	      if(lundid(i).eq.id1)then
		 rm1save = rm(i)
	      endif
	      if(lundid(i).eq.id2)then
		 rm2save = rm(i)
	      endif
c	      write(6,*)i,lundid(i),rm(i),rm1save,rm2save
	   enddo
	   if(rm1save.eq.-1)then
	      write(6,*)'Failed to decode mass values of decay product 1 ',id1
	      call exit()
	   endif
	   if(rm2save.eq.-1)then
	      write(6,*)'Failed to decode mass values of decay product 2 ',id2
	      call exit()
	   endif
c	   write(6,*)'Index, rm1, rm2',index,nbody,ndynamics,
c     1                                 id1,id2,rm1save,rm2save
	endif
	rm1 = rm1save
	rm2 = rm2save
	return
	end
      
c********************************************************************************
c       Do a two-body phase-space weighting of the decaying state
c       P0 4-vector of parent particle
c       P1 4-vector of first daughter
c       P2 4-vector of second daughter
c       PCMIN phase space momentum from initial state
c
c       Detremine the c.m. momentum of the decay.  Set a flag (IPSKEEP)
c       for whether or not to keep the event in proportion to this momentum.
c
c       R.A. Schumacher
c       CMU, October 2009
c
	subroutine phasespaceweight(p0,p1,p2,pcmin,ipskeep)
	real p0(4),p1(4),p2(4)
	data psratiomax/0/
	common /phasespacekeep/psratiomax,threeweightmax,threeweightmin,ifirst
	common /monte/ iseed,xflat,xnormal,xbw
c
	rm0 = p0(4)
	rm1 = p1(4)
	rm2 = p2(4)
	pcmvalue = (rm0**2.-(rm1**2.+rm2**2.))**2.-(2*rm1*rm2)**2.
	pcmvalue = sqrt(pcmvalue)/(2.*rm0)
	psratio  = pcmvalue/pcmin !ratio of output to input momenta
	if(psratio.gt.psratiomax)psratiomax=psratio
	x = ran(iseed)  !random between 0 and 1
	pcmtest = x * psratiomax
	if(pcmtest.le.psratio)then
	   ipskeep = 1   !keep this event
	else
	   ipskeep = 0   !abort this event
	endif
c	write(6,*)'PS>',rm0,rm1,rm2,pcmvalue,pcmin,
c     1                  psratio,psratiomax,pcmtest,ipskeep
	return
	end


c********************************************************************************
c
c       Do a THREE-body phase-space weighting of the decaying state
c       We are computing a 3-body decay as a sequence of two quasi-2-body
c       decays.  That is why we must include the phase space factors by
c       hand.  (See R.A.Sch. notes for the L(1405) lineshape paper and
c       for the L(1520) decay paper.
c
c       Inputs:
c       P0 4-vector of parent particle
c       P1 4-vector of first daughter
c       P2 4-vector of second daughter
c       PCMINTER   cm momentum from intermediate state (K+ or L(1520))
c       PCMVALUEIN cm momentum from initial state (photon on nucleon)
c       WVALUE  value of W for the overall reaction
c
c       Internal variables:
c       PCMVALUE the momentum of the two decay products (pi & Sigma*)
c       THREEWEIGHT phase space momentum for 3-body reaction state
c
c       Detremine the c.m. momentum of the decay.  Set a flag (IPSKEEP)
c       for whether or not to keep the event in proportion to this momentum.
c
c       R.A. Schumacher
c       CMU, November 2013
c
	subroutine phasespaceweight3
     1    (p0,p1,p2,pcminter,pcmvaluein,wvalue,threeweight,ipskeep)
	real p0(4),p1(4),p2(4)
        data ifirst/0/
	common /phasespacekeep3/psratiomax,threeweightmax,threeweightmin,ifirst
	common /monte/ iseed,xflat,xnormal,xbw
c
        if(ifirst.eq.0)then
           ifirst = 1
           threeweightmax = 0.5   !might need to check these value for future work
c           threeweightmin = 1.0
           threeweightmin = 0.0
           write(6,*)'First call to PHASESPACEWEIGHT3' 
           write(6,*)'Maximum phase space ',threeweightmax
           write(6,*)'Minimum phase space ',threeweightmin
        endif
c
	rm0 = p0(4)
	rm1 = p1(4)
	rm2 = p2(4)
	pcmvalue = (rm0**2.-(rm1**2.+rm2**2.))**2.-(2*rm1*rm2)**2.
	pcmvalue = sqrt(pcmvalue)/(2.*rm0)

        threeweight = (pcminter*pcmvalue)/(pcmvaluein*wvalue*wvalue)
        threeweight = threeweight * 10.  !arbitrary scale factor

	if(threeweight.gt.threeweightmax)then
           threeweightmax = threeweight
           write(6,*)'PS3> New maximum phase space ',threeweight
        endif
	if(threeweight.lt.threeweightmin)then
           threeweightmin = threeweight
           write(6,*)'PS3> New minimum phase space ',threeweight
        endif
c
c       This dart-throwing selection method is potentially very inefficient
c       but it has the virtue of being simple and of giving an accuratec
c       weighting of events by the effective 3-body phase space.
c
	x = ran(iseed)  !random number between 0 and 1
	weighttest = x * threeweightmax
	if(weighttest.le.threeweight)then
	   ipskeep = 1   !keep this event
	else
	   ipskeep = 0   !abort this event
	endif
c        write(6,*)' '
c	write(6,*)'PS3>',rm0,rm1,rm2,pcminter,pcmvaluein,wvalue
c	write(6,*)'   >',pcmvalue,threeweight,threeweightmax,weighttest,ipskeep
	return
	end

c*****************************************************************************************
c     Carve out the desired shape from the bremsstrahlung spectrum
c     We need this feature to mock up the distribution of accepted photons
c     when the Monte Carlo model does not include a production cross section
c     for the reaction(s) of interest
c
c     R. A. Schumacher, CMU September 2014
c
c     07-24-2017 - modified to create a "triangle-like" distribution for Spring 2017 data
c     11-08-2019 - reduce minimum energy to ppar treshold
c
      subroutine bremsculpt(ppbeam,ibremkeep)
      data fvaluemax/0.0/
      data ppbeamnot/0.65/  !GeV
      data wparam/0.1/  !dimensionless
      data ifirst/1/
      common /bremsculptcom/fvaluemax
      common /monte/ iseed,xflat,xnormal,xbw
cc
cc     Landau distribution
cc
c      fvalue = exp(-0.5*(ppbeam-ppbeamnot)/wparam + 
c     1           exp(-1*(ppbeam-ppbeamnot)/wparam))

c     
c     Carve out a distribution : sheer hackery to get the shape we want
c     to mimic what we see in GlueX for the Lambda - anti-Lambda
c     and proton-anti-proton
c
c      slope1 = 0.25
c      slope2 = 0.75
      slope1 = 0.15
      slope2 = 0.85
c      slope3 = -0.3333333
c      xlo  = 4.8                !GeV
c      xlo  = 4.0                !GeV
      xlo  = 3.75                !GeV
      xmid = 7.6                !GeV , inflection
      xnot = 8.8                !GeV , maximum of intensity distribution  
      xhi  = 11.8               !GeV
      if(ppbeam .le. xmid)then
         fvalue = slope1 * (ppbeam - xlo)
      elseif(ppbeam .le. xnot)then
         fvalue = slope1 * (xmid- xlo)
         fvalue = fvalue + slope2 * (ppbeam - xmid)
      elseif(ppbeam .le. xhi)then
         fvalue = slope1 * (xmid- xlo) + slope2 * (xnot-xmid)
c         fvalue = fvalue + slope3 * (ppbeam - xnot)
c         fvalue = 0.5 * fvalue
c         fvalue = 0.4 * fvalue
         fvalue = 0.27 * fvalue
      else
         fvalue = 0.0
      endif
      if(fvalue.lt.0.0)fvalue = 0.0
 100  ftest = ran(iseed) * fvaluemax
c      write(6,*)ppbeam,fvalue,ftest
      if(fvalue.gt.fvaluemax)then
         fvaluemax = fvalue
         write(6,*)'BREMSCULPT: New largest value of fvalue',
     1               fvaluemax
      endif
      if(ftest.gt.fvalue)then
         ibremkeep = 0
      else
         ibremkeep = 1
      endif
      return
      end


c*****************************************************************************************
c     Initially read in a flux file of measured beam energy distribution and select at
c     random from this distribution.
c
c     R. A. Schumacher, CMU, June 2019
c
c     This version of the code is intended to be adapted for use with 2017, 2018, 2019 ++
c     GlueX data.
c
      subroutine bremselect(ppbeam,ibremkeep)
      character*200 fluxfilename
      data fvaluemax/0.0/
      data ppbeamnot/0.65/  !GeV
      data wparam/0.1/  !dimensionless
      data ifirst/1/
      real egamma(2000),counts(2000)
      common /bremselectname/fluxfilename,lfluxfilename
      common /bremselectlocal/ifirst,nmin,nmax,fvaluemax,egamma,counts
      common /monte/ iseed,xflat,xnormal,xbw
c
c     Read in the actual experimental distribution of beam photons, a "flux file"
c     The code was written using a specific file for Spring 2017 data, it has 1000 channels, with
c     two columns: bin centers in GeV and counts at that beam energy.
c     The bin width is 10 MeV or .010 GeV.
c     Not clear is whether the given energies are the bin centers or the low edges.
c     There are lots of bins with zero counts, not only below channel 426 (6.255 GeV)
c     and above channel 944 (11.435 GeV), but in between.   That is because the
c     tagger has gaps in it outside the range of the microscope.  Thus, there is a
c     binning mismatch between the tagger energies and the bins in this histogram.
c
c     The bin-center correction used below assumes 10 MeV-wide bins
c     
      if(ifirst.eq.1)then  !Do this only the first time, when setting up the calculation
         ifirst = 0
c         write(6,*)'Reading photon energy distribution from: ',
c     1    fluxfilename(1:lfluxfilename)
         open(unit = 1, file=fluxfilename(1:lfluxfilename),status="old")
         fvaluemax = 0
         nmin = 0
         nmax = 2000  !actually, we think the file has only 1000 bins...
         ifirstn = 1
         ilastn = 1
         do n=1,2000
            read(1,*,end=10)egamma(n),counts(n)
            if(counts(n).gt.fvaluemax)fvaluemax = counts(n)
            egamma(n) = egamma(n) - 0.010/2. !use bin low edge, not the bin center
            if(counts(n).ne.0.0 .and. ifirstn.eq.1)then
               ifirstn = 0
               nmin = n         !first channel with non-zero number of photons in it
            endif
c            write(6,*)n,nmin,nmax,egamma(n),counts(n),fvaluemax
         enddo
 10      ntop = n-1
         close(1)
         do n = ntop,1,-1
            if(counts(n).ne.0 .and. ilastn.eq.1)then
               ilastn = 0
               nmax = n
            endif
         enddo
         write(6,15)'First non-zero channel',nmin,egamma(nmin),counts(nmin)
         write(6,15)'Last  non-zero channel',nmax,egamma(nmax),counts(nmax)
         write(6,20)'Highest flux number',fvaluemax
         write(6,*)'... done reading beam flux file.'
 15      format(1x,a,i5,f10.2,1e10.2)
 20      format(1x,a,1e10.2)
         deltae = abs(egamma(nmin+1) - egamma(nmin))
         if(abs(deltae-0.010).gt.0.001)then
        write(6,*)'The energy bin widths appear not to be 10 MeV',deltae
            call exit()
         endif
      endif
c      call exit()
c
c     Carve out a distribution : we are using a file that has the 'official'
c     GlueX flux distribution for a chosen run period. 
c
c     This code is not elegant, but it is done this way to interface more easily to the
c     style of the calling program.  Ideally, we don't want to return from a subroutine like this
c     without an acceptable candidate photon, but it is the calling program that sets the range of photons.
c
      if(ppbeam.lt.egamma(nmin) .or. ppbeam.gt.egamma(nmax))then
         ibremkeep = 0
         return
      endif
      do n = nmin,nmax               !channel range over which the beam has photons
         if(egamma(n).ge.ppbeam)then !n is the bin ABOVE selected energy
 100        ftest = ran(iseed) * fvaluemax
            if(ftest.gt.counts(n-1))then
               ibremkeep = 0
            else
               ibremkeep = 1
            endif
            goto 200
         endif
      enddo
      ibremkeep = 0
 200  continue
c      if(n.le.500)write(6,*)'Photon select:',
c     1 n,egamma(n),ppbeam,ftest,counts(n-1),ibremkeep
c      call exit()
      return
      end

c**************************************************************************
c     Give the initial-state particle a momentum consistent with the
c     momentum distrbution of deuterium.  Used for g14 analysis
c
c     p(1) - X component of lab momentum
c     p(2) - Y component of lab momentum
c     p(3) - Z component of lab momentum
c     p(4) - mass (input variable)
c
c     R.A.Sch. 12-2014
c
      subroutine getdeuteronmomentum(p)
      real p(4)
      logical ifirst
      data ifirst /.true./
      data phivaluemax /10000.0/
      data alpha / 0.2316 /     !inverse Fermi
      data hbarc / .197337/     !GeV fm
      data pmaxgen/0.40/        !GeV/c
c
      common /monte/ iseed,xflat,xnormal,xbw
      common/deutmomgen/phivaluemax,alphahc,betahc
c
c     Pick Fermi-motion of target particles uniformly within a sphere...
c
      if(ifirst.eqv..true.)then
         ifirst=.false.
         write(6,*)'Max momentum is = ',pmaxgen,' GeV/c'
         alphahc = alpha*hbarc
         betahc  = alphahc * 5.98
      endif
c
 100  pxfg = (2.*ran(iseed)-1.)*pmaxgen
      pyfg = (2.*ran(iseed)-1.)*pmaxgen
      pzfg = (2.*ran(iseed)-1.)*pmaxgen
      prec = sqrt(pxfg*pxfg + pyfg*pyfg + pzfg*pzfg)
c
c     Keep this momentum if the overall recoil momemtum is consistent
c     with the deuteron Hulthen momentum distribution.
c
      phivalue = (1./(prec**2. + alphahc**2.) - 
     1            1./(prec**2. + betahc**2))**2.
      ptest = ran(iseed) * phivaluemax
c     write(6,*)prec,phivalue,ptest
      if(phivalue.gt.phivaluemax)then
         phivaluemax = phivalue
         write(6,*)
     1    'GET DEUTERON MOMENTUM: New largest value of phivalue',
     1               phivaluemax
      endif
      if(ptest.gt.phivalue)goto 100 !keeps accepted values below curve
      p(1) = pxfg
      p(2) = pyfg
      p(3) = pzfg
c      write(6,*)'Get deuteron momentum distribution',p
      return
      end
c

c*******************************************************************************
c       Compute the boost beta that gets us to the rest frame of whatever
c       particle 4-vector we are given
c
c       R.A. Schumacher,  CMU,  5-2015
c
c       PVECTOR (input) is defined as PX, PY, PZ, MASS, ENERGY
c       BETA (output) is the boost 3-vector
c 
        subroutine getbeta(pvector,beta)
c
	real pvector(5),beta(3)
c
c	e0 = etot(pvector)
	e0 = pvector(5)
	do j=1,3
	   beta(j) = pvector(j)/e0
	end do
	return
	end

c************************************************************************
c       Routine to smear the particle's momentum according to a parameter.
c	Smears each mometum component the same way, thereby smearing both
c       energy and angle.
c
c       This routine is about as primitive as it can be:
c       It is just an overall momentum and angle smearing of charged particles
c       without regard to any specific of the detector geometry
c
c       GlueX version
c
c       R. A. Sch., CMU,  ~2016
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

c********************************************************************************
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

c*******************************************************************************
c       Compute Mandelstam t for two given 4-vectors
c       PA is the first ("incident/beam") particle
c       PB is the second ("produced/recoil") particle that couples at a vertex with PA
c
c       tt is the Mandelstam invariant t;  it will be < 0
c       ttprime is the reduced value of t, namely t - tmin, where
c       tmin the the 4-momentum transfer at 0 degrees
c
c       PA - 4-vector of the first  particle (px,py,pz,mass,E)
c       PB - 4-vector of the second particle (px,py,pz,mass,E)
c       P0 - intermediate 4-vector with particle B pointing in the "beam" direction
c       tt - Mandelstam invariant momentum transfer
c       ttmin - minimum value of tt when the particle are moving in the same direction
c       ttprime - reduced tt given be tt - ttmin
c      
c       R.A.Sch. 8-2017
c
	subroutine gett(pa,pb,tt,ttmin,ttprime)
	real pa(5),pb(5),p0(5)
	tt = (pb(1)-pa(1))**2. + (pb(2)-pa(2))**2. + (pb(3)-pa(3))**2.
	tt = (pb(5)-pa(5))**2. - tt
c
	pmagp = sqrt(pb(1)**2.+pb(2)**2.+pb(3)**2.)
c        costh = pa(1)*pb(1)+pa(2)*pb(2)+pa(3)*pb(3) !just need the sign
	p0(1) = 0.0
	p0(2) = 0.0
        if(pa(3).gt.0.0)then
           if(pb(3).gt.0.0)then
              p0(3) = +pmagp
           else
              p0(3) = +pmagp
           endif
        else
           if(pb(3).gt.0.0)then
              p0(3) = -pmagp
           else
              p0(3) = -pmagp
           endif
        endif
        p0(4) = pb(4)
	p0(5) = pb(5)
	ttmin = (p0(1)-pa(1))**2. + (p0(2)-pa(2))**2. + (p0(3)-pa(3))**2.
	ttmin = (p0(5)-pa(5))**2. - ttmin
	ttprime = tt - ttmin
c	write(6,*)'Computing t',tt,ttprime
c	write(6,*)p1
c	write(6,*)p3
c	write(6,*)p0
	return
	end

c********************************************************************************
c     SCULPTMASS     
c     Allow only those invariant masses that are consistent with a pre-selected
c     shape distribution
c     RM23     - input test value for an invariant mass we might accept
c     RMASSMIN - input minimum value of RM23 - the sum of rest masses
c     BSLOPE2  - input value for exponential slope
c     IYESNO   - output verdict: 1 = keep the event;  0 = reject
c     
c     R.A.Sch. 3-23-20
      
      subroutine sculptmass(rm23,rmassmin,rmassmax,bslope2,iyesno)
      common /monte/ iseed,xflat,xnormal,xbw
c
      iyesno = 0
      xx = rm23 - rmassmin
c      write(6,*) 'Sculptmass:',xx,rm23,rmassmin,rmassmax,bslope2
      yy = exp(-xx/bslope2)
      ytest = ran(iseed)
      if(ytest .lt. yy)then
         iyesno = 1
      else
         iyesno = 0
      endif
      return
      end
      
c========================== End of MC_GEN code ==============================
