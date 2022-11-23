#!/usr/bin/python3
#
#   Author:        Curtis A. Meyer
#   Creation Date: March 2022
#
#   Purpose: This code will generate photoproduction events of the form:
#
#                           gamma p -> p X  -> p a b c
#
#            where the intermediate X will decay to three final state mesons.
#            If any of a,b and c are pi0 or eta mesons, those will be decayed
#            to two photons.
#
#            The photproduction reaction chooses a photon energy uniformly distributed
#            between 8 GeV and 9 GeV, and selects the mass of X uniformly between
#            1.02(mass(a) + mass(b) + mass(c)) and 3.0GeV. The photoproduction reaction
#            is then produced with a fixed t-slope, where the default is 2, but this can
#            be changed with the -t or --tslope options.
#
#            The decay to the three photon final state is thrown according to three body
#            phase space, such that for a fixed mass(X), the Dalitz plots will be uniformly
#            populated.
#
#            The program currently has a number of possibiliies for a,b and c which are
#            selected using the -r or --reaction options. The current ones are
#
#            -r pi0pi0pi0 generates pi0 pi0 pi0, where all three pi0s decay to 2 photons.
#            -r pi+pi-pi0 generates pi+ pi- pi0, where the pi0 decays to 2 photons.
#            -r pi+pi-eta generates pi+ pi- eta, where the eta decays to 2 photons.
#            -r pi0pi0eta generates pi0 pi0 eta, where the eta and pi0s all decay to 2 photons.
#            -r etaetapi0 generates eta eta pi0, where the eta's and pi0 all decay to 2 photons.
#            -r etaetaeta generates eta eta eta, where the eta's all decay to 2 photons.
#            -r k+k-pi0 generates K+ K- pi0, where the pi0 decays to 2 photons.
#
#            Other command-line options control the run number, number of events and output
#            file name.
#            -R or --run RunNo sets the run number, the default is 99999.
#            -e or --events NoEvt sets the numbre of events to generate, default is 1,000,000.
#            -f or --file FileName sets the output file name.
#
#   Caveats: This code is written in puython3, so your PYTHONPATH system variable must contain
#            the directory  "  ${HALLS_RECON_HOME}/${BMS_OSNAME}/python3 "
#
#            _________________________________________________________________________________
#
from math import *
import numpy as np
import vector
import hddm_s, fileinput, re
import argparse
#
# Setup the argument parser
#
parser = argparse.ArgumentParser()
parser.add_argument('-R','--run', type=int, help='Run number, default is 99999')
parser.add_argument('-r','--reaction', type=str, help='Reaction tyep: pi+pi-pi0, pi0pi0pi0, pi+pi-eta, pi0pi0eta, etaetapi0, etaetaeta, k+k-pi0',required=True)
parser.add_argument('-e','--events' , type=int, help='Number of evets to generate, default=1,000,000')
parser.add_argument('-f','--file' ,type=str, help='Output file name', required=True)
parser.add_argument('-t','--tslope',type=float,help='t-slope, default=2.0')
#
args = parser.parse_args()
#
#
#
RunNum = 99999
if args.run:
    RunMum = args.run
#
# Number of events to generate
#
NumEvt = 1000000
if args.events:
    NumEvt = args.events
#
# t slope
#
tSlope = 2.0
if args.tslope:
    tSlope = args.tslope
#
if args.file:
    filename = args.file
#
#   Define particle masses
#
Masspi0 = 0.1349770
Masspiq = 0.13957039
Masskq  = 0.493677
Masseta = 0.547862
#
SlopePiqpiqpi0 =  0.56759478
InterPiqpiqpi0 =  0.03055805
#
SlopePi0pi0pi0 =  0.56871802
InterPi0pi0pi0 =  0.02668859
#
SlopePiqpiqeta =  0.62258055
InterPiqpiqeta = -0.24705330
#
SlopePi0pi0eta =  0.62375361
InterPi0pi0eta = -0.25110280
#
SlopeEtaetapi0 =  0.47445698
InterEtaetapi0 =  0.50342619
#
SlopeEtaetaeta =  0.53244812
InterEtaetaeta =  0.21451717
#
SlopeKqkqpi0 = 0.48480714
InterKqkqpi0 = 0.43222451
#
MassP   = 0.938272081
#
MinEgam = 8.0
MaxEgam = 9.0
#
# Set up random number generator
#
rng = np.random.default_rng(12345)
#print(rng)
#print(rng.random())
#
# Set the Geant Particle IDs
#
IDPhoton  = 1
IDPiPlus  = 8
IDPiMinus = 9
IDKplus   = 11
IDKminus  = 12
IDProton  = 14
IDEta     = 17
#
# Set PDG Particle IDs
#
IDPhoPdg  =  22
IDPi0Pdg  =  111
IDPiPPdg  =  211
IDPiMPdg  = -211
IDKpPdg   =  321
IDKmPdg   = -321
IDEtaPdg  =  221
IDProPdg  = 2212
#
if args.reaction=="pi+pi-pi0":
#
# Set up the three daughter particles
#
    m1 = Masspiq
    m2 = Masspiq
    m3 = Masspi0
#
# Set up the geant ID's for the three particles. If we want to have a pi0 or eta decay to
# two photons, and not have geant also track the pi0/eta, we set it geant ID to be zero.
#
    ID1 = IDPiPlus
    ID2 = IDPiMinus
    ID3 = 0
#
# Set up the PDG IDs for the three particles. 
#
    pdgID1 = IDPiPPdg
    pdgID2 = IDPiMPdg
    pdgID3 = IDPi0Pdg
#
#  Choose the slope and intercept for the needed reaction.
#
    Slope = SlopePiqpiqpi0
    Inter = InterPiqpiqpi0
#
elif args.reaction=="pi0pi0pi0":
#
    m1 = Masspi0
    m2 = Masspi0
    m3 = Masspi0
    ID1 = 0
    ID2 = 0
    ID3 = 0
    pdgID1 = IDPi0Pdg
    pdgID2 = IDPi0Pdg
    pdgID3 = IDPi0Pdg
    Slope = SlopePi0pi0pi0
    Inter = InterPi0pi0pi0
#
elif args.reaction=="pi+pi-eta":
#
    m1 = Masspiq
    m2 = Masspiq
    m3 = Masseta
    ID1 = IDPiPlus
    ID2 = IDPiMinus
    ID3 = 0
    pdgID1 = IDPiPPdg
    pdgID2 = IDPiMPdg
    pdgID3 = IDEtaPdg
    Slope = SlopePiqpiqeta
    Inter = InterPiqpiqeta
#
elif args.reaction=="pi0pi0eta":
#
    m1 = Masspi0
    m2 = Masspi0
    m3 = Masseta
    ID1 = 0
    ID2 = 0
    ID3 = 0
    pdgID1 = IDPi0Pdg
    pdgID2 = IDPi0Pdg
    pdgID3 = IDEtaPdg
    Slope = SlopePi0pi0eta
    Inter = InterPi0pi0eta
#
elif args.reaction=="etaetapi0":
#
    m1 = Masseta
    m2 = Masseta
    m3 = Masspi0
    ID1 = 0
    ID2 = 0
    ID3 = 0
    pdgID1 = IDEtaPdg
    pdgID2 = IDEtaPdg
    pdgID3 = IDPi0Pdg
    Slope = SlopeEtaetapi0
    Inter = InterEtaetapi0
#
elif args.reaction=="etaetaeta":
#
    m1 = Masseta
    m2 = Masseta
    m3 = Masseta
    ID1 = 0
    ID2 = 0
    ID3 = 0
    pdgID1 = IDEtaPdg
    pdgID2 = IDEtaPdg
    pdgID3 = IDEtaPdg
    Slope = SlopeEtaetaeta
    Inter = InterEtaetaeta
#
elif args.reaction=="k+k-pi0":
#
    m1 = Masskq
    m2 = Masskq
    m3 = Masspi0
    ID1 = IDKplus
    ID2 = IDKminus
    ID3 = 0
    pdgID1 = IDKpPdg
    pdgID2 = IDKmPdg
    pdgID3 = IDPi0Pdg
    Slope = SlopeKqkqpi0
    Inter = InterKqkqpi0
#
else:
    print("Invalid Reaction")
    exit()
#
pi2 = 2.0 * np.pi
#
# Define functions used in this code.
#
#
# Generate a random mass between MinMass and MaxMass
#
def Uniform(Min,Max):
    a = rng.random()
    return Min + a*(Max-Min)
#
#----------------
#
# Define the Kallen Function
#
def Lambda(a,b,c):
#
    
    lam = a*a - 2.*a*(b+c) + (b-c)*(b-c) 
    if lam < 0.:
        lam = 0.
#
    return lam
#
#    For a given mass M, what is the maximum weight in the Dalitz plot?
#
def MaxWeight(M,M1,M2,M3):
#
#  Use the linear model to get the value of m12 that aximizes the weight.
#
    ma = Slope * M + Inter
#
#  Square needed masses
#
    Msq = M*M
    masq = ma*ma
    mmsq=(M-M3)*(M-M3)
# 
    M1sq = M1*M1
    M2sq = M2*M2
    M12sq = (M1+M2)*(M1+M2)
    M3sq = M3*M3
#
#   Evaluate the four needed Kallen functions
#
    lam1 = Lambda(masq,M1sq,M2sq)
    lam2 = Lambda(Msq,masq,M3sq)
    lam3 = Lambda(Msq,M12sq,M3sq)
    lam4 = Lambda(mmsq,M1sq,M2sq)
#
    wt = np.sqrt((lam1*lam2)/(lam3*lam4)) * (M-m3)/ma
#
    return wt
#
#------------------------
#
# Get the limits on t given Egam, and two masses. This assumes a reaction
#   gamma M1 -> M1 + M2
#
# t.Limits will return a numpy array that contains:
# [0] = s    [1] = tmin  [2] = tmax  [3] = p1cm   [4] = p3cm
# where p1cm is the energy of the beam photon in the center of mass
# fram and p3cm is the momentum of the meson system in the center of
# mass frame.
#
def tLimits(Eg,M1,M2):
#
#  Calculate s
#
    s = M1*(M1+2.0*Eg)
    M1sq = M1*M1
    M2sq = M2*M2
#
#  Calculate the momentum of the photon in the center of mass
#
    p1cm = M1*Eg/np.sqrt(s)
#
# Calculate the center of mass momentum of the outgoing particles.
# 
    p3cm = 0.5*np.sqrt(Lambda(s,M1sq,M2sq)/s)
#
#  Calculate limits on Mandelstam t (which is negative) tmin and tmax
#
    MsFc = M2**4/(4.0*s)
    tmin = MsFc - (p1cm-p3cm)*(p1cm-p3cm)
    tmax = MsFc - (p1cm+p3cm)*(p1cm+p3cm)
#
    return np.array([s,tmin,tmax,p1cm,p3cm]) 
#  
#---------------
#
# Throw t according to an exponential distribution with decay parameter slope
#
def tDist(tLim,slope):
    Rmin = 1.-np.exp(slope*tLim[1])
    Rmax = 1.-np.exp(slope*tLim[2])
    u = Uniform(Rmin,Rmax)
#
    return np.log(1.-u)/slope
#
#----------------
#
#   GenBod will return the three four vectors for m1, m2 and m3 in the rest frame of a
#   particle of mass M. The events are distrubuted according to three body phase space.
#
def GenBod(M,M1,M2,M3):
#
#   Get the maximum weight for the given value of M
#
    WtMx = MaxWeight(M,M1,M2,M3)
    if WtMx <= 0.:
        print("Weight of zero, adjust scale factor. This data may be garbage")
    Msq = M*M
    mmsq = (M-M3)*(M-M3)
    M1sq = M1*M1
    M2sq = M2*M2
    M3sq = M3*M3
    M12sq = (M1+M2)*(M1+M2)
    Pi2 = 2.0*np.pi
#
    lam3 = Lambda(Msq,M12sq,M3sq)
    lam4 = Lambda(mmsq,M1sq,M2sq)
    factor = (M-M3)/np.sqrt(lam3*lam4)
#
#
#
    Prob = 0.0
    Rndm = 1.0
#
    while Prob<=Rndm:
#
#   Choose a value of m12 between (m1+m2) and (M-m3)
#
        ma12 = Uniform((M1+M2),(M-M3))
        ma12sq = ma12 * ma12
#
#    Get the weight asociated with this value of M and m12
#
        lam1 = Lambda(ma12sq,M1sq,M2sq)
        lam2 = Lambda(Msq,ma12sq,M3sq)      
#
        wt = factor * np.sqrt(lam1*lam2)* (M-M3)/ma12
        Prob = wt/WtMx
        Rndm = rng.random()
#
#   We now have a value of ma12 for the event that has been chosen according to phase
#   space, we can now produce the three four vectors.
#
    p12 = np.sqrt(lam1)/(2.0*ma12)
    csth = Uniform(-1.,1.)
    snth = np.sqrt(1-csth*csth)
    phi = Uniform(0.,Pi2)
    p1 = vector.obj(px=p12*snth*np.cos(phi),py=p12*snth*np.sin(phi),pz=p12*csth,M=M1)
    p2 = vector.obj(px=-p1.x,py=-p1.y,pz=-p1.z,M=M2)
#
    pmag3 = np.sqrt(lam2)/(2.0*M)
    csth = Uniform(-1.,1.)
    snth = np.sqrt(1-csth*csth)
    phi = Uniform(0.,Pi2)
    p3 = vector.obj(px=pmag3*snth*np.cos(phi),py=pmag3*snth*np.sin(phi),pz=pmag3*csth,M=M3)
#
#     Get the vector oppoesite ps to boost p1 & p2.
#
    p3m = vector.obj(px=-p3.x,py=-p3.y,pz=-p3.z,M=ma12)
    p1b = p1.boost_p4(p3m)
    p2b = p2.boost_p4(p3m)
#
    p_all = np.zeros((3,4))
#
    p_all[0,0] = p1b.px
    p_all[0,1] = p1b.py
    p_all[0,2] = p1b.pz
    p_all[0,3] = p1b.E
#
    p_all[1,0] = p2b.px
    p_all[1,1] = p2b.py
    p_all[1,2] = p2b.pz
    p_all[1,3] = p2b.E
#
    p_all[2,0] = p3.px
    p_all[2,1] = p3.py
    p_all[2,2] = p3.pz
    p_all[2,3] = p3.E
#
    return p_all
#
#__________________________________________________________
# Start of main code execution
#
fout = hddm_s.ostream(filename)
event = None
#
#  We will generate NumEvt events with Egam between 8 and 9 GeV and MM between
#  1.03*(m1+m2+m3) and 2.8 GeV If we use a scale factor smaller than 3%, we can get divide by zeros.
#
MMmin = 1.03*(m1+m2+m3)
#
# Create an empty array to hold the four vectors from GenBod
#
pcm = np.zeros((3,4))
#
#if TestDal:
#    mm12sq = np.zeros(NumEvt)
#    mm13sq = np.zeros(NumEvt)
#    mm23sq = np.zeros(NumEvt)
#
#
iEvt = 0
while iEvt < NumEvt:
#
#  Get a photon energy bewteen 8.0 and 9.0 GeV. Then s = 2Egam mp + mp^2
#
    Egam  = Uniform(MinEgam,MaxEgam)
    s = MassP*(2.0*Egam+MassP)
#
#  Throw a three-body meson mass uniform between 1.01*(m1+m2+m3) and 3.0 GeV. We will
#  do a sanity check to make sue that (MM+MassP)**2 is smaller than s. This should not 
#  be possible as an 8 GeV photon should be able to produce a meson system with a mass
#  of 3.04GeV.
#
    MM = 30000.
    while (MM+MassP**2) >= s:
        MM = Uniform(MMmin,3.)
#
#    if TestDal:
#        MM = 1.500
#
#   Get s, the limits on t and the center of mass momentum of the incoming photon 
#   and the outgoing meson.
#
    tRange = np.zeros((5))
    tRange = tLimits(Egam,MassP,MM)
#
#   Pick a value of t based on the limits and the tSlope.
#
    t = tDist(tRange,tSlope)
#
#  Get Cos(theta)_cm from t and throw phi uniform in 0 to 2pi
#
    CosTh = 1. - 0.5*(tRange[1]-t)/(tRange[3]*tRange[4])
    SinTh = np.sqrt(1.0-CosTh*CosTh)
    Phi = Uniform(0.,pi2)
#
#  Get the vector to boost the three daughter particles of the meson system from the
#  meson rest frame to the center of mass. Then get the momentum of the final proton 
#  in the center of mass frame.
#
    p3cm = vector.obj(px=tRange[4]*SinTh*np.cos(Phi),py=tRange[4]*SinTh*np.sin(Phi),pz=tRange[4]*CosTh,M=MM)
    ppcm = vector.obj(px=-p3cm.x,py=-p3cm.y,pz=-p3cm.z,M=MassP)
#
# Use genbod to generate the threebody final state according to phase space.
#
    pcm = GenBod(MM,m1,m2,m3)
#
# Copy the results from GenBoD into 4-vectors. These vectors are reported in the
# rest frame of the three meson system, so we will then need to boost them to com.
#
    pa =vector.obj(px=pcm[0,0],py=pcm[0,1],pz=pcm[0,2],E=pcm[0,3])
    pb =vector.obj(px=pcm[1,0],py=pcm[1,1],pz=pcm[1,2],E=pcm[1,3])
    pc =vector.obj(px=pcm[2,0],py=pcm[2,1],pz=pcm[2,2],E=pcm[2,3])
#
# Boost the three daughter mesons into the center of mass frame.
#
    pacm = pa.boost_p4(p3cm)
    pbcm = pb.boost_p4(p3cm)
    pccm = pc.boost_p4(p3cm)
#
# Get the boost from the center of mass into the lab frame. For this, 
# beta = Egam/Egam+mp
#
    Beta = Egam/(Egam+MassP)
#
    pplab = ppcm.boostZ(beta=Beta)
    palab = pacm.boostZ(beta=Beta)
    pblab = pbcm.boostZ(beta=Beta)
    pclab = pccm.boostZ(beta=Beta)
#
#
    event = hddm_s.HDDM()
    pev = event.addPhysicsEvents(1)
    pev[0].runNo = RunNum
    pev[0].eventNo = iEvt
#
# Get incident photon beam
#
    rea = pev[0].addReactions(1)
    bea = rea[0].addBeams(1)
    bea[0].type = IDPhoton
#
    mom = bea[0].addMomenta(1)
    mom[0].px = 0.
    mom[0].py = 0.
    mom[0].pz = Egam
    mom[0].E  = Egam 
# 
    pol = bea[0].addPolarizations(1)
    pol[0].Px = 0.
    pol[0].Py = 0.33
    pol[0].Pz = 0.
#
    tar = rea[0].addTargets(1)
    tar[0].type = IDProton
#
    mom = tar[0].addMomenta(1)
    mom[0].px = 0.
    mom[0].py = 0.
    mom[0].pz = 0.
    mom[0].E  = MassP
#
#  Create a primary vertex
#
    vtx = rea[0].addVertices(1)
    ori = vtx[0].addOrigins(1)
    ori[0].vx = 0.0
    ori[0].vy = 0.0
    ori[0].vz = 0.0
    ori[0].t  = 0.0
#
#  Add the proton, pi+,pi- and eta to the primary vertex. Then add the eta to two photon decay
#
    pro = vtx[0].addProducts(4)
#
    pro[0].decayVertex=0
    pro[0].id = 1             # Particle Number
    pro[0].parentid = 0       # Parent Number
    pro[0].pdgtype = IDProPdg # PDG code
    pro[0].type = IDProton    # Geant Code.
    mom = pro[0].addMomenta(1)
    mom[0].px = pplab.px
    mom[0].py = pplab.py
    mom[0].pz = pplab.pz
    mom[0].E  = pplab.E
#
#  Meson 1 (pi+)
#
    pro[1].decayVertex=0
    pro[1].id = 2             # Particle Number
    pro[1].parentid = 0       # Parent Number
    pro[1].pdgtype = pdgID1   # PDG code
    pro[1].type = ID1         # Geant Code.
    mom = pro[1].addMomenta(1)
    mom[0].px = palab.px
    mom[0].py = palab.py
    mom[0].pz = palab.pz
    mom[0].E  = palab.E
#
# Meson 2 (pi-)
#
    pro[2].decayVertex=0
    pro[2].id = 3             # Particle Number
    pro[2].parentid = 0       # Parent Number
    pro[2].pdgtype = pdgID2   # PDG code
    pro[2].type = ID2         # Geant Code.
    mom = pro[2].addMomenta(1)
    mom[0].px = pblab.px
    mom[0].py = pblab.py
    mom[0].pz = pblab.pz
    mom[0].E  = pblab.E
#
# Meson 3 (eta)
#
    pro[3].decayVertex=0
    pro[3].id = 4             # Particle Number
    pro[3].parentid = 0       # Parent Number
    pro[3].pdgtype = pdgID3   # PDG code
    pro[3].type = ID3         # Geant Code.
    mom = pro[3].addMomenta(1)
    mom[0].px = pclab.px
    mom[0].py = pclab.py
    mom[0].pz = pclab.pz
    mom[0].E  = pclab.E
#
#  Check if any of the mesons decay to two photons
#
    IDmx = 4
    IVtx = 0
#
    if ID1 == 0:
#
        IDmx += 1
#        
#   Decay particle a/1 to two photons.
#   The photon momentum in the a/1 rest frame is m1/2
#
        CosTh = Uniform(-1.,1.)
        SinTh = np.sqrt(1.0-CosTh*CosTh)
        Phi   = Uniform(0.,pi2)
#
        pg1  = vector.obj(px=0.5*m1*SinTh*np.cos(Phi),py=0.5*m1*SinTh*np.sin(Phi),pz=0.5*m1*CosTh,M=0.)
        pg2  = vector.obj(px=-pg1.x,py=-pg1.y,pz=-pg1.z,M=0.)
#
#   Boost the two photons along the meson four momentum to move the four vectors to the lab.
#
        pag1lab = pg1.boost_p4(palab)
        pag2lab = pg2.boost_p4(palab)
#
#  Create a secondary vertex for the eta
#
        vtx = rea[0].addVertices(1)
        ori = vtx[0].addOrigins(1)
        ori[0].vx = 0.0
        ori[0].vy = 0.0
        ori[0].vz = 0.0
        ori[0].t  = 0.0
#
#  Now add the eta to two photon decay
#
        pro = vtx[0].addProducts(2)
#
# We now have the two photons from the deacy of the eta
#
        pro[0].decayVertex=1
        pro[0].id = IDmx          # Particle Number
        pro[0].parentid = 2       # Parent Number
        pro[0].pdgtype = IDPhoPdg # PDG code
        pro[0].type = IDPhoton    # Geant Code.
        mom = pro[0].addMomenta(1)
        mom[0].px = pag1lab.px
        mom[0].py = pag1lab.py
        mom[0].pz = pag1lab.pz
        mom[0].E  = pag1lab.E
#
        IDmx += 1
        pro[1].decayVertex=1
        pro[1].id = IDmx      # Particle Number
        pro[1].parentid = 2       # Parent Number
        pro[1].pdgtype = IDPhoPdg # PDG code
        pro[1].type = IDPhoton    # Geant Code.
        mom = pro[1].addMomenta(1)
        mom[0].px = pag2lab.px
        mom[0].py = pag2lab.py
        mom[0].pz = pag2lab.pz
        mom[0].E  = pag2lab.E
#
    if ID2 == 0:
#
        IDmx += 1
#        
#   Decay particle a/1 to two photons.
#   The photon momentum in the a/1 rest frame is m1/2
#
        CosTh = Uniform(-1.,1.)
        SinTh = np.sqrt(1.0-CosTh*CosTh)
        Phi   = Uniform(0.,pi2)
#
        pg1  = vector.obj(px=0.5*m2*SinTh*np.cos(Phi),py=0.5*m2*SinTh*np.sin(Phi),pz=0.5*m2*CosTh,M=0.)
        pg2  = vector.obj(px=-pg1.x,py=-pg1.y,pz=-pg1.z,M=0.)
#
#   Boost the two photons along the meson four momentum to move the four vectors to the lab.
#
        pbg1lab = pg1.boost_p4(pblab)
        pbg2lab = pg2.boost_p4(pblab)
#
#  Create a secondary vertex for the eta
#
        vtx = rea[0].addVertices(1)
        ori = vtx[0].addOrigins(1)
        ori[0].vx = 0.0
        ori[0].vy = 0.0
        ori[0].vz = 0.0
        ori[0].t  = 0.0
#
#  Now add the eta to two photon decay
#
        pro = vtx[0].addProducts(2)
#
#
# We now have the two photons from the deacy of the eta
#
        pro[0].decayVertex=1
        pro[0].id = IDmx          # Particle Number
        pro[0].parentid = 3       # Parent Number
        pro[0].pdgtype = IDPhoPdg # PDG code
        pro[0].type = IDPhoton    # Geant Code.
        mom = pro[0].addMomenta(1)
        mom[0].px = pbg1lab.px
        mom[0].py = pbg1lab.py
        mom[0].pz = pbg1lab.pz
        mom[0].E  = pbg1lab.E
#
        IDmx += 1
        pro[1].decayVertex=1
        pro[1].id = IDmx      # Particle Number
        pro[1].parentid = 3       # Parent Number
        pro[1].pdgtype = IDPhoPdg # PDG code
        pro[1].type = IDPhoton    # Geant Code.
        mom = pro[1].addMomenta(1)
        mom[0].px = pbg2lab.px
        mom[0].py = pbg2lab.py
        mom[0].pz = pbg2lab.pz
        mom[0].E  = pbg2lab.E
#
#
    if ID3 == 0:
#
        IDmx += 1
#        
#   Decay particle a/1 to two photons.
#   The photon momentum in the a/1 rest frame is m1/2
#
        CosTh = Uniform(-1.,1.)
        SinTh = np.sqrt(1.0-CosTh*CosTh)
        Phi   = Uniform(0.,pi2)
#
        pg1  = vector.obj(px=0.5*m3*SinTh*np.cos(Phi),py=0.5*m3*SinTh*np.sin(Phi),pz=0.5*m3*CosTh,M=0.)
        pg2  = vector.obj(px=-pg1.x,py=-pg1.y,pz=-pg1.z,M=0.)
#
#   Boost the two photons along the meson four momentum to move the four vectors to the lab.
#
        pcg1lab = pg1.boost_p4(pclab)
        pcg2lab = pg2.boost_p4(pclab)
#
#  Create a secondary vertex for the eta
#
        vtx = rea[0].addVertices(1)
        ori = vtx[0].addOrigins(1)
        ori[0].vx = 0.0
        ori[0].vy = 0.0
        ori[0].vz = 0.0
        ori[0].t  = 0.0
#
#  Now add the eta to two photon decay
#
        pro = vtx[0].addProducts(2)
#
#
# We now have the two photons from the deacy of the eta
#
        pro[0].decayVertex=1
        pro[0].id = IDmx          # Particle Number
        pro[0].parentid = 4       # Parent Number
        pro[0].pdgtype = IDPhoPdg # PDG code
        pro[0].type = IDPhoton    # Geant Code.
        mom = pro[0].addMomenta(1)
        mom[0].px = pcg1lab.px
        mom[0].py = pcg1lab.py
        mom[0].pz = pcg1lab.pz
        mom[0].E  = pcg1lab.E
#
        IDmx += 1
        pro[1].decayVertex=1
        pro[1].id = IDmx      # Particle Number
        pro[1].parentid = 4       # Parent Number
        pro[1].pdgtype = IDPhoPdg # PDG code
        pro[1].type = IDPhoton    # Geant Code.
        mom = pro[1].addMomenta(1)
        mom[0].px = pcg2lab.px
        mom[0].py = pcg2lab.py
        mom[0].pz = pcg2lab.pz
        mom[0].E  = pcg2lab.E
#
    fout.write(event)

#
#    if TestDal:
#        p12 = palab + pblab
#        p13 = palab + pclab
#        p23 = pblab + pclab
#
#        mm12sq[iEvt] = p12.M2
#        mm13sq[iEvt] = p13.M2
#        mm23sq[iEvt] = p23.M2
#
    iEvt+=1
#

