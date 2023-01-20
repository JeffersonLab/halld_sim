#include <iostream>

#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"

#include "HDDM/hddm_s.hpp"
#include "particleType.h"
#include "UTILITIES/BeamProperties.h"
#include "UTILITIES/MyReadConfig.h"

#include "carbon_FF.h"

//STRUCTURE TO KEEP THE CONFIGURATION SETTINGS
struct genSettings_t {
  Particle_t target = C12;
  double ma = 0.1349768;
  double width = 7.635e-9;
  double eLower = 5.0;
  double eUpper = 10.8;
  bool verbose = false;
  int runNum = 90504;
  unsigned int nGenFactor = 1000;
  unsigned int nSample = 10000;
  int prescale = 1000000;
  int rSeed = 0;            //seed for random number generator
  string rootFile = "genOut.root";     //name of output ROOT file
  string hddmFile = "genOut.hddm";     //name of output HDDM file
  string beamconfigfile = "beam.cfg"; //beam cfg file
  string genconfigfile = "gen.cfg"; //generator cfg file
};

genSettings_t genSettings;

TH1D * hGvsE;
TH1D * hGvsEout;

// Reservoir variables
vector<double> R;
vector<double> ENERGY;
vector<TLorentzVector> VA;
vector<TLorentzVector> VREC;
double minR;
int minIndex;
double jump;

// Normalization variables
int nGenerated = 0;
double sumW = 0;

// Function templates
void sampleEvent(double weight, double eBeam, TLorentzVector va, TLorentzVector vRec);
void getWeightedEvent(double &weight, double &eBeam, TLorentzVector &va, TLorentzVector &vRec);
void t_scatter(double &weight, double eBeam, TLorentzVector &va, TLorentzVector &vRec);
void isoDecay(TLorentzVector va, TLorentzVector &v1, TLorentzVector &v2);
double dsigmadt(double eBeam, double t);
double getCarbonFF0(double q);
double getGeneralFF(double t, double Anuc, double Znuc);
double getH(double s, double t);
double getSAndRCorr(double q);

const double alpha = 0.0072973525664;
const double GeVfm = 0.1973;
const double cmSqGeVSq = GeVfm*GeVfm*1.E-26;
const double nbGeVSq = cmSqGeVSq*1.E33;

const double mU = 0.9314941024;
const double me = 0.000511;
const double m_12C = 12. * mU - 6*me;

const double mup = 2.79;
const double mp =  0.9382720881;
Bool_t b_generalFF = false;
Bool_t b_atomic_em = false;
TString str_decay = "a>gg";

void Usage()
{

  cerr << "Usage: ./gen_ALP [optional flags]\n\n"
       << "Optional flags:\n"
       << "-v: Verbose\n"
       << "-r: Set random seed\n"
       << "-z: Set simulated run number\n"
       << "-n: Set number of output events\n"
       << "-N: Set number of generated events per output event\n"
       << "-a: Set ALP mass [default pion mass]\n"
       << "-G: Set ALP->2gamma width [default pion width]\n"
       << "-B: Set beam config file\n"
       << "-C: Set generator config file\n"
       << "-R: Set output ROOT file\n"
       << "-H: Set output HDDM file\n"
       << "-h: Print this message and exit\n\n\n";

}

bool init(int argc, char ** argv)
{

  int c;
  while ((c = getopt (argc, &argv[0], "vr:z:n:N:a:G:B:C:R:H:h")) != -1)
    switch(c)
      {

      case 'v':
	genSettings.verbose = true;
	break;
      case 'r':
	genSettings.rSeed = atoi(optarg);
	break;
      case 'z':
	genSettings.runNum = atoi(optarg);
	break;
      case 'n':
	genSettings.nSample = atoi(optarg);
	break;
      case 'N':
	genSettings.nGenFactor = atoi(optarg);
        break;
      case 'a':
	genSettings.ma = atof(optarg);
	break;
      case 'G':
	genSettings.width = atof(optarg);
	break;
      case 'B':
	genSettings.beamconfigfile = optarg;
	break;
      case 'C':
	genSettings.genconfigfile = optarg;
	break;
      case 'R':
	genSettings.rootFile = optarg;
	break;
      case 'H':
	genSettings.hddmFile = optarg;
	break;
      case 'h':
	Usage();
	return false;
      default:
	abort();

      }

  // Parse config file
  if (genSettings.genconfigfile != "") 
    {

      MyReadConfig * ReadFile = new MyReadConfig();
      ReadFile->ReadConfigFile(genSettings.genconfigfile.c_str());

      if (ReadFile->GetConfigName("mass") != "")
	genSettings.ma = ReadFile->GetConfig1Par("mass")[0];
      if (ReadFile->GetConfigName("width") != "")
        genSettings.width = ReadFile->GetConfig1Par("width")[0];
      if (ReadFile->GetConfigName("target") != "") {
	genSettings.target = ParticleEnum((ReadFile->GetConfigName("target")).Data());
	b_generalFF = true;
	if (ReadFile->GetConfigName("atomic") != "")
	  b_atomic_em = true;
      }
      if (ReadFile->GetConfigName("decay") != "")
	str_decay = ReadFile->GetConfigName("decay");
    }

  // Get beam properties from configuration file 
  TH1D * cobrem_vs_E = 0;
  if (genSettings.beamconfigfile != "") 
    {

      BeamProperties beamProp( genSettings.beamconfigfile );
      cobrem_vs_E = (TH1D*)beamProp.GetFlux();
      
    }

  // GET THE HISTOGRAM FOR COHERENT BREMSTRAHLUNG SPECTRUM
  hGvsE=(TH1D*)cobrem_vs_E;
  hGvsEout = (TH1D*)hGvsE->Clone("hGvsEout");
  int eBinLow = hGvsE->GetXaxis()->FindBin(genSettings.eLower);
  int eBinHigh = hGvsE->GetXaxis()->FindBin(genSettings.eUpper);
  hGvsE->GetXaxis()->SetRange(eBinLow,eBinHigh);

  gRandom = new TRandom3(genSettings.rSeed);

  return true;

}

void evnt(int event)
{

  double weight = 0;
  double eBeam;
  TLorentzVector va;
  TLorentzVector vRec;

  while (!(weight > 0.))
    {

      getWeightedEvent(weight,eBeam,va,vRec);
      nGenerated++;

    }

  sampleEvent(weight,eBeam,va,vRec);
  sumW += weight;

}

void fini()
{

  // Initialize output files

  // ROOT
  TFile * outputFile = new TFile(genSettings.rootFile.c_str(),"RECREATE");
  TTree * outputTree = new TTree("genT","ALP MC Tree");
  TLorentzVector vBeam,v1,v2,vRec;
  outputTree->Branch("pBeam",&vBeam);
  outputTree->Branch("pGamma1",&v1);
  outputTree->Branch("pGamma2",&v2);
  outputTree->Branch("pRec",&vRec);
  
  // HDDM
  std::ofstream *outfile = new std::ofstream(genSettings.hddmFile.c_str());
  hddm_s::ostream *outstream = new hddm_s::ostream(*outfile);

  // Loop over all sampled events
  for (unsigned int event = 0; event < genSettings.nSample; event++)
    {

      double eBeam = ENERGY[event];
      TLorentzVector va = VA[event];
      vRec = VREC[event];

      vBeam.SetXYZT(0.,0.,eBeam,eBeam);
      if (str_decay == "a>gg")
	isoDecay(va,v1,v2);
      	
      // Save to ROOT tree
      outputTree->Fill();

      // Format HDDM event
      hddm_s::HDDM record;
      hddm_s::PhysicsEventList pes = record.addPhysicsEvents();
      pes().setRunNo(genSettings.runNum);
      pes().setEventNo(event);
      hddm_s::ReactionList rs = pes().addReactions();

      hddm_s::TargetList ts = rs().addTargets();
      if (!b_atomic_em)
	ts().setType(genSettings.target);
      else
	ts().setType(Electron);
      hddm_s::PropertiesList tpros = ts().addPropertiesList();
      if (!b_atomic_em) {
	tpros().setCharge(ParticleCharge(genSettings.target));
	tpros().setMass(ParticleMass(genSettings.target));
      } else {
	tpros().setCharge(ParticleCharge(Electron));
	tpros().setMass(ParticleMass(Electron));
      }
      hddm_s::MomentumList tmoms = ts().addMomenta();
      tmoms().setPx(0);
      tmoms().setPy(0);
      tmoms().setPz(0);
      if (!b_atomic_em)
	tmoms().setE(ParticleMass(genSettings.target));
      else
	tmoms().setE(ParticleMass(genSettings.target));

      hddm_s::BeamList bs = rs().addBeams();
      bs().setType(Gamma);
      hddm_s::PropertiesList bpros = bs().addPropertiesList();
      bpros().setCharge(ParticleCharge(Gamma));
      bpros().setMass(ParticleMass(Gamma));
      hddm_s::MomentumList bmoms = bs().addMomenta();
      bmoms().setPx(0);
      bmoms().setPy(0);
      bmoms().setPz(eBeam);
      bmoms().setE(eBeam);

      hddm_s::VertexList vs = rs().addVertices();
      hddm_s::OriginList os = vs().addOrigins();
      hddm_s::ProductList ps = vs().addProducts(3);

      os().setT(0);
      os().setVx(0);
      os().setVy(0);
      os().setVz(0);
      if (!b_atomic_em) {
	ps(0).setType(genSettings.target);
	ps(0).setPdgtype(PDGtype(genSettings.target));
      } else {
	ps(0).setType(Electron);
	ps(0).setPdgtype(PDGtype(Electron));
      }
      ps(0).setId(1);
      ps(0).setParentid(0);
      ps(0).setMech(0);
      hddm_s::MomentumList pRec = ps(0).addMomenta();
      pRec().setPx(vRec.X());
      pRec().setPy(vRec.Y());
      pRec().setPz(vRec.Z());
      pRec().setE(vRec.T());
      if (str_decay == "a>gg") {
	ps(1).setType(Gamma);
	ps(1).setPdgtype(PDGtype(Gamma));
	ps(1).setId(2);
	ps(1).setParentid(1);
	ps(1).setMech(1);
	hddm_s::MomentumList p1 = ps(1).addMomenta();
	p1().setPx(v1.X());
	p1().setPy(v1.Y());
	p1().setPz(v1.Z());
	p1().setE(v1.T());
	
	ps(2).setType(Gamma);
	ps(2).setPdgtype(PDGtype(Gamma));
	ps(2).setId(3);
	ps(2).setParentid(2);
	ps(2).setMech(2);
	hddm_s::MomentumList p2 = ps(2).addMomenta();
	p2().setPx(v2.X());
	p2().setPy(v2.Y());
	p2().setPz(v2.Z());
	p2().setE(v2.T());
      } else if (str_decay == "a") {
	ps(1).setType(Unknown);
	ps(1).setPdgtype(PDGtype(Unknown));
	ps(1).setId(2);
	ps(1).setParentid(1);
	ps(1).setMech(1);
	hddm_s::MomentumList p1 = ps(1).addMomenta();
	p1().setPx(va.X());
	p1().setPy(va.Y());
	p1().setPz(va.Z());
	p1().setE(va.T());
      }
      *outstream << record;

    }

  // Write out ROOT tree
  outputFile->cd();
  outputTree->Write();
  outputFile->Close();

  cout << "Total cross section: " << (sumW/nGenerated) << " nb\n";

}

int main(int argc, char ** argv)
{

  if (not init(argc,argv))
    return -1;

  for (unsigned int event = 0; event < genSettings.nGenFactor*genSettings.nSample; event++)
    {

      if ((event % genSettings.prescale == 0) && (genSettings.verbose))
	cout << "Working on reservoir event " << event << "\n";

      evnt(event);

    }

  fini();

  return 0;

}

void sampleEvent(double weight, double eBeam, TLorentzVector va, TLorentzVector vRec)
{

  // Reservoir sampling using the A-ExpJ Algorithm

  if (R.size() < genSettings.nSample)
    {

      double r = pow(gRandom->Uniform(),1./weight);
      
      R.push_back(r);
      ENERGY.push_back(eBeam);
      VA.push_back(va);
      VREC.push_back(vRec);

      minR = *min_element(R.begin(),R.end());
      minIndex = min_element(R.begin(),R.end()) - R.begin();
      jump = log(gRandom->Uniform()) / log(minR);

    }
  else
    {

      jump -= weight;
      if (jump <= 0.)
	{

	  double t = pow(minR,weight);
	  double r = pow(gRandom->Uniform(t,1.),1./weight);

	  R[minIndex] = r;
	  ENERGY[minIndex] = eBeam;
	  VA[minIndex] = va;
	  VREC[minIndex] = vRec;

	  minR = *min_element(R.begin(),R.end());
	  minIndex = min_element(R.begin(),R.end()) - R.begin();
	  jump = log(gRandom->Uniform()) / log(minR);

	}

    }

}

void getWeightedEvent(double &weight, double &eBeam, TLorentzVector &va, TLorentzVector &vRec)
{

  weight = 1; // Start weight at 1 and multiply

  // Sampling photon beam energy
  /*
  eBeam = gRandom->Uniform(genSettings.eLower,genSettings.eUpper);
  eBin = hGvsE->GetXaxis()->FindBin(eBeam);
  double binCounts = hGvsE->GetBinContent(eBin);
  double binWidth = hGvsE->GetBinWidth(eBin);
  weight *= (binCounts/binWidth) / (hGvsE->Integral()/(genSettings.eUpper - genSettings.eLower));

  if (weight <= 0.)
  return;
  */
  eBeam = hGvsE->GetRandom();

  t_scatter(weight,eBeam,va,vRec);

  if (weight <= 0.)
    return;

  TLorentzVector vBeam(0,0,eBeam,eBeam);
  double t = (vBeam - va).M2();

  weight *= dsigmadt(eBeam,t);

}

void t_scatter(double &weight, double eBeam, TLorentzVector &va, TLorentzVector &vRec)
{

  // Getting generator constants
  double ma = genSettings.ma;
  double mA = ParticleMass(genSettings.target);
  if (b_atomic_em) mA = ParticleMass(Electron);
  TLorentzVector vbeam(0,0,eBeam,eBeam);
  TLorentzVector vA(0,0,0,mA);

  TLorentzVector Z_lab = vbeam + vA;
  double s = Z_lab.M2();
  double W_cm = sqrt(s);
  if (W_cm < ma + mA)
    {
      weight=0.;
      return;
    }
    
  // Boost vectors to scattering CM frame
  TVector3 m = Z_lab.BoostVector();
  TLorentzVector vbeam_cm = vbeam;
  TLorentzVector vA_cm = vA;
  vbeam_cm.Boost(-m);
  vA_cm.Boost(-m);

  // Rotate to scattering along z-axis
  double rot_phi = vbeam_cm.Vect().Phi();
  double rot_theta = vbeam_cm.Vect().Theta();

  // Determine scattered energy and momentum
  double Ea_cm = (s + ma*ma - mA*mA)/(2.*W_cm);
  double ERec_cm = W_cm - Ea_cm;
  double p_cm = sqrt(Ea_cm*Ea_cm - ma*ma);

  // Sample t from 1/|t| distribution and randomly assign azimuthal angle
  double k_cm = vbeam_cm.Vect().Mag();
  double A = 2*k_cm*p_cm;
  double B = 2*(mA*mA - sqrt(mA*mA + k_cm*k_cm)*sqrt(mA*mA + p_cm*p_cm));
  double tmin = B - A;
  double tmax = B + A;
  double z = gRandom->Rndm();
  double t = tmin * pow(tmax/tmin,z);
  double PDF = 1/(t*log(tmax/tmin));
  double cosThetaCM = (t-B)/A;
  double theta_cm = acos(cosThetaCM);
  double phi_cm = 2.*TMath::Pi()*gRandom->Rndm();
  weight *= 1/PDF;

  // Set outgoing particle vectors
  TVector3 v_cm;
  v_cm.SetMagThetaPhi(p_cm,theta_cm,phi_cm);

  // Rotate back
  v_cm.RotateY(rot_theta);
  v_cm.RotateZ(rot_phi);
  TLorentzVector va_cm(v_cm,Ea_cm);
  TLorentzVector vRec_cm(-v_cm,ERec_cm);
    
  // Boost to lab system
  va = va_cm;
  vRec = vRec_cm;
  va.Boost(m);
  vRec.Boost(m);
    
}

void isoDecay(TLorentzVector va, TLorentzVector &v1, TLorentzVector &v2)
{

  double k = genSettings.ma/2.;
  
  double x,y,z;
  gRandom->Sphere(x,y,z,k);
  
  v1.SetXYZT(x,y,z,k);
  v2.SetXYZT(-x,-y,-z,k);

  TVector3 b = va.BoostVector();

  v1.Boost(b);
  v2.Boost(b);

}

double dsigmadt(double eBeam, double t)
{

  // Getting generator constants                                                                                                                        
  double width = genSettings.width;
  double Z = ParticleCharge(genSettings.target);
  double mA = ParticleMass(genSettings.target);
  if (b_atomic_em) mA = ParticleMass(Electron);
  double s = 2*eBeam*mA + mA*mA;
  double q = sqrt(t * (t/(4*mA*mA) - 1))/GeVfm;
  double FF0 = 1;
  double H  = getH(s,t);
  if (b_generalFF) 
    FF0 = getGeneralFF(-t, mA, Z) / pow(Z, 2);
  else
    FF0 = getCarbonFF0(q);
  if (b_atomic_em)
    FF0 = sqrt(getSAndRCorr(q));
  
  return nbGeVSq*alpha*Z*Z*FF0*FF0*width*H;
}


double getSAndRCorr(double q) 
{
  //SCREENING AND RADIATIVE CORRECTIONS
  double bohrRadius = (1.0 / alpha) * (1.0 / me);
  double fH = 1.0 / pow(1 + pow(bohrRadius * q / 2.0, 2), 2);
  double sHfactorPair = pow(1.0 - fH,2); //Screening for pair production
  //double sHfactorTrip = 1.0 - pow(fH,2); //Screening for triplet production
  
  double radFactorDelta = 0.0093;
  double radFactorPair =  1.0 + radFactorDelta;
  
  return radFactorPair * sHfactorPair;
}

double getCarbonFF0(double q)
{

  if ((q < 0) or (q > 10))
    return 0;

  double step = 0.01;
  int bin = q/step;
  double x = (q/step) - bin;

  double FF0 = carbon_FF[bin]*(1-x) + carbon_FF[bin+1]*x;
  return FF0;

}

double getGeneralFF(double t, double Anuc, double Znuc)
{

  double a = 111.0 * pow(Znuc, -1. / 3.) / me;
  double d = 0.164/*GeV^2*/ * pow(Anuc, -2. / 3.);
  
  double G2_el = pow(pow(a, 2) * t / (1 + pow(a, 2) * t), 2) * pow(1 / (1. + t / d), 2) * pow(Znuc, 2);
  
  double ap = 773 * pow(Znuc, -2. / 3.) / me;
  
  double G2_in = pow(pow(ap, 2) * t / (1. + pow(ap, 2) * t), 2) * pow((1. + t / (4 * pow(mp, 2)) * (pow(mup, 2) - 1.)) / pow(1 + t / 0.71/*GeV^2*/, 4), 2) * Znuc;
  
  double FF = G2_el + G2_in;//
  
  return FF;
}


double getH(double s, double t)
{

  // Getting generator constants

  double ma = genSettings.ma;
  double mA = ParticleMass(genSettings.target);
  if (b_atomic_em) mA = ParticleMass(Electron);

  double numerator = ma*ma*t*(mA*mA + s) - pow(ma,4)*mA*mA - t*(pow(s - mA*mA,2) + s*t);
  double denominator = t*t*pow(s - mA*mA,2)*pow(t-4*mA*mA,2);

  return 128*TMath::Pi()*pow(mA,4)/pow(ma,3)*numerator/denominator;

}
