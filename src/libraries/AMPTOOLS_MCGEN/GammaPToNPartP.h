#if !defined(GAMMAPTONPARTP)
#define GAMMAPTONPARTP

/*
 *  GammaPToNPartP.h
 *   by Igor Senderovich
 *  structure based on GammaToXYZP
 *  written by Matthew Shepherd
 */

#include "TLorentzVector.h"
#include "TH1.h"
#include "AMPTOOLS_MCGEN/ProductionMechanism.h"

class Kinematics;

class GammaPToNPartP {
  
public:
  GammaPToNPartP();
	
  GammaPToNPartP( float lowMass, float highMass, 
		  vector<double> &ChildMass,
		  ProductionMechanism::Recoil recoil, ProductionMechanism::Type type, float slope = 6.0, double lowT = 0.0, double highT = 12.0, int seed = 0, TString beamConfigFile = "");
  
  //overload constructor for backwards compatibility
  GammaPToNPartP( float lowMass, float highMass, 
		  vector<double> &ChildMass,
		  ProductionMechanism::Type type, float slope = 6.0, double lowT = 0.0, double highT = 12.0, int seed = 0, TString beamConfigFile = "");
  
  Kinematics* generate();
  
  void addResonance( float mass, float width, float bf );
  
private:
  
  ProductionMechanism m_prodMech;
  
  TLorentzVector m_beam;
  TLorentzVector m_target;
  
  //double m_ChildMass[12];
  vector<double> m_ChildMass;
  unsigned int m_Npart;

  TH1D *cobrem_vs_E;
};

#endif
