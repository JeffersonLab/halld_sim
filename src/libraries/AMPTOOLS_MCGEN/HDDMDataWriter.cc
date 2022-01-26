
#include "TLorentzVector.h"
#include "TRandom3.h"

#include "AMPTOOLS_MCGEN/HDDMDataWriter.h"
#include "HDDM/hddm_s.hpp"

HDDMDataWriter::HDDMDataWriter(const string& outFile, int runNumber, int seed)
{
  m_OutputFile = new ofstream(outFile.c_str());
  m_OutputStream = new hddm_s::ostream(*m_OutputFile);
  m_runNumber = runNumber;
  
  m_eventCounter = 1;

  // initialize root's pseudo-random generator
  gRandom->SetSeed(seed);

}

HDDMDataWriter::~HDDMDataWriter()
{
  delete m_OutputStream;
  delete m_OutputFile;
}


void HDDMDataWriter::
writeEvent( const Kinematics& kin, const vector<int>& ptype, bool centeredVertex)
{
  if (centeredVertex)
    // this will trigger hdgeant(4) to generate the vertex distribution
    writeEvent(kin,ptype,0,0,0/*cm*/);
  else
    writeEvent(kin,ptype,0,0,50/*cm*/,80/*cm*/);
}

void HDDMDataWriter::
writeEvent( const Kinematics& kin, const vector<int>& ptype,
	    float vx, float vy, float vz_min, float vz_max)
{
  if (vz_min > vz_max) {
    float tmp=vz_min;
    vz_min=vz_max;
    vz_max=tmp;
  }
  writeEvent(kin, ptype,vx, vy, (vz_max - vz_min) * gRandom->Uniform() + vz_min);
}


void HDDMDataWriter::
writeEvent( const Kinematics& kin, const vector<int>& ptype,
	    float vx, float vy, float vz)
{
  vector< TLorentzVector > particleList = kin.particleList();
  int nParticles=kin.particleList().size();
  
  // Start a new event in the HDDM record
  hddm_s::HDDM record;
  hddm_s::PhysicsEventList pes = record.addPhysicsEvents();
  pes().setRunNo(m_runNumber);
  pes().setEventNo(m_eventCounter);
  hddm_s::ReactionList rs = pes().addReactions();
  hddm_s::VertexList vs = rs().addVertices();
  hddm_s::OriginList os = vs().addOrigins();
  hddm_s::ProductList ps = vs().addProducts(nParticles-1);
  hddm_s::RandomList ranl = rs().addRandoms();

  ranl().setSeed1(gRandom->Integer(std::numeric_limits<int32_t>::max()));
  ranl().setSeed2(gRandom->Integer(std::numeric_limits<int32_t>::max()));
  ranl().setSeed3(gRandom->Integer(std::numeric_limits<int32_t>::max()));
  ranl().setSeed4(gRandom->Integer(std::numeric_limits<int32_t>::max()));

  os().setT(0.0);
  os().setVx(vx);
  os().setVy(vy);
  os().setVz(vz);

  hddm_s::BeamList bs = rs().addBeams();
  bs().setType((Particle_t)1);
  hddm_s::MomentumList bmoms = bs().addMomenta();
  bmoms().setPx(kin.particle(0).Px());
  bmoms().setPy(kin.particle(0).Py());
  bmoms().setPz(kin.particle(0).Pz());
  bmoms().setE(kin.particle(0).E());
  hddm_s::PropertiesList bpros = bs().addPropertiesList();
  bpros().setCharge(0);
  bpros().setMass(0.0);

  hddm_s::TargetList ts = rs().addTargets();
  ts().setType((Particle_t)14);
  hddm_s::MomentumList tmoms = ts().addMomenta();
  tmoms().setPx(0);
  tmoms().setPy(0);
  tmoms().setPz(0);
  tmoms().setE(0.938272);
  hddm_s::PropertiesList tpros = ts().addPropertiesList();
  tpros().setCharge(+1);
  tpros().setMass(0.938272);
  
  for(int i=1; i < nParticles; i++)
  {
      ps(i-1).setType((Particle_t)ptype[i]);
      ps(i-1).setPdgtype(PDGtype((Particle_t)ptype[i]));
      ps(i-1).setId(i);         /* unique value for this particle within the event */
      ps(i-1).setParentid(0);   /* All internally generated particles have no parent */
      ps(i-1).setMech(0);       /* maybe this should be set to something? */
      hddm_s::MomentumList pmoms = ps(i-1).addMomenta();
      pmoms().setPx(kin.particle(i).Px());
      pmoms().setPy(kin.particle(i).Py());
      pmoms().setPz(kin.particle(i).Pz());
      pmoms().setE(kin.particle(i).E());
  }
  
  if (nParticles > 0)
    *m_OutputStream << record;
  m_eventCounter++;
}
