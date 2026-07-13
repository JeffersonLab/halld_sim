
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <bitset>
#include <cctype>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TFile.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/KStarHyperonExtended.h"


static std::bitset<54> parseMask54Hex(std::string h) {
  if (h.size() >= 2 && h[0] == '0' && (h[1] == 'x' || h[1] == 'X')) h = h.substr(2);
  assert(!h.empty());

  unsigned long long x = 0ULL;
  for (char c : h) {
    assert(std::isxdigit(static_cast<unsigned char>(c)));
    x <<= 4;
    if (c >= '0' && c <= '9') x |= (c - '0');
    else if (c >= 'a' && c <= 'f') x |= (c - 'a' + 10);
    else x |= (c - 'A' + 10);
  }

  // 54 bits must fit in 64-bit
  assert((x >> 54) == 0ULL && "mask hex exceeds 54 bits");

  std::bitset<54> m;
  for (int t = 0; t < 54; ++t) if (x & (1ULL << t)) m.set(t);
  return m;
}


static vector<int> parseIndexList(const string& s) {
  vector<int> indices;
  std::stringstream ss(s);
  int val;
  char comma;
  while (ss >> val) {
      indices.push_back(val);
      if (ss.peek() == ',') ss >> comma;
  }
  return indices;
}


KStarHyperonExtended::KStarHyperonExtended( const vector<string>& args ) :
  UserAmplitude<KStarHyperonExtended>( args ),
  polFrac_vs_E( nullptr )
{
  // Expect:
  // 0: kIndices,
  // 1: lambdaIndices,
  // 2: model (0/1/2), 
  //        0     - just x/y/z cross terms; 
  //        1     - cross terms + Schilling; 
  //        2     - factorized model (18 parameters)
  //        3/4/5 - transverse(*, *) + longitudinal(x/y/z) terms + Schilling
  // 3: alpha,
  // 4: N_coeff,
  // 5..5+N-1: coeffs (extended terms only)
  // [if SchillingIncluded==1] 5+N .. 5+N+8: rho000,rho100,rho1m10,rho111,rho001,rho101,rho1m11,rho102,rho1m12
  // 5+N(+9): mask (54 bits hex like 0x3FFFFFFFFFFFFF)
  // 6+N(+9): polAngle
  // 7+N(+9): polMode ("const" or "file")
  // 8+N(+9): polFraction (if const) OR polFile (if file)
  // 9+N(+9): polHist (if file)

  assert(args.size() >= 9); // minimal sanity; stronger checks below

  // --- indices ---
  kIndices      = parseIndexList( args[0] );
  lambdaIndices = parseIndexList( args[1] );
  assert( kIndices.size() == 2 && "K* must have exactly 2 daughter indices" );
  assert( lambdaIndices.size() == 2 && "Lambda must have exactly 2 daughter indices" );
  cout << "K indices: ";
  for(auto i : kIndices) cout << i << " ";
  cout << "Lambda indices: ";
  for(auto i : lambdaIndices) cout << i << " ";
  cout << endl;

  // --- Schilling flag ---
  model = atoi(args[2].c_str());
  bool phaseIncluded = (model >= 3);     
  schillingIncluded = (model != 0);
  int nRho = schillingIncluded ? (model == 2 ? 17 : 9) : 0;

  // --- alpha ---
  alpha = AmpParameter( args[3] );
  registerParameter(alpha);

  // --- N (extended coeff count) ---
  const int N = atoi(args[4].c_str());
  assert(N >= 1);
  if (phaseIncluded) {
    assert(N == 39 && "model>=3 requires exactly 39 extended coefficients (18T + 3phi + 18L)");
  }

  // offsets that depend on SchillingIncluded
  int idxCoeff0  = 5;
  int idxRho0    = idxCoeff0 + N;          // only if schillingIncluded
  int idxMask    = idxCoeff0 + N + nRho;
  int idxPolAng  = idxMask + 1;
  int idxPolMode = idxMask + 2;

  // need at least: ... + mask + polAngle + polMode + polArg1
  assert(static_cast<int>(args.size()) >= idxPolMode + 2);

  // --- extended coeffs (fit parameters) ---
  coeffs.clear();
  coeffs.reserve(N);
  for (int t = 0; t < N; ++t) {
    coeffs.emplace_back( AmpParameter(args[idxCoeff0 + t]) );
    registerParameter( coeffs.back() );
  }

  // --- rhos (fit parameters), if enabled ---
  rhos.clear();
  if (schillingIncluded) {
    const int nRhos = (model == 2 ? 17 : 9);
    rhos.reserve(nRhos);
    for (int q = 0; q < nRhos; ++q) {
      rhos.emplace_back( AmpParameter(args[idxRho0 + q]) );
      registerParameter( rhos.back() );
    }
  }
  // --- mask ---
  if (phaseIncluded) {
    // mask is a dummy for model >= 3; keep arg position but ignore content
    termMask.reset();
    termMask.set(); // all 54 "active" conceptually
    iList.clear(); jList.clear(); kList.clear();  // not used in model >= 3
    std::cout << "[KStarHyperonExtended] "
              << "Schilling=" << (schillingIncluded ? "ON" : "OFF")
              << "  model >= 3 (phase-mixed transverse/longitudinal)  N_ext=" << N
              << std::endl << std::endl;
  } else {
    termMask = parseMask54Hex(args[idxMask]);
    if (model == 2) {
      static const std::bitset<54> kModel2Mask = parseMask54Hex("0x1041041041041");
    
      // ensure no illegal bits are set
      assert((termMask & ~kModel2Mask).none() &&
             "model==2 mask must be a subset of 0x1041041041041");
    
      // optional sanity check
      assert(termMask.count() <= 9 && "model==2 allows at most 9 terms");
    }
    // expand mask -> iList/jList/kList (active terms only)
    iList.clear(); jList.clear(); kList.clear();
    for (int t = 0; t < 54; ++t) {
      if (!termMask.test(t)) continue;
      int ii = t / 18;
      int rem = t % 18;
      int jj = rem / 6;
      int kk = rem % 6;
      iList.push_back(ii);
      jList.push_back(jj);
      kList.push_back(kk);
    }
    assert(static_cast<int>(iList.size()) == N && "N_coeff must equal number of 1-bits in mask");
    
    std::cout << "[KStarHyperonExtended] "
    << "Schilling=" << (schillingIncluded ? "ON" : "OFF")
    << "  N_ext=" << N << "\n";
    for (int t = 0; t < N; ++t)
    std::cout << "  [" << t << "] "
        << coeffs[t].name()
        << "  (" << iList[t] << "," << jList[t] << "," << kList[t] << ")\n";
    std::cout << std::endl;
    // cout << "\n[ARGS] total=" << args.size() << "\n";
    // for (size_t i = 0; i < args.size(); ++i) {
    //   cout << "  [" << i << "] " << args[i];
    //   if (i == idxCoeff0)   cout << "   <-- idxCoeff0";
    //   if (schillingIncluded && i == idxRho0) cout << "   <-- idxRho0";
    //   if (i == idxMask)     cout << "   <-- idxMask";
    //   if (i == idxPolAng)   cout << "   <-- idxPolAng";
    //   if (i == idxPolMode)  cout << "   <-- idxPolMode";
    //   cout << "\n";
    // }
    // cout << endl;
  }



  // --- polAngle ---
  polAngle = AmpParameter(args[idxPolAng]);
  registerParameter(polAngle);

  // --- polarization mode ---
  const std::string polMode = args[idxPolMode];
  if (polMode == "const") {
    assert(static_cast<int>(args.size()) == idxPolMode + 2);
    polFraction = atof(args[idxPolMode + 1].c_str());
    cout << "Fitting with constant polarization " << polFraction << endl;
    polFrac_vs_E = nullptr;
  }
  else if (polMode == "file") {
    assert(static_cast<int>(args.size()) == idxPolMode + 3);
    polFraction = 0.0;
    TFile* f = new TFile(args[idxPolMode + 1].c_str());
    polFrac_vs_E = dynamic_cast<TH1D*>( f->Get(args[idxPolMode + 2].c_str()) );
    assert(polFrac_vs_E != nullptr);
    cout << "Fitting with polarization from " << polFrac_vs_E->GetName() << endl;
  }
  else {
    assert(false && "polMode must be 'const' or 'file'");
  }


  // ---- DIAG: print all args with semantic labels ----
  cout << "\n[KStarHyperonExtended][ARGS] total=" << args.size() << "\n";
  for (size_t i = 0; i < args.size(); ++i) {
    cout << "  [" << i << "] " << args[i];

    // fixed positions
    if (i == 0) cout << "   <-- kIndices";
    if (i == 1) cout << "   <-- lambdaIndices";
    if (i == 2) cout << "   <-- model";
    if (i == 3) cout << "   <-- alpha";
    if (i == 4) cout << "   <-- N_coeff";

    // computed offsets
    if (i == (size_t)idxCoeff0) cout << "   <-- coeffs[0] (first extended coeff)";
    if (i >= (size_t)idxCoeff0 && i < (size_t)(idxCoeff0 + N))
      cout << "   <-- coeffs[" << (i - idxCoeff0) << "]";

    if (schillingIncluded) {
      if (i == (size_t)idxRho0) cout << "   <-- rhos[0] (rho000)";
      if (i >= (size_t)idxRho0 && i < (size_t)(idxRho0 + 9)) {
        static const char* rhoName[9] = {
          "rho000","rho100","rho1m10","rho111","rho001","rho101","rho1m11","rho102","rho1m12"
        };
        cout << "   <-- rhos[" << (i - idxRho0) << "] (" << rhoName[i - idxRho0] << ")";
      }
    }

    if (i == (size_t)idxMask)    cout << "   <-- mask (54-bit hex; dummy if model >= 3)";
    if (i == (size_t)idxPolAng)  cout << "   <-- polAngle (deg)";
    if (i == (size_t)idxPolMode) cout << "   <-- polMode (const/file)";

    // pol mode dependent tail
    if (i == (size_t)(idxPolMode + 1)) {
      if (i < args.size()) {
        if (args[idxPolMode] == "const") cout << "   <-- polFraction";
        else if (args[idxPolMode] == "file") cout << "   <-- polFile";
        else cout << "   <-- polArg1 (unknown polMode)";
      }
    }
    if (i == (size_t)(idxPolMode + 2)) {
      if (i < args.size()) {
        if (args[idxPolMode] == "file") cout << "   <-- polHist";
        else cout << "   <-- polArg2 (unused unless polMode==file)";
      }
    }

    cout << "\n";
  }
  cout << endl;
  // ----- end of DIAG ---- //
}



complex< GDouble >
KStarHyperonExtended::calcAmplitude( GDouble** pKin, GDouble* userVars ) const {

  GDouble Pgamma  = userVars[kPgamma];
  GDouble bigPhi  = polAngle * 0.017453293 + userVars[kBigPhi];

  GDouble g_func[3];
  g_func[0] = 1.0;
  g_func[1] = -Pgamma * cos(2.0 * bigPhi);
  g_func[2] = -Pgamma * sin(2.0 * bigPhi);


  GDouble H[3] = {};
  GDouble FH[3][2] = {
    {1.0, 0.0}, 
    {1.0, 0.0},
    {1.0, 0.0}
  }; // for factorized model

  if (model >= 3) {
    // coeff layout: [0..17]=T(ii,kk), [18..20]=phi(ii), [21..38]=L(ii,kk)
    for (int ii = 0; ii < 3; ++ii) {
      const GDouble phi = coeffs[18 + ii] * 0.017453293; // interpret phase angles in DEGREES to match polAngle convention
      const GDouble cph = cos(phi);
      const GDouble sph = sin(phi);
      
      // pick indices for transverse/longitudinal components among {2,3,4}
      const int jL = model - 1;  // 3/4/5 -> 2/3/4
      assert(jL >= 2 && jL <= 4);
      int jT0 = (jL == 4) ? 2 : jL + 1;
      int jT1 = (jL == 2) ? 4 : jL - 1;
      
      for (int kk = 0; kk < 6; ++kk) {
        const int idxT = ii * 6 + kk;           // 0..17
        const int idxL = 21 + ii * 6 + kk;      // 21..38

        const GDouble T = coeffs[idxT];
        const GDouble L = coeffs[idxL];

        const int kVar = 5 + kk;

        // nx term (jj=0): T*cos(phi)
        H[ii] += (T * cph) * g_func[ii] * userVars[jT0] * userVars[kVar];
        // ny term (jj=1): T*sin(phi)
        H[ii] += (T * sph) * g_func[ii] * userVars[jT1] * userVars[kVar];
        // nz term (jj=2): L
        H[ii] += (L)       * g_func[ii] * userVars[jL] * userVars[kVar];
      }
    }
  } else if (model == 2) {
    int nTerms = static_cast<int>(coeffs.size());
    for (int t = 0; t < nTerms; ++t) {
      int ii = iList[t];
      int jj = 2 + jList[t];
      int kk = 5 + kList[t];
      int signature = (jList[t] == 1 ? 0 : 1);  // y -> 0, x/z -> 1
      FH[ii][signature] += alpha * coeffs[t] * userVars[jj] * userVars[kk];
    }
  } else {
    int nTerms = static_cast<int>(coeffs.size());
    for (int t = 0; t < nTerms; ++t) {
      int ii = iList[t];
      int jj = 2 + jList[t];
      int kk = 5 + kList[t];
      H[ii] += coeffs[t] * g_func[ii] * userVars[jj] * userVars[kk];
    }
  }
 
  // ----- optional Schilling–Wolf SDME -----
  GDouble W[3][2] = {
    {0.0, 0.0},
    {0.0, 0.0},
    {0.0, 0.0}
  };

  if (schillingIncluded) {
    GDouble rho000  = rhos[0];
    GDouble rho100  = rhos[1];
    GDouble rho1m10 = rhos[2];
    GDouble rho111  = rhos[3];
    GDouble rho001  = rhos[4];
    GDouble rho101  = rhos[5];
    GDouble rho1m11 = rhos[6];
    GDouble rho102  = rhos[7];
    GDouble rho1m12 = rhos[8];

    const GDouble rt2 = sqrt(2.0);

    W[0][0] += 0.5*(1. - rho000) + 0.5*(3.*rho000 - 1.)*userVars[6] - rt2*rho100*userVars[7] - rho1m10*userVars[9];
    W[1][0] += rho111*(1. - userVars[6]) + rho001*userVars[6] - rt2*rho101*userVars[7] - rho1m11*userVars[9];
    W[2][1] += rt2*rho102*userVars[8] + rho1m12*userVars[10];

    if (model == 2) {
      GDouble rho100_neg  = rhos[9];
      GDouble rho1m10_neg = rhos[10];

      GDouble rho101_neg  = rhos[11];
      GDouble rho1m11_neg = rhos[12];

      GDouble rho112_neg  = rhos[13];
      GDouble rho002_neg  = rhos[14];
      GDouble rho102_neg  = rhos[15];
      GDouble rho1m12_neg = rhos[16];

      W[0][1] += rt2 * rho100_neg * userVars[8] + rho1m10_neg * userVars[10];

      W[1][1] += rt2 * rho101_neg * userVars[8] + rho1m11_neg * userVars[10];

      W[2][1] += rho112_neg*(1. - userVars[6]) + rho002_neg*userVars[6] - rt2*rho102_neg*userVars[7] - rho1m12_neg*userVars[9];
    }
  }
  else{
    W[0][0] = 1.0;
  }

  GDouble W_normalization = 3.0 / (4.0 * PI);
  GDouble H_normalization = 1.0 / (4.0 * PI);

  GDouble I_total = 0.0;
  if (model != 2) {
    I_total += W_normalization * (W[0][0] + g_func[1] * W[1][0] + g_func[2] * W[2][0])
             + alpha * H_normalization * (H[0] + H[1] + H[2]);
  } else if (model == 2) {
    I_total += W_normalization * H_normalization * (
                              W[0][0] * FH[0][0] + W[0][1] * FH[0][1] +
                 g_func[1] * (W[1][0] * FH[1][0] + W[1][1] * FH[1][1]) +
                 g_func[2] * (W[2][0] * FH[2][1] + W[2][1] * FH[2][0])
               );
  } 

  return complex<GDouble>( sqrt( fabs(I_total) ), 0 );
}



void
KStarHyperonExtended::calcUserVars( GDouble** pKin, GDouble* userVars ) const {
  
  TLorentzVector beam( 
    pKin[0][1], 
    pKin[0][2], 
    pKin[0][3], 
    pKin[0][0] );

  TLorentzVector p1( 
    pKin[kIndices[0]][1],
    pKin[kIndices[0]][2],
    pKin[kIndices[0]][3],
    pKin[kIndices[0]][0] );
    
  TLorentzVector p2( 
    pKin[kIndices[1]][1],
    pKin[kIndices[1]][2],
    pKin[kIndices[1]][3],
    pKin[kIndices[1]][0] );

  TLorentzVector y1( 
    pKin[lambdaIndices[0]][1],
    pKin[lambdaIndices[0]][2],
    pKin[lambdaIndices[0]][3],
    pKin[lambdaIndices[0]][0] );

  TLorentzVector y2( 
    pKin[lambdaIndices[1]][1],
    pKin[lambdaIndices[1]][2],
    pKin[lambdaIndices[1]][3],
    pKin[lambdaIndices[1]][0] );

  TLorentzVector resonance = p1 + p2;
  TLorentzRotation resonanceBoost( -resonance.BoostVector() );
  TLorentzVector hyperon = y1 + y2;
  TLorentzRotation hyperonBoost( -hyperon.BoostVector() );
  
  TLorentzVector beam_resonance    = resonanceBoost * beam;
  TLorentzVector hyperon_resonance = resonanceBoost * hyperon;
  TLorentzVector p1_resonance      = resonanceBoost * p1;
	
  //TLorentzVector beam_hyperon = hyperonBoost * beam;  // beam photon in hyperon rest frame
  TLorentzVector kstar_hyperon = hyperonBoost * resonance; // kstar in hyperon rest frame
  TLorentzVector y1_hyperon = hyperonBoost * y1;      // proton in hyperon rest frame

  // normal to the production plane
  TVector3 y = (beam.Vect().Unit().Cross(-hyperon.Vect().Unit())).Unit();

  // choose helicity frame: z-axis opposite recoil hyperon in resonance rest frame
  TVector3 z = -1. * hyperon_resonance.Vect().Unit();
  TVector3 x = y.Cross(z).Unit();
  TVector3 angles( (p1_resonance.Vect()).Dot(x),
		   (p1_resonance.Vect()).Dot(y),
		   (p1_resonance.Vect()).Dot(z) );

  GDouble theta = angles.Theta();
  GDouble phi   = angles.Phi();
  
  userVars[kB0] = 1.0;                                         // b_0 = 1.0
  userVars[kB1] = cos(theta) * cos(theta);                     // b_00
  userVars[kB2] = sin(2.0 * theta) * cos(phi);                 // b_10c
  userVars[kB3] = sin(2.0 * theta) * sin(phi);                 // b_10s
  userVars[kB4] = sin(theta) * sin(theta) * cos(2.0 * phi);    // b_1−1c
  userVars[kB5] = sin(theta) * sin(theta) * sin(2.0 * phi);    // b_1−1s
	
  // normal to the production plane (formed by beam and kaon)
  TVector3 yH = (beam.Vect().Unit().Cross(-resonance.Vect().Unit())).Unit();

  // choose helicity frame: z-axis opposite kaon in hyperon rest frame
  TVector3 zH = -1. * kstar_hyperon.Vect().Unit();
  TVector3 xH = yH.Cross(zH).Unit();

  userVars[kN0] = cos(y1_hyperon.Vect().Angle(xH));
  userVars[kN1] = cos(y1_hyperon.Vect().Angle(yH));
  userVars[kN2] = cos(y1_hyperon.Vect().Angle(zH));

  // pre-compute and cache information for the linear polarization 
  GDouble pGamma;
  if(polFraction > 0.) { // for fitting with constant polarization 
    pGamma = polFraction;
  }
  else{
    int bin = polFrac_vs_E->GetXaxis()->FindBin(pKin[0][0]);
    if (bin == 0 || bin > polFrac_vs_E->GetXaxis()->GetNbins()){
      pGamma = 0.;
    }
    else pGamma = polFrac_vs_E->GetBinContent(bin);
  }
  userVars[kPgamma] = pGamma;

  TVector3 eps(1.0, 0.0, 0.0); // reference beam polarization vector at 0 degrees
  userVars[kBigPhi] = atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y)));
  

}


#ifdef GPU_ACCELERATION
void
KStarHyperonExtended::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

  const int nTerms = static_cast<int>( coeffs.size() );

  // ---- pack extended coeffs into raw host array ----
  GDouble* host_coeffs = new GDouble[nTerms];
  for( int t = 0; t < nTerms; ++t ){
    host_coeffs[t] = (GDouble)coeffs[t];
  }

  // ---- pack (i,j,k) into raw host struct array ----
  // The GPU kernel ignores `terms` when model==2.
  term_ijk* host_terms = nullptr;                
  int nTerms_kernel = nTerms;                   
  if( model < 3 ){                             
    host_terms = new term_ijk[nTerms];
    for( int t = 0; t < nTerms; ++t ){
      host_terms[t].i = iList[t];
      host_terms[t].j = jList[t];
      host_terms[t].k = kList[t];
    }
  } else {
    host_terms = nullptr;                        
    nTerms_kernel = 0;                           
  }

  // ---- pack Schilling rhos if enabled ----
  GDouble* host_rhos = nullptr;
  if( schillingIncluded ){
    host_rhos = new GDouble[9];
    for( int q = 0; q < 9; ++q ){
      host_rhos[q] = (GDouble)rhos[q];
    }
  }

  if( host_coeffs == nullptr ){
    throw std::runtime_error("Memory allocation failed for host_coeffs in KStarHyperonExtended::launchGPUKernel");
  }
  if( model < 3 && host_terms == nullptr ){     
    throw std::runtime_error("Memory allocation failed for host_terms in KStarHyperonExtended::launchGPUKernel");
  }
  if( schillingIncluded && host_rhos == nullptr ){
    throw std::runtime_error("Memory allocation failed for host_rhos in KStarHyperonExtended::launchGPUKernel");
  }

  // ---- launch the kernel ----
  GPUKStarHyperonExtended_exec(
    dimGrid, dimBlock, GPU_AMP_ARGS,
    (int)model,
    (GDouble)alpha,
    host_coeffs,
    host_terms,
    nTerms_kernel, 
    schillingIncluded,
    host_rhos,          // nullptr if schillingIncluded==false
    (GDouble)polAngle
  );

  // ---- free host memory ----
  delete[] host_coeffs;
  if( host_terms ) delete[] host_terms;         
  if( host_rhos ) delete[] host_rhos;
}
#endif // GPU_ACCELERATION




