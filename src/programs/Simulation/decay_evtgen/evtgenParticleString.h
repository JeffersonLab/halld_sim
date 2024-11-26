
#ifndef __EVTGENPARTICLESTRING_H__
#define __EVTGENPARTICLESTRING_H__

inline static char* EvtGenOutputString(Particle_t p)
{
  //returns string that is exact match to enum name. for auto-generating code
  p = RemapParticleID(p);

  switch (p) {
  case UnknownParticle:
    return (char*)"Unknown";
  case Gamma:
    return (char*)"gamma";
  case Positron:
    return (char*)"e+";
  case Electron:
    return (char*)"e-";
  case Neutrino:
    return (char*)"nu_e";
  case MuonPlus:
    return (char*)"mu+";
  case MuonMinus:
    return (char*)"mu-";
  case Pi0:
    return (char*)"pi0";
  case PiPlus:
    return (char*)"pi+";
  case PiMinus:
    return (char*)"pi-";
  case KLong:
    return (char*)"K_L0";
  case KPlus:
    return (char*)"K+";
  case KMinus:
    return (char*)"K-";
  case Neutron:
    return (char*)"n0";
  case Proton:
    return (char*)"p+";
  case AntiProton:
    return (char*)"anti-p-";
  case KShort:
    return (char*)"K_S0";
  case Eta:
    return (char*)"eta";
  case Lambda:
    return (char*)"Lambda0";
  case SigmaPlus:
    return (char*)"Sigma+";
  case Sigma0:
    return (char*)"Sigma0";
  case SigmaMinus:
    return (char*)"Sigma-";
  case Xi0:
    return (char*)"Xi0";
  case XiMinus:
    return (char*)"Xi-";
  case OmegaMinus:
    return (char*)"Omega-";
  case AntiNeutron:
    return (char*)"anti-n0";
  case AntiLambda:
    return (char*)"anti-Lambda0";
  case AntiSigmaMinus:
    return (char*)"anti-Sigma-";
  case AntiSigma0:
    return (char*)"anti-Sigma0";
  case AntiSigmaPlus:
    return (char*)"anti-Sigma+";
  case AntiXi0:
    return (char*)"anti-Xi0";
  case AntiXiPlus:
    return (char*)"anti-Xi+";
  case AntiOmegaPlus:
    return (char*)"anti-Omega+";
  case Geantino:
    return (char*)"geantino";
  case Rho0:
    return (char*)"rho0";
  case RhoPlus:
    return (char*)"rho+";
  case RhoMinus:
    return (char*)"rho-";
  case omega:
    return (char*)"omega";
  case EtaPrime:
    return (char*)"eta'";
  case phiMeson:
    return (char*)"phi";
  case a0_980:
    return (char*)"a_0";
  case f0_980:
    return (char*)"f_0";
  case KStar_892_0:
    return (char*)"K*0";
  case KStar_892_Plus:
    return (char*)"K*+";
  case KStar_892_Minus:
    return (char*)"K*-";
  case AntiKStar_892_0:
    return (char*)"anti-K*0";
  case K1_1400_Plus:
    return (char*)"K'_1+";
  case K1_1400_Minus:
    return (char*)"K'_1-";
  case b1_1235_Plus:
    return (char*)"b_1+";
  case Sigma_1385_Minus:
    return (char*)"Sigma_1385_Minus";
  case Sigma_1385_0:
    return (char*)"Sigma_1385_0";
  case Sigma_1385_Plus:
    return (char*)"Sigma_1385_Plus";
  case Deuteron:
    return (char*)"deuteron";
  case Triton:
    return (char*)"Triton";  // FIX
  case Helium:
    return (char*)"Helium";  // FIX
  case He3:
    return (char*)"He3";
  case Pb208:
    return (char*)"Pb208";   // FIX??
  case DeltaPlusPlus:
    return (char*)"Delta++";
  case Jpsi:
    return (char*)"J/psi";
  case Eta_c:
    return (char*)"eta_c";
  case Chi_c0:
    return (char*)"chi_c0";
  case Chi_c1:
    return (char*)"chi_c1";
  case Chi_c2:
    return (char*)"chi_c2";
  case Psi2s:
    return (char*)"psi(2S)";
  case D0:
    return (char*)"D0";
  case AntiD0:
    return (char*)"anti-D0";
  case DPlus:
    return (char*)"D+";
  case Dstar0:
    return (char*)"D*0";
  case DstarPlus:
    return (char*)"D*+";
  case Lambda_c:
    return (char*)"Lambda_c0";
  default:
    return (char*)"Unknown";
  }
}


#endif  // __EVTGENPARTICLESTRING_H__
