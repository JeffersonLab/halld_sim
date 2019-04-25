//
// Initial revision 2019/4/24 11:59:01 jonesrt
//

#include <stdio.h>
#include <particleType.h>

void gspart_(int *ipart, char *chnpar, int *itrtyp,
             float *amass, float *charge, float *tlife, float *ub, int *nwb);
void gsdk_(int *ipart, float *bratio, int *mode);

const double hbar = 0.197 / 2.998e+23; // GeV.s
const int decay1=1, decay2=100, decay3=10000;

void gupart_()
{
   // Defines unstable particle types not in the standard library 

   int ipart;
   char chnpar[20];
   int itrtyp;
   float amass;
   float charge;
   float tlife;
   float ub[9]={0};
   int zero=0;
   float bratio[12] = {0};
   int mode[12] = {0};

   ipart = Rho0;
   snprintf(chnpar, 20, "%20s", ParticleType(Rho0));
   itrtyp = 4;
   amass = 0.770;
   charge = 0;
   tlife = hbar / 0.149;
   gspart_(&ipart,chnpar,&itrtyp,&amass,&charge,&tlife,ub,&zero);
   bratio[5] = 100;
   mode[5] = PiPlus*decay1 + PiMinus*decay2;
   gsdk_(&ipart, bratio+5, mode+5);

   ipart = RhoPlus;
   snprintf(chnpar, 20, "%20s", ParticleType(RhoPlus));
   itrtyp = 4;
   amass = 0.770;
   charge = 1;
   tlife = hbar / 0.149;
   gspart_(&ipart,chnpar,&itrtyp,&amass,&charge,&tlife,ub,&zero);
   bratio[5] = 100;
   mode[5] = Pi0*decay1 + PiPlus*decay2;
   gsdk_(&ipart, bratio+5, mode+5);

   ipart = RhoMinus;
   snprintf(chnpar, 20, "%20s", ParticleType(RhoMinus));
   itrtyp = 4;
   amass = 0.770;
   charge = -1;
   tlife = hbar / 0.149;
   gspart_(&ipart,chnpar,&itrtyp,&amass,&charge,&tlife,ub,&zero);
   bratio[5] = 100;
   mode[5] = Pi0*decay1 + PiMinus*decay2;
   gsdk_(&ipart, bratio+5, mode+5);

   ipart = omega;
   snprintf(chnpar, 20, "%20s", ParticleType(omega));
   itrtyp = 4;
   amass = 0.783;
   charge = 0;
   tlife = hbar / 0.00849;
   gspart_(&ipart,chnpar,&itrtyp,&amass,&charge,&tlife,ub,&zero);
   bratio[5] = 89.2;
   mode[5] = Pi0*decay1 + PiMinus*decay2 + PiPlus*decay3;
   bratio[4] = 8.3;
   mode[4] = Pi0*decay1 + Gamma*decay2;
   bratio[3] = 1.5;
   mode[3] = PiPlus*decay1 + PiMinus*decay2;
   bratio[2] = 0.87;
   mode[2] = Pi0*decay1 + Pi0*decay2 + Gamma*decay3; // other all-neutral decays
   bratio[1] = 0.05;
   mode[1] = Eta*decay1 + Gamma*decay2;
   bratio[0] = 0.08;
   mode[0] = Pi0*decay1 + Electron*decay2 + Positron*decay3;
   gsdk_(&ipart, bratio, mode);

   ipart = EtaPrime;
   snprintf(chnpar, 20, "%20s", ParticleType(EtaPrime));
   itrtyp = 4;
   amass = 0.958;
   charge = 0;
   tlife = hbar / 0.197e-3;
   gspart_(&ipart,chnpar,&itrtyp,&amass,&charge,&tlife,ub,&zero);
   bratio[5] = 42.9;
   mode[5] = PiPlus*decay1 + PiMinus*decay2 + Eta*decay3;
   bratio[4] = 29.1;
   mode[4] = Rho0*decay1 + Gamma*decay2;
   bratio[3] = 22.3;
   mode[3] = Pi0*decay1 + Pi0*decay2 + Eta*decay3;
   bratio[2] = 2.6;
   mode[2] = omega*decay1 + Gamma*decay2;
   bratio[1] = 2.1;
   mode[1] = Gamma*decay1 + Gamma*decay2;
   bratio[0] = 1.0;
   mode[0] = Pi0*decay1 + Rho0*decay2;
   gsdk_(&ipart, bratio, mode);

   ipart = phiMeson;
   snprintf(chnpar, 20, "%20s", ParticleType(phiMeson));
   itrtyp = 4;
   amass = 1.019;
   charge = 0;
   tlife = hbar / 4.27e-3;
   gspart_(&ipart,chnpar,&itrtyp,&amass,&charge,&tlife,ub,&zero);
   bratio[5] = 48.9;
   mode[5] = KPlus*decay1 + KMinus*decay2;
   bratio[4] = 34.2;
   mode[4] = KLong*decay1 + KShort*decay2;
   bratio[3] = 15.46;
   mode[3] = Pi0*decay1 + PiPlus*decay2 + PiMinus*decay3;
   bratio[2] = 1.31;
   mode[2] = Eta*decay1 + Gamma*decay2;
   bratio[1] = 0.13;
   mode[1] = Pi0*decay1 + Gamma*decay2;
   gsdk_(&ipart, bratio+1, mode+1);

   ipart = a0_980; // the neutral a0(980)
   snprintf(chnpar, 20, "%20s", ParticleType(a0_980));
   itrtyp = 4;
   amass = 0.980;
   charge = 0;
   tlife = hbar / 0.050;
   gspart_(&ipart,chnpar,&itrtyp,&amass,&charge,&tlife,ub,&zero);
   bratio[5] = 90;
   mode[5] = Eta*decay1 + Pi0*decay2;
   bratio[4] = 5;
   mode[4] = KLong*decay1 + KLong*decay2;
   bratio[3] = 5;
   mode[3] = KShort*decay1 + KShort*decay2;
   gsdk_(&ipart, bratio+3, mode+3);

   ipart = f0_980;
   snprintf(chnpar, 20, "%20s", ParticleType(f0_980));
   itrtyp = 4;
   amass = 0.990;
   charge = 0;
   tlife = hbar / 0.050;
   gspart_(&ipart,chnpar,&itrtyp,&amass,&charge,&tlife,ub,&zero);
   bratio[5] = 90 * 1./ 3;
   mode[5] = Pi0*decay1 + Pi0*decay2;
   bratio[5] = 90 * 2./3;
   mode[5] = PiPlus*decay1 + PiMinus*decay2;
   bratio[4] = 5;
   mode[4] = KLong*decay1 + KLong*decay2;
   bratio[3] = 5;
   mode[3] = KShort*decay1 + KShort*decay2;
   gsdk_(&ipart, bratio+3, mode+3);

   ipart = KStar_892_0;
   snprintf(chnpar, 20, "%20s", ParticleType(KStar_892_0));
   itrtyp = 4;
   amass = 0.896;
   charge = 0;
   tlife = hbar / 0.047;
   gspart_(&ipart,chnpar,&itrtyp,&amass,&charge,&tlife,ub,&zero);
   bratio[5] = 67;
   mode[5] = KPlus*decay1 + PiMinus*decay2;
   bratio[4] = 16.5;
   mode[4] = KShort*decay1 + Pi0*decay2;
   bratio[3] = 16.5;
   mode[3] = KLong*decay1 + Pi0*decay2;
   gsdk_(&ipart, bratio+3, mode+3);

   ipart = AntiKStar_892_0;
   snprintf(chnpar, 20, "%20s", ParticleType(AntiKStar_892_0));
   itrtyp = 4;
   amass = 0.896;
   charge = 0;
   tlife = hbar / 0.047;
   gspart_(&ipart,chnpar,&itrtyp,&amass,&charge,&tlife,ub,&zero);
   bratio[5] = 67;
   mode[5] = KMinus*decay1 + PiPlus*decay2;
   bratio[4] = 16.5;
   mode[4] = KShort*decay1 + Pi0*decay2;
   bratio[3] = 16.5;
   mode[3] = KLong*decay1 + Pi0*decay2;
   gsdk_(&ipart, bratio+3, mode+3);

   ipart = KStar_892_Plus;
   snprintf(chnpar, 20, "%20s", ParticleType(KStar_892_Plus));
   itrtyp = 4;
   amass = 0.892;
   charge = +1;
   tlife = hbar / 0.046;
   gspart_(&ipart,chnpar,&itrtyp,&amass,&charge,&tlife,ub,&zero);
   bratio[5] = 33;
   mode[5] = KPlus*decay1 + Pi0*decay2;
   bratio[4] = 33.5;
   mode[4] = KShort*decay1 + PiPlus*decay2;
   bratio[3] = 33.5;
   mode[3] = KLong*decay1 + PiPlus*decay2;
   gsdk_(&ipart, bratio+3, mode+3);

   ipart = K1_1400_Plus;
   snprintf(chnpar, 20, "%20s", ParticleType(K1_1400_Plus));
   itrtyp = 4;
   amass = 1.403;
   charge = +1;
   tlife = hbar / 0.174;
   gspart_(&ipart,chnpar,&itrtyp,&amass,&charge,&tlife,ub,&zero);
   bratio[5] = 67;
   mode[5] = KStar_892_0*decay1 + PiPlus*decay2;
   bratio[4] = 33;
   mode[4] = KStar_892_Plus*decay1 + Pi0*decay2;
   gsdk_(&ipart, bratio+4, mode+4);

   ipart = K1_1400_Minus;
   snprintf(chnpar, 20, "%20s", ParticleType(K1_1400_Minus));
   itrtyp = 4;
   amass = 1.403;
   charge = -1;
   tlife = hbar / 0.174;
   gspart_(&ipart,chnpar,&itrtyp,&amass,&charge,&tlife,ub,&zero);
   bratio[5] = 67;
   mode[5] = KStar_892_0*decay1 + PiMinus*decay2;
   bratio[4] = 33;
   mode[4] = KStar_892_Minus*decay1 + Pi0*decay2;
   gsdk_(&ipart, bratio+4, mode+4);

   ipart = b1_1235_Plus;
   snprintf(chnpar, 20, "%20s", ParticleType(b1_1235_Plus));
   itrtyp = 4;
   amass = 1.230;
   charge = +1;
   tlife = hbar / 0.142;
   gspart_(&ipart,chnpar,&itrtyp,&amass,&charge,&tlife,ub,&zero);
   bratio[5] = 100;
   mode[5] = omega*decay1 + PiPlus*decay2;
   gsdk_(&ipart, bratio+5, mode+5);

   ipart = Sigma_1385_Minus;
   snprintf(chnpar, 20, "%20s", ParticleType(Sigma_1385_Minus));
   itrtyp = 4;
   amass = 1.383;
   charge = -1;
   tlife = hbar / 0.036;
   gspart_(&ipart,chnpar,&itrtyp,&amass,&charge,&tlife,ub,&zero);
   bratio[5] = 87;
   mode[5] = Lambda*decay1 + PiMinus*decay2;
   bratio[4] = 6.5;
   mode[4] = SigmaMinus*decay1 + Pi0*decay2;
   bratio[3] = 6.5;
   mode[3] = Sigma0*decay1 + PiMinus*decay2;
   gsdk_(&ipart, bratio+3, mode+3);

   ipart = Sigma_1385_Plus;
   snprintf(chnpar, 20, "%20s", ParticleType(Sigma_1385_Plus));
   itrtyp = 4;
   amass = 1.383;
   charge = +1;
   tlife = hbar / 0.036;
   gspart_(&ipart,chnpar,&itrtyp,&amass,&charge,&tlife,ub,&zero);
   bratio[5] = 87;
   mode[5] = Lambda*decay1 + PiPlus*decay2;
   bratio[4] = 6.5;
   mode[4] = SigmaPlus*decay1 + Pi0*decay2;
   bratio[3] = 6.5;
   mode[3] = Sigma0*decay1 + PiPlus*decay2;
   gsdk_(&ipart, bratio+3, mode+3);

   ipart = Sigma_1385_0;
   snprintf(chnpar, 20, "%20s", ParticleType(Sigma_1385_0));
   itrtyp = 4;
   amass = 1.384;
   charge = 0;
   tlife = hbar / 0.036;
   gspart_(&ipart,chnpar,&itrtyp,&amass,&charge,&tlife,ub,&zero);
   bratio[5] = 87;
   mode[5] = Lambda*decay1 + Pi0*decay2;
   bratio[4] = 6;
   mode[4] = SigmaPlus*decay1 + PiMinus*decay2;
   bratio[3] = 6;
   mode[3] = SigmaMinus*decay1 + PiPlus*decay2;
   bratio[2] = 1;
   mode[2] = Lambda*decay1 + Gamma*decay2;
   gsdk_(&ipart, bratio+2, mode+2);

   ipart = DeltaPlusPlus;
   snprintf(chnpar, 20, "%20s", ParticleType(DeltaPlusPlus));
   itrtyp = 4;
   amass = 1.232;
   charge = 0;
   tlife = hbar / 0.115;
   gspart_(&ipart,chnpar,&itrtyp,&amass,&charge,&tlife,ub,&zero);
   bratio[5] = 100;
   mode[5] = Proton*decay1 + PiPlus*decay2;
   gsdk_(&ipart, bratio+5, mode+5);
}
