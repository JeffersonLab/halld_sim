#include <cassert>
#include <complex>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <math.h>
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "dblRegge.h"

dblRegge::dblRegge( const vector< string >& args ) :
        UserAmplitude< dblRegge >( args )
{
        assert( args.size() == 3 );
        j = atoi( args[0].c_str() );
        fast = atoi( args[1].c_str() );
        b = atof( args[2].c_str() );
}

complex< GDouble >
dblRegge::calcAmplitude( GDouble** pKin ) const {
        TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] );
        TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] );
        TLorentzVector p1     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] );
        TLorentzVector p2     ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] );

        TLorentzVector resonance = p1 + p2;

        double aPrime = 0.9;
///////////
        double tau1 = -1; //for Neutral case. Will be -1^J otherwise
        double tau2 = -1;
///////////
        complex< GDouble > Amp1;
        complex< GDouble > Amp2;
        complex<GDouble> TotalAmp;
        complex<GDouble> coeff1 = 0;
        complex<GDouble> coeff2 = 0;
        double s12 = resonance.M2();
        double s13 = (p1 + recoil).M2();
        double s23 =  (p2 + recoil).M2();
        double t1 = (beam - p1).M2();
        double t2 = (beam - p2).M2();
        double s = (recoil + p1 + p2).M2();
        double u3 = t1 + t2 +s12 - (beam.M2() + p1.M2() + p2.M2());

        complex<GDouble> ui (0,1);


//fast eta:
        double a1 = aPrime*t1 + 0.5;
        double a2 = aPrime*u3 + 0.5;

        complex <GDouble> Xi1 =0.5* (tau1 + exp(-ui*M_PI*a1)) ;
        complex <GDouble> Xi21 = 0.5*(tau2*tau1 + exp(-ui*M_PI*(a2 - a1)));
        double V1, V2;
        complex <GDouble> Xi2 =0.5* (tau2 + exp(-ui*M_PI*a2));

        double lambda = 0.5;
        double lambdaP = 0.5;
        if(a1==a2 || fast == 2){
                Amp1 = 0;
        }
        else{
                V1 = exp(b*t1) / (a1-a2);
                V2 = exp(b*u3) / (a2 - a1);

                Amp1 = -(TMath::Power(aPrime*s,a1)*TMath::Power(aPrime*s23, a2-a1)*Xi1*Xi21*V1 +  TMath::Power(aPrime*s,a2)*TMath::Power(aPrime*s12, a1-a2)*Xi2*Xi21*V2);


                coeff1 = (TMath::Power(-t1, 0.5) / p2.M()) * TMath::Power( (-u3/(4*recoil.M2())), 0.5*abs(lambda - lambdaP) ) ;
        }

//fast pion:

        a1 = aPrime*t2 + 0.5;
        a2 = aPrime*u3 + 0.5;

        Xi1 =0.5* (tau1 + exp(-ui*M_PI*a1)) ;
        Xi21 = 0.5*(tau2*tau1 + exp(-ui*M_PI*(a2 - a1)));
        Xi2 =0.5* (tau2 + exp(-ui*M_PI*a2));

        if(a1==a2 || fast == 1){
                Amp2 = 0;
        }
        else{
                V1 = exp(b*t2)/ (a1-a2);
                V2 = exp(b*u3) / (a2 - a1);
                Amp2 = -(TMath::Power(aPrime*s,a1)*TMath::Power(aPrime*s13, a2-a1)*Xi1*Xi21*V1 +  TMath::Power(aPrime*s,a2)*TMath::Power(aPrime*s12, a1-a2)*Xi2*Xi21*V2);


        coeff2 =  (TMath::Power(-t2, 0.5) / p1.M()) * TMath::Power( (-u3/(4*recoil.M2())), 0.5*abs(lambda - lambdaP) );
        }

        TotalAmp = coeff1*Amp1 + coeff2*Amp2;
        return TotalAmp;
}
