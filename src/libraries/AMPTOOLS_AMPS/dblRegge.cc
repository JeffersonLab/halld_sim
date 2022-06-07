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

        GDouble aPrime = 0.9;
///////////
        GDouble tau1 = -1; //for Neutral case. Will be -1^J otherwise
        GDouble tau2 = -1;
///////////
        complex< GDouble > Amp1;
        complex< GDouble > Amp2;
        complex<GDouble> TotalAmp;
        complex<GDouble> coeff1 = 0;
        complex<GDouble> coeff2 = 0;
        GDouble s12 = resonance.M2();
        GDouble s13 = (p1 + recoil).M2();
        GDouble s23 =  (p2 + recoil).M2();
        GDouble t1 = (beam - p1).M2();
        GDouble t2 = (beam - p2).M2();
        GDouble s = (recoil + p1 + p2).M2();
        GDouble u3 = t1 + t2 +s12 - (beam.M2() + p1.M2() + p2.M2());

        complex<GDouble> ui (0,1);


//fast eta:
        GDouble a1 = aPrime*t1 + 0.5;
        GDouble a2 = aPrime*u3 + 0.5;

        complex <GDouble> Xi1 = GDouble(0.5)*(tau1 + exp(-ui*GDouble(M_PI)*a1)) ;
        complex <GDouble> Xi21 = GDouble(0.5)*(tau2*tau1 + exp(-ui*GDouble(M_PI)*(a2 - a1)));
        GDouble V1, V2;
        complex <GDouble> Xi2 = GDouble(0.5)*(tau2 + exp(-ui*GDouble(M_PI)*a2));
        complex <GDouble> Xi12 = GDouble(0.5)*(tau2*tau1 + exp(-ui*GDouble(M_PI)*(a1 - a2)));

        GDouble lambda = 0.5;
        GDouble lambdaP = 0.5;
        if(a1==a2 || fast == 2){
                Amp1 = 0;
        }
        else{
                V1 = exp(b*t1) / (a1-a2);
                V2 = exp(b*u3) / (a2 - a1);

                Amp1 = -(G_POW(aPrime*s,a1)*G_POW(aPrime*s23, a2-a1)*Xi1*Xi21*V1 +  G_POW(aPrime*s,a2)*G_POW(aPrime*s12, a1-a2)*Xi2*Xi12*V2);


                coeff1 = (G_POW(-t1, 0.5) / p2.M()) * G_POW( (-u3/(4*recoil.M2())), 0.5*abs(lambda - lambdaP) ) ;
        }

//fast pion:

        a1 = aPrime*t2 + 0.5;
        a2 = aPrime*u3 + 0.5;

        Xi1 = GDouble(0.5)*(tau1 + exp(-ui*GDouble(M_PI)*a1)) ;
        Xi21 = GDouble(0.5)*(tau2*tau1 + exp(-ui*GDouble(M_PI)*(a2 - a1)));
        Xi2 = GDouble(0.5)*(tau2 + exp(-ui*GDouble(M_PI)*a2));
        Xi12 = GDouble(0.5)*(tau2*tau1 + exp(-ui*GDouble(M_PI)*(a1 - a2)));
        
        if(a1==a2 || fast == 1){
                Amp2 = 0;
        }
        else{
                V1 = exp(b*t2)/ (a1-a2);
                V2 = exp(b*u3) / (a2 - a1);
                Amp2 = -(G_POW(aPrime*s,a1)*G_POW(aPrime*s13, a2-a1)*Xi1*Xi21*V1 +  G_POW(aPrime*s,a2)*G_POW(aPrime*s12, a1-a2)*Xi2*Xi12*V2);


        coeff2 =  (G_POW(-t2, 0.5) / p1.M()) * G_POW( (-u3/(4*recoil.M2())), 0.5*abs(lambda - lambdaP) );
        }

        TotalAmp = coeff1*Amp1 + coeff2*Amp2;
         
        GDouble m1 = p1.M();
        GDouble m2 = p2.M();
        GDouble mP = recoil.M();

        GDouble breakupP = G_POW( (s*s + G_POW(m1, 4) + G_POW(m2,4)- 2*(m1*m1*(s + m2*m2) + m2*m2*s )), 0.5) / (2*G_POW(s12, 0.5));

        GDouble numericCoeff = G_POW(0.125 *G_POW((1/(4*GDouble(M_PI))), 4)* (breakupP /(s*s +mP*mP*mP*mP - 2*s*mP*mP)  ), 0.5);


        TotalAmp = TotalAmp*numericCoeff;
        
        return TotalAmp;
}
