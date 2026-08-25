#define FormFactors_cxx
#include "FormFactors.h"
#include <limits>
#include <unistd.h>
using namespace std;

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

float G_E(int nucleon, float QQ){
        float result=0;
        float MagM = 2.7928;
        if (nucleon==1) result = G_M(nucleon,QQ)/MagM * (1. - 0.1306 * QQ + 0.004174 * pow(QQ,2.) - 0.000752 * pow(QQ,3.));
        return result;

}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

float G_M (int nucleon, float QQ){
        float result=0;
        float MagM = 2.7928;
        if (nucleon==1) result = MagM /( 1. + 0.116 * sqrt(QQ) + 2.874 * QQ
                        + 0.241 * pow(QQ,1.5) + 1.006 * pow(QQ,2.) + 0.345 * pow(QQ,2.5));
        return result;
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

float FormFactors(int FF,int nucleon, float QQ){

        float result=0;
        float G_E = 0, G_M=0;
        float MagM = 2.7928;

        if (nucleon==1){
                G_M = MagM /( 1. + 0.116 * sqrt(QQ) + 2.874 * QQ
                        + 0.241 * pow(QQ,1.5) + 1.006 * pow(QQ,2.) + 0.345 * pow(QQ,2.5));

                G_E = G_M/MagM * (1. - 0.1306 * QQ + 0.004174 * pow(QQ,2.) - 0.000752 * pow(QQ,3.));

        } else if (nucleon==2){
                cout<<"not implemented yet"<<endl;
                G_E=0;
                G_M=0;

        } else {
                G_E=0;
                G_M=0;
        }

        if (FF==1){
                result = ( QQ / pow(2. * M_Nucleon, 2.) * G_M + G_E ) /( 1. + QQ / pow(2. * M_Nucleon, 2.) );
        } else if (FF ==2){
                result =  ( G_M - G_E ) / ( 1. + QQ / pow(2. * M_Nucleon, 2.) );
        }
        return result;
}

