#ifndef TABUTILS_H
#define TABUTILS_H

#include <iostream>


const int NB_VAR_CIN_tcs = 5;
const int NB_VAR_SEF_tcs = 6;
const int NB_VAR_CIN_tcs_pol = 6;
const int NB_VAR_SEF_tcs_pol = 8;


struct Vars_Cin_tcs {
    bool operator<(const Vars_Cin_tcs& vc) const;
    int vars[NB_VAR_CIN_tcs];
};

struct Vars_Cin_tcs_pol {
    bool operator<(const Vars_Cin_tcs_pol& vc) const;
    int vars[NB_VAR_CIN_tcs_pol];
};


struct Vars_Sef_tcs {
    double vars[NB_VAR_SEF_tcs];
    enum {BH, TCS, SUMPP,SUMPM, SUMMP, SUMMM};
};

struct Vars_Sef_tcs_pol {
    double vars[NB_VAR_SEF_tcs_pol];
    enum {BH, TCS, SUMPP,SUMPM, SUMMP, SUMMM, BHpara, BHperp};
};


std::ostream& operator<<(std::ostream& out, const Vars_Cin_tcs& vc_tcs);
std::ostream& operator<<(std::ostream& out, const Vars_Sef_tcs& vs_tcs);
std::ostream& operator<<(std::ostream& out, const Vars_Cin_tcs_pol& vc_tcs_pol);
std::ostream& operator<<(std::ostream& out, const Vars_Sef_tcs_pol& vs_tcs_pol);

#endif // TABUTILS_H

