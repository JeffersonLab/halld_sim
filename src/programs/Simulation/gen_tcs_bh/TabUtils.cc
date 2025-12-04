#include "TabUtils.h"

bool Vars_Cin_tcs::operator<(const Vars_Cin_tcs& vc_tcs) const {
    int i = 0;
    bool ret = false;

    for(; i < NB_VAR_CIN_tcs; i++) {
        if(vars[i] > vc_tcs.vars[i]) {
            break;
        } else if(vars[i] < vc_tcs.vars[i]) {
            ret = true;
            break;
        }
    }

    return ret;
}

std::ostream& operator<<(std::ostream& out, const Vars_Cin_tcs& vc_tcs) {
    int i = 0;

    out << "{";

    for(; i < NB_VAR_CIN_tcs; i++) {
        out << vc_tcs.vars[i] << (i == NB_VAR_CIN_tcs-1 ? "}" : ",");
    }

    return out;
}

std::ostream& operator<<(std::ostream& out, const Vars_Sef_tcs& vs_tcs) {
    int i = 0;

    out << "{";

    for(; i < NB_VAR_SEF_tcs; i++) {
        out << vs_tcs.vars[i] << (i == NB_VAR_SEF_tcs-1 ? "}" : ",");
    }

    return out;
}

bool Vars_Cin_tcs_pol::operator<(const Vars_Cin_tcs_pol& vc_tcs_pol) const {
    int i = 0;
    bool ret = false;

    for(; i < NB_VAR_CIN_tcs_pol; i++) {
        if(vars[i] > vc_tcs_pol.vars[i]) {
            break;
        } else if(vars[i] < vc_tcs_pol.vars[i]) {
            ret = true;
            break;
        }
    }

    return ret;
}

std::ostream& operator<<(std::ostream& out, const Vars_Cin_tcs_pol& vc_tcs_pol) {
    int i = 0;

    out << "{";

    for(; i < NB_VAR_CIN_tcs_pol; i++) {
        out << vc_tcs_pol.vars[i] << (i == NB_VAR_CIN_tcs_pol-1 ? "}" : ",");
    }

    return out;
}

std::ostream& operator<<(std::ostream& out, const Vars_Sef_tcs_pol& vs_tcs_pol) {
    int i = 0;

    out << "{";

    for(; i < NB_VAR_SEF_tcs_pol; i++) {
        out << vs_tcs_pol.vars[i] << (i == NB_VAR_SEF_tcs_pol-1 ? "}" : ",");
    }

    return out;
}