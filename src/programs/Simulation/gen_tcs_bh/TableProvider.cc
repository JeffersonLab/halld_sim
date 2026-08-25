#include "TableProvider.h"
#include <iostream>

bool TableProvider::insert_tcs(const Vars_Cin_tcs& key, const Vars_Sef_tcs& value) {
    return (__sef_table_tcs.insert(Table_Element_tcs({key, value}))).second;
}

bool TableProvider::insert_tcs_pol(const Vars_Cin_tcs_pol& key, const Vars_Sef_tcs_pol& value) {
    return (__sef_table_tcs_pol.insert(Table_Element_tcs_pol({key, value}))).second;
}


Vars_Sef_tcs TableProvider::get_sect_eff_tcs(const Vars_Cin_tcs& vc) const {
    Vars_Sef_tcs ret = {0,0,0,0,0,0};
    auto it = __sef_table_tcs.find(vc);

    if(it != __sef_table_tcs.end()) {
        ret = it->second;
    }

    return ret;
}

Vars_Sef_tcs_pol TableProvider::get_sect_eff_tcs_pol(const Vars_Cin_tcs_pol& vc) const {
    Vars_Sef_tcs_pol ret = {0,0,0,0,0,0, 0, 0};
    auto it = __sef_table_tcs_pol.find(vc);

    if(it != __sef_table_tcs_pol.end()) {
        ret = it->second;
    }

    return ret;
}


unsigned long long int TableProvider::size() const {
    return __sef_table_tcs.size();
    //return __sef_table_.size();
}
