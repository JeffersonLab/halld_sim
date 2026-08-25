#ifndef TABLEPROVIDER_H
#define TABLEPROVIDER_H

#include <map>
#include <string>
#include "TabUtils.h"

typedef std::pair<Vars_Cin_tcs, Vars_Sef_tcs> Table_Element_tcs;
typedef std::map<Vars_Cin_tcs, Vars_Sef_tcs> Cross_Section_Table_tcs;

typedef std::pair<Vars_Cin_tcs_pol, Vars_Sef_tcs_pol> Table_Element_tcs_pol;
typedef std::map<Vars_Cin_tcs_pol, Vars_Sef_tcs_pol> Cross_Section_Table_tcs_pol;

class TableProvider {
    public:
	bool insert_tcs(const Vars_Cin_tcs& key, const Vars_Sef_tcs& value);
	bool insert_tcs_pol(const Vars_Cin_tcs_pol& key, const Vars_Sef_tcs_pol& value);
	Vars_Sef_tcs get_sect_eff_tcs(const Vars_Cin_tcs& vc) const;
	Vars_Sef_tcs_pol get_sect_eff_tcs_pol(const Vars_Cin_tcs_pol& vc) const;
        unsigned long long int size() const;

    private:
//        Cross_Section_Table __sef_table;
	Cross_Section_Table_tcs __sef_table_tcs;
	Cross_Section_Table_tcs_pol __sef_table_tcs_pol;
};

#endif /* TABLEPROVIDER_H */
