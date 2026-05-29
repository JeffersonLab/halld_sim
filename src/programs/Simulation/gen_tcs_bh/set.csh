#!/bin/csh

rm -rf Data
mkdir Data
cd Data
# create data folder with low Q'2 (LQ2) region cross sections  tables
ln -f -s /work/halla/solid/mboer/public/analysis_tools/cut_tables/scanBHsing_t02_allQ_fullrange.dat scanBHsing_fullrange_LQ2.dat

ln -f -s /work/halla/solid/mboer/grid/Jobout/merged_tcs/v26/grid_table_tcs_circlong.dat grid_table_tcs_circlong_LQ2.dat
ln -f -s /work/halla/solid/mboer/grid/Jobout/merged_tcs/v26/grid_table_tcs_circperp.dat grid_table_tcs_circperp_LQ2.dat
ln -f -s /work/halla/solid/mboer/grid/Jobout/merged_tcs/v26/grid_table_tcs_linlong.dat grid_table_tcs_linlong_LQ2.dat

ln -f -s /work/halla/solid/mboer/grid/Jobout/merged_tcs/v26/grid_table_tcs_neutron_circlong.dat grid_table_tcs_neutron_circlong_LQ2.dat
ln -f -s /work/halla/solid/mboer/grid/Jobout/merged_tcs/v26/grid_table_tcs_neutron_circperp.dat grid_table_tcs_neutron_circperp_LQ2.dat
ln -f -s /work/halla/solid/mboer/grid/Jobout/merged_tcs/v26/grid_table_tcs_neutron_linlong.dat grid_table_tcs_neutron_linlong_LQ2.dat

# create data folder with high Q'2 (HQ2) region cross sections  tables
ln -f -s /work/halla/solid/mboer/public/analysis_tools/cut_tables/scanBHsing_t02_fullrange.dat scanBHsing_fullrange_HQ2.dat

ln -f -s /work/halla/solid/mboer/grid/Jobout/merged_tcs/v11/grid_table_tcs_circlong.dat grid_table_tcs_circlong_HQ2.dat
ln -f -s /work/halla/solid/mboer/grid/Jobout/merged_tcs/v11/grid_table_tcs_circperp.dat grid_table_tcs_circperp_HQ2.dat
ln -f -s /work/halla/solid/mboer/grid/Jobout/merged_tcs/v11/grid_table_tcs_linlong.dat grid_table_tcs_linlong_HQ2.dat

ln -f -s  /work/halla/solid/mboer/grid/Jobout/merged_tcs/v11/grid_table_tcs_circlong_neutron.dat grid_table_tcs_neutron_circlong_HQ2.dat
ln -f -s  /work/halla/solid/mboer/grid/Jobout/merged_tcs/v11/grid_table_tcs_circperp_neutron.dat grid_table_tcs_neutron_circperp_HQ2.dat
ln -f -s  /work/halla/solid/mboer/grid/Jobout/merged_tcs/v11/grid_table_tcs_linlong_neutron.dat grid_table_tcs_neutron_linlong_HQ2.dat

# create data folder with full Q'2 (FQ2) region cross sections tables: Names to be use when table becomes available
# scanBHsing_fullrange_FQ2.dat

# grid_table_tcs_circlong_FQ2.dat
# grid_table_tcs_circperp_FQ2.dat
# grid_table_tcs_linlong_FQ2.dat

# grid_table_tcs_neutron_circlong_FQ2.dat
# grid_table_tcs_neutron_circperp_FQ2.dat
# grid_table_tcs_neutron_linlong_FQ2.dat

cd ..
