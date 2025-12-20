#!/bin/csh

rm -rf Data
mkdir Data
cd Data

ln -s /work/halla/solid/mboer/public/analysis_tools/cut_tables/scanBHsing_t02_fullrange.dat scanBHsing_fullrange.dat

ln -s /work/halla/solid/mboer/grid/Jobout/merged_tcs/v11/grid_table_tcs_circlong.dat grid_table_tcs_circlong.dat
ln -s /work/halla/solid/mboer/grid/Jobout/merged_tcs/v11/grid_table_tcs_circperp.dat grid_table_tcs_circperp.dat
ln -s /work/halla/solid/mboer/grid/Jobout/merged_tcs/v11/grid_table_tcs_linlong.dat grid_table_tcs_linlong.dat

ln -s  /work/halla/solid/mboer/grid/Jobout/merged_tcs/v11/grid_table_tcs_circlong_neutron.dat grid_table_tcs_neutron_circlong.dat
ln -s  /work/halla/solid/mboer/grid/Jobout/merged_tcs/v11/grid_table_tcs_circperp_neutron.dat grid_table_tcs_neutron_circperp.dat
ln -s  /work/halla/solid/mboer/grid/Jobout/merged_tcs/v11/grid_table_tcs_linlong_neutron.dat grid_table_tcs_neutron_linlong.dat

cd ..
