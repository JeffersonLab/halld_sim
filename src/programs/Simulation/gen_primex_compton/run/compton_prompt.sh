#
# HallD software                                                       
# Copyright(C) 2020       HallD Group
#                                                                      
# Author: The HallD Group                        
# Contributors: Igal Jaegle                                            
#                                                                      
# This software is provided "as is" without any warranty.              
#

#source /work/halld/home/ijaegle/Env/custom_GlueX_dev.sh

path=$1
egam=$2
evtnb=$3
runnb_min=$4
runnb_max=$5
filename=$7
linenumber=$8

store=$path/run_${runnb_min}_${runnb_max}_egam_${egam}GeV

file=sdc

errdir=$store/err
mkdir -p $errdir
outdir=$store/out
mkdir -p $outdir

cp $HALLD_SIM_HOME/src/programs/Simulation/gen_primex_compton/sd_compton/bases-init.dbf $store/
config=bases-init.dbf

run=$store/$file.sh
echo '#source /work/halld/home/ijaegle/Env/custom_GlueX_dev.sh' > $run
echo "cd ${store}" >> $run
echo "sd_compton ${config} ${file} ${egam} ${evtnb} ${filename} ${linenumber} > ${file}.log" >> $run
echo "h2root ${file}.hbook" >> $run
echo "rm ${file}.bin ${file}.dat ${file}.hbook" >> $run
source $run

echo "completed"
