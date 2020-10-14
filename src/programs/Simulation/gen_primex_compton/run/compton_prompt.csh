#!/bin/csh -f
# HallD software                                                       
# Copyright(C) 2020       HallD Group
#                                                                      
# Author: The HallD Group                        
# Contributors: Igal Jaegle                                            
#                                                                      
# This software is provided "as is" without any warranty.              
#

#source /work/halld/home/ijaegle/Env/custom_GlueX_dev.sh

set path=$1
set egam=$2
set evtnb=$3
set runnb_min=$4
set runnb_max=$5
set wf=$6
set filename=$7
set linenumber=$8

set store=$path/run_${runnb_min}_${runnb_max}_egam_${egam}GeV

set file=sdc

set errdir=$store/err
mkdir -p $errdir
set outdir=$store/out
mkdir -p $outdir

cp $HALLD_SIM_HOME/src/programs/Simulation/gen_primex_compton/sd_compton/bases-init.dbf $store/
config=bases-init.dbf

set run=$store/$file.sh
echo '#source /work/halld/home/ijaegle/Env/custom_GlueX_dev.sh' > $run
echo "cd ${store}" >> $run
echo "sd_compton ${config} ${file} ${egam} ${evtnb} ${filename} ${linenumber}> ${file}.log" >> $run
echo "h2root ${file}.hbook" >> $run
echo "rm ${file}.bin ${file}.dat ${file}.hbook" >> $run
source $run

echo "completed"
