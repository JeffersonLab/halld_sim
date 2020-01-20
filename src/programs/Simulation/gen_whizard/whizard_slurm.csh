#!/bin/csh -f
# HallD software                                                       
# Copyright(C) 2020       GlueX and PrimEX-D Collaborations            
#                                                                      
# Author: The GlueX and PrimEX-D Collaborations                        
# Contributors: Igal Jaegle                                            
#                                                                      
# This software is provided "as is" without any warranty.              
#
set type=$1
set evtnb=$2
set egam=$3
set runnb_min=$4
set runnb_max=$5
set seed=$6
set pathdir=$7

if ( "$type" == "ae_to_ae" ) then
    set vis='\"gamma\"\,\"e\-\"'
else if ( "$type" == "ae_to_aae" ) then
    set vis='\"gamma\"\,\"gamma\"\,\"e\-\"'
else if ( "$type" == "ae_to_eee" ) then
    set vis='\"e\+\"\,\"e\-\"\,\"e\-\"'
else if ( "$type" == "ae_to_aeee" ) then
    set vis='\"gamma\"\,\"e\+\"\,\"e\-\"\,\"e\-\"'
endif

set file=${type}_${runnb_min}_${runnb_max}_${egam}
set store=$path/$file
set errdir=$store/err
mkdir -p $errdir
set outdir=$store/out
mkdir -p $outdir

sed 's,SEED,'${seed}',g; s,VIS,'${vis}',g; s,EVTNB,'${evtnb}',g; s,EGAM,'${egam}',g' temp.sin > $store/run.sin

echo 'source /work/halld/home/ijaegle/Env/custom_wo.csh' > $store/run.csh
echo "cd ${store}" >> $store/run.csh 
echo 'whizard run.sin' >> $store/run.csh
echo "mv NAME.lhe ../${file}.lhe" >> $store/run.csh 
echo 'rm *.vg *.evx *.log *.la *.lo *.phs *.o *.mod *.f90 *.makefile' >> $store/run.csh
echo 'cd ..'  >> $store/run.csh
echo "rm -rf ${file}" >> $store/run.csh
chmod +x $store/run.csh

set errfile=$errdir/stderr_wo.err
set outfile=$outdir/stdout_wo.out
swif add-job -workflow $wf -project gluex -track analysis -stdout $outfile -stderr $errfile $store/./run.sh 

echo "completed"
