#
# HallD software                                                       
# Copyright(C) 2020       GlueX and PrimEX-D Collaborations            
#                                                                      
# Author: The GlueX and PrimEX-D Collaborations                        
# Contributors: Igal Jaegle                                            
#                                                                      
# This software is provided "as is" without any warranty.              
#
# NB: only work gcc 7.2, do not forget to source
source /work/halld/home/ijaegle/Env/custom_wo.sh
# Usage: ./go_all.sh nbofevt ROOTFluxFile ROOTFluxName RunNbMin RunNbMax output_dir workflow_name
#

nbofevt=$1
lowe=$2
highe=$3
fluxfile=$4
fluxname=$5
runmin=$6
runmax=$7
out_dir=$8
wf=$9
myshell=$10

swif create $wf 

for fs in ae aae eee aeee; do
    genconfig=whizard_$fs.cfg
    sed 's,FS,'${fs}',g; s,LOWE,'${lowe}',g; s,HIGHE,'${highe}',g; s,FLUXFILE,'${fluxfile}',g; s,FLUXNAME,'${fluxname}',g; s,RUNMIN,'${runmin}',g; s,RUNMAX,'${runmax}',g' temp.cfg > $genconfig
    ./run_whizard $nbofevt $genconfig $runmin $runmax $wf $out_dir $myshell 
done

swif run $wf
