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
source /work/halld/home/ijaegle/Env/custom_wo.csh
# Usage: ./go_all.sh nbofevt ROOTFluxFile ROOTFluxName RunNbMin RunNbMax output_dir workflow_name
#

set nbofevt=$1
set lowe=$2
set highe=$3
set fluxfile=$4
set fluxname=$5
set runmin=$6
set runmax=$7
set out_dir=$8
set wf=$9
set myshell=$10

swif create $wf 

foreach fs (ae aae eee aeee)
    set genconfig=whizard_$fs.cfg
    sed 's,FS,'${fs}',g; s,LOWE,'${lowe}',g; s,HIGHE,'${highe}',g; s,FLUXFILE,'${fluxfile}',g; s,FLUXNAME,'${fluxname}',g; s,RUNMIN,'${runmin}',g; s,RUNMAX,'${runmax}',g' temp.cfg > $genconfig
    ./run_whizard $nbofevt $genconfig $runmin $runmax $wf $out_dir $myshell 
done

swif run $wf
