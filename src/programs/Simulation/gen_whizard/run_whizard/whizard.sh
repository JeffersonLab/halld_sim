#
# HallD software                                                       
# Copyright(C) 2020       GlueX and PrimEX-D Collaborations            
#                                                                      
# Author: The GlueX and PrimEX-D Collaborations                        
# Contributors: Igal Jaegle                                            
#                                                                      
# This software is provided "as is" without any warranty.              
#

source /work/halld/home/ijaegle/Env/custom_wo.sh

type=$1
evtnb=$2
egam=$3
runnb_min=$4
runnb_max=$5
seed=$6
wf=$7
path=$8

echo 'type ' $type ' evtnb ' $evtnb ' egam ' $egam ' seed ' $seed
if [ "$type" == "ae_to_ae" ]; then
    vis="\"gamma\"\,\"e\-\""
fi
if [ "$type" == "ae_to_aae" ]; then
    vis="\"gamma\"\,\"gamma\"\,\"e\-\""
fi
if [ "$type" == "ae_to_eee" ]; then
    vis="\"e\+\"\,\"e\-\"\,\"e\-\""
fi
if [ "$type" == "ae_to_aeee" ]; then
    vis="\"gamma\"\,\"e\+\"\,\"e\-\"\,\"e\-\""
fi
#seed=`od -An -N3 -l /dev/random |  sed 's/^ *\(.*\) *$/\1/'`
file=${type}_${runnb_min}_${runnb_max}_$egam

store=$path/$file
errdir=$store/err
mkdir -p $errdir
outdir=$store/out
mkdir -p $outdir

sed 's,SEED,'${seed}',g; s,VIS,'${vis}',g; s,EVTNB,'${evtnb}',g; s,EGAM,'${egam}',g' temp.sin > $store/run.sin

echo 'source /work/halld/home/ijaegle/Env/custom_wo.sh' > $store/run.sh
echo "cd ${store}" >> $store/run.sh 
echo 'whizard run.sin' >> $store/run.sh
echo "mv NAME.lhe ../${file}.lhe" >> $store/run.sh 
echo 'rm *.vg *.evx *.log *.la *.lo *.phs *.o *.mod *.f90 *.makefile .libs' >> $store/run.sh
chmod +x $store/run.sh

errfile=$errdir/stderr_wo.err
outfile=$outdir/stdout_wo.out
swif add-job -workflow $wf -project gluex -track analysis -stdout $outfile -stderr $errfile $store/./run.sh 
