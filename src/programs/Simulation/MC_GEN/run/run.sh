 #########################################################
 #  Usage: ./run.sh 
 ########################################################
 #                      # grab beam config files        
 #  flux.cc             # invoke beamProperties class        
 #                      # translate the beam energy distribution into ASCII file
 ########################################################
 # 
 # created by:  Hao Li
 #              Carnegie Mellon University
 #              15-Feb-2020
 ######################################################## 

# Run number
RUN_NUMBER=31057

# Reaction Mechanisms --> pointing to the specific generator's
#                         definition template and the correspo
#						  -nding translator
# MECH CODE:
#            PPBAR_MECH_0 -- ppbar 3body phase space
#            PPBAR_MECH_1 -- ppbar t-channel (single Regge)
#            PPBAR_MECH_2 -- ppbar u-channel (single Regge)
#            PPBAR_MECH_3 -- ppbar (double Regge)
#            LAMLAMBAR_MECH_0 -- lamlambar 3body phase space
#            LAMLAMBAR_MECH_1 -- lamlambar mechanism 1
#            LAMLAMBAR_MECH_2 -- lamlambar mechanism 2
#            LAMLAMBAR_MECH_3 -- lamlambar mechanism 3
#            LAMLAMBAR_MECH_4 -- lamlambar mechanism 4
#            LAMLAMBAR_MECH_5 -- lamlambar (double Regge)
REACTION_CHANNEL=lamlambar
MECH_CODE=LAMLAMBAR_MECH_5 

# Model-related Parameters
DYNCODE0=0.53
DYNCODE1=0.24
DYNCODE2=0.7

# Simulation Parameters
TOTALSIZE=10000
SAMPLESIZE=5000
MOMENTUM_MIN=6.4
MOMENTUM_MAX=11.4

#FLUX_DIR=`printf '/home/haoli/build/test/MC_GEN_sdobbs/halld_sim/src/programs/Simulation/MC_GEN/run/flux_%d_%d.ascii' "${RUN_NUMBER}" "${RUN_NUMBER}"`
#FLUX_DIR=/home/haoli/build/test/MC_GEN_sdobbs/halld_sim/src/programs/Simulation/MC_GEN/run/flux_30274_31057.ascii
FLUX_DIR=./flux_30274_31057.ascii
GEN_DIR=/home/haoli/build/test/MC_GEN_sdobbs/halld_sim/src/programs/Simulation/MC_GEN/run/build/test/MC_GEN_sdobbs/halld_sim/src/programs/Simulation/MC_GEN/run/test/MC_GEN_sdobbs/halld_sim/src/programs/Simulation/MC_GEN/run/test/MC_GEN_sdobbs/halld_sim/src/programs/Simulation/MC_GEN/run/test/MC_GEN_sdobbs/halld_sim/src/programs/Simulation/MC_GEN/run

echo $FLUX_DIR
echo $GEN_DIR

# Get tagged flux hist from ccdb
python $HD_UTILITIES_HOME/psflux/plot_flux_ccdb.py -b ${RUN_NUMBER} -e ${RUN_NUMBER}

# Translate the hist into ASCII format
ROOTSCRIPT=`printf 'Flux_to_Ascii.C("flux_%s_%s.root")' "$RUN_NUMBER" "$RUN_NUMBER" `
root -l -b -q $ROOTSCRIPT

# Edit the generator's definition file
cp template/gen_template.def gen.def
EDITLIST=DYNCODE0,DYNCODE1,DYNCODE2,TOTALSIZE,SAMPLESIZE,MOMENTUM_MIN,MOMENTUM_MAX,FLUX_DIR,GEN_DIR
for KEYWORD in ${EDITLIST//,/ }
	do
		REGEXSTRING=`printf 's=%s=%s=g' "${KEYWORD}" "${!KEYWORD}"`
		echo ${REGEXSTRING}
		sed -i ${REGEXSTRING} gen.def
	done


# Initiate the generator with its definition file
#         Note: out put will be an ASCII file
mc_gen gen.def

# Translate the output into HDDM format
END=$(($TOTALSIZE/SAMPLESIZE))
for i in $(seq 1 $END)
do
	FILENAME=`printf '%sgen_%04d.ascii' "${GEN_DIR}" "$i"`
	GEN2HDDM_lamlambar -r${RUN_NUMBER} $FILENAME
done










