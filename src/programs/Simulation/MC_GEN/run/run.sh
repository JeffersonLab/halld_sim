 #!/bin/sh 
 #########################################################
 #  Usage: ./run.sh /path/to/config
 ########################################################
 #  grab generator parameters from config file        
 #  invoke plot_flux_ccdb.py to get tagged_flux.root 
 #  translate flux into ASCII file as mc_gen's beam profile
 #  create new mc_gen definition file from given parameters
 #  run mc_gen with ASCII output
 #  translate ASCII into HDDM
 ########################################################
 #              Hao Li
 #  created by: Carnegie Mellon University
 #              15-Feb-2020
 ######################################################## 



CONFIGFILE=$@
echo "CONFIGFILE = " $CONFIGFILE
eval $(sed '/:/!d;/^ *#/d;s/:/ /;' < "$CONFIGFILE" | while read -r key val
do
    #verify here
    #...
    str="$key='$val'"
    echo "$str"
done)


# Prepare
echo "Print list of control parameters:"
echo "RUN_NUMBER = " $RUN_NUMBER
echo "REACTION_CHANNEL = " $REACTION_CHANNEL
echo "MECH = " $MECH
echo "DYNCODE0 = " $DYNCODE0
echo "DYNCODE1 = " $DYNCODE1
echo "DYNCODE2 = " $DYNCODE2
echo "TOTALSIZE = " $TOTALSIZE
echo "SAMPLESIZE = " $SAMPLESIZE
echo "MOMENTUM_MIN = " $MOMENTUM_MIN
echo "MOMENTUM_MAX = " $MOMENTUM_MAX
echo "GEN_DIR = " ${GEN_DIR}
mkdir -p ${GEN_DIR}
echo $FILE


# Get tagged flux hist from ccdb
python $HD_UTILITIES_HOME/psflux/plot_flux_ccdb.py -b ${RUN_NUMBER} -e ${RUN_NUMBER}
FLUX_DIR=`printf './flux_%d_%d.ascii' "${RUN_NUMBER}" "${RUN_NUMBER}"`
echo "FLUX_DIR = " $FLUX_DIR 

# Translate the hist into ASCII format
ROOTSCRIPT=`printf '$HALLD_SIM_HOME/src/programs/Simulation/MC_GEN/run/Flux_to_Ascii.C("flux_%s_%s.root")' "$RUN_NUMBER" "$RUN_NUMBER" `
root -l -b -q $ROOTSCRIPT

# Edit the generator's definition file
cp $HALLD_SIM_HOME/src/programs/Simulation/MC_GEN/run/template/gen_${MECH}.def $GEN_DIR/standard_name.configuration
EDITLIST=DYNCODE0,DYNCODE1,DYNCODE2,TOTALSIZE,SAMPLESIZE,MOMENTUM_MIN,MOMENTUM_MAX,FLUX_DIR,GEN_DIR
for KEYWORD in ${EDITLIST//,/ }
	do
		REGEXSTRING=`printf 's=%s=%s=g' "${KEYWORD}" "${!KEYWORD}"`
		#echo ${REGEXSTRING}
		sed -i ${REGEXSTRING} $GEN_DIR/standard_name.configuration
	done


# Initiate the generator with its definition file
#         Note: out put will be an ASCII file
mc_gen $GEN_DIR/standard_name.configuration

# Translate the output into HDDM format
END=$(($TOTALSIZE/SAMPLESIZE))
for i in $(seq 1 $END)
do
	FILENAME=`printf '%sstandard_name_%04d.ascii' "${GEN_DIR}" "$i"`
	if [ "$REACTION_CHANNEL" == "lamlambar" ]; then
		GEN2HDDM_lamlambar -r${RUN_NUMBER} $FILENAME
	fi

	if [ "$REACTION_CHANNEL" == "ppbar" ]; then
		GEN2HDDM_ppbar -r${RUN_NUMBER} $FILENAME
	fi


done


#clear files
mv flux_*.ascii ${GEN_DIR}.
rm flux_*.root
rm randomseed.num









