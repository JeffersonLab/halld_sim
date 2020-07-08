To run mc_gen with run.sh and config files:

#### setup the environments
source $/path/to/halld_sim

#### copy run.sh to local dir
cp $HALLD_SIM_HOME/src/programs/Simulation/MC_GEN/run/run.sh .

#### execute run.sh with example.config
./run.sh $HALLD_SIM_HOME/src/programs/Simulation/MC_GEN/run/template/example.config
