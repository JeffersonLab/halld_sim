Compton event generator works in two steps

Step 1: produce the input root files
eg: gen_primex_compton -c beam.cfg -e compton.cfg -hd test.hddm -n 1000000 -r 61321

beam.cfg
PhotonBeamLowEnergy LOWE
PhotonBeamHighEnergy HIGHE
ROOTFluxFile FLUXFILE/Run_RUNMIN_RUNMAX_flux.root
ROOTFluxName FLUXNAME

compton.cfg
#Run Compton PRIMEX true/false
run: true

#What is your shell
#bash/tcsh
shell: bash

#Path to the lhe file location
#If you want to produce one single HDDM file for all lhe files contained in one directory
#leave it blank/give directory path of the lhe file
dir:

#Path to the lhe file location
#leave it blank/give the full path name of one single lhe file
file: 

#Path to output root file
out_dir: PATH

#work flow name
#leave it blank for prompt terminal usage/give a name for slurm usage
workflow: FLOW



Step 2: convert root file into hddm file 
gen_primex_compton -c beam.cfg -e compton.cfg -hd output.hddm
compton.cfg
#Run Compton PRIMEX true/false
run: false

#What is your shell
#bash/tcsh
shell: bash

#Path to the lhe file location
#If you want to produce one single HDDM file for all lhe files contained in one directory
#leave it blank/give directory path of the lhe file
dir:

#Path to the lhe file location
#leave it blank/give the full path name of one single lhe file
file: input.root

#Path to output root file
out_dir: 

#work flow name
#leave it blank for prompt terminal usage/give a name for slurm usage
workflow: 

The Fortran code is here": sd_compton and works as followed:
sd_compton [cal_param.cfg] [output] [photon-beam] [event number]
