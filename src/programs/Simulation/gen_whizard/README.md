18.01.2020 Igal Jaegle

To use gen_whizard:

Step 1: produce the lhe files with WHIZARD

    a/ cd run_whizard

    b/ source /work/halld/home/ijaegle/Env/custom_wo (WHIZARD only works with gcc 5 >=)

    c/ make clean; make

    d/ There are 2 methods "in prompt terminal" and "with slurm" for more details read run_whizard/README.md

    e/ it is producing the lhe files in the your chosen directory

Step 2: convert the lhe files into root files

    a/ gluex env. variables

    b/ cd lhe_to_root

    c/ make clean,make

    d/ ./lhe_to_root path to the lhe files directory/

Step 3: 

     1/ read the root file(s) produced in "Step 2" and

     2/ correct the free electron cross-section for your target choice by applying screening and radiative corrections

     3/ produce the hddm files

gen_whizard usage:

    a/ gluex env. variables

    b/ edit the beam.cfg and change 
     PhotonBeamLowEnergy 

     PhotonBeamHighEnergy

     ROOTFluxFile 

     ROOTFluxName 

     to your own chosing

    c/ edit the whizard.cfg

    d/ there 3 usages: 

        1/ "in prompt terminal" (NOT YET POSSIBLE) 

    	2/ "with slurm" (NOT YET POSSIBLE) do first 

    	swif create workflowname 

    	and the last thing do once all jobs are send

    	swif run workflowname 
    


