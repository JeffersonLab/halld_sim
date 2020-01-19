18.01.2020 Igal Jaegle

Usage: in prompt terminal 
      Method 1: prompt terminal
      	       ./go_all_prompt.csh nbofevt LowEnergy HighEnbergy ROOTFluxFile ROOTFluxName RunNbMin RunNbMax output_dir your_shell
      Method 2: with slurm	       
              ./go_all_slurm.csh nbofevt LowEnergy HighEnbergy ROOTFluxFile ROOTFluxName RUNMIN RUNMAX workflow_name output_dir your_shell


NB: a/ ROOTFluxFile: is the path the directory containing the root file with the flux histogram
    b/ The rootfile name MUST BE: Run_RUNMIN_RUNMAX_flux.root
    c/ ROOTFluxName: histogram name, it can be whatever the user choose
    d/ You must precise your $SHELL ie either "bash" or "tcsh"
