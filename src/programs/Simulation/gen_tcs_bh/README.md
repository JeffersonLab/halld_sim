# 20251124
Time Like Compton Scattering Event Generator

Based on Marie Boer Code:
https://hallaweb.jlab.org/wiki/index.php/DDVCS_and_TCS_event_generator
https://hallaweb.jlab.org/wiki/index.php/DEEPGen_event_generator

----------------------------------------------
            * Configuration file
----------------------------------------------

tcs_bh.cfg


..........................................................
            * To Run code
..........................................................

Zero Step: if Data directory is missing or symbolic links to data are missing: run
./set.csh

First Step: Set your desire events generator settings in the configuration file "tcs_bh.cfg"

Second Step: Be in the directory: "src" (go 3 directories back from "gen_tcs_bh" directory) and run the command
scons -j20 install
It compile the project with scons

Third Step: be in the directory: "gen_tcs_bh"  then run the command below in a terminal
../../../../Linux_Alma9-x86_64-gcc11.5.0/bin/gen_tcs_bh

copy the command; paste in terminal; press "Enter";


----------------------------------------------
            * Output files
----------------------------------------------

Output are store in: tcs_bh_output.hddm i.e name use in tcs_bh.cfg file,  tcs_bh_output_#.root and the logs files are stored in tcs_bh_output_#.log

Note:
"#" is any number set as index in the configuration file to distinguish different runs

Note:
tcs_bh_output_#.root = ROOT file. 4-vectors, kinematics, flags, event info
tcs_bh_output_#.log = text file. contains run number, total events, events in file, phase space in first line. track input file and all options
tcs_bh_output_hep_#.dat = HEP or text file for simc


----------------------------------------------
            * Models
----------------------------------------------
The TCS generator uses GPDs from VGG model.
Polarization options: unpolarized cross section, beam and/or target polarized


-------------------------------------------
            * Kinematics
-------------------------------------------

All processes are generated as a function of invariants + angular variables
Except for elastic scattering
See below in Options for the dimensions.


-------------------------------------------
            * Options
-------------------------------------------

Flag to call the reactions (set reaction flag in "tcs_bh.cfg"):
(process)

Hard exclusive reactions:
1. tcs = generation of TCS+BH 5 or 6-fold (un)polarized differential cross section for p and n weights in dt.dQ'2.dphi.dtheta

Phase-space:
11. ps_eephoto_fix = phase space generation for dilepton pair photoproduction (real photon) in fix target experiment


