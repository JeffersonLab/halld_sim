Basic scripts in order to simulate events using gen_amp_V2. 
# Getting Started 

## Item needed
  * config file (see 'genPS.cfg'). 

## config file 
   The config file uses the same structure as the gen_amp config file. The basic feature needed to simulate events are the reaction name: 
   * EX: 'reaction omegapi Beam Proton Pi0 Pi+ Pi-'

## How to use
   * In this example we generate 10k four-body phase phase where the π^+^ π^-^ π^0^ are in the upper vertex and the proton is in the lower vertex.  
   * run 'gen_amp_V2 -ac genPS.cfg -o test.root -uv 234 -lv 1 -f'.
 Things to note:
   * Unless a beam config file is used, a local file will be produced within the coherent peak.
   * The values from the "uv & lv' flags indicate the particle index from the reaction name.
   * The '-f' flag skips the amplitude accept/reject process and generates phase space.
   * Two output files will be created, 'test.root' and 'gen_amp_diagnostic.root'.

## Documentation
    For further options and descriptions see 'GlueX_Fixed_Target_Generator.pdf'.

## Using in MCWrapper
 
  When using this generator in MCWrapper one needs to place the flags in the first line of the config file.
  Using the same example on the first line '# -uv 234 -lv 1 -f'.
  The '-ac -bc -o -n -r -m -p -a -b' flags will be defined using the MCWrapper config file.


