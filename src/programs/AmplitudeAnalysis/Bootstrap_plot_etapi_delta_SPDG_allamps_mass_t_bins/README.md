This code can be used for calculating the uncertainties from bootstrapping on the moments calculated 
using project_moments_polarized code.
The process of doing fits for N bootstrapping samples by submitting jobs is described in details in 
hd_utilities/PWA_scripts/Bootstrapping_M_t_bins .
Once the jobs are succesfully completed one can use Bootstrap_plot_etapi_delta_SPDG_allamps_mass_t_bins
and run it doing "Bootstrap_plot_etapi_delta_SPDG_allamps_mass_t_bins -o etapi_fit.txt" in the directory 
that contains EtaPi_fit directory with fit results. This will calculate the moments for given bootstrapping sample
and givent Ma nd t bin and write it to a file "etapi_fit.txt", where the first line will include the name of 
the variable in each column.
After this one can plot the moments using a python code that I will add in hd_utilities.
