%README 
%
%      All the examples have been run with Matlab R2016b. 

%      The main script is main.m, where the user specifies a gap width in 
%      nanometers, an incident wavelength in microns and a value for the b 
%      factor that multiplies beta^2 (b=0 is LRA, and b=1,1.5 have been used 
%      for the paper). This is the only routine that needs to be executed
%      by the user. All the other routines are in a subfolder that is
%      automatically added to the path
%      
%      run_example.m loads the mesh used for the simulations and calls
%      the maxwell solver. Runtimes of each subpart are shown in screen. If
%      b>0 the hydrodynamic model is solved, which requires more RAM than 
%      the b = 0 case (all the examples can be run in a machine with 128GB 
%      of RAM).
%
%      The solver outputs the values of the electric (EDG), magnetic (HDG),
%      current density (JDG) and electron density (RDG) at the high-order
%      nodes of the discretization as complex numbers. In run_example.m we
%      then compute the transmission and the field enhancement as specified
%      in the SI.
%
%      Back to main.m, we also provide a routine for plotting any of the
%      above solution fields or modifications thereof (e.g. real part,
%      absolute magnitude) on any z-constant surface. The user must specify
%      the component of the solution field to plot as well as the z-value 
%      where to plot (if it does not coincide with a z-value that belongs to the
%      discretization, we instead show the solution field at the z that is
%      closest to z-value and is in the discretization). A couple of
%      examples are provided. 
%
  