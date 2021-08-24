# HawkesCausalInfection
Using Hawkes Processes to identify cases in infectious disease outbreaks with no known causal infection

Cluster based simulation of a Hawkes process done using Simulation_Cluster.R.

Inputs:
Events - any events at time 0 (vector)
T_max - maximum time of simulation (float)
N_max - maximum number of events used as a safety in case of especially long runs. Should be specified higher than anticipated need. (integer)
mu_fn - the background intensity to be simulated from. Use one of Background_Intensities.R (function)
phi_fn - the triggering function to use. Use one of Kernel_Times.R (function)
branching_factor - the branching factor of the chosen triggering function. Use one of Branching_Factors.R (function)
parameters - the parameters of the process to be simulated from. To be user specified (list)

EM for Rayleigh kernel and Exponential kernel can be done using their specific code options.

Inputs:
starting_values -  initial values of parameters of process (list)
times - time points of point process to be estimated (vector)
T_max - maximum time point of process (float)
thresh - threshold of sensitivity of process. 1e-4 often used. (float)
