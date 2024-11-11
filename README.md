Run Eigen_Solver.m first to find the linear solution for the atmospheric equations. 
The eigenvectors outputted from this code are used in the other codes to help reconstruct the MJO signal. 
Run SST_bivar.m to produce the bivariate regression coefficients for SSTAs. 
main_modeling_code.m calls Eigen_Solver and SST_bivar before producing the coupled model simulation. 
main_modeling_code.m calls ENSO_statistics_results_code.m at the end which is used to produce all of the statistics for the ENSO state variables.
After running main_modeling_code.m, running MJO_spectrum will produce the power spectra of the MJO variables from the model simulation and observation.
wind_stats.m  produces the statistic for observation and model simulation.

Before running MJO_std.m, plot_MJO_SST_obs.m, or conditional_correlation.m run calc_MJO_obs.m to construct the data set for MJO from observational data.
MJO_std.m plots the standard deviation of MJO conditioned on different ENSO events as a function of longitude. 
plot_MJO_SST.m plots the SST and MJO Hovmoller diagrams  produced by the conceptual model. plot_MJO_SST_obs.m plots the same but for observations. 
conditional_correlation.m plots the lagged correlation between MJO and SST conditioned on different ENSO events for both the model simulation and observation.
