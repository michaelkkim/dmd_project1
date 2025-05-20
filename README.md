# dmd_project1
A new algorithm for sampling parameters in a structured correlation matrix with application to estimating optimal combinations of muscles to quantify progression in Duchenne muscular dystrophy

### Inside _code_1/..._ folder:
- _data_mcmc_analy_1.R_ executes MCMC and posterior analysis for Duchenne muscular dystrophy (DMD) data (data not shared)
- _functions_1.R_ contains various helper functions
- _sim_datasets_1st_1.R_ creates datasets needed for simulation
- _sim_mcmc_2nd_1.R_ executes MCMC for simulation dataset(s) (same MCMC details as _data_mcmc_analy_1.R_)
- _sim_analy_3rd_1.R_ executes posterior analysis for simulation dataset(s)

### Inside _y_sim_1/..._ folder:
- _target/..._ folder contains posterior estimates computed from _data_mcmc_analy_1.R_, which is needed to run 'sim_datasets_1st_1.R'
