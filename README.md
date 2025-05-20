# dmd_project1
A new algorithm for sampling parameters in a structured correlation matrix with application to estimating optimal combinations of muscles to quantify progression in Duchenne muscular dystrophy (DMD)

### Inside _code_1/..._ folder:
- _data_mcmc_analy_1.R_ executes MCMC and posterior analysis for DMD data
- _functions_1.R_ contains various helper functions
- _sim_datasets_1st_1.R_ creates simulation data
- _sim_mcmc_2nd_1.R_ executes MCMC for simulation data
- _sim_analy_3rd_1.R_ executes posterior analysis for simulation data

### Inside _y_sim_1/target/..._ folder:
- Posterior estimates of model parameters computed from _data_mcmc_analy_1.R_ -- needed for _sim_datasets_1st_1.R_.
