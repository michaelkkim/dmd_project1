# dmd_project1
A new algorithm for sampling parameters in a structured correlation matrix with application to estimating optimal combinations of muscles to quantify progression in Duchenne muscular dystrophy (DMD)

### Inside `code_1/...` folder:
- `data_mcmc_analy_1.R` executes MCMC and posterior analysis for DMD data
- `functions_1.R` contains various helper functions
- `sim_datasets_1st_1.R` creates simulation data
- `sim_mcmc_2nd_1.R` executes MCMC for simulation data
- `sim_analy_3rd_1.R` executes posterior analysis for simulation data

### Inside `y_sim_1/target/...` folder:
- Posterior estimates of model parameters computed from `data_mcmc_analy_1.R` -- needed to run `sim_datasets_1st_1.R`.

### `pd_intervals_example`
- Detailed toy example to help others compute approximate positive definite (PD) bounds for correlation elements inside of their own structured correlation matrix. The key idea is to offload any computational burden from the actual MCMC by applying numerous helper functions before the MCMC when dealing with a particular structured correlation matrix.
