### cntrl+find "MODIFY AS NEEDED"
set.seed(1)

library(mvtnorm) # 'dmvnorm'
library(invgamma) # 'dinvgamma'
library(abind) # 'abind'

library(parallel) # 'mclapply'
library(bayesplot) # things like 'mcmc_trace'
library(rstan) # thing(s) like 'Rhat'
library(bayestestR) # 'ci' for credible intervals (use 'hdi' in most cases)

library(ggplot2) # 'ggplot'
library(reshape) # 'melt' (needed for 'ggplot')
library(shiny) # I think this is to save html

thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript (shell script)
    return(dirname(normalizePath(sub(needle, "", cmdArgs[match]))))
  } else if(Sys.getenv("RSTUDIO")=="1") {
    # Rstudio
    return(dirname(rstudioapi::getSourceEditorContext()$path))
  }
}

amb_names <- c("Early","Late","Non")

code_dir <- thisFile()
setwd(code_dir)

parent_source <- dirname(code_dir)
base_name <- basename(code_dir)
sub_base_length <- nchar(base_name)-nchar("code")
main_base_name <- substr(base_name,
                         1, nchar(base_name)-sub_base_length)
sub_base_name <- substr(base_name,
                        nchar(base_name)-sub_base_length+1, nchar(base_name))
output_dir <- paste0(parent_source,"/","Soutputs",sub_base_name)
y_sim_dir <- paste0(parent_source,"/","y_sim",sub_base_name)

source(paste0("functions",sub_base_name,".R"))



### MODIFY AS NEEDED
amb_indx <- 1 # 1, 2, 3 (L=4); # 1, 2 (L=5)
muscle_names_simple <- c("SOL","VL","BB","DEL") # "6MWT"
muscle_names_simple_subset <- muscle_names_simple
L <- length(muscle_names_simple)
n_sim_data <- 500
nchains <- 1
post_center <- "post_med" # post_med, post_mean
optim_method <- "L-BFGS-B"
grad_bool <- TRUE
upper_bool <- TRUE
amb_name <- amb_names[amb_indx]

if ((amb_indx==1) && (L==4)) {
  output_sub_name <- paste0("Soutputs",sub_base_name,"sa")
} else if ((amb_indx==2) && (L==4)) {
  output_sub_name <- paste0("Soutputs",sub_base_name,"sb")
} else if ((amb_indx==3) && (L==4)) {
  output_sub_name <- paste0("Soutputs",sub_base_name,"sc")
} else if ((amb_indx==1) && (L==5)) {
  output_sub_name <- paste0("Soutputs",sub_base_name,"sd")
} else if ((amb_indx==2) && (L==5)) {
  output_sub_name <- paste0("Soutputs",sub_base_name,"se")
}



output_sub_dir <- paste0(output_dir,"/",output_sub_name)
for (sim_dat_indx in 1:n_sim_data) {
  if (!file.exists(paste0(output_sub_dir,"/","RDS","/",
                          "mcmc",sim_dat_indx,".RDS"))) {
    print(sim_dat_indx)
  }
}

mult_chains <- readRDS(paste0(output_sub_dir,"/","RDS","/",
                              "mcmc",1,".RDS"))
N <- mult_chains[[1]]$N
iters <- mult_chains[[1]]$iters
burn_in_prop <- mult_chains[[1]]$burn_in_prop
burn_in_number <- iters*burn_in_prop
r_jump_dist_name <- mult_chains[[1]]$r_jump_dist_name
r_updates <- mult_chains[[1]]$r_updates
Js <- mult_chains[[1]]$Js
ps <- mult_chains[[1]]$ps
L <- mult_chains[[1]]$L
Q <- mult_chains[[1]]$Q
J_max <- mult_chains[[1]]$J_max
p_max <- mult_chains[[1]]$p_max
j_order0_max <- mult_chains[[1]]$j_order0_max
l_order0_max <- mult_chains[[1]]$l_order0_max
nu <- mult_chains[[1]]$nu
kappa <- mult_chains[[1]]$kappa
mu_names_detailed <- mult_chains[[1]]$mu_names_detailed
s_names_detailed <- mult_chains[[1]]$s_names_detailed
r_names_detailed <- mult_chains[[1]]$r_names_detailed
uniqueR_locs_upper <- mult_chains[[1]]$uniqueR_locs_upper
uniqueR_locs_all <- mult_chains[[1]]$uniqueR_locs_all



y_sim_dir_target <- paste0(y_sim_dir,"/","target")
target_params <- readRDS(paste0(y_sim_dir_target,"/",
                                amb_name,"_",L,"L","_",post_center,".RDS"))
target_short_mu <- target_params$target_short_mu
target_short_s <- target_params$target_short_s
target_r <- target_params$target_r
target_w <- target_params$target_w
target_SRM <- target_params$target_SRM
w_names_detailed <- paste0("w",1:L,"[",muscle_names_simple,"]")










setwd(output_sub_dir)
output_w_and_SRM_RDS_dir <- paste0(output_sub_dir,"/","w&SRM_RDS")
dir.create(output_w_and_SRM_RDS_dir)
setwd(output_w_and_SRM_RDS_dir)

n_etas <- choose(L,2)
muscle_combns <- combn(1:L,2)
eta_names <- sapply(1:n_etas, function(combn_indx)
  paste0("eta_{",muscle_combns[1,combn_indx],muscle_combns[2,combn_indx],"}"))
rho_names <- paste0("rho_{(",1:L,")}")
r_names <- c(eta_names, rho_names, "gamma")

save_sim_w_and_SRM_RDS_func(n_sim_data, output_sub_dir,
                            nchains, burn_in_prop,
                            muscle_names_simple, muscle_names_simple_subset,
                            r_names,
                            optim_method, grad_bool, upper_bool,
                            w_names_detailed)

setwd(code_dir)



# Get coverage for weights and SRM
w_cover_count <- rep(0,L)
SRM_cover_count <- 0
cred_int_lengths <- matrix(nrow=n_sim_data,ncol=L+1)
for (sim_dat_indx in 1:n_sim_data) {
  w_and_SRM_temp <- readRDS(paste0(output_w_and_SRM_RDS_dir,"/",
                                   "w&SRM",sim_dat_indx,".RDS"))
  weights_temp <- w_and_SRM_temp$SRMw
  SRM_temp <- w_and_SRM_temp$SRMvalues
  
  cred_ints <- matrix(nrow=L,ncol=2)
  for (l in 1:L) {
    cred_ints[l,] <- quantile(weights_temp[,l], probs=c(0.025,0.975), names=FALSE, type=8)
    
    lower_bound <- cred_ints[l,1]
    upper_bound <- cred_ints[l,2]
    if ((target_w[l] >= lower_bound) && (target_w[l] <= upper_bound)) {
      w_cover_count[l] <- w_cover_count[l] + 1
    }
    cred_int_lengths[sim_dat_indx,l] <- upper_bound-lower_bound
  }
  
  cred_int_SRM <- quantile(SRM_temp, probs=c(0.025,0.975), names=FALSE, type=8)
  lower_bound_SRM <- cred_int_SRM[1]
  upper_bound_SRM <- cred_int_SRM[2]
  if ((target_SRM >= lower_bound_SRM) && (target_SRM <= upper_bound_SRM)) {
    SRM_cover_count <- SRM_cover_count + 1
  }
  cred_int_lengths[sim_dat_indx,L+1] <- upper_bound_SRM-lower_bound_SRM
}
w_and_SRM_coverages <- 
  matrix(c(w_cover_count/n_sim_data, SRM_cover_count/n_sim_data), ncol=L+1)
colnames(w_and_SRM_coverages) <- c(paste0(w_names_detailed," coverage"),
                                   "SRM coverage")
cred_int_lengths_avg <- colMeans(cred_int_lengths)
cred_length_avg <- matrix(cred_int_lengths_avg, ncol=L+1, byrow=TRUE)
colnames(cred_length_avg) <- c(paste0(w_names_detailed),
                               "SRM")
setwd(output_sub_dir)
write.csv(w_and_SRM_coverages,"w&SRM_cover.csv")
write.csv(cred_length_avg,"cred_length_avg.csv")
setwd(code_dir)










n_etas <- choose(L,2)
setwd(output_sub_dir)
output_w_and_SRM_from_post_med_RDS_dir <- paste0(output_sub_dir,"/","w&SRM_from_post_med_RDS")
dir.create(output_w_and_SRM_from_post_med_RDS_dir)
setwd(output_w_and_SRM_from_post_med_RDS_dir)

save_sim_w_and_SRM_from_post_med_RDS_func(n_sim_data, output_sub_dir,
                                          burn_in_number, iters, n_etas,
                                          optim_method, grad_bool, upper_bool)

setwd(code_dir)



# Get bias and MSE for weights and SRM
w_and_SRMval_from_post_med_sims <- matrix(nrow=n_sim_data,ncol=L+1)
for (sim_dat_indx in 1:n_sim_data) {
  w_and_SRMval_from_post_med_temp <- 
    readRDS(paste0(output_w_and_SRM_from_post_med_RDS_dir,"/",
                   "w&SRM_from_post_med",sim_dat_indx,".RDS"))
  w_from_post_med_temp <- w_and_SRMval_from_post_med_temp$w_from_post_med_temp
  SRMval_from_post_med_temp <- w_and_SRMval_from_post_med_temp$SRMval_from_post_med_temp
  
  w_and_SRMval_from_post_med_sims[sim_dat_indx,] <- c(w_from_post_med_temp, SRMval_from_post_med_temp)
}
target_w_and_SRM <- c(target_w,target_SRM)
bias <- colMeans(w_and_SRMval_from_post_med_sims)-target_w_and_SRM
relative_bias <- bias/target_w_and_SRM
mse <- sapply(1:(L+1), function(l)
  mean((w_and_SRMval_from_post_med_sims[,l]-target_w_and_SRM[l])^2) )
root_mse <- sqrt(mse)
bias_and_mse <- rbind(bias, relative_bias, mse, root_mse)
rownames(bias_and_mse) <- c("bias", "rel_bias", "mse", "root_mse")
colnames(bias_and_mse) <- c(paste0(w_names_detailed),
                            "SRM")
setwd(output_sub_dir)
write.csv(bias_and_mse,"bias&mse.csv")
setwd(code_dir)










setwd(output_sub_dir)
output_sub_r_accept_dir <- paste0(output_sub_dir,"/","r_accept%")
dir.create(output_sub_r_accept_dir)
setwd(output_sub_r_accept_dir)
save_r_accept_RDS_func(n_sim_data, output_sub_dir,
                       nchains, burn_in_prop)
setwd(code_dir)

setwd(output_sub_dir)
output_sub_R_PD_dir <- paste0(output_sub_dir,"/","%times_R_PD")
dir.create(output_sub_R_PD_dir)
setwd(output_sub_R_PD_dir)
save_R_PD_RDS_func(n_sim_data, output_sub_dir,
                   nchains, burn_in_prop)
setwd(code_dir)



# Get averages across simulations of r_accept% and %times_R_PD
r_accept_all_sim <- matrix(nrow=n_sim_data, ncol=Q)
R_PD_all_sim <- matrix(nrow=n_sim_data, ncol=Q)
for (sim_dat_indx in 1:n_sim_data) {
  r_accept_temp_sim <- readRDS(paste0(output_sub_r_accept_dir,"/",
                                  "r_accept%",sim_dat_indx,".RDS"))
  r_accept_temp_sim <- r_accept_temp_sim[,"average"] # average of chain(s)
  r_accept_all_sim[sim_dat_indx,] <- r_accept_temp_sim
  
  R_PD_temp_sim <- readRDS(paste0(output_sub_R_PD_dir,"/",
                              "%times_R_PD",sim_dat_indx,".RDS"))
  R_PD_temp_sim <- R_PD_temp_sim[,"average"]
  R_PD_all_sim[sim_dat_indx,] <- R_PD_temp_sim
}
r_accept_avg_sim <- colMeans(r_accept_all_sim)
R_PD_avg_sim <- colMeans(R_PD_all_sim)
setwd(output_sub_dir)
write.csv(r_accept_avg_sim,"r_accept%_avg.csv")
write.csv(R_PD_avg_sim, "%times_R_PD_avg.csv")
setwd(code_dir)