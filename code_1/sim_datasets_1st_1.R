set.seed(1)
library(mvtnorm) # 'rmvnorm', 'mvtnorm::rmvt' (more official)
library(sn) # 'sn::rmsn' (more official)
library(MomTrunc) # 'MomTrunc::rmvSN' (should be same as 'sn::rmsn')

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

code_dir <- thisFile()
setwd(code_dir)

parent_source <- dirname(code_dir)
base_name <- basename(code_dir)
sub_base_length <- nchar(base_name)-nchar("code")
main_base_name <- substr(base_name,
                         1, nchar(base_name)-sub_base_length)
sub_base_name <- substr(base_name,
                        nchar(base_name)-sub_base_length+1, nchar(base_name))
output_dir <- paste0(parent_source,"/","outputs",sub_base_name)
y_sim_dir <- paste0(parent_source,"/","y_sim",sub_base_name)
dir.create(y_sim_dir)

source(paste0("functions",sub_base_name,".R"))



### MODIFY AS NEEDED
actual_data <- read.csv(paste0(parent_source,"/2023_Lockdown_Wide_20230201.csv"))
muscle_names <- c("MRS_FF_SOL","MRS_FF_VL","MRS_FF_BB","MRS_FF_DEL") # "X6MWT"
L <- length(muscle_names)
Q <- choose(L,2)+L+1
if (L==4) {
  n_amb_stages <- 3
  amb_names <- c("Early", "Late", "Non")
} else if (L==5) {
  n_amb_stages <- 2
  amb_names <- c("Early", "Late")
}

### MODIFY AS NEEDED
y_actuals <- prepare_dat_groups(actual_data, muscle_names, "annualized_delta")

### MODIFY AS NEEDED
J_max_unique <- 4
N_sim_unique <- 100
post_center <- "post_med"
type_of_data <- "bal1" # bal1, unbal1, unbal2
type_of_model <- 'N'  # "N", "t10", "t3", "SN"

test_bool <- FALSE
n_sim_rep <- 1 # 0, 1, 5
n_sim_data <- 500

### MODIFY AS NEEDED
if (L==4) {
  missing_perc_by_muscle_vec_unique <- c(0.05, 0.05, 0.75, 0.75)
  # last element should be 0 (don't want ALL MRS_FF's = NA's in a row)
  missing_perc_by_rowNAs_vec_unique <- c(0.1, 0.6, 0.1, 0)
} else if (L==5) {
  missing_perc_by_muscle_vec_unique <- c(0.05, 0.05, 0.75, 0.75, 0.2)
  # last 2 elems should be 0 (don't want ALL MRS_FF's = NA's in a row)
  missing_perc_by_rowNAs_vec_unique <- c(0.1, 0.65, 0.15, 0, 0)
}
missing_perc_by_muscle_mat_unique <- 
  matrix(rep(missing_perc_by_muscle_vec_unique,J_max_unique), ncol=L, byrow=TRUE)
missing_perc_by_rowNAs_mat_unique <- 
  matrix(rep(missing_perc_by_rowNAs_vec_unique,J_max_unique), ncol=L, byrow=TRUE)
### MODIFY AS NEEDED
# these should sum to 1
mat_Js_dist_unique <- matrix(c(1:J_max_unique, rep(1/J_max_unique,J_max_unique)), ncol=2) # note: byrow=FALSE here
colnames(mat_Js_dist_unique) <- c("Js","%") # note: no 'Freq' column here
perc_Js_dist_unique <- mat_Js_dist_unique[,"%"]



post_center_short <- substr(post_center,
                            nchar("post_")+1, nchar(post_center))
sim_data_info <- paste0(L,"L_",J_max_unique,"J_",N_sim_unique,"N_",post_center_short)

y_sim_dir_target <- paste0(y_sim_dir,"/","target")
dir.create(y_sim_dir_target)
y_sim_dir_target_L <- paste0(y_sim_dir_target,"/",L,"L","_",post_center)
dir.create(y_sim_dir_target_L)

y_sim_dir_test <- paste0(y_sim_dir,"/","test")
dir.create(y_sim_dir_test)
y_sim_dir_test_L <- paste0(y_sim_dir_test,"/",sim_data_info)
dir.create(y_sim_dir_test_L)
y_sim_dir_test_L_child <- paste0(y_sim_dir_test_L,"/",type_of_data,"_",type_of_model)
dir.create(y_sim_dir_test_L_child)

y_sim_dir_data <- paste0(y_sim_dir,"/","data")
dir.create(y_sim_dir_data)
y_sim_dir_data_L <- paste0(y_sim_dir_data,"/",sim_data_info)
dir.create(y_sim_dir_data_L)
y_sim_dir_data_L_full <- paste0(y_sim_dir_data_L,"/","full","_",type_of_model)
dir.create(y_sim_dir_data_L_full)
y_sim_dir_data_L_full_amb_csv_vec <- sapply(1:n_amb_stages, function(amb_indx)
  paste0(y_sim_dir_data_L_full,"/",amb_names[amb_indx],"_csv","_",n_sim_data,"nsim"))
sapply(1:n_amb_stages, function(amb_indx) dir.create(y_sim_dir_data_L_full_amb_csv_vec[amb_indx]))
y_sim_dir_data_L_full_amb_vec <- sapply(1:n_amb_stages, function(amb_indx)
  paste0(y_sim_dir_data_L_full,"/",amb_names[amb_indx],"_",n_sim_data,"nsim"))
sapply(1:n_amb_stages, function(amb_indx) dir.create(y_sim_dir_data_L_full_amb_vec[amb_indx]))

y_sim_dir_data_L_child <- paste0(y_sim_dir_data_L,"/",type_of_data,"_",type_of_model)



y_mat_Lcols <- lapply(1:n_amb_stages, function(amb_indx)
  matrix(unlist(y_actuals[[amb_indx]]), ncol=L, byrow=TRUE))
Ns <- sapply(1:n_amb_stages, function(amb_indx) length(y_actuals[[amb_indx]])) # subjects
N.Js <- sapply(1:n_amb_stages, function(amb_indx) nrow(y_mat_Lcols[[amb_indx]])) # visits

ps_list <- lapply(1:n_amb_stages, function(amb_indx)
  sapply(1:Ns[amb_indx], function(i) length(y_actuals[[amb_indx]][[i]])))
Js_list <- lapply(1:n_amb_stages, function(amb_indx)
  ps_list[[amb_indx]]/L) # 8max for ff1, 6max for ff2, 7max for ff3; use 'max(Js)'
j_order0s_list <- lapply(1:n_amb_stages, function(amb_indx)
  lapply(1:Ns[amb_indx], function(i) rep(1:Js_list[[amb_indx]][i], each=L)))
l_order0s_list <- lapply(1:n_amb_stages, function(amb_indx)
  lapply(1:Ns[amb_indx], function(i) rep(1:L,Js_list[[amb_indx]][i])))
p_max_vec <- sapply(1:n_amb_stages, function(amb_indx)
  ps_list[[amb_indx]][which(Js_list[[amb_indx]]==max(Js_list[[amb_indx]]))[1]])
J_max_vec <- p_max_vec/L
j_order0_max_list <- lapply(1:n_amb_stages, function(amb_indx)
  j_order0s_list[[amb_indx]][[which(Js_list[[amb_indx]]==max(Js_list[[amb_indx]]))[1]]])
l_order0_max_list <- lapply(1:n_amb_stages, function(amb_indx)
  l_order0s_list[[amb_indx]][[which(Js_list[[amb_indx]]==max(Js_list[[amb_indx]]))[1]]])

### names for some functions
n_etas <- choose(L,2)
muscle_combns <- combn(1:L,2)
eta_names <- sapply(1:n_etas, function(combn_indx)
  paste0("eta_{",muscle_combns[1,combn_indx],muscle_combns[2,combn_indx],"}"))
rho_names <- paste0("rho_{(",1:L,")}")
r_names <- c(eta_names, rho_names, "gamma")



target_short_mu_list <- vector("list",length=n_amb_stages)
target_short_s_list <- vector("list",length=n_amb_stages)
target_r_list <- vector("list",length=n_amb_stages)
target_w_list <- vector("list",length=n_amb_stages)
target_SRM_list <- vector("list",length=n_amb_stages)
target_mu_list <- vector("list",length=n_amb_stages)
target_Sigma_list <- vector("list",length=n_amb_stages)
for (amb_indx in 1:n_amb_stages) {
  amb_name <- amb_names[amb_indx]
  target_params <- readRDS(paste0(y_sim_dir_target,"/",
                                  amb_name,"_",L,"L","_",post_center,".RDS"))
  target_short_mu <- target_params$target_short_mu
  target_short_s <- target_params$target_short_s
  target_r <- target_params$target_r
  target_w <- target_params$target_w
  target_SRM <- target_params$target_SRM
  
  target_mu <- target_short_mu[l_order0_max_list[[amb_indx]]]
  target_s <- target_short_s[l_order0_max_list[[amb_indx]]]
  target_S <- diag(target_s)
  target_R <- rtoR_old(target_r, r_names, p_max_vec[amb_indx],
                       j_order0_max_list[[amb_indx]], l_order0_max_list[[amb_indx]])
  target_Sigma <- target_S %*% target_R %*% target_S
  
  target_short_mu_list[[amb_indx]] <- target_short_mu
  target_short_s_list[[amb_indx]] <- target_short_s
  target_r_list[[amb_indx]] <- target_r
  target_w_list[[amb_indx]] <- target_w
  target_SRM_list[[amb_indx]] <- target_SRM
  
  target_mu_list[[amb_indx]] <- target_mu
  target_Sigma_list[[amb_indx]] <- target_Sigma
}










if (test_bool==FALSE) {
  for (amb_indx in 1:n_amb_stages) {
    
    target_mu_amb <- target_mu_list[[amb_indx]]
    target_Sigma_amb <- target_Sigma_list[[amb_indx]]
    p_max_unique <- J_max_unique*L
    
    target_mu_temp <- target_mu_amb[1:p_max_unique]
    target_Sigma_temp <- target_Sigma_amb[1:p_max_unique,1:p_max_unique]
    
    for (sim_dat_indx in 1:n_sim_data) {

      if (type_of_model=="N") {
        y_sim_full_mat <- 
          mvtnorm::rmvnorm(N_sim_unique, target_mu_temp, target_Sigma_temp)
      } else if ((type_of_model=="t10") || (type_of_model=="t3")) {
        if (type_of_model=="t10") {
          deg_f <- 10
        } else if (type_of_model=="t3") {
          deg_f <- 3
        }
        y_sim_full_mat <- 
          mvtnorm::rmvt(n=N_sim_unique, sigma=target_Sigma_temp*(deg_f-2)/(deg_f),
                        df=deg_f, delta=target_mu_temp)
      } else if (type_of_model=="SN") {
        lambda <- rep(0.1, p_max_unique)
        y_sim_full_mat <- 
          MomTrunc::rmvSN(n=N_sim_unique, mu=target_mu_temp,
                          Sigma=target_Sigma_temp, lambda=lambda)
      }
      
      y_sim_full_list <- as.list(data.frame(t(y_sim_full_mat))) # convert to list form
      y_sim_full_list <- unname(y_sim_full_list) # unname 'X1', 'X2', 'X3', ...
      
      setwd(y_sim_dir_data_L_full_amb_csv_vec[amb_indx])
      if (!file.exists(paste0("y_sim_full",sim_dat_indx,".csv"))) {
        write.csv(y_sim_full_mat, paste0("y_sim_full",sim_dat_indx,".csv"))
      }
      setwd(code_dir)
      
      setwd(y_sim_dir_data_L_full_amb_vec[amb_indx])
      if (!file.exists(paste0("y_sim",sim_dat_indx,".RDS"))) {
        saveRDS(y_sim_full_list, file=paste0("y_sim",sim_dat_indx,".RDS"))
      }
      setwd(code_dir)
      
    }
    
  }
  
  
  
  
  
  # only target parameters are determined by ambulatory data
  for (sim_rep_indx in 1:n_sim_rep) {
    y_sim_dir_data_L_child_rep <- paste0(y_sim_dir_data_L_child,"_",sim_rep_indx) 
    dir.create(y_sim_dir_data_L_child_rep)
    for (amb_indx in 1:n_amb_stages) {
      y_sim_dir_data_L_child_rep_amb <- paste0(y_sim_dir_data_L_child_rep,"/",amb_names[amb_indx],"_",n_sim_data,"nsim")
      dir.create(y_sim_dir_data_L_child_rep_amb)
      for (sim_dat_indx in 1:n_sim_data) {
        # 'y_sim_full' same as 'y_sim_full_mat'
        y_sim_full <- read.csv(paste0(y_sim_dir_data_L_full_amb_csv_vec[amb_indx],"/","y_sim_full",sim_dat_indx,".csv"))
        y_sim_full <- as.matrix(y_sim_full)[,-1]
        y_sim_full <- unname(y_sim_full)
        
        if (type_of_data=="bal1") {
          more_sim_results_temp <-
            balanced1_sim_func(y_sim_full,
                               test_bool, N_sim_unique, L,
                               missing_perc_by_muscle_vec_unique, # vector
                               missing_perc_by_rowNAs_vec_unique, # vector
                               J_max_unique
                               )
        } else if (type_of_data=="unbal1") {
          more_sim_results_temp <-
            unbalanced1_sim_func(y_sim_full,
                                 test_bool, N_sim_unique, L,
                                 missing_perc_by_muscle_vec_unique, # vector
                                 missing_perc_by_rowNAs_vec_unique, # vector
                                 mat_Js_dist_unique, perc_Js_dist_unique
                                 )
        } else if (type_of_data=="unbal2") {
          more_sim_results_temp <-
            unbalanced2_sim_func(y_sim_full,
                                 test_bool, N_sim_unique, L,
                                 missing_perc_by_muscle_mat_unique, # matrix
                                 missing_perc_by_rowNAs_mat_unique, # matrix
                                 mat_Js_dist_unique, perc_Js_dist_unique
                                 )
        }
        
        y_sim_wNA_list_temp <- more_sim_results_temp$y_sim_wNA_list
        
        setwd(y_sim_dir_data_L_child_rep_amb)
        if (!file.exists(paste0("y_sim",sim_dat_indx,".RDS"))) {
          saveRDS(y_sim_wNA_list_temp, file=paste0("y_sim",sim_dat_indx,".RDS"))
        }
        setwd(code_dir)
        
      }
      
    }
    
  }
  
  
  
  
  
} else if (test_bool==TRUE) {
  # only target parameters are determined by ambulatory data
  MSE_names <- c("MSE_sim_full_mu","MSE_sim_wNA_mu",
                 "MSE_sim_full_s","MSE_sim_wNA_s",
                 "MSE_sim_full_etas","MSE_sim_wNA_etas",
                 "MSE_sim_full_w","MSE_sim_wNA_w",
                 "MSE_sim_full_SRM","MSE_sim_wNA_SRM",
                 "MSE_avg")
  more_MSE_mat <- matrix(nrow=n_amb_stages,ncol=length(MSE_names),
                         dimnames=list(amb_names,MSE_names))
  more_sim_results <- vector("list",length=n_amb_stages)
  for (amb_indx in 1:n_amb_stages) {
    y_sim_dir_test_L_child_amb <- paste0(y_sim_dir_test_L_child,"/",amb_names[amb_indx])
    dir.create(y_sim_dir_test_L_child_amb)
    
    # 'y_sim_full' same as 'y_sim_full_mat'
    y_sim_full <- read.csv(paste0(y_sim_dir_data_L_full_amb_csv_vec[amb_indx],"/","y_sim_full",1,".csv"))
    y_sim_full <- as.matrix(y_sim_full)[,-1]
    y_sim_full <- unname(y_sim_full)
    
    if (type_of_data=="bal1") {
      more_sim_results[[amb_indx]] <-
        balanced1_sim_func(y_sim_full,
                           test_bool, N_sim_unique, L,
                           missing_perc_by_muscle_vec_unique, # vector
                           missing_perc_by_rowNAs_vec_unique, # vector
                           J_max_unique,
                           amb_indx,
                           target_short_mu_list, target_short_s_list, target_r_list,
                           target_w_list, target_SRM_list
                           )
    } else if (type_of_data=="unbal1") {
      more_sim_results[[amb_indx]] <-
        unbalanced1_sim_func(y_sim_full,
                             test_bool, N_sim_unique, L,
                             missing_perc_by_muscle_vec_unique, # vector
                             missing_perc_by_rowNAs_vec_unique, # vector
                             mat_Js_dist_unique, perc_Js_dist_unique,
                             amb_indx,
                             target_short_mu_list, target_short_s_list, target_r_list,
                             target_w_list, target_SRM_list
                             )
    } else if (type_of_data=="unbal2") {
      more_sim_results[[amb_indx]] <-
        unbalanced2_sim_func(y_sim_full,
                             test_bool, N_sim_unique, L,
                             missing_perc_by_muscle_mat_unique, # matrix
                             missing_perc_by_rowNAs_mat_unique, # matrix
                             mat_Js_dist_unique, perc_Js_dist_unique,
                             amb_indx,
                             target_short_mu_list, target_short_s_list, target_r_list,
                             target_w_list, target_SRM_list)
    }
    
    ### Evaluate simulated data with missingness
    # Note that we used posterior median of mu, s, r (compare this - overall)
    # and used our own set of column/row missingness and J_dist
    more_sim_results[[amb_indx]]
    
    the_MSEs <- more_sim_results[[amb_indx]]$the_MSEs
    more_MSE_mat[amb_indx,] <- the_MSEs
    
    ps_sim <- more_sim_results[[amb_indx]]$ps_sim
    dat_sim_full <- more_sim_results[[amb_indx]]$dat_sim_full
    dat_sim_wNA <- more_sim_results[[amb_indx]]$dat_sim_wNA
    setwd(y_sim_dir_test_L_child_amb)
    if (!file.exists("ps_sim.csv")) {
      write.csv(ps_sim,"ps_sim.csv")
    }
    if (!file.exists("dat_sim_full.csv")) {
      write.csv(dat_sim_full, "dat_sim_full.csv")
    }
    if (!file.exists("dat_sim_miss.csv")) {
      write.csv(dat_sim_wNA, "dat_sim_miss.csv")
    }
    setwd(code_dir)
    
  }
  
  setwd(y_sim_dir_test_L_child)
  if (!file.exists("more_MSE.csv")) {
    write.csv(more_MSE_mat, "more_MSE.csv")
  }
  
  
}
