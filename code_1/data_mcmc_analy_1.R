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

library(devtools) # didn't restart R (no), didn't compile from source code (no)
library(rgl)

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
ff_dat_names <- c("Early Ambulatory Data", "Late Ambulatory Data", "Non-Ambulatory Data")

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

source(paste0("functions",sub_base_name,".R"))

### MODIFY AS NEEDED
if (detectCores()==12) {
  ncores <- 12-1
} else {
  ncores <- 32 
}
amb_indx <- 1 # 1, 2, 3 (L=4); # 1, 2 (L=5)
muscle_names_simple <- c("SOL","VL","BB","DEL") # "6MWT"
muscle_names_simple_subset <- muscle_names_simple # muscle_names_simple | c("SOL","VL")
muscle_names <- c("MRS_FF_SOL","MRS_FF_VL","MRS_FF_BB","MRS_FF_DEL") # "X6MWT"

actual_data <- read.csv(paste0(parent_source,"/2023_Lockdown_Wide_20230201.csv"))
L <- length(muscle_names_simple)
measure_name <- "annualized_delta"
y_actuals <- prepare_dat_groups(actual_data, muscle_names, measure_name)

y <- y_actuals[[amb_indx]]
nchains <- 4
iters <- 2000
burn_in_prop <- 1000/iters
r_jump_dist_name <- "sbeta" # sbeta, sbeta_one, sbeta11, unif, unif_one, unif11
optim_method <- "L-BFGS-B"
grad_bool <- TRUE
upper_bool <- TRUE



amb_name <- amb_names[amb_indx]
ff_dat_name <- ff_dat_names[amb_indx]
N <- length(y) # 140 for ff1, 74 for ff2, 51 for ff3
Q <- choose(L,2)+L+1
ps <- sapply(1:N, function(i) length(y[[i]]))
Js <- ps/L
Js_sum <- sum(Js)
J_max <- max(Js) # 8max for ff1, 6max for ff2, 7max for ff3
j_order0s <- lapply(1:N, function(i) rep(1:Js[i], each=L))
l_order0s <- lapply(1:N, function(i) rep(1:L, Js[i]))
p_max <- ps[which(Js==J_max)[1]]
j_order0_max <- j_order0s[[which(Js==J_max)[1]]]
l_order0_max <- l_order0s[[which(Js==J_max)[1]]]

target_short_mu <- rep("",L)
target_short_s <- rep("",L)
target_r <- rep("",Q)
target_w <- rep("",L)
target_SRM <- ""



# HAVE DIFFUSE/UNINFORMATIVE PRIOR > 2 FOR JUMPING DIST
# (TO HAVE DEFINED MEAN AND VARIANCE)
# KEEP IT SMALL SO LIKELIHOOD DOMINATES PRIOR FOR JUMPING DIST (POSTERIOR DIST)
nu <- rep(2.1,L)



### MODIFY AS NEEDED
sbeta_strings <- c("sbeta", "sbeta_one", "sbeta11")
unif_strings <- c("unif", "unif_one", "unif11")
if ((r_jump_dist_name %in% sbeta_strings) && (N==140) && (L==4)) { # ff1, J_max=8
  kappa <- c(26.5, 17.2, 19.4,   16.9, 18.5,   23.0,   33.3, 27.6, 25.9, 30.5,   320)
} else if ((r_jump_dist_name %in% sbeta_strings) && (N==74) && (L==4)) { # ff2, J_max=6
  kappa <- c(7.2, 4.7, 5.07,   4.4, 4.81,   8.5,   10.02, 52.00, 6.16, 5.17,   209)
} else if ((r_jump_dist_name %in% sbeta_strings) && (N==51) && (L==4)) { # ff3, J_max=7
  kappa <- c(8.53, 5.0, 7.82,   4.72, 4.52,   7.3,   12.09, 28.05, 6.23, 5.02,   114)
} else if ((r_jump_dist_name %in% sbeta_strings) && (N==140) && (L==5)) { # ff1, J_max=8
  kappa <- c(26.7, 18.2, 20.1, 18.2,   17.0, 18.2, 18.8,   21.8, 25.65,   19.9,   26.2, 26.3, 25.0, 29.5, 46.6,   480)
} else if ((r_jump_dist_name %in% sbeta_strings) && (N==74) && (L==5)) { # ff2, J_max=6
  kappa <- c(6.37, 4.72, 4.95, 5.8,   4.32, 4.83, 4.59,   8.7, 5.6,   5.2,   17.1, 56.3, 6.2, 5.1, 12.0,   410)

} else if (r_jump_dist_name %in% unif_strings) {
  kappa <- rep(NA, Q)
}



### MODIFY AS NEEDED
if (iters < 50000) {
  r_updates <- rep(1,Q)
} else {
  if ((r_jump_dist_name %in% sbeta_strings) && (N==140) && (L==4)) { # ff1, J_max=8
    r_updates <- c(3,6,5,6,4,5,3,3,11,4,3)
  } else if ((r_jump_dist_name %in% sbeta_strings) && (N==74) && (L==4)) { # ff2, J_max=6
    r_updates <- c(3,7,5,6,4,7,4,3,4,4,3)
  } else if ((r_jump_dist_name %in% sbeta_strings) && (N==51) && (L==4)) { # ff3, J_max=7
    r_updates <- c(3,4,3,4,4,4,3,3,3,3,3)
  } else if ((r_jump_dist_name %in% sbeta_strings) && (N==140) && (L==5)) { # ff1, J_max=8
    r_updates <- c(2,4,5,2,5,5,2,5,20,15,2,2,8,6,3,3)
  } else if ((r_jump_dist_name %in% sbeta_strings) && (N==74) && (L==5)) { # ff2, J_max=6
    r_updates <- c(3,7,7,4,4,4,4,6,6,6,4,4,4,5,4,4)
    
    
    
  } else {
    r_updates <- rep(3,Q)
  }
}



### for sweep_once
Hs <- lapply(1:N, function(i)
  matrix(nrow=ps[i]+1,ncol=ps[i]+1))
one_to_p1s <- lapply(1:N, function(i) 1:(ps[i]+1))

### for mu
Xs <- list()
X.Ts <- list()
for (i in 1:N) {
  X <- matrix(nrow=ps[i], ncol=L)
  for (m in 1:ps[i]) {
    X[m,] <- 1*(l_order0s[[i]][m]==c(1:L))
  }
  X.T <- t(X)

  Xs[[i]] <- X
  X.Ts[[i]] <- X.T
}

### for s
uniqueS_locs <- uniqueS_locs_func(l_order0_max)

### names for some functions
n_etas <- choose(L,2)
muscle_combns <- combn(1:L,2)
eta_names <- sapply(1:n_etas, function(combn_indx)
  paste0("eta_{",muscle_combns[1,combn_indx],muscle_combns[2,combn_indx],"}"))
rho_names <- paste0("rho_{(",1:L,")}")
r_names <- c(eta_names, rho_names, "gamma")

### names for plots
mu_names_detailed <- paste0("mu",1:L,"[",muscle_names_simple,"]")
s_names_detailed <- paste0("s",1:L,"[", muscle_names_simple,"]")
muscle_names_simple_combn <- combn(muscle_names_simple,2)
eta_names_detailed <- sapply(1:n_etas, function(combn_indx)
  paste0("eta",muscle_combns[1,combn_indx],muscle_combns[2,combn_indx],
         "[",muscle_names_simple_combn[1,combn_indx],",",muscle_names_simple_combn[2,combn_indx],"]"))
rho_names_detailed <- paste0("rho(",1:L,")","[",muscle_names_simple,"]")
r_names_detailed <- c(eta_names_detailed, rho_names_detailed, "gamma")

### Locations of unique elements
R_latex_upper <- R_order_latex_upper(p_max,j_order0_max,l_order0_max)
R_latex_all <- R_latex_upper
R_latex_all[lower.tri(R_latex_all)] <- t(R_latex_all)[lower.tri(R_latex_all)]
uniqueR_locs_all <- uniqueR_locs_func(R_latex_all, r_names)

### For Barnard stuff later
pivot_implied_locs <- pivot_and_implied_locs(R_latex_upper, r_names)
pivot_locs <- pivot_implied_locs[[1]]
goal_indxs <- pivot_locs[-Q,]

subblock_indxs_list <- vector("list",length=Q-1)
for (k in 1:(Q-1)) {
  goal_indx1 <- goal_indxs[k,1]
  goal_indx2 <- goal_indxs[k,2]
  if (k %in% 1:n_etas) {
    subblock_indxs_list[[k]] <- 
      eta_biggest_submatrix_indxs(goal_indx1, goal_indx2, J_max, L,
                                  j_order0_max, l_order0_max)
  } else {
    subblock_indxs_list[[k]] <- 
      rho_biggest_submatrix_indxs(goal_indx1, goal_indx2, p_max,
                                  j_order0_max, l_order0_max)
  }
}
n_eta_subs <- nrow(subblock_indxs_list[[1]])

eta_names_twice <- rep(eta_names,each=2)
rho_names_twice <- paste0("rho_{(",as.vector(muscle_combns),")}")
n_etas2 <- 2*n_etas
gamma_sub_ncol <- 3
gamma_subs_latex <- gamma_subs_latex_func(n_etas2, gamma_sub_ncol,
                                          eta_names_twice, rho_names_twice)
unique_subb_locs <- unique_subb_locs_func(r_names, gamma_subs_latex)
empty_gamma_subs <- NULL
for (m_indx in 1:n_etas2) {
  empty_gamma_sub_temp <- diag(gamma_sub_ncol)
  empty_gamma_subs <- cbind(empty_gamma_subs, empty_gamma_sub_temp)
}

### r updates stuff
itr <- rep(1,Q)
itr_v2 <- rep(1,Q)

itr_names <- vector("list", length=Q)
for (k in 1:Q) {
  if (r_updates[k]==1) {
    itr_names[[k]] <- paste0("iter",1:iters)
  } else {
    itr_names[[k]] <- numeric(length=1+(iters-1)*r_updates[k])
    itr_names[[k]][1] <- c("iter1")
    for (it in 2:iters) {
      for (ru in 1:r_updates[k]) {
        itr_v2[k] <- itr_v2[k] + 1
        itr_names[[k]][itr_v2[k]] <- paste0("iter",it,"_upd",ru)
      }
    }
  }
}



y_etc <- list(
  y=y, N=N, iters=iters, burn_in_prop=burn_in_prop,
  r_jump_dist_name=r_jump_dist_name,
  r_updates=r_updates, itr=itr, itr_names=itr_names,
  Js=Js, Js_sum=Js_sum, L=L, Q=Q, ps=ps,
  J_max=J_max, p_max=p_max, j_order0_max=j_order0_max, l_order0_max=l_order0_max,
  target_short_mu=target_short_mu, target_short_s=target_short_s, target_r=target_r,
  target_w=target_w, target_SRM=target_SRM,
  nu=nu, kappa=kappa,
  Hs=Hs, one_to_p1s=one_to_p1s,
  Xs=Xs, X.Ts=X.Ts,
  uniqueS_locs=uniqueS_locs,
  mu_names_detailed=mu_names_detailed, s_names_detailed=s_names_detailed,
  r_names_detailed=r_names_detailed,
  uniqueR_locs_all=uniqueR_locs_all,
  subblock_indxs_list=subblock_indxs_list, n_etas=n_etas, n_eta_subs=n_eta_subs,
  n_etas2=n_etas2, gamma_sub_ncol=gamma_sub_ncol,
  unique_subb_locs=unique_subb_locs, empty_gamma_subs=empty_gamma_subs)





mvn_disc <- function(y_etc) {
  y <- y_etc$y
  N <- y_etc$N
  iters <- y_etc$iters
  burn_in_prop <- y_etc$burn_in_prop
  
  r_jump_dist_name <- y_etc$r_jump_dist_name
  
  r_updates <- y_etc$r_updates
  itr <- y_etc$itr
  itr_names <- y_etc$itr_names
  
  Js <- y_etc$Js # 8max for ff1, 6max for ff2, 7max for ff3
  Js_sum <- y_etc$Js_sum
  L <- y_etc$L
  Q <- y_etc$Q
  ps <- y_etc$ps
  J_max <- y_etc$J_max
  p_max <- y_etc$p_max
  j_order0_max <- y_etc$j_order0_max
  l_order0_max <- y_etc$l_order0_max

  target_short_mu <- y_etc$target_short_mu
  target_short_s <- y_etc$target_short_s
  target_r <- y_etc$target_r
  target_w <- y_etc$target_w
  target_SRM <- y_etc$target_SRM
  
  nu <- y_etc$nu
  kappa <- y_etc$kappa
  
  Hs <- y_etc$Hs
  one_to_p1s <- y_etc$one_to_p1s
  
  Xs <- y_etc$Xs
  X.Ts <- y_etc$X.Ts

  uniqueS_locs <- y_etc$uniqueS_locs
  
  mu_names_detailed <- y_etc$mu_names_detailed
  s_names_detailed <- y_etc$s_names_detailed
  r_names_detailed <- y_etc$r_names_detailed

  uniqueR_locs_all <- y_etc$uniqueR_locs_all
  
  ### For Barnard stuff later
  subblock_indxs_list <- y_etc$subblock_indxs_list
  n_etas <- y_etc$n_etas
  n_eta_subs <- y_etc$n_eta_subs

  ### For gamma + Barnard
  n_etas2 <- y_etc$n_etas2
  gamma_sub_ncol <- y_etc$gamma_sub_ncol
  unique_subb_locs <- y_etc$unique_subb_locs
  empty_gamma_subs <- y_etc$empty_gamma_subs
  
  
  
  
  
  success_count_s <- matrix(c(rep(NA,L),rep(0,(iters-1)*L)),nrow=iters,ncol=L,byrow=TRUE)
  
  success_count_r <- lapply(1:Q, function(k) c(NA,rep(0,(iters-1)*r_updates[k])))
  log_prior_stars_r <- lapply(1:Q, function(k) c(NA,rep(NA,(iters-1)*r_updates[k])))
  
  
  
  
  
  y_mat_Lcol <- matrix(unlist(y), ncol=L, byrow=TRUE)
  sample_var <- apply(y_mat_Lcol, 2, var, na.rm=TRUE)
  mu_0 <- colMeans(y_mat_Lcol, na.rm=TRUE)
  column_mins_maxs <- apply(y_mat_Lcol, 2, range, na.rm=TRUE)
  column_ranges <- column_mins_maxs[2,]-column_mins_maxs[1,]
  column_ranges_transformed <- (0.25*column_ranges)^2
  Sigma_0 <- diag(column_ranges_transformed)
  Sigma_0.inv <- diag(1/column_ranges_transformed)
  mu_0.Sigma_0.inv <- Sigma_0.inv%*%mu_0
  
  if (!any(is.na(y_mat_Lcol))) {
    y_imputed <- y
    matrix_imputed <- matrix(unlist(y_imputed), ncol=L, byrow=TRUE)
  } else {
    y_wNA <- y
    y_imputed <- y_wNA
    mis_vars <- list()
    obs_vars <- list()
    obs_mis <- list()
    ps_tilde <- c()
    for (i in 1:N) {
      mis_bool <- is.na(y_wNA[[i]])
      mis_vars[[i]] <- which(mis_bool==T)
      obs_vars[[i]] <- which(mis_bool==F)
      obs_mis[[i]] <- c(obs_vars[[i]],mis_vars[[i]])
      ps_tilde[i] <- length(obs_vars[[i]])
    }
  }
  
  
  
  ### Create vectors, arrays, lists
  short_mu <- matrix(nrow=iters,ncol=L)
  short_mu[1,] <- colMeans(y_mat_Lcol, na.rm=TRUE)
  mu_special <- short_mu[1,l_order0_max]
  
  short_s <- matrix(nrow=iters,ncol=L)
  short_s[1,] <- apply(y_mat_Lcol, 2, sd, na.rm=TRUE)
  S_special <- diag(short_s[1,l_order0_max])
  
  r <- matrix(nrow=iters,ncol=Q)
  r[1,] <- rep(0,Q)
  R_special <- diag(p_max)
  
  r_all <- lapply(1:Q, function(k) c(0,rep(NA,(iters-1)*r_updates[k])))
  
  Sigma_special <- S_special%*%R_special%*%S_special
  Sigma_special.inv <- solve(Sigma_special)
  
  
  
  for (it in 2:iters) {
    if (any(is.na(y_mat_Lcol))) {
      for (i in 1:N) {
        p_i <- ps[i]
        p_tilde_i <- ps_tilde[i]
        
        bottom_aug_cov <- cbind(mu_special[1:p_i], Sigma_special[1:p_i,1:p_i])
        top_row <- cbind(c(-1), t(mu_special[1:p_i]))
        aug_cov <- rbind(top_row, bottom_aug_cov)
        if (p_tilde_i < p_i) {
          G <- aug_cov[c(1,obs_mis[[i]]+1), c(1,obs_mis[[i]]+1)]
          
          for (m in 1:p_tilde_i) {
            G <- sweep_once(G, Hs[[i]], one_to_p1s[[i]], m+1)
          }
          
          E_mat <- as.matrix(G[1:(p_tilde_i+1),(p_tilde_i+2):(p_i+1)])
          DA_mean <- numeric(p_i-p_tilde_i)
          for (n in 1:(p_i-p_tilde_i)) {
            DA_mean[n] <- E_mat[1,n]+t(E_mat[2:(p_tilde_i+1),n])%*%y_wNA[[i]][obs_vars[[i]]]
          }
          
          resid_cov <- G[(p_tilde_i+2):(p_i+1),(p_tilde_i+2):(p_i+1)]
          
          y_imputed[[i]] <- y_wNA[[i]]
          if (p_tilde_i != p_i-1) {
            y_imputed[[i]][mis_vars[[i]]] <- rmvnorm(1,DA_mean,resid_cov)
          } else {
            y_imputed[[i]][mis_vars[[i]]] <- rnorm(1,DA_mean,resid_cov)
          }
        }
      }
      matrix_imputed <- matrix(unlist(y_imputed), ncol=L, byrow=TRUE)
    }
    
    first_sum <- matrix(0, nrow=L, ncol=L)
    second_sum <- matrix(0, nrow=L, ncol=1)
    for (i in 1:N) {
      p_i <- ps[i]
      first_sum <- first_sum + X.Ts[[i]]%*%Sigma_special.inv[1:p_i,1:p_i]%*%Xs[[i]]
      second_sum <- second_sum + X.Ts[[i]]%*%Sigma_special.inv[1:p_i,1:p_i]%*%y_imputed[[i]]
    }
    first_sum <- solve(first_sum + Sigma_0.inv)
    second_sum <- second_sum + mu_0.Sigma_0.inv
    
    short_mu[it,] <- rmvnorm(1,first_sum%*%second_sum,first_sum)
    mu_special <- short_mu[it,l_order0_max]

    ### MH for s
    for (l in 1:L) {
      # old
      s_l_old <- short_s[it-1,l]
      S_old <- S_special
      S_old[uniqueS_locs[[l]]] <- s_l_old # faster than 'diag' method
      Sigma_old <-  S_old%*%R_special%*%S_old
      
      # new
      y_l <- matrix_imputed[,l]
      the_sum <- sum((y_l-short_mu[it,l])^2)
      the_shape <- nu[l]+0.5*Js_sum
      the_rate <- (nu[l]+1)*sample_var[l] + 0.5*the_sum
      s2_l_jump <- rinvgamma(n=1, shape=the_shape, rate=the_rate)
      s_l_jump <- sqrt(s2_l_jump)
      S_star <- S_special
      S_star[uniqueS_locs[[l]]] <- s_l_jump
      Sigma_star <- S_star%*%R_special%*%S_star
      
      log_lik_old <- 0
      log_lik_star <- 0
      for (i in 1:N) {
        p_i <- ps[i]
        log_lik_old <- log_lik_old + dmvnorm(y_imputed[[i]],
                                             mu_special[1:p_i],
                                             Sigma_old[1:p_i,1:p_i],
                                             log=TRUE)
        log_lik_star <- log_lik_star + dmvnorm(y_imputed[[i]],
                                               mu_special[1:p_i],
                                               Sigma_star[1:p_i,1:p_i],
                                               log=TRUE)
      }

      log_prior_old <- dinvgamma(s_l_old^2,shape=3,rate=4*sample_var[l],log=TRUE)
      log_prior_star <- dinvgamma(s2_l_jump,shape=3,rate=4*sample_var[l],log=TRUE)
      
      log_gc_old <- dinvgamma(s_l_old^2,shape=the_shape,rate=the_rate,log=TRUE)
      log_gc_star <- dinvgamma(s2_l_jump,shape=the_shape,rate=the_rate,log=TRUE)
      
      log_ratio_s <- (log_lik_star+log_prior_star-log_gc_star)-
        (log_lik_old+log_prior_old-log_gc_old)

      log_alpha_s <- min(log(1),log_ratio_s)
      u <- runif(1,min=0,max=1)
      if (log_alpha_s>log(u)) {
        short_s[it,l] <- s_l_jump
        success_count_s[it,l] <- success_count_s[it,l]+1 # max value should be 1 for matrix element
      } else if (log_alpha_s<=log(u)) {
        short_s[it,l] <- s_l_old
      }
      S_special[uniqueS_locs[[l]]] <- short_s[it,l]
    }

    ### MH for r
    for (k in 1:Q) {
      for (ru in 1:r_updates[k]) {
        itr[k] <- itr[k] + 1
        # old
        r_k_old <- r[it-1,k]
        R_old <- R_special
        R_old[uniqueR_locs_all[[k]]] <- r_k_old
        Sigma_old <- S_special%*%R_old%*%S_special
        
        #new
        if ( (r_jump_dist_name=="sbeta") || (r_jump_dist_name=="unif") ) {
          # intersecting PD intervals
          if (k %in% 1:n_etas) {
            indxs_temp_mat <- subblock_indxs_list[[k]]
            r_k_interval <- c(-1,1)
            for (eta_sub_indx in 1:n_eta_subs) {
              indxs_temp <- indxs_temp_mat[eta_sub_indx,]
              R_old_temp <- R_old[indxs_temp,indxs_temp]
              barnard_interval <- barnard(1,2,R_old_temp) # not same 1,2 as below
              r_k_interval <- c(max(r_k_interval[1],barnard_interval[1]),
                                min(r_k_interval[2],barnard_interval[2]))
            }
            
          } else if (k %in% (n_etas+1):(Q-1)) {
            indxs_temp <- subblock_indxs_list[[k]]
            R_old_temp <- R_old[indxs_temp,indxs_temp]
            barnard_interval <- barnard(1,2,R_old_temp) # not same 1,2 as below
            r_k_interval <- barnard_interval
            
          } else if (k==Q) {
            r_old_no_gamma <- r[it,1:(Q-1)]
            gamma_subs <- gamma_subs_func(Q, r_old_no_gamma, empty_gamma_subs,
                                          unique_subb_locs, gamma_sub_ncol)
            r_k_interval <- c(-1,1)
            for (indx in 1:n_etas2) {
              indx_m1 <- indx-1
              subb_indx <- (gamma_sub_ncol*indx_m1+1):(gamma_sub_ncol*indx_m1+gamma_sub_ncol)
              barnard_interval <- barnard(1,2,gamma_subs[,subb_indx]) # not same 1,2 as below
              r_k_interval <- c(max(r_k_interval[1],barnard_interval[1]),
                                min(r_k_interval[2],barnard_interval[2]))
            }
          }
          
        } else if ( (r_jump_dist_name=="sbeta_one") || (r_jump_dist_name=="unif_one") ) {
          # just use one largest possible submatrix
          if (k %in% 1:n_etas) {
            indxs_temp_mat <- subblock_indxs_list[[k]]
            eta_sub_indx <- sample(1:n_eta_subs, size=1)
            indxs_temp <- indxs_temp_mat[eta_sub_indx,]
            R_old_temp <- R_old[indxs_temp,indxs_temp]
            r_k_interval <- barnard(1,2,R_old_temp)
          } else if (k %in% (n_etas+1):(Q-1)) {
            indxs_temp <- subblock_indxs_list[[k]]
            R_old_temp <- R_old[indxs_temp,indxs_temp]
            r_k_interval <- barnard(1,2,R_old_temp)
          } else if (k==Q) {
            r_old_no_gamma <- r[it,1:(Q-1)]
            gamma_subs <- gamma_subs_func(Q, r_old_no_gamma, empty_gamma_subs,
                                          unique_subb_locs, gamma_sub_ncol)
            indx <- sample(1:n_etas2, size=1)
            indx_m1 <- indx-1
            subb_indx <- (gamma_sub_ncol*indx_m1+1):(gamma_sub_ncol*indx_m1+gamma_sub_ncol)
            r_k_interval <- barnard(1,2,gamma_subs[,subb_indx])
          }
          
        } else if ( (r_jump_dist_name=="sbeta11") || (r_jump_dist_name=="unif11") ) {
          r_k_interval <- c(-1,1)
        }
        
        
        
        lower_bound <- min(r_k_interval)
        upper_bound <- max(r_k_interval)
        
        if ( (r_jump_dist_name=="sbeta") || (r_jump_dist_name=="sbeta_one") || (r_jump_dist_name=="sbeta11") ) {
          alphaB_old <- ( (kappa[k]-1)*lower_bound+(2-kappa[k])*r_k_old
                          -upper_bound )/(lower_bound-upper_bound)
          betaB_old <- kappa[k]-alphaB_old
          
          r_k_jump <- rbeta(1,alphaB_old,betaB_old)*(upper_bound-lower_bound)+lower_bound
        } else if ( (r_jump_dist_name=="unif") || (r_jump_dist_name=="unif_one") || (r_jump_dist_name=="unif11") ) {
          r_k_jump <- runif(1,lower_bound,upper_bound)
        }
        
        
        
        R_star <- R_special
        R_star[uniqueR_locs_all[[k]]] <- r_k_jump
        Sigma_star <- S_special%*%R_star%*%S_special
        
        log_lik_old <- 0
        log_lik_star <- 0
        for (i in 1:N) {
          p_i <- ps[i]
          log_lik_old <- log_lik_old + dmvnorm(y_imputed[[i]],
                                               mu_special[1:p_i],
                                               Sigma_old[1:p_i,1:p_i],
                                               log=TRUE)
          log_lik_star <- log_lik_star + dmvnorm(y_imputed[[i]],
                                                 mu_special[1:p_i],
                                                 Sigma_star[1:p_i,1:p_i],
                                                 log=TRUE)
        }
        
        log_prior_old <- log(1*(all(eigen(R_old)$values > 0)))
        log_prior_star <- log(1*(all(eigen(R_star)$values > 0)))
        log_prior_stars_r[[k]][itr[k]] <- log_prior_star
        
        if ( (r_jump_dist_name=="sbeta") || (r_jump_dist_name=="sbeta_one") || (r_jump_dist_name=="sbeta11") ) {
          alphaB_jump <- ( (kappa[k]-1)*lower_bound+(2-kappa[k])*r_k_jump
                           -upper_bound )/(lower_bound-upper_bound)
          betaB_jump <- kappa[k]-alphaB_jump
          
          log_fc_old <- dlogSHIFTEDbeta(r_k_old,alphaB_jump,betaB_jump,lower_bound,upper_bound)
          log_fc_star <- dlogSHIFTEDbeta(r_k_jump,alphaB_old,betaB_old,lower_bound,upper_bound)
        } else if ( (r_jump_dist_name=="unif") || (r_jump_dist_name=="unif_one") || (r_jump_dist_name=="unif11") ) {
          log_fc_old <- log(dunif(r_k_old,lower_bound,upper_bound))
          log_fc_star <- log(dunif(r_k_jump,lower_bound,upper_bound))
        }
        
        log_ratio_r <- (log_lik_star+log_prior_star-log_fc_star)-
          (log_lik_old+log_prior_old-log_fc_old)
        
        log_alpha_r <- min(log(1),log_ratio_r)
        
        u <- runif(1,min=0,max=1)
        if (ru < r_updates[k]) {
          if (log_alpha_r>log(u)) {
            r[it-1,k] <- r_k_jump
            success_count_r[[k]][itr[k]] <- 1
            r_all[[k]][itr[k]] <- r_k_jump
          } else if (log_alpha_r<=log(u)) {
            r_all[[k]][itr[k]] <- r_k_old
          }
        } else if (ru == r_updates[k]) {
          if (log_alpha_r>log(u)) {
            r[it,k] <- r_k_jump
            success_count_r[[k]][itr[k]] <- 1
            r_all[[k]][itr[k]] <- r_k_jump
          } else if (log_alpha_r<=log(u)) {
            r[it,k] <- r_k_old
            r_all[[k]][itr[k]] <- r_k_old
          }
        }
      }
      
      R_special[uniqueR_locs_all[[k]]] <- r_all[[k]][itr[k]]
      
    }
    
    Sigma_special <- S_special%*%R_special%*%S_special
    Sigma_special.inv <- solve(Sigma_special)
  }
  
  r_actual <- matrix(nrow=iters,ncol=Q)
  for (k in 1:Q) {
    r_actual[,k] <- r_all[[k]][seq(1,length(r_all[[k]]),r_updates[k])]
  }
  
  return(list(
    y=y, N=N, iters=iters, burn_in_prop=burn_in_prop,
    r_jump_dist_name=r_jump_dist_name,
    r_updates=r_updates, itr_names=itr_names,
    Js=Js, L=L, Q=Q, ps=ps,
    J_max=J_max, p_max=p_max, j_order0_max=j_order0_max, l_order0_max=l_order0_max,
    target_short_mu=target_short_mu, target_short_s=target_short_s, target_r=target_r,
    target_w=target_w, target_SRM=target_SRM,
    nu=nu, kappa=kappa,
    mu_names_detailed=mu_names_detailed, s_names_detailed=s_names_detailed,
    r_names_detailed=r_names_detailed,
    short_mu=short_mu, short_s=short_s, 
    r_actual=r_actual,

    success_count_s=success_count_s,
    success_count_r=success_count_r, log_prior_stars_r=log_prior_stars_r
    ))
}

start_time <- Sys.time()
list_same_dataset <- lapply(1:nchains, function(x) y_etc)
mult_chains <- mclapply(list_same_dataset, mvn_disc, mc.cores=ncores)
end_time <- Sys.time()

setwd(parent_source)
dir.create(output_dir)
setwd(output_dir)
time_txt(mult_chains,ff_dat_name,burn_in_prop,nchains,ncores,start_time,end_time)
setwd(code_dir)










setwd(output_dir)
saveRDS(mult_chains, file=paste0("mcmc_chains.RDS"))
#mult_chains <- readRDS(paste0(output_dir,"/mcmc_chains.RDS"))
setwd(code_dir)



y <- mult_chains[[1]]$y
N <- mult_chains[[1]]$N
iters <- mult_chains[[1]]$iters
r_jump_dist_name <- mult_chains[[1]]$r_jump_dist_name
r_updates <- mult_chains[[1]]$r_updates
itr_names <- mult_chains[[1]]$itr_names
Js <- mult_chains[[1]]$Js
L <- mult_chains[[1]]$L
Q <- mult_chains[[1]]$Q
ps <- mult_chains[[1]]$ps
J_max <- mult_chains[[1]]$J_max
p_max <- mult_chains[[1]]$p_max
j_order0_max <- mult_chains[[1]]$j_order0_max
l_order0_max <- mult_chains[[1]]$l_order0_max
nu <- mult_chains[[1]]$nu
kappa <- mult_chains[[1]]$kappa
mu_names_detailed <- mult_chains[[1]]$mu_names_detailed
s_names_detailed <- mult_chains[[1]]$s_names_detailed
r_names_detailed <- mult_chains[[1]]$r_names_detailed





setwd(output_dir)
main_diagnostics_stuff <- main_diagnostics(mult_chains, nchains, burn_in_prop)
MH_accept_percentage_s <- main_diagnostics_stuff$MH_accept_percentage_s
MH_accept_percentage_r <- main_diagnostics_stuff$MH_accept_percentage_r
R_PD_percentage <- main_diagnostics_stuff$R_PD_percentage
write.csv(MH_accept_percentage_s, paste0("s_accept%.csv"))
write.csv(kappa, paste0("kappa.csv"))
write.csv(MH_accept_percentage_r, paste0("r_accept%.csv"))
write.csv(R_PD_percentage, paste0("%times_R_PD_for_r_jump.csv"))
setwd(code_dir)










# Prepare objects for plotting using functions from 'bayesplot' library
r_list_actual <- list()
for (c_indx in 1:nchains) {
  r_actual <- mult_chains[[c_indx]]$r_actual
  colnames(r_actual) <- r_names_detailed
  r_list_actual[[c_indx]] <- r_actual
}
s_list <- list()
for (c_indx in 1:nchains) {
  short_s <- mult_chains[[c_indx]]$short_s
  colnames(short_s) <- s_names_detailed
  s_list[[c_indx]] <- short_s
}
mu_list <- list()
for (c_indx in 1:nchains) {
  short_mu <- mult_chains[[c_indx]]$short_mu
  colnames(short_mu) <- mu_names_detailed
  mu_list[[c_indx]] <- short_mu
}

setwd(output_dir)
dir.create(paste0(output_dir,"/","target_params"))
setwd(paste0(output_dir,"/","target_params"))
short_mu_post_centers <- target_param_func(mu_list, mu_names_detailed, burn_in_prop)
short_s_post_centers <- target_param_func(s_list, s_names_detailed, burn_in_prop)
r_post_centers <- target_param_func(r_list_actual, r_names_detailed, burn_in_prop)
short_mu_post_med <- short_mu_post_centers$post_med
short_s_post_med <- short_s_post_centers$post_med
r_post_med <- r_post_centers$post_med
muscle_Sigma_post_med <- get_muscle_Sigma_func(short_s_post_med, r_post_med,
                                               L, n_etas)
w_from_post_med_stuffs <- muscle_optim_func(short_mu_post_med, muscle_Sigma_post_med,
                                            optim_method, grad_bool, upper_bool)
w_from_post_med <- w_from_post_med_stuffs$w
SRMval_from_post_med <- SRM(w_from_post_med, short_mu_post_med, muscle_Sigma_post_med)
params_post_med <- list(target_short_mu = short_mu_post_med,
                        target_short_s = short_s_post_med,
                        target_r = r_post_med,
                        target_w = w_from_post_med,
                        target_SRM = SRMval_from_post_med)
saveRDS(params_post_med, file=paste0(amb_name,"_",L,"L","_","post_med",".RDS"))
setwd(code_dir)



setwd(output_dir)
dir.create(paste0(output_dir,"/","mcmc_plots"))
setwd(paste0(output_dir,"/","mcmc_plots"))
dir.create(paste0(output_dir,"/","mcmc_plots","/","r"))
setwd(paste0(output_dir,"/","mcmc_plots","/","r"))
mcmc_plots(r_list_actual, "r", r_names_detailed, ff_dat_name, rep("",Q), burn_in_prop)
setwd(code_dir)

setwd(output_dir)
setwd(paste0(output_dir,"/","mcmc_plots"))
dir.create(paste0(output_dir,"/","mcmc_plots","/","s"))
setwd(paste0(output_dir,"/","mcmc_plots","/","s"))
mcmc_plots(s_list, "s", s_names_detailed, ff_dat_name, rep("",L), burn_in_prop)
setwd(code_dir)

setwd(output_dir)
setwd(paste0(output_dir,"/","mcmc_plots"))
dir.create(paste0(output_dir,"/","mcmc_plots","/","mu"))
setwd(paste0(output_dir,"/","mcmc_plots","/","mu"))
mcmc_plots(mu_list, "mu", mu_names_detailed, ff_dat_name, rep("",L), burn_in_prop)
setwd(code_dir)










SRM_opt_stuff <- 
  SRM_opt_subset(mult_chains, nchains, burn_in_prop,
                 muscle_names_simple, muscle_names_simple_subset, r_names,
                 optim_method, grad_bool, upper_bool)
SRMw_orig <- SRM_opt_stuff$SRMw_orig
SRMw <- SRM_opt_stuff$SRMw
SRMw[which(SRMw<0,arr.ind=TRUE)] <- 0
SRMvalues <- SRM_opt_stuff$SRMvalues

setwd(output_dir)
dir.create(paste0(output_dir,"/","SRM"))
setwd(paste0(output_dir,"/","SRM"))
if (all(muscle_names_simple==muscle_names_simple_subset)) {
  opt_out_by_data(SRMw, "SRM", iters, burn_in_prop, ff_dat_name, NA)
  #opt_out_by_data(SRMw, "SRM", iters, burn_in_prop, ff_dat_name, target_w)
}
w_and_muscleSRMresults <- 
  w_and_muscleSRM_func_subset(mult_chains, burn_in_prop, SRMw, SRMvalues,
                              muscle_names_simple, muscle_names_simple_subset, 
                              ff_dat_name,
                              optim_method, grad_bool, upper_bool)
setwd(code_dir)










# MIXTURE DIAGNOSTIC
setwd(output_dir)
dir.create(paste0(output_dir,"/","mixture_diagnostic"))
setwd(paste0(output_dir,"/","mixture_diagnostic"))
mix_stuff <- mixture_diagnostic(mult_chains,chunk_length=100)
the_rounding <- 3

dir.create(paste0(output_dir,"/","mixture_diagnostic","/","mu"))
dir.create(paste0(output_dir,"/","mixture_diagnostic","/","mu","/","avg"))
setwd(paste0(output_dir,"/","mixture_diagnostic","/","mu","/","avg"))
lapply(1:L, function(l) 
  write.csv(round(mix_stuff$short_mu_avgs_all_chains[[l]], the_rounding),
            paste0(mu_names_detailed[l],".csv")))
dir.create(paste0(output_dir,"/","mixture_diagnostic","/","mu","/","sd"))
setwd(paste0(output_dir,"/","mixture_diagnostic","/","mu","/","sd"))
lapply(1:L, function(l) 
  write.csv(round(mix_stuff$short_mu_sds_all_chains[[l]], the_rounding),
            paste0(mu_names_detailed[l],".csv")))

dir.create(paste0(output_dir,"/","mixture_diagnostic","/","s"))
dir.create(paste0(output_dir,"/","mixture_diagnostic","/","s","/","avg"))
setwd(paste0(output_dir,"/","mixture_diagnostic","/","s","/","avg"))
lapply(1:L, function(l) 
  write.csv(round(mix_stuff$short_s_avgs_all_chains[[l]], the_rounding),
            paste0(s_names_detailed[l],".csv")))
dir.create(paste0(output_dir,"/","mixture_diagnostic","/","s","/","sd"))
setwd(paste0(output_dir,"/","mixture_diagnostic","/","s","/","sd"))
lapply(1:L, function(l) 
  write.csv(round(mix_stuff$short_s_sds_all_chains[[l]], the_rounding),
            paste0(s_names_detailed[l],".csv")))

dir.create(paste0(output_dir,"/","mixture_diagnostic","/","r"))
dir.create(paste0(output_dir,"/","mixture_diagnostic","/","r","/","avg"))
setwd(paste0(output_dir,"/","mixture_diagnostic","/","r","/","avg"))
lapply(1:Q, function(k) 
  write.csv(round(mix_stuff$r_actual_avgs_all_chains[[k]], the_rounding),
            paste0(r_names_detailed[k],".csv")))
dir.create(paste0(output_dir,"/","mixture_diagnostic","/","r","/","sd"))
setwd(paste0(output_dir,"/","mixture_diagnostic","/","r","/","sd"))
lapply(1:Q, function(k) 
  write.csv(round(mix_stuff$r_actual_sds_all_chains[[k]], the_rounding),
            paste0(r_names_detailed[k],".csv")))

setwd(code_dir)










# Ran gemplot stuff on local computer and for shifted-beta jumping distribution:
if( (Sys.getenv("RSTUDIO")=="1") && (detectCores()==12) && (r_jump_dist_name %in% sbeta_strings) ) {
  if (L==3) {
    cc_names <- c("x","y")
    vertices <- matrix(c(0,0,
                         1,0,
                         1/2,sqrt(3)/2,nrow=L,byrow=TRUE))
  } else if (L==4) {
    cc_names <- c("x","y","z")
    vertices <- matrix(c(0,0,0,
                         1,0,0,
                         1/2,sqrt(3)/2,0,
                         1/2,sqrt(3)/6,sqrt(6)/3),nrow=L,byrow=TRUE)
  } else if (L==5) {
    vertices <- matrix(c(0,0,0,4/sqrt(5),
                         1,1,1,-1/sqrt(5),
                         1,-1,-1,-1/sqrt(5),
                         -1,1,-1,-1/sqrt(5),
                         -1,-1,1,-1/sqrt(5)),nrow=L,byrow=TRUE)
    cc_names <- c("w","x","y","z")
  }
  
  setwd(output_dir)
  output_rgl_html_dir <- paste0(output_dir,"/","rgl_html")
  dir.create(output_rgl_html_dir)
  setwd(code_dir)
  
  start_time_hldepth <- Sys.time()
  if (L==4) {
    html_name <- paste0(output_rgl_html_dir,"/gem_klaR.html")
    
    w_names_detailed <- paste0("w",1:L,"[",muscle_names_simple,"]")
    SRMw_chain <- SRMw[,3,]
    colnames(SRMw_chain) <- w_names_detailed
    
    cc <- SRMw_chain %*% vertices
    colnames(cc) <- cc_names
    D_cc <- as.data.frame(cc)
    D_cc <- D_cc[(nrow(SRMw_chain)-1000+1):nrow(SRMw_chain),]
    
    ### Grid size and graphic parameters
    grid.size <- 20
    red <- rgb(200, 100, 100, alpha = 100, maxColorValue = 255)
    blue <- rgb(100, 100, 200, alpha = 100, maxColorValue = 255)
    yel <- rgb(255, 255, 102, alpha = 100, maxColorValue = 255)
    white <- rgb(255, 255, 255, alpha = 100, maxColorValue = 255)
    material3d(color=c(red, blue, yel, white), alpha=c(0.5, 0.5, 0.5, 0.5), smooth=FALSE, specular="black")
    
    ### Calucation and visualization of gemplot for D_cc
    G <- gridfun(D_cc, grid.size=20)
    G$H <- hldepth(D_cc, G, verbose=TRUE)
    dm <- depmed(G)
    B <- bag(D_cc, G)
    the_loop <- loop(D_cc, B, inflation = 3, dm = dm)
    
    rgl.open()
    bg3d(color="white")
    points3d(D_cc[the_loop$outliers==0,1], D_cc[the_loop$outliers==0,2], D_cc[the_loop$outliers==0,3],
             col=blue) ### red
    lines3d(vertices[c(1,2,3,4,1,3,1,2,4),1],
            vertices[c(1,2,3,4,1,3,1,2,4),2],
            vertices[c(1,2,3,4,1,3,1,2,4),3],
            col="black")
    text3d(vertices[,1], vertices[,2], vertices[,3],
           w_names_detailed, col="black")
    text3d(D_cc[the_loop$outliers==1,1], D_cc[the_loop$outliers==1,2], D_cc[the_loop$outliers==1,3],
           rep("o",length(which(the_loop$outliers==1))), col="black")
    gem(B$coords, B$hull, yel)
    gem(the_loop$coords.loop, the_loop$hull.loop, red)
    axes3d(col="black")
    
    # Indicate the point of equal weights in gemPlot
    cc_equal_w <- matrix(rep(1/L,L) %*% vertices, nrow=1, ncol=length(cc_names))
    colnames(cc_equal_w) <- cc_names
    D_cc_equal_w <- as.data.frame(cc_equal_w)
    text3d(D_cc_equal_w[,1], D_cc_equal_w[,2], D_cc_equal_w[,3],
           "E", col="black")
    
    
    setwd(output_rgl_html_dir)
    widgets <- rglwidget()
    htmltools::save_html(widgets,html_name)
    setwd(code_dir)
  }
  end_time_hldepth <- Sys.time()
  difftime(end_time_hldepth, start_time_hldepth, units='hours')
  
}