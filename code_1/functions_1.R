prepare_dat_groups <- function(actual_data, muscle_names, measure_name) {
  # Prepare data by function/ambulatory groups
  
  # annualized_delta = measure_value_diff/age_diff, where
  # measure_value_diff = measure_value_xfer at current visit - measure_value_xfer at previous visit
  # age_diff = age at current visit - age at previous visit
  # It is only reported where visits are >6 and <24 months apart.
  
  extra_variables <- c("person_id","age_at_visit","function_group")
  n_extra_vars <- length(extra_variables)
  column_names <- c(extra_variables, paste0(muscle_names,".",measure_name))
  
  # use (0,1) scaled 'measure_value_xfer' to compute a scaled 'annualized_delta' for non-MRS_FF variables
  L <- length(muscle_names)
  MRS_FF_nchar <- nchar("MRS_FF")
  weird_column_indxs <- c()
  new_column_names <- column_names
  MRS_FF_col_indxs <- (n_extra_vars+1):length(column_names)
  for (l in 1:L) { # L should == length(MRS_FF_col_indxs)
    muscle_name <- muscle_names[l]
    if ((substr(muscle_name,1,MRS_FF_nchar)!="MRS_FF") && (measure_name=="annualized_delta")) {
      # for computation
      weird_column_indxs <- append(weird_column_indxs,n_extra_vars+l) # think I have to keep this 'append' function
      # for accessing 'actual_data' to create 'ff_all'
      column_names[n_extra_vars+l] <- paste0(muscle_name,".","measure_value_xfer")
      # for later naming edited data.frame of 'ff_all_new'
      new_column_names[n_extra_vars+l] <- paste0(muscle_name,".","annualized_delta2")
      # for later deleting rows where all FF's are NA
      MRS_FF_col_indxs <- MRS_FF_col_indxs[-l]
    }
  }
  
  # have to calculate 'annualized_delta' across the ENTIRE dataset for the 'weird_column_indxs'
  # (i.e. think of 2 rows where first row is Early Ambulatory and next row is Late/Non Ambulatory)
  ff_all <- actual_data[,column_names] # keep as data.frame because variables like "function_group" would make all elems of matrix strings
  
  # delete original row indexing
  row.names(ff_all) <- NULL
  
  # order by 'person_id' then  'age_at_visit' so we can compute 'annualized_delta' from 'measure_value_xfer' for non-MRS_FF variables
  # [This is also just good to do]
  ff_all <- ff_all[order(ff_all[,"person_id"],ff_all[,"age_at_visit"],decreasing=FALSE),]
  
  if (length(weird_column_indxs)>=1) {
    # scale 'measure_value_xfer' for non-MRS_FF variables (currently 'X6MWT') so they're between 0 and 1
    for (weird_col_indx in weird_column_indxs) {
      ff_all[,weird_col_indx] <- ff_all[,weird_col_indx]/max(ff_all[,weird_col_indx],na.rm=TRUE)
    }
    
    # compute 'annualized_delta' for non-MRS_FF variables (currently 'X6MWT') AFTER scaling all the 'measure_value_xfer' values
    uniq_ids_before <- unique(ff_all[,"person_id"])
    ff_all_new <- NULL
    for (i in 1:length(uniq_ids_before)) {
      ff_all_i <- ff_all[ff_all[,"person_id"]==uniq_ids_before[i],]
      visit_times_i <- nrow(ff_all_i)
      if (visit_times_i>1) {
        for (weird_col_indx in weird_column_indxs) {
          AD_i_temp <- c()
          for (row_indx in 1:(visit_times_i-1)) {
            AD_i_temp[row_indx] <- 
              (ff_all_i[row_indx+1,weird_col_indx]-ff_all_i[row_indx,weird_col_indx])/
              (ff_all_i[row_indx+1,"age_at_visit"]-ff_all_i[row_indx,"age_at_visit"])
          }
          AD_i_temp[visit_times_i] <- NA
          ff_all_i[,weird_col_indx] <- AD_i_temp
        }
      } else if (visit_times_i==1) {
        for (weird_col_indx in weird_column_indxs) {
          ff_all_i[,weird_col_indx] <- NA
        }
      }
      
      ff_all_new <- rbind(ff_all_new,ff_all_i)
    }
    colnames(ff_all_new) <- new_column_names
  } else {
    ff_all_new <- ff_all
  }
  
  # delete rows where all FF's are NA 
  # (do this after scaling 'measure_value_xfer' and calculating new 'annualized_delta2' 
  # so that you have all the necessary 'measure_value_xfer' values)
  # this also automatically ensures you don't have a full row of NA's
  ff_all_new <- ff_all_new[rowSums(is.na(ff_all_new[,MRS_FF_col_indxs]))!= length(MRS_FF_col_indxs), ]
  
  
  
  the_column_names <- colnames(ff_all_new)
  # remove categorical variable, 'function_group' to not have matrix of strings
  # (need a numeric matrix or data.frame to work with functions like 'colMeans' or 'sweep')
  ff1 <- as.matrix(ff_all_new[ff_all_new[,"function_group"] %in% "Early Ambulatory", 
                              the_column_names[! the_column_names %in% c("function_group")]])
  ff2 <- as.matrix(ff_all_new[ff_all_new[,"function_group"] %in% "Late Ambulatory", 
                              the_column_names[! the_column_names %in% c("function_group")]])
  ff3 <- as.matrix(ff_all_new[ff_all_new[,"function_group"] %in% "Nonambulatory", 
                              the_column_names[! the_column_names %in% c("function_group")]])
  
  # delete original row indexing
  row.names(ff1) <- NULL; row.names(ff2) <- NULL; row.names(ff3) <- NULL
  
  ff_dats <- list(ff1, ff2, ff3)
  n_amb_stages <- length(ff_dats)
  for (amb_indx in 1:n_amb_stages) {
    ff_dat <- ff_dats[[amb_indx]]
    
    cols_w_neg_mean <- as.numeric(which(colMeans(ff_dat,na.rm=TRUE)<0))
    ff_dat[,cols_w_neg_mean] <- -1*ff_dat[,cols_w_neg_mean]
    
    ff_dats[[amb_indx]] <- ff_dat
  }
  
  y_actuals <- list()
  for (amb_indx in 1:n_amb_stages) {
    dat <- ff_dats[[amb_indx]] # keep naming as 'dat' rather than 'ff_dat' if future testing needed
    uniq_ids <- unique(dat[,"person_id"])
    N_dat <- length(uniq_ids)
    
    y_vecs <- list()
    for (i in 1:N_dat) {
      y_mat_temp <- dat[dat[,"person_id"]==uniq_ids[i],
                        the_column_names[! the_column_names %in% extra_variables]] # remove the unnecessary variables
      y_vec_temp <- c(t(y_mat_temp))
      if ((!all(is.na(y_vec_temp))) && (length(y_vec_temp)!=0))  { # don't add if all NA or if y_vec_temp is empty (PROBABLY NOT NEEDED)
        y_vecs[[i]] <- y_vec_temp
      } else { # FOR TESTING (HOPEFULLY THIS DOESN'T HAPPEN)
        print(paste0("ff",amb_indx,"'s subject i=",i," has all p_i outcomes = NA"))
      }
    }
    y_actuals[[amb_indx]] <- y_vecs
  }
  
  # mostly for testing purposes
  N_actuals_before <- sapply(1:n_amb_stages, function(amb_indx) length(unique(ff_dats[[amb_indx]][,"person_id"])))
  N_actuals_after <- sapply(1:n_amb_stages, function(amb_indx) length(y_actuals[[amb_indx]]))
  
  return(y_actuals)
}

more_neg_than_pos <- function(some_vec) {
  length(which(some_vec<=0)) > length(which(some_vec>0)) 
}










R_order_latex <- function(p, js, ls, matrix_format="upper") {
  matrix_structure <- diag(nrow=p,ncol=p)
  for (m in 1:p) {
    for (n in 1:p) {
      row_j <- js[m]
      row_l <- ls[m]
      col_j <- js[n]
      col_l <- ls[n]
      if ((row_j==col_j) && (row_l!=col_l)) {
        smaller_l <- min(row_l, col_l)
        bigger_l <- max(row_l, col_l)
        matrix_structure[m,n] <- paste0("eta_{",smaller_l,bigger_l,"}")
      } else if ((row_j!=col_j) && (row_l==col_l)) {
        the_l <- row_l
        matrix_structure[m,n] <- paste0("rho_{(",the_l,")}")
      } else if ((row_j!=col_j) && (row_l!=col_l)) {
        matrix_structure[m,n] <- paste0("gamma")
      }
    }
  }
  
  if (matrix_format=="all") {
    return(matrix_structure)
  } else if (matrix_format=="upper") {
    upper_matrix_structure <- matrix_structure
    upper_matrix_structure[lower.tri(matrix_structure)] <- NA
    return(upper_matrix_structure)
  } else if (matrix_format=="lower") {
    lower_matrix_structure <- matrix_structure
    lower_matrix_structure[upper.tri(matrix_structure)] <- NA
    return(lower_matrix_structure)
  }
}



R_order_latex_upper <- function(p, js, ls) {
  upper_matrix_structure <- diag(nrow=p,ncol=p)
  for (m in 1:p) {
    for (n in 1:p) {
      if (m < n) {
        row_j <- js[m]
        row_l <- ls[m]
        col_j <- js[n]
        col_l <- ls[n]
        if ((row_j==col_j) && (row_l!=col_l)) {
          smaller_l <- min(row_l, col_l)
          bigger_l <- max(row_l, col_l)
          upper_matrix_structure[m,n] <- paste0("eta_{",smaller_l,bigger_l,"}")
        } else if ((row_j!=col_j) && (row_l==col_l)) {
          the_l <- row_l
          upper_matrix_structure[m,n] <- paste0("rho_{(",the_l,")}")
        } else if ((row_j!=col_j) && (row_l!=col_l)) {
          upper_matrix_structure[m,n] <- paste0("gamma")
        } 
      }
    }
  }
  upper_matrix_structure[lower.tri(upper_matrix_structure)] <- NA
  return(upper_matrix_structure)
}



uniqueR_locs_func <- function(R_latex_something, r_names) {
  unique_elems_R <- r_names
  uniqueR_locs <- list()
  for (u in 1:length(unique_elems_R)) {
    elem_locs <- which(R_latex_something==unique_elems_R[u],arr.ind=TRUE)
    if (nrow(elem_locs)>1) { 
      # sort
      elem_locs <- elem_locs[order(elem_locs[,1],elem_locs[,2],decreasing=FALSE),]
    }
    uniqueR_locs[[u]] <- elem_locs
  }
  return(uniqueR_locs)
}



uniqueS_locs_func <- function(l_order0) {
  S_mat <- diag(l_order0)
  
  unique_elems_S <- unique(l_order0)
  uniqueS_locs <- list()
  for (u in 1:length(unique_elems_S)) {
    elem_locs <- which(S_mat==unique_elems_S[u],arr.ind=TRUE)
    if (nrow(elem_locs)>1) { 
      # sort
      elem_locs <- elem_locs[order(elem_locs[,1],elem_locs[,2],decreasing=FALSE),]
    }
    uniqueS_locs[[u]] <- elem_locs
  }
  return(uniqueS_locs)
}



pivot_and_implied_locs <- function(R_latex_upper, r_names) {
  unique_elems_R <- unique(c(t(R_latex_upper))) # transpose lets us read unique/pivot elems by row
  unique_elems_R <- unique_elems_R[!unique_elems_R %in% NA]
  unique_elems_R <- unique_elems_R[!unique_elems_R %in% "1"]
  r_order <- match(r_names, unique_elems_R)
  unique_elems_R <- unique_elems_R[r_order]
  #r_names == unique_elems_R
  # As seen by the commented out code '#r_names == unique_elems_R',
  # the 'unique_elems_R' is the same as 'r_names' because R_latex_upper is the correlation matrix represented in terms of strings.
  # However note that the above logic is useful if the correlation matrix wasn't a matrix of strings
  # but rather a matrix of actual numbers, and you needed the 'unique_elems_R' (numbers) from the correlation matrix
  # to be listed in the same ordering as in 'r_names'
  
  # 'pivot_locs' listed in same order as in 'r_names'
  pivot_locs <- matrix(nrow=length(unique_elems_R),ncol=2)
  implied_locs <- list()
  for (u in 1:length(unique_elems_R)) {
    elem_locs <- which(R_latex_upper==unique_elems_R[u],arr.ind=TRUE)
    if (nrow(elem_locs)>1) { 
      # sort s.t. first row of matrix = pivot location
      elem_locs <- elem_locs[order(elem_locs[,1],elem_locs[,2],decreasing=FALSE),]
    }
    # automatically stores as row in 'pivot_locs' matrix
    pivot_locs[u,] <- elem_locs[1,]
    # coerce that every list's element = matrix form
    if (is.vector(elem_locs[-1,])) {
      # case 1: only one implied location => turn vector into matrix form
      implied_locs[[u]] <- t(as.matrix(elem_locs[-1,]))
    } else {
      # case 2: 0 implied location => already matrix form
      # case 3: >2 implied locations => already matrix form
      implied_locs[[u]] <- elem_locs[-1,]
    }
  }
  
  return(list(pivot_locs, implied_locs))
}



eta_biggest_submatrix_indxs <- function(goal_indx1, goal_indx2, J_max, L,
                                        j_order0_max, l_order0_max) {
  goal_row_l <- l_order0_max[goal_indx1]
  goal_col_l <- l_order0_max[goal_indx2]
  smaller_l <- min(goal_row_l, goal_col_l)
  bigger_l <- max(goal_row_l, goal_col_l)
  
  omit_l_mat <- gtools::combinations(n=2, r=J_max-1, v=c(smaller_l,bigger_l),
                                     repeats=TRUE)
  subblock_length <- L*J_max-J_max+1
  
  p_max <- J_max*L
  jl0_mat <- matrix(c(j_order0_max, l_order0_max),ncol=2)
  jl0_keys <- apply(jl0_mat, 1, paste0, collapse='')
  jl0_values <- 1:p_max
  jl0_hash <- as.list(jl0_values)
  names(jl0_hash) <- jl0_keys
  
  j_order1_max <- c(rep(1,L), rep(2:J_max, each=L-1), 2:J_max)
  subblock_indxs_mat <- matrix(nrow=nrow(omit_l_mat), ncol=subblock_length)
  the_muscles <- 1:L
  for (row_indx in 1:nrow(omit_l_mat)) {
    l_order1_max <- c(smaller_l,bigger_l,the_muscles[-c(smaller_l,bigger_l)])
    for (col_indx in 1:ncol(omit_l_mat)) {
      the_muscles_m1 <- the_muscles[-omit_l_mat[row_indx,col_indx]]
      l_order1_max <- c(l_order1_max, the_muscles_m1)
    }
    l_order1_max <- c(l_order1_max, omit_l_mat[row_indx,])
    jl1_mat <- matrix(c(j_order1_max, l_order1_max),ncol=2)
    jl1_keys <- apply(jl1_mat, 1, paste0, collapse='')
    jl1_values <- as.numeric(jl0_hash[jl1_keys])
    final_indxs <- jl1_values
    subblock_indxs <- final_indxs[1:subblock_length]
    subblock_indxs_mat[row_indx,] <- subblock_indxs
    
  }
  
  return(subblock_indxs_mat)
}



rho_biggest_submatrix_indxs <- function(goal_indx1, goal_indx2, p_max,
                                        j_order0_max, l_order0_max) {
  goal_row_l <- l_order0_max[goal_indx1]
  the_l <- goal_row_l
  goal_corr <- paste0("rho_{(",the_l,")}")
  
  subblock_indxs <- c(goal_indx1, goal_indx2)
  all_indxs <- 1:p_max
  indxs_to_go_through <- all_indxs[-subblock_indxs]
  for (m in indxs_to_go_through) {
    j_poss_fill <- j_order0_max[c(subblock_indxs,m)]
    l_poss_fill <- l_order0_max[c(subblock_indxs,m)]
    
    R_poss <- R_order_latex(length(j_poss_fill), j_poss_fill, l_poss_fill, "upper")
    goal_corr_freq <- nrow(which(R_poss==goal_corr, arr.ind=TRUE))
    if (goal_corr_freq==1) {
      subblock_indxs <- append(subblock_indxs,m)
    }
  }
  return(subblock_indxs)
}



gamma_subs_latex_func <- function(n_etas2, gamma_sub_ncol, eta_names_twice, rho_names_twice) {
  gamma_subs_latex <- NULL
  for (m_indx in 1:n_etas2) {
    gamma_sub_latex_temp <- matrix(c("1","gamma",eta_names_twice[m_indx],
                                     "gamma","1",rho_names_twice[m_indx],
                                     eta_names_twice[m_indx],rho_names_twice[m_indx],"1"),
                                   nrow=gamma_sub_ncol,ncol=gamma_sub_ncol,byrow=TRUE)
    gamma_subs_latex <- cbind(gamma_subs_latex, gamma_sub_latex_temp)
  }
  return(gamma_subs_latex)
}



unique_subb_locs_func <- function(r_names, gamma_subs_latex) {
  unique_elems_R <- r_names
  unique_subb_locs <- list()
  for (u in 1:length(unique_elems_R)) {
    elem_locs <- which(gamma_subs_latex==unique_elems_R[u],arr.ind=TRUE)
    if (nrow(elem_locs)>1) { 
      # sort
      elem_locs <- elem_locs[order(elem_locs[,1],elem_locs[,2],decreasing=FALSE),]
    }
    unique_subb_locs[[u]] <- elem_locs
  }
  return(unique_subb_locs)
}



gamma_subs_func <- function(Q, r_vec, empty_gamma_subs, unique_subb_locs, gamma_sub_ncol) {
  for (k in 1:(Q-1)) {
    empty_gamma_subs[unique_subb_locs[[k]]] <- r_vec[k]
  }
  
  return(empty_gamma_subs)
}



barnard <- function(i,j,R) {
  R_1 <- R; R_m1 <- R; R_0 <- R
  R_1[i,j] <- 1; R_1[j,i] <- 1
  R_m1[i,j] <- -1; R_m1[j,i] <- -1
  R_0[i,j] <- 0; R_0[j,i] <- 0
  det.R_1 <- det(R_1)
  det.R_m1 <- det(R_m1)
  det.R_0 <- det(R_0)
  
  a <- 0.5*(det.R_1 + det.R_m1 - 2*det.R_0)
  b <- 0.5*(det.R_1 - det.R_m1)
  c <- det.R_0
  sqrt.b24ac <- sqrt(b^2-4*a*c)
  two_times_a <- 2*a
  
  r1 <- (-b - sqrt.b24ac)/(two_times_a)
  r2 <- (-b + sqrt.b24ac)/(two_times_a)
  root1 <- min(r1,r2)
  root2 <- max(r1,r2)
  
  interval <- c(root1,root2)
  return(interval)
}



sweep_once <- function(G, H, one_to_p1, k) { # G must be symmetric
  H[k,k] <- -1/G[k,k]
  
  indxs_to_go_through <- one_to_p1[-k]
  
  for (j in indxs_to_go_through) {
    H[j,k] <- -G[j,k]*H[k,k]
    H[k,j] <- -G[j,k]*H[k,k]
  }
  for (j in indxs_to_go_through) {
    for (l in indxs_to_go_through) {
      if (j <= l) {
        H[j,l] <- G[j,l]-H[j,k]*G[k,l]
      }
    }
  }
  H[lower.tri(H)] <- t(H)[lower.tri(H)]
  
  return(H)
}



# DON'T DELETE
rtoR_old <- function(r_vec, r_names, p, j_order, l_order) {
  #r_names <- c("eta_{12}", "eta_{13}", "eta_{14}",
  #             "eta_{23}", "eta_{24}", "eta_{34}",
  #             "rho_{(1)}", "rho_{(2)}", "rho_{(3)}", "rho_{(4)}",
  #             "gamma")
  #names(r_vec) <- r_names
  R_mat <- diag(nrow=p,ncol=p)
  for (m in 1:p) {
    for (n in 1:p) {
      row_j <- j_order[m]
      row_l <- l_order[m]
      col_j <- j_order[n]
      col_l <- l_order[n]
      
      if ((row_j==col_j) && (row_l!=col_l)) {
        smaller_l <- min(row_l, col_l)
        bigger_l <- max(row_l, col_l)
        r_name <- paste0("eta_{",smaller_l,bigger_l,"}")
        r_indx <- which(r_names %in% r_name)
        R_mat[m,n] <- r_vec[r_indx]
      } else if ((row_j!=col_j) && (row_l==col_l)) {
        the_l <- row_l
        r_name <- paste0("rho_{(",the_l,")}")
        r_indx <- which(r_names %in% r_name)
        R_mat[m,n] <- r_vec[r_indx]
      } else if ((row_j!=col_j) && (row_l!=col_l)) {
        r_name <- paste0("gamma")
        r_indx <- which(r_names %in% r_name)
        R_mat[m,n] <- r_vec[r_indx]
      }
    }
  }
  
  return(R_mat)
}



rtoR <- function(r_vec, Q, I_mat, uniqueR_locs) {
  for (k in 1:Q) {
    I_mat[uniqueR_locs[[k]]] <- r_vec[k]
  }
  return(I_mat) # really returning the R_mat
}



dlogSHIFTEDbeta <- function(y,alphaB,betaB,a,b) {
  #the_density <- (alphaB-1)*log(y-a)+(betaB-1)*log(b-y)-(alphaB+betaB-1)*log(b-a)-lbeta(alphaB,betaB)
  the_density <- (alphaB-1)*log(y-a)+(betaB-1)*log(b-y)-(alphaB+betaB-1)*log(b-a)-(lgamma(alphaB)+lgamma(betaB)-lgamma(alphaB+betaB))
  return(the_density)
}










time_txt <- function(mult_chains, ff_dat_name, burn_in_prop, nchains, ncores,
                     start_time, end_time, id="") {
  if (id!="") {
    id <- mult_chains[[1]]$id
  }
  orig_iters <- mult_chains[[1]]$iters
  burn_in_number <- orig_iters*burn_in_prop
  r_updates <- mult_chains[[1]]$r_updates
  
  r_jump_dist_name <- mult_chains[[1]]$r_jump_dist_name
  L <- mult_chains[[1]]$L
  N <- mult_chains[[1]]$N
  J_max <- mult_chains[[1]]$J_max
  
  nu <- mult_chains[[1]]$nu
  kappa <- mult_chains[[1]]$kappa
  
  runtime <- paste0(difftime(end_time, start_time, units='hours'))
  
  DA_MH_time <- paste0("id=",id,", ",ff_dat_name," > L=",L, " > r_JUMP_DIST=",r_jump_dist_name,
                       " (N=",N,", J_max=",J_max,", Q=",Q,")",
                       "; chains=",nchains,
                       ", cores=", ncores,
                       "; orig_iters=",orig_iters,
                       ", burn-in=",burn_in_number,
                       ", r_updates=[",paste0(r_updates,collapse="|"),"]",
                       "; nu=[",paste0(nu,collapse="|"),"]",
                       "; kappa=[",paste0(kappa,collapse="|"),"]",
                       "; ",runtime, " hours")
  
  write(DA_MH_time, file=paste0("time.txt"), append=TRUE)
  write("", file=paste0("time.txt"), append=TRUE)
}



main_diagnostics <- function(mult_chains, nchains, burn_in_prop) {
  orig_iters <- mult_chains[[1]]$iters
  burn_in_number <- orig_iters*burn_in_prop
  r_updates <- mult_chains[[1]]$r_updates
  L <- mult_chains[[1]]$L
  Q <- mult_chains[[1]]$Q
  r_names_simple <- mult_chains[[1]]$r_names_simple
  
  MH_accept_percentage_s <- matrix(nrow=L, ncol=nchains)
  for (c_indx in 1:nchains) {
    success_count_s <- mult_chains[[c_indx]]$success_count_s
    for (l in 1:L) {
      # to test nu_l
      # Calculate after burn-in:
      MH_accept_percentage_s[l,c_indx] <- sum(success_count_s[(burn_in_number+1):orig_iters,l])/(orig_iters-burn_in_number)
    }
  }
  MH_accept_percentage_s <- cbind(MH_accept_percentage_s, rowMeans(MH_accept_percentage_s))
  MH_accept_percentage_s <- round(MH_accept_percentage_s,3)
  rownames(MH_accept_percentage_s) <- paste0("l=",1:L)
  colnames(MH_accept_percentage_s) <- c(paste0("chain",1:nchains),"average")
  
  r_update_burnin_indxs <- lapply(1:Q, function(k)
    #(burn_in_number+1):(1+(orig_iters-1)*r_updates[k]))
    # accounts for r_updates when doing burn-in
    # first index below same as 'which(itr_names[[k]]=="iter1001_upd1")'
    (1+(burn_in_number-1)*r_updates[k]+1):(1+(orig_iters-1)*r_updates[k]))
  
  
  MH_accept_percentage_r <- matrix(nrow=Q, ncol=nchains)
  for (c_indx in 1:nchains) {
    success_count_r <- mult_chains[[c_indx]]$success_count_r
    for (k in 1:Q) {
      # to test kappa_q
      # Calculate after burn-in:
      MH_accept_percentage_r[k,c_indx] <-
        sum(success_count_r[[k]][r_update_burnin_indxs[[k]]])/length(r_update_burnin_indxs[[k]])
    }
  }
  MH_accept_percentage_r <- cbind(MH_accept_percentage_r, rowMeans(MH_accept_percentage_r))
  MH_accept_percentage_r <- round(MH_accept_percentage_r,3)
  rownames(MH_accept_percentage_r) <- paste0("k=",1:Q,": ",r_names_simple)
  colnames(MH_accept_percentage_r) <- c(paste0("chain",1:nchains),"average")

  R_PD_percentage <- matrix(nrow=Q, ncol=nchains)
  for (c_indx in 1:nchains) {
    log_prior_stars_r <- mult_chains[[c_indx]]$log_prior_stars_r
    for (k in 1:Q) {
      # Calculate after burn-in:
      R_PD_percentage[k,c_indx] <-
        sum(exp(log_prior_stars_r[[k]])[r_update_burnin_indxs[[k]]])/length(r_update_burnin_indxs[[k]])
    }
  }
  R_PD_percentage <- cbind(R_PD_percentage, rowMeans(R_PD_percentage))
  R_PD_percentage <- round(R_PD_percentage,3)
  rownames(R_PD_percentage) <- paste0("k=",1:Q,": ",r_names_simple)
  colnames(R_PD_percentage) <- c(paste0("chain",1:nchains),"average")

  return(list(MH_accept_percentage_s=MH_accept_percentage_s,
              MH_accept_percentage_r=MH_accept_percentage_r,
              R_PD_percentage=R_PD_percentage))
}



target_param_func <- function(mcmc_list_orig, param_names_detailed, burn_in_prop) {
  orig_iters <- nrow(mcmc_list_orig[[1]])
  burn_in_number <- orig_iters*burn_in_prop
  n_params <- ncol(mcmc_list_orig[[1]])
  nchains <- length(mcmc_list_orig)
  
  # Account for burn-in:
  mcmc_list_burnin <- vector("list", length=nchains)
  for (c_indx in 1:nchains) {
    mcmc_list_burnin[[c_indx]] <- mcmc_list_orig[[c_indx]][(burn_in_number+1):orig_iters,]
  }
  # redefine mcmc_list to only iterations after burn-in
  mcmc_list <- mcmc_list_burnin
  
  mcmc_array <- array(dim=c(orig_iters-burn_in_number, nchains, n_params),
                      dimnames=list(paste0("iter",(burn_in_number+1):orig_iters),
                                    paste0("chain",1:nchains),
                                    param_names_detailed))
  for (c_indx in 1:nchains) {
    for (p_indx in 1:n_params) {
      mcmc_array[,c_indx,p_indx] <- mcmc_list[[c_indx]][,p_indx]
    }
  }
  
  inner_prob <- 0.8; outer_prob <- 0.95
  
  param_names_for_UI <- param_names_detailed
  mcmc_array_for_UI <- mcmc_array
  dimnames(mcmc_array_for_UI)[[3]] <- param_names_for_UI
  
  post_med <- as.numeric(as.matrix(as.data.frame(
    mcmc_intervals_data(mcmc_array_for_UI,
                        prob=inner_prob, prob_outer=outer_prob,
                        point_est="median")[,"m"])))
  post_mean <- as.numeric(as.matrix(as.data.frame(
    mcmc_intervals_data(mcmc_array_for_UI,
                        prob=inner_prob, prob_outer=outer_prob,
                        point_est="mean")[,"m"])))
  return(list(post_med=post_med, post_mean=post_mean))
}



mcmc_plots <- function(mcmc_list_orig, mcmc_string, param_names_detailed, ff_dat_name, target_mcmc, burn_in_prop, id="") {
  if (id!="") {
    id <- mult_chains[[1]]$id
  }
  orig_iters <- nrow(mcmc_list_orig[[1]])
  burn_in_number <- orig_iters*burn_in_prop
  n_params <- ncol(mcmc_list_orig[[1]])
  nchains <- length(mcmc_list_orig)
  
  # Part of title for plots:
  if (mcmc_string=="r") {
    mcmc_type <- "Posterior Draws of Correlations"
  } else if (mcmc_string=="s") {
    mcmc_type <- "Posterior Draws of Standard Deviations"
  } else if (mcmc_string=="mu") {
    mcmc_type <- "Posterior Draws of Means"
  }
  
  # Possible part of title for plots:
  if (is.numeric(target_mcmc)) {
    target_mcmc_rounded <- round(target_mcmc,5)
    target_mcmc_rounded <- paste0(", Target=",target_mcmc_rounded)
  } else { #if ((is.na(target_mcmc)) || (is.null(target_mcmc))) {
    target_mcmc_rounded <- rep("",length(param_names_detailed))
  }
  
  # set color scheme for plots
  color_scheme_set("brewer-Dark2")
  # darkgray, pink>red,
  #RColorBrewer::brewer.pal.info
  # colorblind FALSE: brewer-Accent, Set1 (meh...), Set3
  # colorblind TRUE: brewer-Dark2, brewer-Paired
  
  # Account for burn-in:
  mcmc_list_burnin <- vector("list", length=nchains)
  for (c_indx in 1:nchains) {
    mcmc_list_burnin[[c_indx]] <- mcmc_list_orig[[c_indx]][(burn_in_number+1):orig_iters,]
  }
  # redefine mcmc_list to only iterations after burn-in
  mcmc_list <- mcmc_list_burnin
  
  mcmc_array <- array(dim=c(orig_iters-burn_in_number, nchains, n_params),
                      dimnames=list(paste0("iter",(burn_in_number+1):orig_iters),
                                    paste0("chain",1:nchains),
                                    param_names_detailed))
  for (c_indx in 1:nchains) {
    for (p_indx in 1:n_params) {
      mcmc_array[,c_indx,p_indx] <- mcmc_list[[c_indx]][,p_indx]
    }
  }
  
  mcmc_Rhats <- numeric(n_params)
  for (p_indx in 1:n_params) {
    mcmc_Rhats[p_indx] <- Rhat(mcmc_array[,,p_indx])
  }
  mcmc_Rhats_rounded <- round(mcmc_Rhats,3)
  param_names_for_trace <- paste0(param_names_detailed,", Rhat=",mcmc_Rhats_rounded,target_mcmc_rounded)
  mcmc_array_for_trace <- mcmc_array
  dimnames(mcmc_array_for_trace)[[3]] <- param_names_for_trace
  png(file=paste0(mcmc_string,"_trace",id,".png"), width=1920, height=1027, pointsize=12)
  par(mfrow=c(1,1))
  #g <- mcmc_trace(mcmc_list) + theme_gray(base_size=20) # same results with 'mcmc_array'
  g <- mcmc_trace(mcmc_array_for_trace) + 
    ggtitle(paste0("Trace Plots for ",mcmc_type," for ",ff_dat_name)) +
    theme_gray(base_size=20)
  print(g)
  dev.off()
  
  
  png(file=paste0(mcmc_string,"_acf",id,".png"), width=1920, height=1027, pointsize=12)
  par(mfrow=c(1,1))
  # Grid of autocorrelation plots by chain and parameter. The lags argument gives the maximum number
  # of lags at which to calculate the autocorrelation function.
  # line plot:
  g <- mcmc_acf(mcmc_list) +
    ggtitle(paste0("Autocorrelation Plots for ",mcmc_type," for ",ff_dat_name)) +
    theme_gray(base_size=20)
  print(g)
  dev.off()
  
  
  param_names_for_dens <- paste0(param_names_detailed,target_mcmc_rounded)
  mcmc_array_for_dens <- mcmc_array
  dimnames(mcmc_array_for_dens)[[3]] <- param_names_for_dens
  png(file=paste0(mcmc_string,"_dens_overlay",id,".png"), width=1920, height=1027, pointsize=12)
  par(mfrow=c(1,1))
  # Kernel density plots of posterior draws with chains separated but overlaid on a single plot.
  g <- mcmc_dens_overlay(mcmc_array_for_dens) +
    ggtitle(paste0("Kernel Density Plots for ",mcmc_type," for ",ff_dat_name)) +
    theme_gray(base_size=20)
  print(g)
  dev.off()
  
  png(file=paste0(mcmc_string,"_dens",id,".png"), width=1920, height=1027, pointsize=12)
  par(mfrow=c(1,1))
  # Kernel density plots of posterior draws with chains separated but overlaid on a single plot.
  g <- mcmc_dens(mcmc_array_for_dens) +
    ggtitle(paste0("Kernel Density Plots for ",mcmc_type," for ",ff_dat_name)) +
    theme_gray(base_size=20)
  print(g)
  dev.off()
  
  
  inner_prob <- 0.8; outer_prob <- 0.95; pt_est <- "median"
  param_names_for_UI <- paste0(param_names_detailed,target_mcmc_rounded)
  mcmc_array_for_UI <- mcmc_array
  dimnames(mcmc_array_for_UI)[[3]] <- param_names_for_UI
  png(file=paste0(mcmc_string,"_UI_",pt_est,id,".png"), width=1920, height=1027, pointsize=12)
  par(mfrow=c(1,1))
  # Plots of uncertainty intervals computed from posterior draws with all chains merged.
  g <- mcmc_intervals(mcmc_array_for_UI, prob=inner_prob, prob_outer=outer_prob, point_est=pt_est) +
    ggtitle(paste0("Credible Intervals for ",mcmc_type," for ",ff_dat_name)) +
    labs(subtitle=paste0("(combined the chains)\n",
                         "(point estimate = ",pt_est, ", inner probability = ",inner_prob,", outer probability = ",outer_prob,")")) +
    theme_gray(base_size=20)
  print(g)
  dev.off()
  write.csv(as.data.frame(mcmc_intervals_data(mcmc_array_for_UI, prob=inner_prob, prob_outer=outer_prob, point_est=pt_est)),
            file=paste0(mcmc_string,"_UI_",pt_est,id,".csv"))
}













SRM <- function(w,mu,Sigma) {
  return( c(t(w)%*%mu) / c(sqrt(t(w)%*%Sigma%*%w)) )
}

SRM_grad <- function(w,mu,Sigma) {
  return( c((t(w)%*%Sigma%*%w)^{-0.5})*mu
          -c(t(w)%*%mu)*c((t(w)%*%Sigma%*%w)^{-1.5})*(Sigma%*%w) )
}

softmax <- function(w) {
  w1 <- w - max(w)
  exp(w1)/sum(exp(w1))
}

get_muscle_Sigma_func <- function(short_s, r,
                                  L, n_etas) {
  # proved in microbenchmark (this is faster than rtoR or rtoR_old method)
  short_S <- diag(short_s)
  etas <- r[1:n_etas]
  muscle_R <- diag(1,L)
  # assume l_order0 (original) ordering for R
  muscle_R[lower.tri(muscle_R,diag=FALSE)] <- etas
  muscle_R <- t(muscle_R)
  muscle_R[lower.tri(muscle_R,diag=FALSE)] <-
    t(muscle_R)[lower.tri(muscle_R)]
  muscle_Sigma <- short_S%*%muscle_R%*%short_S
  
  if (all(eigen(muscle_Sigma)$values>0)) {
    return(muscle_Sigma)
  } else {
    print("muscle_Sigma is not PD")
  }
}

muscle_optim_func <- function(muscle_mu, muscle_Sigma,
                              optim_method="L-BFGS-B", grad_bool=TRUE, upper_bool=TRUE) {
  L <- nrow(muscle_Sigma)
  if (optim_method=="L-BFGS-B") {
    lower_optim_bound <- rep(0,L)
    upper_optim_bound <- rep(1,L)
  } else {
    lower_optim_bound <- -Inf
    upper_optim_bound <- Inf
  }
  
  if (grad_bool==TRUE) {
    if (upper_bool==TRUE) {
      w_summ <- optim(par=rep(1/L,L), fn=SRM, gr=SRM_grad,
                               mu=muscle_mu, Sigma=muscle_Sigma,
                               method=optim_method,
                               lower=lower_optim_bound, upper=upper_optim_bound,
                               control=list(fnscale=-1))
    } else {
      w_summ <- optim(par=rep(1/L,L), fn=SRM, gr=SRM_grad,
                               mu=muscle_mu, Sigma=muscle_Sigma,
                               method=optim_method,
                               lower=lower_optim_bound, 
                               control=list(fnscale=-1))
    }
  } else {
    if (upper_bool==TRUE) {
      w_summ <- optim(par=rep(1/L,L), fn=SRM,
                               mu=muscle_mu, Sigma=muscle_Sigma,
                               method=optim_method,
                               lower=lower_optim_bound, upper=upper_optim_bound,
                               control=list(fnscale=-1))
    } else {
      w_summ <- optim(par=rep(1/L,L), fn=SRM,
                               mu=muscle_mu, Sigma=muscle_Sigma,
                               method=optim_method,
                               lower=lower_optim_bound, 
                               control=list(fnscale=-1))
    }
  }
  
  w_orig <- w_summ$par
  if (optim_method=="L-BFGS-B") {
    w <- abs(w_orig)/sum(abs(w_orig))
  } else {
    w <- softmax(w_orig)
  }
  #SRMval <- SRM(w, muscle_mu, muscle_Sigma)
  
  return(list(w_orig=w_orig,
              w=w#,
              #SRMval=SRMval
              ))
}



SRM_opt_subset <- function(mult_chains, nchains, burn_in_prop,
                           muscle_names_simple, muscle_names_simple_subset, r_names,
                           optim_method, grad_bool, upper_bool) {
  
  subset_indxs <- which(muscle_names_simple %in% muscle_names_simple_subset)
  L_subset <- length(subset_indxs)
  
  n_etas_subset <- choose(L_subset,2)
  muscle_combns_subset <- combn(subset_indxs,2)
  eta_names_subset <- sapply(1:ncol(muscle_combns_subset), function(combn_indx) # note: ncol(muscle_combns_subset) == n_etas_subset
    paste0("eta_{",muscle_combns_subset[1,combn_indx],muscle_combns_subset[2,combn_indx],"}"))
  rho_names_subset <- paste0("rho_{(",subset_indxs,")}")
  r_names_subset <- c(eta_names_subset, rho_names_subset, "gamma")
  
  r_subset_indxs <- which(r_names %in% r_names_subset)
  
  #w_names_simple <- paste0("w",1:L_subset)
  w_names_detailed_subset <- paste0("w",1:L_subset,"-",muscle_names_simple_subset)
  
  orig_iters <- mult_chains[[1]]$iters # rest of code inside this function accounts for BURN-IN !!!
  burn_in_number <- orig_iters*burn_in_prop
  
  SRMw_orig <- array(dim=c(orig_iters-burn_in_number,nchains,L_subset),
                     dimnames=list(paste0("iter",(burn_in_number+1):orig_iters),
                                   paste0("chain",1:nchains),
                                   w_names_detailed_subset))
  SRMw <- array(dim=c(orig_iters-burn_in_number,nchains,L_subset),
                dimnames=list(paste0("iter",(burn_in_number+1):orig_iters),
                              paste0("chain",1:nchains),
                              w_names_detailed_subset))
  SRMvalues <- matrix(nrow=orig_iters-burn_in_number, ncol=nchains, 
                      dimnames=list(paste0("iter",(burn_in_number+1):orig_iters),
                                    paste0("chain",1:nchains)))
  for (c_indx in 1:nchains) {
    short_mu <- mult_chains[[c_indx]]$short_mu
    short_mu_subset <- short_mu[,subset_indxs]
    short_s <- mult_chains[[c_indx]]$short_s
    short_s_subset <- short_s[,subset_indxs]
    r_actual <- mult_chains[[c_indx]]$r_actual
    r_actual_subset <- r_actual[,r_subset_indxs]
    for (it in 1:(orig_iters-burn_in_number)) {
      short_mu_iter_subset <- short_mu_subset[burn_in_number+it,]
      beta_tilde_iter_subset <- matrix(short_mu_iter_subset, nrow=L_subset, ncol=1, byrow=TRUE)
      short_s_iter_subset <- short_s_subset[burn_in_number+it,]
      r_iter_subset <- r_actual_subset[burn_in_number+it,]
      
      muscle_Sigma_iter_subset <- get_muscle_Sigma_func(short_s_iter_subset, r_iter_subset,
                                                        L_subset, n_etas_subset)
      w_iter_stuffs <- muscle_optim_func(beta_tilde_iter_subset, muscle_Sigma_iter_subset,
                                         optim_method, grad_bool, upper_bool)
      w_iter_orig <- w_iter_stuffs$w_orig
      w_iter <- w_iter_stuffs$w
      
      SRMw_orig[it,c_indx,] <- w_iter_orig
      SRMw[it,c_indx,] <- w_iter
      SRMvalues[it,c_indx] <- SRM(SRMw[it,c_indx,], beta_tilde_iter_subset, muscle_Sigma_iter_subset)
    }
  }
  return(list(SRMw_orig=SRMw_orig, SRMw=SRMw, SRMvalues=SRMvalues))
  
}



merge_chains <- function(x) { # from 'bayesplot' package's source code
  #parameter_names <- dimnames(x)$parameters
  parameter_names <- dimnames(x)[[3]]
  xdim <- dim(x)
  mat <- array(x, dim = c(prod(xdim[1:2]), xdim[3]))
  colnames(mat) <- parameter_names
  mat
}

opt_out_by_data <- function(SRMw, obj_func_string, orig_iters, burn_in_prop, ff_dat_name, target_w, id="") {
  if (id!="") {
    id <- mult_chains[[1]]$id
  }
  burn_in_number <- orig_iters*burn_in_prop
  nchains <- dim(SRMw)[2]
  
  SRMw_wOrigNames <- SRMw
  # Part of title for plots:
  if (obj_func_string=="SRM") {
    w_type <- "SRM Weights"
    w_type_partial <- "SRM"
  }
  # Possible part of title for plots:
  if (is.numeric(target_w)) {
    target_w_rounded <- round(target_w,5)
    target_w_rounded <- paste0(", Target=",target_w_rounded)
  } else { #if ((is.na(target_w)) || (is.null(target_w)))  {
    target_w_rounded <- rep("",dim(SRMw)[3]) # dim(SRMw)[3]==L
  }
  w_names_old <- dimnames(SRMw)[[3]]
  w_names_new <- paste0(w_names_old,target_w_rounded)
  
  # use "reflection" on density and then density plot only from 0 to 1 ('ggplot' seems to work well)
  # note that the 'from' and 'to' parameters for the 'density' function doesn't account for area under density = 1 between 0 and 1
  dimnames(SRMw)[[3]] <- w_names_new
  neg_SRMw <- -1*SRMw
  SRMw2 <- abind(SRMw,neg_SRMw,along=1)
  
  
  ### get posterior distribution of SRMw (deterministic function of beta_tilde, muscle_Sigma) ###
  png(file=paste0("w_hist",id,".png"), width = 1920, height = 1027, pointsize=12)
  par(mfrow=c(1,1))
  g <- mcmc_hist(SRMw2) +
    scale_x_continuous(breaks=seq(0,1,0.2), labels=seq(0,1,0.2), limits=c(-0.05,1.05), expand=c(0,0)) + 
    ggtitle(paste0("Histograms for ",w_type," for ", ff_dat_name)) +
    labs(subtitle=paste0("(combined the chains)")) +
    theme_gray(base_size=20)
  print(g)
  dev.off()
  
  png(file=paste0("w_dens_overlay",id,".png"), width=1920, height=1027, pointsize=12)
  par(mfrow=c(1,1))
  # Kernel density plots of posterior draws with chains separated but overlaid on a single plot.
  g <- mcmc_dens_overlay(SRMw2) + expand_limits(x=-0.1) +
    scale_x_continuous(breaks=seq(0,1,0.2), labels=seq(0,1,0.2), limits=c(0,1), expand=c(0,0)) + 
    ggtitle(paste0("Kernel Density Plots for ",w_type," for ", ff_dat_name)) + 
    theme_gray(base_size=20)
  print(g)
  dev.off()
  
  png(file=paste0("w_dens",id,".png"), width=1920, height=1027, pointsize=12)
  par(mfrow=c(1,1))
  # Kernel density plots of posterior draws with chains separated but overlaid on a single plot.
  g <- mcmc_dens(SRMw2) + expand_limits(x=-0.1) +
    scale_x_continuous(breaks=seq(0,1,0.2), labels=seq(0,1,0.2), limits=c(0,1), expand=c(0,0)) + 
    ggtitle(paste0("Kernel Density Plots for ",w_type," for ", ff_dat_name)) + 
    theme_gray(base_size=20)
  print(g)
  dev.off()
  
  
  
  ### posterior mean/median and respective credible interval: data and plot ###
  inner_prob <- 0.8; outer_prob <- 0.95; pt_est <- "median"
  png(file=paste0("w_UI_",pt_est,id,".png"), width=1000, height=1000, pointsize=12)
  par(mfrow=c(1,1))
  # function automatically combines the chains, so "SRMw_comb" gives same results
  # 'prob=0.8, prob_outer=0.95' to match w/ 'stan_plot' function arguments
  g <- mcmc_intervals(SRMw, prob=inner_prob, prob_outer=outer_prob, point_est=pt_est)
  g <- g +
    ggtitle(paste0("Credible Intervals for ",w_type," for ",ff_dat_name)) + 
    labs(subtitle=paste0("(combined the chains)\n",
                         "(point estimate = ",pt_est, ", inner probability = ",inner_prob,", outer probability = ",outer_prob,")")) +
    theme_gray(base_size=20)
  print(g)
  dev.off()
  write.csv(as.data.frame(mcmc_intervals_data(SRMw, prob=inner_prob, prob_outer=outer_prob, point_est=pt_est)),
            file=paste0("w_UI_",pt_est,id,".csv"))
}



w_and_muscleSRM_func_subset <- function(mult_chains, burn_in_prop, SRMw, SRMvalues,
                                        muscle_names_simple, muscle_names_simple_subset, 
                                        ff_dat_name,
                                        optim_method, grad_bool, upper_bool,
                                        id="") {
  if (id!="") {
    id <- mult_chains[[1]]$id
  }
  
  the_rounding <- 3
  inner_prob <- 0.8; outer_prob <- 0.95
  nchains <- length(mult_chains)
  orig_iters <- mult_chains[[1]]$iters
  burn_in_number <- orig_iters*burn_in_prop
  
  SRMvalues_combined <- c(SRMvalues)
  SRMw_combined <- merge_chains(SRMw)
  
  SRMvalues_densities <- density(SRMvalues_combined)
  SRMval_mode_theory <- SRMvalues_densities$x[which.max(SRMvalues_densities$y)]
  SRMval_mode_indx <- which.min(abs(SRMvalues_combined-SRMval_mode_theory))
  SRMval_mode_actual <- SRMvalues_combined[SRMval_mode_indx]
  w_from_postSRMmode <- SRMw_combined[SRMval_mode_indx,]
  
  
  
  subset_indxs <- which(muscle_names_simple %in% muscle_names_simple_subset)
  L_subset <- length(subset_indxs)
  
  n_etas_subset <- choose(L_subset,2)
  muscle_combns_subset <- combn(subset_indxs,2)
  eta_names_subset <- sapply(1:ncol(muscle_combns_subset), function(combn_indx)
    paste0("eta",muscle_combns_subset[1,combn_indx],muscle_combns_subset[2,combn_indx]))
  rho_names_subset <- paste0("rho_{(",subset_indxs,")}")
  r_names_subset <- c(eta_names_subset, rho_names_subset, "gamma")
  
  r_subset_indxs <- which(r_names %in% r_names_subset)
  
  
  
  initial_array1 <- array(dim=c(orig_iters-burn_in_number,nchains,L_subset),
                          dimnames=list(paste0("iter",(burn_in_number+1):orig_iters),
                                        paste0("chain",1:nchains),
                                        muscle_names_simple_subset))
  muscleSRM_post <- initial_array1
  short_mu_post <- initial_array1
  short_s_post <- initial_array1
  
  etas_post <- array(dim=c(orig_iters-burn_in_number,nchains,n_etas_subset),
                     dimnames=list(paste0("iter",(burn_in_number+1):orig_iters),
                                   paste0("chain",1:nchains),
                                   eta_names_subset))
  initial_matrix1 <- matrix(nrow=nchains,ncol=L_subset,
                            dimnames=list(paste0("chain",1:nchains),
                                          muscle_names_simple_subset))
  muscleSRM_post_avg_by_chain <- initial_matrix1
  muscleSRM_post_med_by_chain <- initial_matrix1
  
  for (c_indx in 1:nchains) {
    short_mu <- mult_chains[[c_indx]]$short_mu
    short_mu_subset <- short_mu[,subset_indxs]
    short_s <- mult_chains[[c_indx]]$short_s
    short_s_subset <- short_s[,subset_indxs]
    r_actual <- mult_chains[[c_indx]]$r_actual
    r_actual_subset <- r_actual[,r_subset_indxs]
    for (it in 1:(orig_iters-burn_in_number)) {
      short_mu_iter_subset <- short_mu_subset[burn_in_number+it,]
      short_s_iter_subset <- short_s_subset[burn_in_number+it,]
      r_iter_subset <- r_actual[burn_in_number+it,]
      etas_iter_subset <- r_iter_subset[1:n_etas_subset]
      
      muscleSRM_post[it,c_indx,] <- short_mu_iter_subset/short_s_iter_subset
      short_mu_post[it,c_indx,] <- short_mu_iter_subset
      short_s_post[it,c_indx,] <- short_s_iter_subset
      etas_post[it,c_indx,] <- etas_iter_subset
    }
    muscleSRM_post_avg_by_chain[c_indx,] <- colMeans(muscleSRM_post[,c_indx,])
    muscleSRM_post_med_by_chain[c_indx,] <- matrixStats::colMedians(muscleSRM_post[,c_indx,])
  }
  SRMval_from_postSRMmode <- c(SRMval_mode_actual, rep(NA,L_subset-1))
  muscleSRM_post_combined <- merge_chains(muscleSRM_post)
  muscleSRM_from_postSRMmode <- muscleSRM_post_combined[SRMval_mode_indx,]
  
  
  
  
  
  muscleSRM_post_avg <- colMeans(muscleSRM_post_avg_by_chain)
  muscleSRM_post_med <- matrixStats::colMedians(muscleSRM_post_med_by_chain)
  
  w_post_med <- as.numeric(as.matrix(as.data.frame(mcmc_intervals_data(
    SRMw, prob=inner_prob, prob_outer=outer_prob, point_est="median")[,"m"])))
  w_post_med_sumto1 <- abs(w_post_med)/sum(abs(w_post_med))
  w_post_mean <- as.numeric(as.matrix(as.data.frame(mcmc_intervals_data(
    SRMw, prob=inner_prob, prob_outer=outer_prob, point_est="mean")[,"m"])))
  w_post_mean_sumto1 <- abs(w_post_mean)/sum(abs(w_post_mean))
  
  short_mu_post_med <- as.numeric(as.matrix(as.data.frame(mcmc_intervals_data(
    short_mu_post, prob=inner_prob, prob_outer=outer_prob, point_est="median")[,"m"])))
  beta_tilde_post_med <- short_mu_post_med
  short_mu_post_mean <- as.numeric(as.matrix(as.data.frame(mcmc_intervals_data(
    short_mu_post, prob=inner_prob, prob_outer=outer_prob, point_est="mean")[,"m"])))
  beta_tilde_post_mean <- short_mu_post_mean
  
  short_s_post_med <- as.numeric(as.matrix(as.data.frame(mcmc_intervals_data(
    short_s_post, prob=inner_prob, prob_outer=outer_prob, point_est="median")[,"m"])))
  short_S_post_med <- diag(short_s_post_med)
  short_s_post_mean <- as.numeric(as.matrix(as.data.frame(mcmc_intervals_data(
    short_s_post, prob=inner_prob, prob_outer=outer_prob, point_est="mean")[,"m"])))
  short_S_post_mean <- diag(short_s_post_mean)
  
  etas_post_med <- as.numeric(as.matrix(as.data.frame(mcmc_intervals_data(
    etas_post, prob=inner_prob, prob_outer=outer_prob, point_est="median")[,"m"])))
  muscle_R_post_med <- diag(1,L_subset)
  muscle_R_post_med[lower.tri(muscle_R_post_med,diag=FALSE)] <-
    etas_post_med
  muscle_R_post_med <- t(muscle_R_post_med)
  muscle_R_post_med[lower.tri(muscle_R_post_med,diag=FALSE)] <- 
    t(muscle_R_post_med)[lower.tri(muscle_R_post_med)]
  etas_post_mean <- as.numeric(as.matrix(as.data.frame(mcmc_intervals_data(
    etas_post, prob=inner_prob, prob_outer=outer_prob, point_est="mean")[,"m"])))
  muscle_R_post_mean <- diag(1,L_subset)
  muscle_R_post_mean[lower.tri(muscle_R_post_mean,diag=FALSE)] <-
    etas_post_mean
  muscle_R_post_mean <- t(muscle_R_post_mean)
  muscle_R_post_mean[lower.tri(muscle_R_post_mean,diag=FALSE)] <- 
    t(muscle_R_post_mean)[lower.tri(muscle_R_post_mean)]
  
  muscle_Sigma_post_med <- 
    short_S_post_med %*% muscle_R_post_med %*% short_S_post_med
  muscle_Sigma_post_mean <- 
    short_S_post_mean %*% muscle_R_post_mean %*% short_S_post_mean
  
  w_from_post_med_stuffs <- muscle_optim_func(beta_tilde_post_med, muscle_Sigma_post_med,
                                              optim_method, grad_bool, upper_bool)
  w_from_post_med <- w_from_post_med_stuffs$w
  SRMval_from_post_med <- SRM(w_from_post_med, beta_tilde_post_med, muscle_Sigma_post_med)
  SRMval_from_post_med <- c(SRMval_from_post_med, rep(NA,L_subset-1))
  
  w_from_post_mean_stuffs <- muscle_optim_func(beta_tilde_post_mean, muscle_Sigma_post_mean,
                                               optim_method, grad_bool, upper_bool)
  w_from_post_mean <- w_from_post_mean_stuffs$w
  SRMval_from_post_mean <- SRM(w_from_post_mean, beta_tilde_post_mean, muscle_Sigma_post_mean)
  SRMval_from_post_mean <- c(SRMval_from_post_mean, rep(NA,L_subset-1))
  
  SRMval_post_med <- 
    SRM(w_post_med_sumto1, beta_tilde_post_med, muscle_Sigma_post_med)
  SRMval_post_med <- c(SRMval_post_med, rep(NA,L_subset-1))
  
  SRMval_equal_w <- 
    SRM(rep(1/L_subset,L_subset), beta_tilde_post_med, muscle_Sigma_post_med)
  SRMval_equal_w <- c(SRMval_equal_w, rep(NA,L_subset-1))
  
  
  
  
  
  y_mat_Lcol <- matrix(unlist(mult_chains[[1]]$y), ncol=length(muscle_names_simple), byrow=TRUE)
  y_mat_Lcol <- y_mat_Lcol[,subset_indxs]
  short_mu_dat <- colMeans(y_mat_Lcol, na.rm=TRUE)
  short_s_dat <- apply(y_mat_Lcol, 2, sd, na.rm=TRUE)
  muscleSRM_dat <- short_mu_dat/short_s_dat
  
  muscle_R_dat <- cor(y_mat_Lcol,use="pairwise.complete.obs")
  print(all(eigen(muscle_R_dat)$values > 0))
  
  beta_tilde_dat <- matrix(short_mu_dat, nrow=L_subset, ncol=1, byrow=TRUE)
  short_S_dat <- diag(short_s_dat)
  muscle_Sigma_dat <- short_S_dat %*% muscle_R_dat %*% short_S_dat
  
  w_dat_stuffs <- muscle_optim_func(beta_tilde_dat, muscle_Sigma_dat,
                                    optim_method, grad_bool, upper_bool)
  w_dat <- w_dat_stuffs$w
  
  SRMval_dat <- SRM(w_dat, beta_tilde_post_med, muscle_Sigma_post_med)
  SRMval_dat <- c(SRMval_dat, rep(NA,L_subset-1))
  
  
  
  
  
  SRMval_median <- median(SRMvalues_combined)
  SRMval_mean <- mean(SRMvalues_combined)
  SRMval_min <- min(SRMvalues_combined)
  SRMval_max <- max(SRMvalues_combined)

  ci_hdi <- ci(SRMvalues_combined, method="HDI")
  cred_int_lower <- ci_hdi$CI_low
  cred_int_upper <- ci_hdi$CI_high
  SRM_stats_vec <- c(SRMval_mode_actual, SRMval_median, SRMval_mean,
                     cred_int_lower, cred_int_upper,
                     SRMval_min, SRMval_max)
  SRM_stats_mat <- matrix(SRM_stats_vec, nrow=1,
                          dimnames=list("postSRM",
                                        c("mode", "median", "mean",
                                          "95%lower", "95%higher",
                                          "min", "max")))
  if (all(muscle_names_simple==muscle_names_simple_subset)) {
    write.csv(round(SRM_stats_mat,the_rounding), paste0("SRM_stats",id,".csv"))
  } else {
    write.csv(round(SRM_stats_mat,the_rounding), paste0("SRM_stats_subset",id,".csv"))
  }
  
  
  
  
  
  # added 'group' param to get facet label at top rather than x-axis label at bottom
  SRMvalues_df <- data.frame(SRMvalues_combined=SRMvalues_combined, group="SRM")
  if (all(muscle_names_simple==muscle_names_simple_subset)) {
    png(file=paste0("SRM_dens",id,".png"), width=1920, height=1027, pointsize=12)
  } else {
    png(file=paste0("SRM_dens_subset",id,".png"), width=1920, height=1027, pointsize=12)
  }
  par(mfrow=c(1,1))
  g <- ggplot(SRMvalues_df, aes(x=SRMvalues_combined)) + geom_density() +
    ggtitle(paste0("Kernel Density Plot for SRM", " for ", ff_dat_name)) +
    labs(subtitle=paste0("(combined the chains)")) +
    ylab(NULL) + xlab(NULL) + facet_grid(. ~ group) +
    theme_gray(base_size=20)
  print(g)
  dev.off()
  
  
  
  
  
  w_and_muscleSRM_list <- 
    list(w_from_post_mean, w_from_post_med, w_post_med, w_post_med_sumto1, w_dat,
         rep(NA,L_subset),
         muscleSRM_from_postSRMmode, muscleSRM_post_med, muscleSRM_dat,
         rep(NA,L_subset),
         
         SRMval_from_postSRMmode, SRMval_from_post_med, SRMval_post_med, SRMval_dat, SRMval_equal_w
    )
  w_and_muscleSRM <- 
    matrix(unlist(w_and_muscleSRM_list), nrow=length(w_and_muscleSRM_list),
           dimnames=list(c(
             "*1. weights (posterior [w] from posterior SRM mode)",
             "**2a weights (target [w] using median of posterior [mu], [s], [r])",
             "2b. weights (median of posterior [w])",
             "2b. weights (median of posterior [w] sumto1)",
             "3. weights (data [w])",
             
             "",
             
             "1. muscle SRMs (posterior [mu]/[s] from posterior SRM mode)",
             "**2. muscle SRMs (median of posterior [mu]/[s])",
             "3. muscle SRMs (data [mu]/[s])",
             
             "",
             
             "*1. SRM (posterior [w] from posterior SRM mode)",
             "**2a SRM (target [w] using median of posterior [mu], [s], [r])",
             "2b. SRM (median of posterior [w] sumto1)",
             "3. SRM (data [w])",
             "4. SRM (equal weights)"
           ),
           muscle_names_simple_subset),
           byrow=TRUE)
  
  if (all(muscle_names_simple==muscle_names_simple_subset)) {
    write.csv(round(w_and_muscleSRM,the_rounding),paste0("w_and_muscleSRM",id,".csv"))
    write.csv(w_and_muscleSRM,paste0("!w_and_muscleSRM",id,".csv"))
  } else {
    write.csv(round(w_and_muscleSRM,the_rounding),paste0("w_and_muscleSRM_subset",id,".csv"))
    write.csv(w_and_muscleSRM,paste0("!w_and_muscleSRM_subset",id,".csv"))
  }
  
  
  
  
  
  return(list(w_and_muscleSRM=w_and_muscleSRM))
}



mixture_diagnostic <- function(mult_chains,chunk_length) {
  orig_iters <- mult_chains[[1]]$iters
  burn_in_number <- orig_iters*burn_in_prop
  
  main_iters <- (burn_in_number+1):orig_iters
  chunk_indxs_list <- split(main_iters, ceiling(seq_along(main_iters)/chunk_length))
  n_chunks <- length(chunk_indxs_list)
  chunk_indxs_names <- sapply(1:n_chunks, function(chunk) paste0(chunk_indxs_list[[chunk]][1],":",chunk_indxs_list[[chunk]][chunk_length]))
  
  mu_names_detailed <- mult_chains[[1]]$mu_names_detailed
  s_names_detailed <- mult_chains[[1]]$s_names_detailed
  r_names_detailed <- mult_chains[[1]]$r_names_detailed
  L <- mult_chains[[1]]$L
  Q <- mult_chains[[1]]$Q
  stats_names <- c("mean","0.025","0.25","0.5","0.75","0.975")
  quantile_probs <- c(0.025,0.25,0.5,0.75,0.975)
  n_stats <- length(stats_names)
  
  short_mu_stats_all_chunks_all_chains <- vector("list",length=nchains)
  short_s_stats_all_chunks_all_chains <- vector("list",length=nchains)
  r_actual_stats_all_chunks_all_chains <- vector("list",length=nchains)
  for (c_indx in 1:nchains) {
    short_mu <- mult_chains[[c_indx]]$short_mu
    short_s <- mult_chains[[c_indx]]$short_s
    r_actual <- mult_chains[[c_indx]]$r_actual
    
    short_mu_stats_all_chunks <- lapply(1:L, function(l) matrix(nrow=n_chunks,ncol=n_stats,dimnames=list(chunk_indxs_names,stats_names)))
    names(short_mu_stats_all_chunks) <- mu_names_detailed
    short_s_stats_all_chunks <- lapply(1:L, function(l) matrix(nrow=n_chunks,ncol=n_stats,dimnames=list(chunk_indxs_names,stats_names)))
    names(short_s_stats_all_chunks) <- s_names_detailed
    r_actual_stats_all_chunks <- lapply(1:Q, function(k) matrix(nrow=n_chunks,ncol=n_stats,dimnames=list(chunk_indxs_names,stats_names)))
    names(r_actual_stats_all_chunks) <- r_names_detailed
    
    for (chunk in 1:n_chunks) {
      chunk_indxs <- chunk_indxs_list[[chunk]]
      
      short_mu_chunk <- short_mu[chunk_indxs,]
      short_mu_avgs_chunk <- colMeans(short_mu_chunk)
      short_mu_quantiles_chunk <- apply(short_mu_chunk,2,quantile,probs=quantile_probs)
      short_mu_stats_chunk <- rbind(short_mu_avgs_chunk,short_mu_quantiles_chunk)
      
      short_s_chunk <- short_s[chunk_indxs,]
      short_s_avgs_chunk <- colMeans(short_s_chunk)
      short_s_quantiles_chunk <- apply(short_s_chunk,2,quantile,probs=quantile_probs)
      short_s_stats_chunk <- rbind(short_s_avgs_chunk, short_s_quantiles_chunk)
      
      r_actual_chunk <- r_actual[chunk_indxs,]
      r_actual_avgs_chunk <- colMeans(r_actual_chunk)
      r_actual_quantiles_chunk <- apply(r_actual_chunk,2,quantile,probs=quantile_probs)
      r_actual_stats_chunk <- rbind(r_actual_avgs_chunk, r_actual_quantiles_chunk)
      for (l in 1:L) {
        short_mu_stats_all_chunks[[l]][chunk,] <- short_mu_stats_chunk[,l]
        short_s_stats_all_chunks[[l]][chunk,] <- short_s_stats_chunk[,l]
      }
      for (k in 1:Q) {
        r_actual_stats_all_chunks[[k]][chunk,] <- r_actual_stats_chunk[,k]
      }
    }
    short_mu_stats_all_chunks_all_chains[[c_indx]] <- short_mu_stats_all_chunks
    short_s_stats_all_chunks_all_chains[[c_indx]] <- short_s_stats_all_chunks
    r_actual_stats_all_chunks_all_chains[[c_indx]] <- r_actual_stats_all_chunks
  }
  
  
  
  short_mu_avgs_all_chains <- vector("list",length=L)
  names(short_mu_avgs_all_chains) <- mu_names_detailed
  short_s_avgs_all_chains <- vector("list",length=L)
  names(short_s_avgs_all_chains) <- s_names_detailed
  r_actual_avgs_all_chains <- vector("list",length=Q)
  names(r_actual_avgs_all_chains) <- r_names_detailed
  chain_names1 <-  c(paste0("chain",1:nchains," (mean across chunks of iterations)"), "range")
  for (l in 1:L) {
    mu_l_avgs_all_chains <- matrix(nrow=nchains+1,ncol=n_stats,dimnames=list(chain_names1,stats_names))
    s_l_avgs_all_chains <- matrix(nrow=nchains+1,ncol=n_stats,dimnames=list(chain_names1,stats_names))
    for (c_indx in 1:nchains) {
      mu_l_avgs_all_chains[c_indx,] <- colMeans(short_mu_stats_all_chunks_all_chains[[c_indx]][[l]])
      s_l_avgs_all_chains[c_indx,] <- colMeans(short_s_stats_all_chunks_all_chains[[c_indx]][[l]])
    }
    mu_l_avgs_range <- apply(mu_l_avgs_all_chains[1:nchains,],2,range)
    mu_l_avgs_all_chains[c_indx+1,] <- mu_l_avgs_range[2,]-mu_l_avgs_range[1,]
    s_l_avgs_range <- apply(s_l_avgs_all_chains[1:nchains,],2,range)
    s_l_avgs_all_chains[c_indx+1,] <- s_l_avgs_range[2,]-s_l_avgs_range[1,]
    
    short_mu_avgs_all_chains[[l]] <- mu_l_avgs_all_chains
    short_s_avgs_all_chains[[l]] <- s_l_avgs_all_chains
  }
  for (k in 1:Q) {
    r_k_avgs_all_chains <- matrix(nrow=nchains+1,ncol=n_stats,dimnames=list(chain_names1,stats_names))
    for (c_indx in 1:nchains) {
      r_k_avgs_all_chains[c_indx,] <- colMeans(r_actual_stats_all_chunks_all_chains[[c_indx]][[k]])
    }
    r_k_avgs_range <- apply(r_k_avgs_all_chains[1:nchains,],2,range)
    r_k_avgs_all_chains[c_indx+1,] <- r_k_avgs_range[2,]-r_k_avgs_range[1,]
    
    r_actual_avgs_all_chains[[k]] <- r_k_avgs_all_chains
  }
  
  
  
  short_mu_sds_all_chains <- vector("list",length=L)
  names(short_mu_sds_all_chains) <- mu_names_detailed
  short_s_sds_all_chains <- vector("list",length=L)
  names(short_s_sds_all_chains) <- s_names_detailed
  r_actual_sds_all_chains <- vector("list",length=Q)
  names(r_actual_sds_all_chains) <- r_names_detailed
  chain_names2 <-  c(paste0("chain",1:nchains," (sd across chunks of iterations)"), "range")
  for (l in 1:L) {
    mu_l_sds_all_chains <- matrix(nrow=nchains+1,ncol=n_stats,dimnames=list(chain_names2,stats_names))
    s_l_sds_all_chains <- matrix(nrow=nchains+1,ncol=n_stats,dimnames=list(chain_names2,stats_names))
    for (c_indx in 1:nchains) {
      mu_l_sds_all_chains[c_indx,] <- apply(short_mu_stats_all_chunks_all_chains[[c_indx]][[l]], 2, sd)
      s_l_sds_all_chains[c_indx,] <- apply(short_s_stats_all_chunks_all_chains[[c_indx]][[l]], 2, sd)
    }
    mu_l_sds_range <- apply(mu_l_sds_all_chains[1:nchains,],2,range)
    mu_l_sds_all_chains[c_indx+1,] <- mu_l_sds_range[2,]-mu_l_sds_range[1,]
    s_l_sds_range <- apply(s_l_sds_all_chains[1:nchains,],2,range)
    s_l_sds_all_chains[c_indx+1,] <- s_l_sds_range[2,]-s_l_sds_range[1,]
    
    short_mu_sds_all_chains[[l]] <- mu_l_sds_all_chains
    short_s_sds_all_chains[[l]] <- s_l_sds_all_chains
  }
  for (k in 1:Q) {
    r_k_sds_all_chains <- matrix(nrow=nchains+1,ncol=n_stats,dimnames=list(chain_names2,stats_names))
    for (c_indx in 1:nchains) {
      r_k_sds_all_chains[c_indx,] <- apply(r_actual_stats_all_chunks_all_chains[[c_indx]][[k]], 2, sd)
    }
    r_k_sds_range <- apply(r_k_sds_all_chains[1:nchains,],2,range)
    r_k_sds_all_chains[c_indx+1,] <- r_k_sds_range[2,]-r_k_sds_range[1,]
    
    r_actual_sds_all_chains[[k]] <- r_k_sds_all_chains
  }
  
  return(list(short_mu_stats_all_chunks_all_chains=short_mu_stats_all_chunks_all_chains, # just to look at
              short_s_stats_all_chunks_all_chains=short_s_stats_all_chunks_all_chains, # just to look at 
              r_actual_stats_all_chunks_all_chains=r_actual_stats_all_chunks_all_chains, # just to look at
              
              short_mu_avgs_all_chains=short_mu_avgs_all_chains, # to create csv files
              short_s_avgs_all_chains=short_s_avgs_all_chains, # to create csv files
              r_actual_avgs_all_chains=r_actual_avgs_all_chains, # to create csv files
              
              short_mu_sds_all_chains=short_mu_sds_all_chains, # to create csv files
              short_s_sds_all_chains=short_s_sds_all_chains, # to create csv files
              r_actual_sds_all_chains=r_actual_sds_all_chains # to create csv files
              )
         )
}





missing_perc1_func <- function(y_mat_Lcols_amb) { # for: balanced1_sim_func, unbalanced1_sim_func
  N.J <- nrow(y_mat_Lcols_amb)
  L <- ncol(y_mat_Lcols_amb)
  
  missing_perc_by_muscle <- colSums(is.na(y_mat_Lcols_amb))/N.J # vector
  
  mat_NAs_for_each_row <- as.matrix(as.data.frame(
    table(rowSums(is.na(y_mat_Lcols_amb)))))
  class(mat_NAs_for_each_row) <- "numeric"
  for (l_indx in 0:L) {
    if (!l_indx %in% mat_NAs_for_each_row[,"Var1"]) {
      mat_NAs_for_each_row <- rbind(mat_NAs_for_each_row, c(l_indx,0))
    }
  }
  mat_NAs_for_each_row <- mat_NAs_for_each_row[
    order(mat_NAs_for_each_row[,"Var1"]),]
  mat_NAs_for_each_row <- cbind(mat_NAs_for_each_row,
                                mat_NAs_for_each_row[,"Freq"]/N.J)
  colnames(mat_NAs_for_each_row) <- c("NAs_for_each_row","Freq","%")
  missing_perc_by_rowNAs <- mat_NAs_for_each_row[
    mat_NAs_for_each_row[,"NAs_for_each_row"]!=0,"%"] # vector
  
  return(list(missing_perc_by_muscle=missing_perc_by_muscle,
              missing_perc_by_rowNAs=missing_perc_by_rowNAs))
}

missing_perc2_func <- function(y_list, L) { # for: unbalanced2_sim_func
  N <- length(y_list) # N == nrow(y_sim_full) == length(y_sim_full_unbal) == length(y_sim_wNA_unbal)
  ps <- sapply(1:N, function(i) length(y_list[[i]]))
  Js <- ps/L
  p_max <- max(ps)
  J_max <- max(Js)

  y_mat_Lcols_amb <- matrix(unlist(y_list), ncol=L, byrow=TRUE)
  
  unique_Js <- sort(unique(Js))
  
  missing_perc_by_muscle <- matrix(nrow=J_max, ncol=L)
  missing_perc_by_rowNAs <- matrix(nrow=J_max, ncol=L)
  for (j_indx in unique_Js) {
    p_indx <- j_indx*L
    y_list_temp <- Filter(function(x) length(x)==p_indx,
                          y_list)
    y_mat_Lcols_temp <- matrix(unlist(y_list_temp), ncol=L, byrow=TRUE)
    N.J_temp <- nrow(y_mat_Lcols_temp)
    
    missing_perc_by_muscle[j_indx,] <- colSums(is.na(y_mat_Lcols_temp)/N.J_temp)
    
    mat_NAs_for_each_row_temp <- as.matrix(as.data.frame(
      table(rowSums(is.na(y_mat_Lcols_temp)))))
    class(mat_NAs_for_each_row_temp) <- "numeric"
    for (l_indx in 0:L) {
      if (!l_indx %in% mat_NAs_for_each_row_temp[,"Var1"]) {
        mat_NAs_for_each_row_temp <- rbind(mat_NAs_for_each_row_temp, c(l_indx,0))
      }
    }
    mat_NAs_for_each_row_temp <- mat_NAs_for_each_row_temp[
      order(mat_NAs_for_each_row_temp[,"Var1"]),]
    mat_NAs_for_each_row_temp <- cbind(mat_NAs_for_each_row_temp,
                                       mat_NAs_for_each_row_temp[,"Freq"]/N.J_temp)
    colnames(mat_NAs_for_each_row_temp) <- c("NAs_for_each_row","Freq","%")
    missing_perc_by_rowNAs_temp <- mat_NAs_for_each_row_temp[
      mat_NAs_for_each_row_temp[,"NAs_for_each_row"]!=0,"%"]
    missing_perc_by_rowNAs[j_indx,] <- missing_perc_by_rowNAs_temp
  }
  
  return(list(missing_perc_by_muscle=missing_perc_by_muscle,
              missing_perc_by_rowNAs=missing_perc_by_rowNAs))
}

Js_dist_func <- function(Js, N) { # for: unbalanced1_sim_func, unbalanced2_sim_func
  mat_Js_dist <- as.matrix(as.data.frame(table(Js)))
  class(mat_Js_dist) <- "numeric"
  mat_Js_dist <- cbind(mat_Js_dist, mat_Js_dist[,"Freq"]/N)
  colnames(mat_Js_dist) <- c("Js", "Freq", "%")
  perc_Js_dist <- mat_Js_dist[,"%"]
  
  return(list(mat_Js_dist=mat_Js_dist,
              perc_Js_dist=perc_Js_dist))
}

MSE_func <- function(dat_sim_full, dat_sim_wNA,
                     target_short_mu, target_short_s, target_r,
                     target_w, target_SRM) { # for: balanced1_sim_func, unbalanced1_sim_func, unbalanced2_sim_func
  L <- ncol(dat_sim_full)
  
  #target_short_mu
  sim_full_short_mu <- colMeans(dat_sim_full, na.rm=TRUE)
  sim_wNA_short_mu <- colMeans(dat_sim_wNA, na.rm=TRUE)
  MSE_sim_full_short_mu <- mean((target_short_mu - sim_full_short_mu)^2)
  MSE_sim_wNA_short_mu <- mean((target_short_mu - sim_wNA_short_mu)^2)
  
  #target_short_s
  sim_full_short_s <- apply(dat_sim_full, 2, sd, na.rm=TRUE)
  sim_wNA_short_s <- apply(dat_sim_wNA, 2, sd, na.rm=TRUE)
  MSE_sim_full_short_s <- mean((target_short_s - sim_full_short_s)^2)
  MSE_sim_wNA_short_s <- mean((target_short_s - sim_wNA_short_s)^2)
  
  target_etas <- target_r[1:choose(L,2)]
  sim_full_etas_mat <- cor(dat_sim_full)
  sim_full_etas <- unique(c(t(sim_full_etas_mat)))[-1]
  sim_wNA_etas_mat <- cor(dat_sim_wNA,use="complete.obs")
  sim_wNA_etas <- unique(c(t(sim_wNA_etas_mat)))[-1]
  MSE_sim_full_etas <- mean((target_etas - sim_full_etas)^2)
  MSE_sim_wNA_etas <- mean((target_etas - sim_wNA_etas)^2)
  
  sim_full_muscle_Sigma <- diag(sim_full_short_s)%*%sim_full_etas_mat%*%diag(sim_full_short_s)
  sim_wNA_muscle_Sigma <- diag(sim_wNA_short_s)%*%sim_wNA_etas_mat%*%diag(sim_wNA_short_s)
  
  # optim_method="L-BFGS-B", grad_bool==TRUE, upper_bool==TRUE:
  sim_full_w_stuffs <- muscle_optim_func(sim_full_short_mu, sim_full_muscle_Sigma)
  sim_full_w <- sim_full_w_stuffs$w
  MSE_sim_full_w <- mean((target_w - sim_full_w)^2)
  
  sim_full_SRM <- SRM(sim_full_w, sim_full_short_mu, sim_full_muscle_Sigma)
  MSE_sim_full_SRM <- mean((target_SRM - sim_full_SRM)^2)
  
  # optim_method="L-BFGS-B", grad_bool==TRUE, upper_bool==TRUE:
  sim_wNA_w_stuffs <- muscle_optim_func(sim_wNA_short_mu, sim_wNA_muscle_Sigma)
  sim_wNA_w <- sim_wNA_w_stuffs$w
  MSE_sim_wNA_w <- mean((target_w - sim_wNA_w)^2)
  
  sim_wNA_SRM <- SRM(sim_wNA_w, sim_wNA_short_mu, sim_wNA_muscle_Sigma)
  MSE_sim_wNA_SRM <- mean((target_SRM - sim_wNA_SRM)^2)
  
  MSEs <- c(MSE_sim_full_short_mu, MSE_sim_wNA_short_mu,
            MSE_sim_full_short_s, MSE_sim_wNA_short_s,
            MSE_sim_full_etas, MSE_sim_wNA_etas,
            MSE_sim_full_w, MSE_sim_wNA_w,
            MSE_sim_full_SRM, MSE_sim_wNA_SRM)
  avg_MSEs <- mean(MSEs)
  the_MSEs <- c(MSEs, avg_MSEs)
  
  return(the_MSEs)
}

balanced1_sim_func <- function(y_sim_full,
                               test_bool, N_sim, L,
                               missing_perc_by_muscle, missing_perc_by_rowNAs, # vecs
                               J_max,
                               amb_indx=NULL,
                               target_short_mu_list=NULL, target_short_s_list=NULL, target_r_list=NULL,
                               target_w_list=NULL, target_SRM_list=NULL) {
  
  p_max <- J_max*L
  
  target_short_mu <- target_short_mu_list[[amb_indx]]
  target_short_s <- target_short_s_list[[amb_indx]]
  target_r <- target_r_list[[amb_indx]]
  target_w <- target_w_list[[amb_indx]]
  target_SRM <- target_SRM_list[[amb_indx]]
  
  N.J_sim <- N_sim*J_max
  dat_sim_full <- matrix(t(y_sim_full), nrow=N.J_sim, ncol=L, byrow=TRUE)
  
  dat_sim_row_indxs_temp <- 1:(N.J_sim)
  missing_row_indxs_list <- vector("list",length=L)
  for (l in 1:L) {
    # this will produce NULL in 'missing_row_indxs_list[[l]]' if FALSE
    if (missing_perc_by_rowNAs[l]!=0) {
      if (length(dat_sim_row_indxs_temp)==1) {
        # reason: run 'sample(x=5,size=1)' bunch of times
        missing_row_indxs_list[[l]] <- dat_sim_row_indxs_temp
      } else {
        missing_row_indxs_list[[l]] <-
          sample(dat_sim_row_indxs_temp,
                 size=round(missing_perc_by_rowNAs[l]*N.J_sim,0))
        dat_sim_row_indxs_temp <- dat_sim_row_indxs_temp[
          -which(dat_sim_row_indxs_temp %in% missing_row_indxs_list[[l]])]
      }
    }
  }
  
  dat_sim_wNA <- dat_sim_full
  for (l in 1:L) {
    # need this conditional because
    # 'size' has to be <= non-zero probabilities in 'prob' for the 'sample' function
    if (missing_perc_by_rowNAs[l]!=0) {
      if (l < L) {
        for (indx in 1:length(missing_row_indxs_list[[l]])) {
          col_indxs_temp <- sample(1:L,size=l,prob=missing_perc_by_muscle)
          row_indx_temp <- missing_row_indxs_list[[l]][indx]
          dat_sim_wNA[row_indx_temp,col_indxs_temp] <- NA
        }
      } else {
        dat_sim_wNA[missing_row_indxs_list[[l]],] <- NA
      }
    }
  }
  y_sim_wNA_mat <- matrix(t(dat_sim_wNA), nrow=N_sim, ncol=J_max*L, byrow=TRUE)
  rows_w_no_obs_vars <- which(rowSums(is.na(y_sim_wNA_mat)) == p_max)
  for (indx in rows_w_no_obs_vars) {
    y_sim_wNA_mat[indx,] <- y_sim_full[indx,]
  }
  y_sim_wNA_list <- as.list(data.frame(t(y_sim_wNA_mat))) # convert to list form
  y_sim_wNA_list <- unname(y_sim_wNA_list) # unname 'X1', 'X2', 'X3', ...
  
  
  
  if (test_bool==TRUE) {
    ### Evaluate simulated data with missingness
    # Note that we used posterior median of mu, s, r (compare this - overall)
    # and used missingness of actual data (compare this - by Js)
    missing_perc_vecs_sim <- missing_perc1_func(dat_sim_wNA)
    missing_perc_by_muscle_sim <-  missing_perc_vecs_sim$missing_perc_by_muscle
    missing_perc_by_rowNAs_sim <- missing_perc_vecs_sim$missing_perc_by_rowNAs
    
    the_MSEs <- MSE_func(dat_sim_full, dat_sim_wNA,
                         target_short_mu, target_short_s, target_r,
                         target_w, target_SRM)
  }
  
  
  
  if (test_bool==TRUE) {
    return(list(missing_perc_by_muscle=missing_perc_by_muscle,
                missing_perc_by_muscle_sim=missing_perc_by_muscle_sim,
                missing_perc_by_rowNAs=missing_perc_by_rowNAs,
                missing_perc_by_rowNAs_sim=missing_perc_by_rowNAs_sim,
                the_MSEs=the_MSEs,
                
                ps_sim=rep(p_max,N_sim),
                dat_sim_full=dat_sim_full,
                dat_sim_wNA=dat_sim_wNA
    )
    )
  } else if (test_bool==FALSE) {
    return(list(
                y_sim_wNA_mat=y_sim_wNA_mat,
                y_sim_wNA_list=y_sim_wNA_list
    )
    )
  }
  
  
  
  ### SIMPLE MISSINGNESS
  ###missing_perc_by_y_col <- rep(0.2, p_max)
  
  ###y_sim_wNA_simple <- y_sim_full
  ###for (m in 1:p_max) {
  ###row_indxs_temp_simple <- 
  ###sample(1:N_sim, 
  ###size=round(missing_perc_by_y_col[m]*N_sim, 0))
  ###y_sim_wNA_simple[row_indxs_temp_simple,m] <- NA
  ###}
  ###dat_sim_wNA_simple <- matrix(t(y_sim_wNA_simple),
  ###nrow=N.J_sim, ncol=L, byrow=TRUE)
}

unbalanced1_sim_func <- function(y_sim_full,
                                 test_bool, N_sim, L,
                                 missing_perc_by_muscle, missing_perc_by_rowNAs, # vecs
                                 mat_Js_dist, perc_Js_dist,
                                 amb_indx=NULL,
                                 target_short_mu_list=NULL, target_short_s_list=NULL, target_r_list=NULL,
                                 target_w_list=NULL, target_SRM_list=NULL) {
  
  J_max <- max(mat_Js_dist[,"Js"]) # 'J_max' changes based on 'mat_Js_dist' input
  p_max <- J_max*L
  unique_Js <- mat_Js_dist[,"Js"]
  
  target_short_mu <- target_short_mu_list[[amb_indx]]
  target_short_s <- target_short_s_list[[amb_indx]]
  target_r <- target_r_list[[amb_indx]]
  target_w <- target_w_list[[amb_indx]]
  target_SRM <- target_SRM_list[[amb_indx]]
  
  subject_indxs <- 1:N_sim
  subjects_for_each_Js_list <- vector("list", length=J_max)
  for (j_indx in unique_Js) {
    if (length(subject_indxs)==1) {
      # reason: run 'sample(x=5,size=1)' bunch of times
      subjects_for_each_Js_list[[j_indx]] <- subject_indxs
    } else {
      subjects_for_each_Js_list[[j_indx]] <-
        sample(subject_indxs,
               size=round(mat_Js_dist[mat_Js_dist[,"Js"]==j_indx,"%"]*N_sim,0))
      subject_indxs <- subject_indxs[-which(subject_indxs %in%
                                              subjects_for_each_Js_list[[j_indx]])]
    }
  }
  missing_subjects <- 1:N_sim
  missing_subjects <- missing_subjects[
    !missing_subjects %in% unlist(subjects_for_each_Js_list)]
  
  y_sim_full_unbal <- vector("list", length=N_sim)
  muscle_group_indxs <- matrix(1:p_max,ncol=L,byrow=TRUE)
  for (j_indx in unique_Js) {
    for (indx in 1:length(subjects_for_each_Js_list[[j_indx]])) {
      i_indx <- subjects_for_each_Js_list[[j_indx]][indx]
      y_sim_full_row_temp <- y_sim_full[i_indx,]
      if (j_indx < J_max) {
        # you sample on '1:J_max' and not 'unique_Js', b/c
        # these are meant for 'muscle_group_indxs'
        remove_time_point_indxs_temp <- sample(1:J_max, size=J_max-j_indx)
        # technically don't need to sort
        remove_muscle_group_indxs_temp <- sort(c(muscle_group_indxs[remove_time_point_indxs_temp,]))
        y_sim_full_unbal_row_temp <- y_sim_full_row_temp[-remove_muscle_group_indxs_temp]
        y_sim_full_unbal[[i_indx]] <- y_sim_full_unbal_row_temp
      } else {
        y_sim_full_unbal[[i_indx]] <- y_sim_full_row_temp
      }
    }
  }
  if (length(missing_subjects)!=0) {
    for (i_indx in missing_subjects) {
      y_sim_full_unbal[[i_indx]] <- y_sim_full[i_indx,]
    }
  }
  
  # to later convert 'dat_sim_wNA_unbal' back to 'y_sim_wNA_unbal'
  ps_sim <- sapply(1:N_sim, function(i) length(y_sim_full_unbal[[i]]))
  #table(ps_sim); table(ps) # to test 'y_sim_full_unbal' has same 'ps' dist
  Js_sim <- ps_sim/L
  N.J_sim <- sum(Js_sim) # nrow(dat_sim_full_unbal) == nrow(dat_sim_wNA_unbal)
  
  # same idea as y_mat_Lcols
  dat_sim_full_unbal <- matrix(unlist(y_sim_full_unbal), ncol=L, byrow=TRUE)
  
  dat_sim_row_indxs_temp <- 1:N.J_sim
  missing_row_indxs_list <- vector("list",length=L)
  for (l in 1:L) {
    # this will produce NULL in 'missing_row_indxs_list[[l]]' if FALSE
    if (missing_perc_by_rowNAs[l]!=0) {
      if (length(dat_sim_row_indxs_temp)==1) {
        # reason: run 'sample(x=5,size=1)' bunch of times
        missing_row_indxs_list[[l]] <- dat_sim_row_indxs_temp
      } else {
        missing_row_indxs_list[[l]] <- 
          sample(dat_sim_row_indxs_temp,
                 size=round(missing_perc_by_rowNAs[l]*N.J_sim,0))
        dat_sim_row_indxs_temp <- dat_sim_row_indxs_temp[
          -which(dat_sim_row_indxs_temp %in% missing_row_indxs_list[[l]])]
      }
    }
  }
  
  dat_sim_wNA_unbal <- dat_sim_full_unbal
  for (l in 1:L) {
    # need this conditional because
    # 'size' has to be <= non-zero probabilities in 'prob' for the 'sample' function
    if (missing_perc_by_rowNAs[l]!=0) {
      if (l < L) {
        for (indx in 1:length(missing_row_indxs_list[[l]])) {
          col_indxs_temp <- sample(1:L,size=l,prob=missing_perc_by_muscle)
          row_indx_temp <- missing_row_indxs_list[[l]][indx]
          dat_sim_wNA_unbal[row_indx_temp,col_indxs_temp] <- NA
        }
      } else if (l==L) {
        dat_sim_wNA_unbal[missing_row_indxs_list[[l]],] <- NA
      }
    }
  }
  
  y_sim_wNA_unbal <- vector("list", length=N_sim)
  vec_sim_wNA_unbal <- c(t(dat_sim_wNA_unbal))
  for (i in 1:N_sim) {
    num_outcomes_subject_i <- ps_sim[i]
    y_sim_wNA_unbal[[i]] <- vec_sim_wNA_unbal[1:num_outcomes_subject_i]
    vec_sim_wNA_unbal <- vec_sim_wNA_unbal[-c(1:num_outcomes_subject_i)]
  }
  
  
  
  if (test_bool==TRUE) {
    ### Evaluate simulated data with missingness
    # Note that we used posterior median of mu, s, r (compare this - overall)
    # and used missingness of actual data (compare this - by Js)
    missing_perc_vecs_sim <- missing_perc1_func(dat_sim_wNA_unbal)
    missing_perc_by_muscle_sim <-  missing_perc_vecs_sim$missing_perc_by_muscle
    missing_perc_by_rowNAs_sim <- missing_perc_vecs_sim$missing_perc_by_rowNAs
    
    Js_dist_stuff_sim <- Js_dist_func(Js_sim, N_sim)
    mat_Js_dist_sim <- Js_dist_stuff_sim$mat_Js_dist
    perc_Js_dist_sim <- Js_dist_stuff_sim$perc_Js_dist
    
    the_MSEs <- MSE_func(dat_sim_full_unbal, dat_sim_wNA_unbal,
                         target_short_mu, target_short_s, target_r,
                         target_w, target_SRM)
  }
  
  
  
  if (test_bool==TRUE) {
    return(list(missing_perc_by_muscle=missing_perc_by_muscle,
                missing_perc_by_muscle_sim=missing_perc_by_muscle_sim,
                missing_perc_by_rowNAs=missing_perc_by_rowNAs,
                missing_perc_by_rowNAs_sim=missing_perc_by_rowNAs_sim,
                perc_Js_dist=perc_Js_dist,
                perc_Js_dist_sim=perc_Js_dist_sim,
                the_MSEs=the_MSEs,
                
                ps_sim=ps_sim,
                dat_sim_full=dat_sim_full_unbal, # 'dat_sim_full' b/c want same naming as other functions
                dat_sim_wNA=dat_sim_wNA_unbal # 'dat_sim_wNA' b/c want same naming as other functions
    )
    )
  } else if (test_bool==FALSE) {
    return(list(y_sim_full_unbal=y_sim_full_unbal,
                y_sim_wNA_list=y_sim_wNA_unbal # 'y_sim_wNA_list' b/c want same naming as other functions
    )
    )
  }
}

unbalanced2_sim_func <- function(y_sim_full,
                                 test_bool, N_sim, L,
                                 missing_perc_by_muscle, missing_perc_by_rowNAs, # mats
                                 mat_Js_dist, perc_Js_dist,
                                 amb_indx=NULL,
                                 target_short_mu_list=NULL, target_short_s_list=NULL, target_r_list=NULL,
                                 target_w_list=NULL, target_SRM_list=NULL) {
  
  J_max <- max(mat_Js_dist[,"Js"]) # 'J_max' changes based on 'mat_Js_dist' input
  p_max <- J_max*L
  unique_Js <- mat_Js_dist[,"Js"]
  
  target_short_mu <- target_short_mu_list[[amb_indx]]
  target_short_s <- target_short_s_list[[amb_indx]]
  target_r <- target_r_list[[amb_indx]]
  target_w <- target_w_list[[amb_indx]]
  target_SRM <- target_SRM_list[[amb_indx]]
  
  subject_indxs <- 1:N_sim
  subjects_for_each_Js_list <- vector("list", length=J_max)
  for (j_indx in unique_Js) {
    if (length(subject_indxs)==1) {
      # reason: run 'sample(x=5,size=1)' bunch of times
      subjects_for_each_Js_list[[j_indx]] <- subject_indxs
    } else {
      subjects_for_each_Js_list[[j_indx]] <-
        sample(subject_indxs,
               size=round(mat_Js_dist[mat_Js_dist[,"Js"]==j_indx,"%"]*N_sim,0))
      subject_indxs <- subject_indxs[-which(subject_indxs %in%
                                              subjects_for_each_Js_list[[j_indx]])]
    }
  }
  missing_subjects <- 1:N_sim
  missing_subjects <- missing_subjects[
    !missing_subjects %in% unlist(subjects_for_each_Js_list)]
  
  y_sim_full_unbal <- vector("list", length=N_sim)
  muscle_group_indxs <- matrix(1:p_max,ncol=L,byrow=TRUE)
  for (j_indx in unique_Js) {
    for (indx in 1:length(subjects_for_each_Js_list[[j_indx]])) {
      i_indx <- subjects_for_each_Js_list[[j_indx]][indx]
      y_sim_full_row_temp <- y_sim_full[i_indx,]
      if (j_indx < J_max) {
        # you sample on '1:J_max' and not 'unique_Js', b/c
        # these are meant for 'muscle_group_indxs'
        remove_time_point_indxs_temp <- sample(1:J_max, size=J_max-j_indx)
        # technically don't need to sort
        remove_muscle_group_indxs_temp <- sort(c(muscle_group_indxs[remove_time_point_indxs_temp,]))
        y_sim_full_unbal_row_temp <- y_sim_full_row_temp[-remove_muscle_group_indxs_temp]
        y_sim_full_unbal[[i_indx]] <- y_sim_full_unbal_row_temp
      } else {
        y_sim_full_unbal[[i_indx]] <- y_sim_full_row_temp
      }
    }
  }
  if (length(missing_subjects)!=0) {
    for (i_indx in missing_subjects) {
      y_sim_full_unbal[[i_indx]] <- y_sim_full[i_indx,]
    }
  }
  
  # reorder 'y_sim_full_unbal' so it matches 'y_sim_wNA_unbal'
  y_sim_full_unbal2 <- list()
  for (j_indx in unique_Js) {
    p_indx <- j_indx*L
    y_sim_full_unbal_temp <- Filter(function(x) length(x)==p_indx,
                                    y_sim_full_unbal)
    # I think keep 'append' for now (unless better way later) b/c 'y_sim_full_unbal_temp' has varying lengths
    y_sim_full_unbal2 <- append(y_sim_full_unbal2, y_sim_full_unbal_temp)
  }
  
  ps_sim <- sapply(1:N_sim, function(i) length(y_sim_full_unbal2[[i]]))
  #table(ps_sim); table(ps) # to test 'y_sim_full_unbal' has same 'ps' dist
  Js_sim <- ps_sim/L
  N.J_sim <- sum(Js_sim) # nrow(dat_sim_full_unbal) == nrow(dat_sim_wNA_unbal); not needed here I think
  
  # same idea as y_mat_Lcols
  dat_sim_full_unbal <- matrix(unlist(y_sim_full_unbal2), ncol=L, byrow=TRUE)
  
  y_sim_wNA_unbal <- list()
  for (j_indx in unique_Js) {
    p_indx <- j_indx*L
    y_sim_full_unbal_temp <- Filter(function(x) length(x)==p_indx,
                                    y_sim_full_unbal2)
    dat_sim_full_unbal_temp <- matrix(unlist(y_sim_full_unbal_temp), ncol=L, byrow=TRUE)
    N.J_sim_temp <- nrow(dat_sim_full_unbal_temp)
    
    missing_perc_by_muscle_temp <- missing_perc_by_muscle[j_indx,]
    missing_perc_by_rowNAs_temp <- missing_perc_by_rowNAs[j_indx,]
    
    dat_sim_row_indxs_temp <- 1:N.J_sim_temp
    missing_row_indxs_list <- vector("list",length=L)
    for (l in 1:L) {
      # this will produce NULL in 'missing_row_indxs_list[[l]]' if FALSE
      if (missing_perc_by_rowNAs_temp[l]!=0) {
        if (length(dat_sim_row_indxs_temp)==1) {
          # reason: run 'sample(x=5,size=1)' bunch of times
          missing_row_indxs_list[[l]] <- dat_sim_row_indxs_temp
        } else {
          missing_row_indxs_list[[l]] <- 
            sample(dat_sim_row_indxs_temp,
                   size=round(missing_perc_by_rowNAs_temp[l]*N.J_sim_temp,0))
          dat_sim_row_indxs_temp <- dat_sim_row_indxs_temp[
            -which(dat_sim_row_indxs_temp %in% missing_row_indxs_list[[l]])]
        }
      }
    }
    
    dat_sim_wNA_unbal_temp <- dat_sim_full_unbal_temp
    for (l in 1:L) {
      # need this conditional because
      # 'size' has to be <= non-zero probabilities in 'prob' for the 'sample' function
      if (missing_perc_by_rowNAs_temp[l]!=0) {
        if (l < L) {
          for (indx in 1:length(missing_row_indxs_list[[l]])) {
            col_indxs_temp <- sample(1:L,size=l,prob=missing_perc_by_muscle_temp)
            row_indx_temp <- missing_row_indxs_list[[l]][indx]
            dat_sim_wNA_unbal_temp[row_indx_temp,col_indxs_temp] <- NA
          }
        } else if (l==L) {
          dat_sim_wNA_unbal_temp[missing_row_indxs_list[[l]],] <- NA
        }
      }
    }
    
    N_temp <- length(y_sim_full_unbal_temp)
    y_sim_wNA_unbal_temp <- vector("list", length=N_temp)
    vec_sim_wNA_unbal_temp <- c(t(dat_sim_wNA_unbal_temp))
    for (i in 1:N_temp) {
      num_outcomes_subject_i <- p_indx
      y_sim_wNA_unbal_temp[[i]] <- vec_sim_wNA_unbal_temp[1:num_outcomes_subject_i]
      vec_sim_wNA_unbal_temp <- vec_sim_wNA_unbal_temp[-c(1:num_outcomes_subject_i)]
    }
    
    # I think keep 'append' for now (unless better way later) b/c 'y_sim_wNA_unbal_temp' has varying lengths
    y_sim_wNA_unbal <- append(y_sim_wNA_unbal, y_sim_wNA_unbal_temp)
    
  }
  
  # same idea as y_mat_Lcols
  dat_sim_wNA_unbal <- matrix(unlist(y_sim_wNA_unbal), ncol=L, byrow=TRUE)
  
  
  
  if (test_bool==TRUE) {
    ### Evaluate simulated data with missingness
    # Note that we used posterior median of mu, s, r (compare this - overall)
    # and used missingness of actual data (compare this - by Js)
    missing_perc_mats_sim <- missing_perc2_func(y_sim_wNA_unbal, L)
    missing_perc_by_muscle_sim <-  missing_perc_mats_sim$missing_perc_by_muscle
    missing_perc_by_rowNAs_sim <- missing_perc_mats_sim$missing_perc_by_rowNAs
    
    Js_dist_stuff_sim <- Js_dist_func(Js_sim, N_sim)
    mat_Js_dist_sim <- Js_dist_stuff_sim$mat_Js_dist
    perc_Js_dist_sim <- Js_dist_stuff_sim$perc_Js_dist
    
    the_MSEs <- MSE_func(dat_sim_full_unbal, dat_sim_wNA_unbal,
                         target_short_mu, target_short_s, target_r,
                         target_w, target_SRM)
  }
  
  
  
  if (test_bool==TRUE) {
    return(list(missing_perc_by_muscle=missing_perc_by_muscle,
                missing_perc_by_muscle_sim=missing_perc_by_muscle_sim,
                missing_perc_by_rowNAs=missing_perc_by_rowNAs,
                missing_perc_by_rowNAs_sim=missing_perc_by_rowNAs_sim,
                perc_Js_dist=perc_Js_dist,
                perc_Js_dist_sim=perc_Js_dist_sim,
                the_MSEs=the_MSEs,
                
                ps_sim=ps_sim,
                dat_sim_full=dat_sim_full_unbal, # 'dat_sim_full' b/c want same naming as other functions
                dat_sim_wNA=dat_sim_wNA_unbal # 'dat_sim_wNA' b/c want same naming as other functions
    )
    )
  } else if (test_bool==FALSE) {
    return(list(y_sim_full_unbal2=y_sim_full_unbal2,
                y_sim_wNA_list=y_sim_wNA_unbal # 'y_sim_wNA_list' b/c want same naming as other functions
    )
    )
  }
  
}





simplex_vertices_func <- function(L) {
  qr.Q(qr(matrix(1, nrow=L)) ,complete = TRUE)[,-1]
}



# Below are gemPlot functions (source code from a paper)
gridfun <- function (D, grid.size, k=4) 
{
  if (dim(D)[2]<3) message("Error: number of dimensions must be greater or equal 3.")
  if (dim(D)[2]==3) {
    grid.x <- seq(min(D$x)[1], max(D$x)[1], length.out = grid.size)
    grid.y <- seq(min(D$y)[1], max(D$y)[1], length.out = grid.size)
    grid.z <- seq(min(D$z)[1], max(D$z)[1], length.out = grid.size)
    H <- array(NA, c(grid.size, grid.size, grid.size))
    return(list(grid.x = grid.x, grid.y = grid.y, grid.z = grid.z, H = H))
  }              
  if (dim(D)[2]>3) {
    grid.k <- matrix(NA, k, grid.size)
    for (i in 1:k) grid.k[i,] <- seq(min(D[,i])[1], max(D[,i])[1], length.out = grid.size)
    H <- array(NA, rep(grid.size, k))
    return(list(grid.k = grid.k, H = H))
  }
}

hldepth = function (D, G, verbose = TRUE) 
{
  if (dim(D)[2]==3) {
    n <- dim(D)[1]
    grid.size <- length(G$grid.x)
    perc <- 10
    H <- G$H
    for (i in 1:grid.size) {
      for (j in 1:grid.size) {
        for (k in 1:grid.size) {
          u <- c(G$grid.x[i], G$grid.y[j], G$grid.z[k])
          H[i, j, k] <- n * as.numeric(DepthProc::depth(u, as.matrix(D),
                                                        method="Tukey"))
          
          
        }
      }
      if (100 * i/grid.size >= perc && verbose == TRUE) {
        message(paste("Calculation of halfspace location depths: ", 
                      round(100 * i/grid.size, 0), " % of grid points done", 
                      sep = ""))
        perc <- perc + 10
      }
    }
  }
  if (dim(D)[2]>3) {
    n <- dim(D)[1]
    k = dim(G$grid.k)[1]
    grid.size <- dim(G$grid.k)[2]
    grid.kk <- grid.size^k
    H <- G$H
    H2 <- H
    H2[1:grid.kk] <- array(1:grid.kk, rep(grid.size, k))
    perc <- 10
    for (j in 1:grid.kk) {
      index <- which(H2==j, arr.ind=TRUE)
      u <- rep(NA, k)
      for (i in 1:k) u[i] <- G$grid.k[i,index[i]]
      H[j] <- n * as.numeric(DepthProc::depth(u, as.matrix(D),
                                              method="Tukey"))
      if (100 * j/grid.kk >= perc && verbose == TRUE) {
        message(paste("Calculation of halfspace location depths: ", 
                      round(100 * j/grid.kk, 0), " % of grid points done", 
                      sep = ""))
        perc <- perc + 10
      }
    }
  }
  return(H)            
}
depmed = function (G) 
{
  if (length(G$grid.k)==0) {
    grid.size <- length(G$grid.x)
    maxi <- max(G$H)[1]
    H2 <- (G$H == maxi)
    i2 <- rep(0, grid.size)
    j2 <- rep(0, grid.size)
    k2 <- rep(0, grid.size)
    for (i in 1:grid.size) {
      for (j in 1:grid.size) {
        for (k in 1:grid.size) {
          if (H2[i, j, k] == TRUE) {
            i2[i] <- i
            j2[j] <- j
            k2[k] <- k
          }
        }
      }
    }
    med.i <- median(G$grid.x[which(i2 != 0)])
    med.j <- median(G$grid.y[which(j2 != 0)])
    med.k <- median(G$grid.z[which(k2 != 0)])
    return(c(med.i, med.j, med.k))
  }
  if (length(G$grid.k)>0) {
    k = dim(G$grid.k)[1]
    maxi = max(G$H)[1]
    index = which(G$H==maxi, arr.ind=TRUE)
    med = matrix(NA, k, dim(index)[1])
    for (i in 1:k) med[i,] = G$grid.k[i,index[i]]
    med = apply(med, 1, median)
    return(med)
  }
}

bag <- function (D, G) 
{
  require(geometry)
  if (dim(D)[2]==3) {
    grid.size = dim(G$H)[1]
    n <- dim(D)[1]
    D.k <- rep(NA, n)
    for (i in 1:n) {
      I <- matrix(NA, 8, 3)
      I[1, ] <- c(G$grid.x[max(which(G$grid.x <= D[i, 1]))], 
                  G$grid.y[max(which(G$grid.y <= D[i, 2]))],
                  G$grid.z[max(which(G$grid.z <= D[i, 3]))])
      I[2, ] <- c(G$grid.x[max(which(G$grid.x <= D[i, 1]))], 
                  G$grid.y[max(which(G$grid.y <= D[i, 2]))],
                  G$grid.z[min(which(G$grid.z >= D[i, 3]))])
      I[3, ] <- c(G$grid.x[max(which(G$grid.x <= D[i, 1]))], 
                  G$grid.y[min(which(G$grid.y >= D[i, 2]))],
                  G$grid.z[max(which(G$grid.z <= D[i, 3]))])
      I[4, ] <- c(G$grid.x[max(which(G$grid.x <= D[i, 1]))], 
                  G$grid.y[min(which(G$grid.y >= D[i, 2]))],
                  G$grid.z[min(which(G$grid.z >= D[i, 3]))])
      I[5, ] <- c(G$grid.x[min(which(G$grid.x >= D[i, 1]))], 
                  G$grid.y[max(which(G$grid.y <= D[i, 2]))],
                  G$grid.z[max(which(G$grid.z <= D[i, 3]))])
      I[6, ] <- c(G$grid.x[min(which(G$grid.x >= D[i, 1]))], 
                  G$grid.y[max(which(G$grid.y <= D[i, 2]))],
                  G$grid.z[min(which(G$grid.z >= D[i, 3]))])
      I[7, ] <- c(G$grid.x[min(which(G$grid.x >= D[i, 1]))], 
                  G$grid.y[min(which(G$grid.y >= D[i, 2]))],
                  G$grid.z[max(which(G$grid.z <= D[i, 3]))])
      I[8, ] <- c(G$grid.x[min(which(G$grid.x >= D[i, 1]))], 
                  G$grid.y[min(which(G$grid.y >= D[i, 2]))],
                  G$grid.z[min(which(G$grid.z >= D[i, 3]))])
      I <- cbind(I, NA)
      for (t in 1:8) {
        index1 <- match(I[t, 1], G$grid.x)
        index2 <- match(I[t, 2], G$grid.y)
        index3 <- match(I[t, 3], G$grid.z)
        I[, 4] <- G$H[index1, index2, index3]
      }
      D.k[i] <- min(I[, 4])
    }
    H2 <- (G$H >= max(which(cumsum(table(D.k)) <= (n/2))))
    BAG <- matrix(NA, 0, 3)
    for (i in 1:grid.size) {
      for (j in 1:grid.size) {
        for (k in 1:grid.size) {
          if (H2[i, j, k] == TRUE) {
            BAG <- rbind(BAG, c(G$grid.x[i], G$grid.y[j], 
                                G$grid.z[k]))
          }
        }
      }
    }
    convH <- convhulln(BAG)
  }
  if (dim(D)[2]>3) {
    n = dim(D)[1]
    d = dim(D)[2]
    k = dim(G$grid.k)[2]
    I = matrix(NA, 2^d, d)
    D.k = rep(NA, n)
    for (i in 1:n) {
      U = cbind(G$grid.k, as.numeric(D[i,]))
      U2 = U
      for (t in 1:d) {
        U2[t,] = rank(U[t,], ties.method="first")
        dimnames(G$H)[[t]] = 1:k
      }
      I = t(matrix(U2[,k+1]-1, nrow=d, ncol=2^d))
      for (t in 1:d) I[,t] = I[,t] + as.numeric(gl(2, 2^(d-t), 2^d)) - 1
      I[which(I>k)] = k
      I = cbind(I, NA)
      for (s in 1:(2^d)) I[s,d+1] = G$H[t(matrix(as.character(I[s,1:d])))]
      D.k[i] = min(I[,d+1])
    }
    H2 <- (G$H >= max(which(cumsum(table(D.k)) <= (n/2))))
    H3 <- H2
    H3[1:(k^d)] <- array(1:(k^d), rep(k, d))
    BAG = matrix(NA, table(H2)[2], d)
    index = which(H2==TRUE)
    for (t in 1:length(index)) {
      index2 = which(H3==index[t], arr.ind=TRUE)
      for (j in 1:d) BAG[t,j] = G$grid.k[j,index2[j]]
    }
    convH <- convhulln(BAG)
  }
  return(list(coords = BAG, hull = convH))
}

loop <- function (D, B, inflation = 3, dm) 
{
  n <- dim(D)[1]
  d = dim(D)[2]
  if (d==3) dm = matrix(dm, 1, 3)
  index.F <- sort(intersect(as.vector(B$hull), as.vector(B$hull)))
  FENCE <- B$coords[index.F, ]
  MED.MAT <- t(matrix(dm, d, dim(FENCE)[1]))
  FENCE <- MED.MAT + inflation * (FENCE - MED.MAT)
  colnames(FENCE) <- colnames(D)
  convH <- convhulln(FENCE)
  outliers <- rep(0, n)
  for (i in 1:n) {
    Z <- rbind(FENCE, D[i, ])
    convH.Z <- convhulln(Z)
    if (!is.na(match(dim(FENCE)[1] + 1, convH.Z))) {
      outliers[i] <- 1
    }
  }
  LOOP <- D[which(outliers == 0), ]
  convH2 <- convhulln(LOOP)
  return(list(coords.loop = LOOP, hull.loop = convH2, coords.fence = FENCE, 
              hull.fence = convH, outliers = outliers))
}

gem <- function (coords, hull, clr) 
{
  require(geometry)
  for (i in 1:dim(hull)[1]) {
    x <- coords[hull[i, ], 1]
    y <- coords[hull[i, ], 2]
    z <- coords[hull[i, ], 3]
    triangles3d(x, y, z, col = clr)
  }
}





save_sim_w_and_SRM_RDS_func <- function(n_sim_data, output_sub_dir,
                                        nchains, burn_in_prop,
                                        muscle_names_simple, muscle_names_simple_subset,
                                        r_names,
                                        optim_method, grad_bool, upper_bool,
                                        w_names_detailed) {
  for (sim_dat_indx in 1:n_sim_data) {
    if (!file.exists(paste0("w&SRM",sim_dat_indx,".RDS"))) {
      print(paste0("Saving simulation #",sim_dat_indx,"'s weights and SRM values."))
      mult_chains_temp <- readRDS(paste0(output_sub_dir,"/","RDS","/",
                                         "mcmc",sim_dat_indx,".RDS"))
      grad_bool <- TRUE
      upper_bool <- TRUE
      SRM_opt_stuff <- SRM_opt_subset(mult_chains_temp, nchains, burn_in_prop,
                                      muscle_names_simple, muscle_names_simple_subset, r_names,
                                      optim_method, grad_bool, upper_bool)
      
      SRMw <- SRM_opt_stuff$SRMw # 3D-array, accounting for burn-in
      SRMw[which(SRMw<0,arr.ind=TRUE)] <- 0
      SRMw <- SRMw[,1,] # since only 1 chain for simulation
      colnames(SRMw) <- w_names_detailed
      
      SRMw_orig <- SRM_opt_stuff$SRMw_orig
      SRMw_orig[which(SRMw_orig<0,arr.ind=TRUE)] <- 0
      SRMw_orig <- SRMw_orig[,1,] # since only 1 chain for simulation
      colnames(SRMw_orig) <- w_names_detailed
      
      SRMvalues <- SRM_opt_stuff$SRMvalues
      colnames(SRMvalues) <- "SRM"
      
      saveRDS(list(SRMw=SRMw,
                   SRMw_orig=SRMw_orig,
                   SRMvalues=SRMvalues),
              file=paste0("w&SRM",sim_dat_indx,".RDS")) 
    }
    
  }
}

save_sim_w_and_SRM_from_post_med_RDS_func <- function(n_sim_data, output_sub_dir,
                                                      burn_in_number, iters, n_etas,
                                                      optim_method, grad_bool, upper_bool) {
  for (sim_dat_indx in 1:n_sim_data) {
    if (!file.exists(paste0("w&SRM_from_post_med",sim_dat_indx,".RDS"))) {
      print(paste0("Saving simulation #",sim_dat_indx,"'s optimal weights"))
      mult_chains_temp <- readRDS(paste0(output_sub_dir,"/","RDS","/",
                                         "mcmc",sim_dat_indx,".RDS"))
      # I can do this because 1 chain in simulation (different from 'target_param_func')
      short_mu_post_med_temp <- apply(mult_chains_temp[[1]]$short_mu[(burn_in_number+1):iters,], 2, median)
      short_s_post_med_temp <- apply(mult_chains_temp[[1]]$short_s[(burn_in_number+1):iters,], 2, median)
      r_post_med_temp <- apply(mult_chains_temp[[1]]$r_actual[(burn_in_number+1):iters,], 2, median)
      
      muscle_Sigma_post_med_temp <- 
        get_muscle_Sigma_func(short_s_post_med_temp, r_post_med_temp,
                              L, n_etas)
      
      w_from_post_med_stuffs_temp <-
        muscle_optim_func(short_mu_post_med_temp, muscle_Sigma_post_med_temp,
                          optim_method, grad_bool, upper_bool)
      w_from_post_med_temp <- w_from_post_med_stuffs_temp$w
      
      SRMval_from_post_med_temp <- 
        SRM(w_from_post_med_temp,
            short_mu_post_med_temp, muscle_Sigma_post_med_temp)
      
      saveRDS(list(w_from_post_med_temp=w_from_post_med_temp,
                   SRMval_from_post_med_temp=SRMval_from_post_med_temp),
              file=paste0("w&SRM_from_post_med",sim_dat_indx,".RDS"))
    }
  } 
}

save_r_accept_RDS_func <- function(n_sim_data, output_sub_dir,
                                   nchains, burn_in_prop) {
  for (sim_dat_indx in 1:n_sim_data) {
    if (!file.exists(paste0("r_accept%",sim_dat_indx,".RDS"))) {
      print(paste0("Saving simulation #",sim_dat_indx,"'s r acceptance rates."))
      mult_chains_temp <- readRDS(paste0(output_sub_dir,"/","RDS","/",
                                         "mcmc",sim_dat_indx,".RDS"))
      main_diagnostics_stuff_temp <- main_diagnostics(mult_chains, nchains, burn_in_prop)
      MH_accept_percentage_r_temp <- main_diagnostics_stuff_temp$MH_accept_percentage_r
      
      saveRDS(MH_accept_percentage_r_temp, paste0("r_accept%",sim_dat_indx,".RDS"))
    }
  }
}

save_R_PD_RDS_func <- function(n_sim_data, output_sub_dir,
                               nchains, burn_in_prop) {
  for (sim_dat_indx in 1:n_sim_data) {
    if (!file.exists(paste0("%times_R_PD",sim_dat_indx,".RDS"))) {
      print(paste0("Saving simulation #",sim_dat_indx,"'s % times R is PD."))
      mult_chains_temp <- readRDS(paste0(output_sub_dir,"/","RDS","/",
                                         "mcmc",sim_dat_indx,".RDS"))
      main_diagnostics_stuff_temp <- main_diagnostics(mult_chains, nchains, burn_in_prop)
      R_PD_percentage_temp <- main_diagnostics_stuff_temp$R_PD_percentage
      
      saveRDS(R_PD_percentage_temp, paste0("%times_R_PD",sim_dat_indx,".RDS"))
    }
  }
}