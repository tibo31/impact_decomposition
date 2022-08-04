#############################################################################
###################################

compute_effect <- function(Z, rho, OW, DW, ind_o, ind_d, change_z = NULL,
                           site = "all", Y_in_log = F, hat_y = NULL) {
  
  # verification
  if (Y_in_log & is.null(hat_y)) {
    stop("If Y is in log, hat_y must be given to compute decomposition")
  }
  
  # Z is a DF of size p x 5 (beta_o, beta_d, delta_o, delta_d, beta_i)
  p <- nrow(Z)
  names_p <- row.names(Z)
  
  # determine the type of each variable
  type_z <- numeric(p)
  for (k in 1:p) {
    if (Z[k, "beta_d"] == 0 & Z[k, "delta_d"] == 0 & Z[k, "beta_i"] == 0) {
      type_z[k] <- 1
    } else {
      if (Z[k, "beta_o"] == 0 & Z[k, "delta_o"] == 0 & Z[k, "beta_i"] == 0) {
        type_z[k] <- 2 
      } else {
        type_z[k] <- 3
      }
    }
  }
  
  # how much do we increase Z
  if(is.null(change_z) | length(change_z) != p) {
    change_z <- rep(1, p)
  }
  
  # initialization 
  O <- unique(ind_o)
  n_o <- length(O)
  D <- unique(ind_d)
  n_d <- length(D)
  OD <- intersect(O, D)
  n_od <- length(OD)
  N <- length(ind_o)
  rho_o <- rho[1] 
  rho_d <- rho[2] 
  rho_w <- rho[3]
  
  # compute the number of ongoing resp. arriving flows per site
  n_outgoing <- table(ind_o)
  n_ingoing <- table(ind_d)
  
  # compute the effect site by site for all sites
  if (site[1] == "all") {
    S <- c(setdiff(O, OD), setdiff(D, OD), OD)  
  } else {
    S <- site
  }
  
  # preparation of the results 
  n_s <- length(S)
  effect_o <- matrix(0, n_s, 5)
  rownames(effect_o) <- S
  colnames(effect_o) <- c("IE", "OE", "DE", "NE", "TE")
  effect_d <- effect_delta_o <- effect_delta_d <- effect_o
  
  # we detect the flows that have intra flows 
  intra_flows_index <- intersect(OD, intersect(ind_o[ind_o == ind_d], S))
  
  #####################################
  # computation of the number of local impact
  n_effect <- effect_o
  
  # Total
  n_effect[S, "TE"] <- N
  # intra 
  n_effect[intra_flows_index, "IE"] <- 1
  # origin
  n_effect[setdiff(S[S %in% O], D), "OE"] <- n_outgoing[setdiff(S[S %in% O], D)] 
  n_effect[intersect(OD, S[S %in% O]), "OE"] <- n_outgoing[intersect(OD, S[S %in% O])] - 1
  # destination
  n_effect[setdiff(S[S %in% D], O), "DE"] <- n_ingoing[setdiff(S[S %in% D], O)] 
  n_effect[intersect(OD, S[S %in% D]), "DE"] <- n_ingoing[intersect(OD, S[S %in% D])] - 1
  # network 
  n_effect[setdiff(S[S %in% O], D), "NE"] <- N - n_outgoing[setdiff(S[S %in% O], D)] 
  n_effect[setdiff(S[S %in% D], O), "NE"] <- N - n_ingoing[setdiff(S[S %in% D], O)] 
  n_effect[intersect(OD, S), "NE"] <- N - n_outgoing[intersect(OD, S)] - n_ingoing[intersect(OD, S)] + 1
  
  ############################
  # we compute lambda_prime
  lambda_prime <- matrix(0, n_o, n_o)
  for(i in 1:n_o) {
    for (j in 1:n_o) {
      if (i == j) {
        lambda_prime[i, j] <- rho_d
      }
      lambda_prime[i, j] <- lambda_prime[i, j] + (rho_o + rho_w) * OW[O[i], O[j]] 
    }
  }
  lambda <- solve(diag(n_o) - lambda_prime) 
  rownames(lambda) <- colnames(lambda) <- O
  
  # compute the gamma_prime function
  gamma_prime <- matrix(0, n_d, n_d)
  for(i in 1:n_d) {
    for (j in 1:n_d) {
      if (i == j) {
        gamma_prime[i, j] <- rho_o
      }
      # gamma_prime[i, j] <- sum(Aw[ind_i, ind_j])
      gamma_prime[i, j] <- gamma_prime[i, j] + (rho_d + rho_w) * DW[D[i], D[j]]
    }
  }
  gamma <- solve(diag(n_d) - gamma_prime) 
  rownames(gamma) <- colnames(gamma) <- D
  
  # Take into account the beta_i parameter if necessary
  # If s has an intra characteristic, we need to add the effect due to beta_i 
  if (any(Z[, "beta_i"] != 0) & length(intra_flows_index) > 0) {
    A <- diag(N) - rho_o * kronecker(OW, diag(n_d)) - rho_d * kronecker(diag(n_o), DW) - 
      rho_w * kronecker(OW, DW)
    b_mat <- matrix(0, N, length(intra_flows_index))
    for(k in 1:length(intra_flows_index)) {
      b_mat[which(ind_o == intra_flows_index[k] & ind_d == intra_flows_index[k]), k] <- 1
    }
    A_w_col <- solve(A, b_mat)
    colnames(A_w_col) <- intra_flows_index
  }
  
  ################################
  # computation of the local impact
  
  # Option 1 : Y is not in log 
  # in that case, some simplifications can be done and we do not require \hat Y 

  if (!Y_in_log) {
    # loop on the intersections
    for(s in intra_flows_index) {
      effect_o[s, "IE"] <- lambda[s, s] 
      effect_d[s, "IE"] <- gamma[s, s] 
      
      # for the lag  part 
      effect_delta_o[s, "IE"] <- sum(lambda[s, ] * OW[, s])  
      effect_delta_d[s, "IE"] <- sum(gamma[s, ] * DW[, s])  
    }
    
    # loop on the origin sites
    for(s in S[S %in% O]) {
      effect_o[s, "OE"] <- lambda[s, s] * n_outgoing[s] 
      if (s %in% D) {
        effect_d[s, "OE"] <-  sum(gamma[ind_d[ind_o == s], s]) 
      }
      
      # for the lagged variables 
      effect_delta_o[s, "OE"] <- sum(lambda[s, ] * OW[, s]) * n_outgoing[s] 
      if (s %in% D) {
        effect_delta_d[s, "OE"] <-  sum(
          sapply(D, function(y) sum(gamma[ind_d[ind_o == s], y])) * DW[, s]
          ) 
      }
    }
    
    # loop on the destination sites  
    for(s in S[S %in% D]) {
      if (s %in% O) {
        effect_o[s, "DE"] <- sum(lambda[ind_o[ind_d == s], s]) 
      }
      
      effect_d[s, "DE"] <-  gamma[s, s] * n_ingoing[s]
      
      # for the lagged variables 
      if (s %in% O) {
        effect_delta_o[s, "DE"] <-  sum(
          sapply(O, function(y) sum(lambda[ind_o[ind_d == s], y])) * OW[, s]
          ) 
      } 
      effect_delta_d[s, "DE"] <-  sum(gamma[s, ] * DW[, s]) * n_ingoing[s]
    }
    
    # loop on the total effects
    for(s in S) {
      if (s %in% O)
        effect_o[s, "TE"] <- sum(lambda[ind_o, s]) 
      if (s %in% D)
        effect_d[s, "TE"] <-  sum(gamma[ind_d, s]) 
      if (s %in% O)
        effect_delta_o[s, "TE"] <- sum(lambda[ind_o, ] %*% OW[, s]) 
      if (s %in% D)
        effect_delta_d[s, "TE"] <- sum(gamma[ind_d, ] %*% DW[, s]) 
    }
    
    ###############
    res <- vector("list", p)
    for(k in 1:p) {
      # definition of S
      if (type_z[k] == 1) {
        S <- O 
        } else {
          if (type_z[k] == 2) {
            S <- D
            } else {
              if (type_z[k] == 3) {
                # ordering S with respect to O, D and O\cap D
                S <- c(setdiff(O, OD), setdiff(D, OD), OD)  
              }
            }
        }
      
      if (site[1] != "all") {
        S <- intersect(S, site)
      }
      n_s <- length(S)
      effect <- matrix(0, n_s, 5)
      rownames(effect) <- S
      colnames(effect) <- c("IE", "OE", "DE", "NE", "TE")

      effect[, "IE"] <- Z[k, "beta_o"] * effect_o[S, "IE"]  + Z[k, "beta_d"] * effect_d[S, "IE"] + 
        Z[k, "delta_o"] * effect_delta_o[S, "IE"] + Z[k, "delta_d"] * effect_delta_d[S, "IE"]  
  
      effect[, "OE"] <- Z[k, "beta_o"] * effect_o[S, "OE"] + Z[k, "beta_d"] * effect_d[S, "OE"] +
        Z[k, "delta_o"] * effect_delta_o[S, "OE"] + Z[k, "delta_d"] * effect_delta_d[S, "OE"] - effect[, "IE"]
  
      effect[, "DE"] <- Z[k, "beta_o"] * effect_o[S, "DE"] + Z[k, "beta_d"] * effect_d[S, "DE"] +
        Z[k, "delta_o"] * effect_delta_o[S, "DE"] + Z[k, "delta_d"] * effect_delta_d[S, "DE"] - effect[, "IE"]
  
      effect[, "TE"] <- Z[k, "beta_o"] * effect_o[S, "TE"] + Z[k, "beta_d"] * effect_d[S, "TE"] +
        Z[k, "delta_o"] * effect_delta_o[S, "TE"] + Z[k, "delta_d"] * effect_delta_d[S, "TE"] 
  
      effect[, "NE"] <- effect[, "TE"] - effect[, "IE"] - effect[, "OE"] - effect[, "DE"]
  
      # case beta_i is non null  
      if (Z[, "beta_i"] != 0 & length(intra_flows_index) > 0) {
        for(j in 1:length(intra_flows_index)) {
          s <- intra_flows_index[j]
          effect[s, "IE"] <- effect[s, "IE"] + Z[, "beta_i"] * A_w_col[ind_o == s & ind_d == s, s]
          effect[s, "OE"] <- effect[s, "OE"] + Z[, "beta_i"] * sum(A_w_col[ind_o == s & ind_d != s, s])
          effect[s, "DE"] <- effect[s, "DE"] + Z[, "beta_i"] * sum(A_w_col[ind_d == s & ind_o != s, s]) 
          effect[s, "NE"] <- effect[s, "NE"] + Z[, "beta_i"] * sum(A_w_col[ind_d != s & ind_o != s, s]) 
          effect[s, "TE"] <- effect[s, "TE"] + Z[, "beta_i"] * sum(A_w_col[, j])
        }  
      }
      res[[k]] <- list(local_effect = change_z[k] * effect, n_effect = n_effect[S, ])
    }
  } else {
    # Option 2 : Y is in log 
    # we do require \hat Y 

    res <- vector("list", p)
    for(k in 1:p) {
      names_k <- names_p[k]
      
      # definition of S
      if (type_z[k] == 1) {
        S <- O 
      } else {
        if (type_z[k] == 2) {
          S <- D
        } else {
          if (type_z[k] == 3) {
            # ordering S with respect to O, D and O\cap D
            S <- c(setdiff(O, OD), setdiff(D, OD), OD)  
          }
        }
      }
      
      if (site[1] != "all") {
        S <- intersect(S, site)
      }
      n_s <- length(S)
      effect <- matrix(0, n_s, 5)
      rownames(effect) <- S
      colnames(effect) <- c("IE", "OE", "DE", "NE", "TE")

      for (j in 1:n_s){
        s <- S[j]
        temp <- rep(0, N)
        if (s %in% O)
          temp <- temp + lambda[ind_o, s] * Z[names_k, "beta_o"] 
        
        if (s %in% D)
          temp <- temp + gamma[ind_d, s] * Z[names_k, "beta_d"] 
        
        if (s %in% O)
          temp <- temp + OW[, s] %*% t(lambda[ind_o, ]) * Z[names_k, "delta_o"] 
        
        if (s %in% D)
          temp <- temp + DW[, s] %*% t(gamma[ind_d, ]) * Z[names_k, "delta_d"]
        
        if (Z[names_k, "beta_i"] != 0 & s %in% OD) {
          temp <- temp +  Z[names_k, "beta_i"] * A_w_col[, s]
        }
        
        temp <- temp * change_z[k]
        temp <- hat_y * (exp(temp) - 1)
        
        # index
        ind_intra <- which(ind_o == s & ind_d == s)
        ind_d_eff <- which(ind_d == s & ind_o != s)
        ind_o_eff <- which(ind_o == s & ind_d != s)
        ind_net <- which(ind_d != s & ind_o != s)
      
        if(length(ind_intra) > 0) {
          effect[s, "IE"] <- temp[ind_intra]
        }
        
        if(length(ind_o_eff) > 0) {
          effect[s, "OE"] <- sum(temp[ind_o_eff])
        }
        
        if(length(ind_d_eff) > 0) {
          effect[s, "DE"] <- sum(temp[ind_d_eff])
        }    
        
        if(length(ind_net) > 0) {
          effect[s, "NE"] <- sum(temp[ind_net])
        }   
        
        effect[s, "TE"] <- sum(temp)
      }
      res[[k]] <- list(local_effect = effect, n_effect = n_effect[S, ])
    }
  }  
  names(res) <- names_p
  return(res)
}

# ######################################
# Z <- data.frame(row.names = c("electric", "population", "island", "inflation", "mult_for_exchange", "for_exc_volatility", "transfer_cost", "political_sta", "landlocked", "natural_disaster"),
#                 beta_o = c(0.432, 0.984, 2.244, -1.537, -1.089, -0.019, -0.084, 0, 0, 0),
#                 beta_d = c(-1.431, 0.612, 0, -0.532, 0, 0, 0, -0.552, 0.592, -0.482),
#                 delta_o = c(-0.474, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#                 delta_d = c(0.107, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#                 beta_i = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
# rho <- c(rho_o, rho_d, rho_w)
# 
# compute_effect(Z, rho, OW = w, DW = w, ind_o = index_o, ind_d = index_d)
# 
# #######################################
# my_x <- data.frame(names = "x",
#                    beta_o = 0.5,
#                    beta_d = 1,
#                    beta_i = 2,
#                    delta_o = 0.1,
#                    delta_d = 0.3)
# compute_effect(my_x, rho, OW = w, DW = w, ind_o = index_o, ind_d = index_d)
# 
# #######################################
# my_x <- data.frame(names = c("x_type_1", "x_type_2", "x_type_3"),
#                    beta_o = c(0.5, 0, 0.5),
#                    beta_d = c(0, 1, 1),
#                    beta_i = c(0, 0, 0),
#                    delta_o = c(0.1, 0, 0.1),
#                    delta_d = c(0, 0.3, 0.3))
# compute_effect(my_x, rho, OW = ow, DW = dw, ind_o = index_o_2, ind_d = index_d_2)
