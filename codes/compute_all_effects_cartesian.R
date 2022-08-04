
#########################################################
#####  Compute all effects

compute_all_effect <- function(Z, rho, OW, DW, ind_o, ind_d, change_z = NULL,
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
  
  # computation of the local impact
  res_local <- data.frame(site = character(), type = character(), 
                          var = character(),
                          values = numeric())

  for(k in 1:p) {
    beta_o <- Z[k, "beta_o"] * change_z[k]
    beta_d <- Z[k, "beta_d"] * change_z[k]    
    delta_o <- Z[k, "delta_o"] * change_z[k]
    delta_d <- Z[k, "delta_d"] * change_z[k]
    beta_i <- Z[k, "beta_i"] * change_z[k]
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
    
    # compute the effect for variable k 
    for(j in 1:length(S)) {
      s <- S[j]
      if (s %in% O) {
        temp <- beta_o * lambda[ind_o, s]
        for(j in 1:N) {
          temp[j] <- temp[j] + delta_o * sum(OW[, s] * lambda[ind_o[j], ]) 
        }
        } else {
          temp <- rep(0, N)
        }
      
      if (s %in% D) {
        temp <- temp + beta_d * gamma[ind_d, s]  
        for(j in 1:N) {
          temp[j] <- temp[j] + delta_d * sum(DW[, s] * gamma[ind_d[j], ])
        }
      }
      
      if (beta_i != 0 & s %in% intra_flows_index) {
        temp <- temp + beta_i * A_w_col[, s]
      }
      if (Y_in_log) {
        temp <- hat_y * (exp(temp) - 1)
      }
      ind_intra <- which(ind_o == s & ind_d == s)
      ind_origin <- which(ind_o == s & ind_d != s)
      ind_dest <- which(ind_d == s & ind_o != s)
      ind_net <- which(ind_d != s & ind_o != s)
      
      if (length(ind_intra) > 0) {
        res_local <- rbind(res_local, 
                           data.frame(site = s, type = "intra", 
                                      var = names_p[k], values = temp[ind_intra])
                           )
      }
      
      if (length(ind_origin) > 0) {
        res_local <- rbind(res_local, 
                           data.frame(site = s, type = "origin", 
                                      var = names_p[k], values = temp[ind_origin])
      )
      }
      
      if (length(ind_dest) > 0) {
        res_local <- rbind(res_local, 
                           data.frame(site = s, type = "destination", 
                                      var = names_p[k], values = temp[ind_dest])
      )
      }
      res_local <- rbind(res_local, 
                         data.frame(site = s, type = "network", 
                                    var = names_p[k], 
                                    values = temp[-c(ind_intra, ind_origin, ind_dest)])
      )
      }
    }
  return(res_local = res_local)
}

####################################################################
