#######################################################################
####  Version with simplified formula

# Input parameters:
# beta_o, beta_d, beta_i, delta_o, delta_d: coefficients in model (3)
# ind_o, ind_d: indexes of the flows
# OW, DW: spatial weight matrices of size no and nd
# rho_o, rho_d, rho_w: spatial autocorrelation parameters
# type_x and change_x: TBC

compute_effect <- function(beta_o, beta_d, beta_i, delta_o, delta_d, 
                           ind_o, ind_d,
                           OW, DW, 
                           rho_o, rho_d, rho_w,
                           type_x = 3, change_x = 1,
                           log_Y = F) {
  
  # initialization 
  O <- unique(ind_o)
  n_o <- length(O)
  D <- unique(ind_d)
  n_d <- length(D)
  OD <- intersect(O, D)
  n_od <- length(OD)
  N <- length(ind_o)
  p <- length(beta_o)

  # definition of S
  if (type_x == 1) {
    S <- O 
  } else {
    if (type_x == 2) {
      S <- D
    } else {
      if (type_x == 3) {
        # ordering S with respect to O, D and O\cap D
        S <- c(setdiff(O, OD), setdiff(D, OD), OD)  
      }
    }
  }
  n_s <- length(S)
  
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
  lambda <- solve(diag(n_o) - lambda_prime) * change_x
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
  gamma <- solve(diag(n_d) - gamma_prime) * change_x
  rownames(gamma) <- colnames(gamma) <- D
  
  # computation of the local impact
  effect <- matrix(0, n_s, 5)
  rownames(effect) <- S
  colnames(effect) <- c("IE", "OE", "DE", "NE", "TE")
  effect_o <- effect_d <- effect_delta_o <- effect_delta_d <- effect
  
  # compute the number of ongoing resp. arriving flows per site
  n_outgoing <- table(ind_o)
  n_ingoing <- table(ind_d)
  
  # loop on the intersections
  # we detect the flows that have intra flows 
  intra_flows_index <- ind_o[ind_o == ind_d]
  for(s in intersect(OD, intra_flows_index)) {
    effect_o[s, "IE"] <- lambda[s, s] * beta_o
    effect_d[s, "IE"] <- gamma[s, s] * beta_d 
    
    # for the lag  part 
    effect_delta_o[s, "IE"] <- sum(lambda[s, ] * OW[, s]) * delta_o 
    effect_delta_d[s, "IE"] <- sum(gamma[s, ] * DW[, s]) * delta_d    
  }
  
  # loop on the origin sites
  for(s in S[S %in% O]) {
    effect_o[s, "OE"] <- lambda[s, s] * n_outgoing[s] * beta_o
    if (s %in% D) {
      effect_d[s, "OE"] <-  sum(gamma[ind_d[ind_o == s], s]) * beta_d 
    }
    
    # for the lagged variables 
    effect_delta_o[s, "OE"] <- sum(lambda[s, ] * OW[, s]) * n_outgoing[s] * delta_o 
    if (s %in% D) {
      effect_delta_d[s, "OE"] <-  sum(
        sapply(D, function(y) sum(gamma[ind_d[ind_o == s], y])) * DW[, s]
      ) * delta_d
    }
  }
  
  for(s in S[S %in% D]) {
    if (s %in% O) {
      effect_o[s, "DE"] <- sum(lambda[ind_o[ind_d == s], s]) * beta_o
    }
    effect_d[s, "DE"] <-  gamma[s, s] * n_ingoing[s] * beta_d 
    
    # for the lagged variables 
    if (s %in% O) {
      effect_delta_o[s, "DE"] <-  sum(
        sapply(O, function(y) sum(lambda[ind_o[ind_d == s], y])) * OW[, s]
      ) * delta_o 
    } 
    effect_delta_d[s, "DE"] <-  sum(gamma[s, ] * DW[, s]) * n_ingoing[s] * delta_d
  }
  
  # loop on the total effects
  for(s in S) {
    if (s %in% O)
      effect_o[s, "TE"] <- sum(lambda[ind_o, s]) * beta_o
    if (s %in% D)
      effect_d[s, "TE"] <-  sum(gamma[ind_d, s]) * beta_d 
    
    if (s %in% O)
      effect_delta_o[s, "TE"] <- sum(lambda[ind_o, ] %*% OW[, s]) * delta_o  
    
    if (s %in% D)
      effect_delta_d[s, "TE"] <- sum(gamma[ind_d, ] %*% DW[, s]) * delta_d
  }
  
  
  effect[, "IE"] <- effect_o[, "IE"] + effect_d[, "IE"] +   
    effect_delta_o[, "IE"] + effect_delta_d[, "IE"]  
  
  effect[, "OE"] <- effect_o[, "OE"] + effect_d[, "OE"] +
    effect_delta_o[, "OE"] + effect_delta_d[, "OE"] - effect[, "IE"]
  
  effect[, "DE"] <- effect_o[, "DE"] + effect_d[, "DE"] +
    effect_delta_o[, "DE"] + effect_delta_d[, "DE"] - effect[, "IE"]
  
  effect[, "TE"] <- effect_o[, "TE"] + effect_d[, "TE"] +
    effect_delta_o[, "TE"] + effect_delta_d[, "TE"] 
  
  effect[, "NE"] <- effect[, "TE"] - effect[, "IE"] - effect[, "OE"] - effect[, "DE"]
  
  
  # If s has an intra characteristic, we need to add the effect due to beta_i 
  if (beta_i != 0 & length(OD) > 0) {
    A <- diag(N) - rho_o * kronecker(OW, diag(n_d)) - rho_d * kronecker(diag(n_o), DW) - 
      rho_w * kronecker(OW, DW)
    b_mat <- matrix(0, N, length(OD))
    for(k in 1:length(OD)) {
      b_mat[which(ind_o == OD[k] & ind_d == OD[k]), k] <- 1
    }
    A_w_col <- solve(A, b_mat)
      
    # if (delta_i == 0) {
    #  
    #  A_w_col <- solve(A, )
    #}
    #index_A <- paste(rep(O, times = n_d), rep(D, each = n_o), sep = "_")
    #index_od <- paste(ind_o, ind_d, sep = "_")
    #new_index <- sapply(index_od, function(x) grep(x, index_A))
    #A <- A[new_index, new_index]
    
    # order the matrix A with respect to A 
    
    for(k in 1:length(OD)) {
      s <- OD[k]
      effect[s, "IE"] <- effect[s, "IE"] + beta_i * A_w_col[ind_o == s & ind_d == s, k]
      effect[s, "OE"] <- effect[s, "OE"] + beta_i * sum(A_w_col[ind_o == s & ind_d != s, k])
      effect[s, "DE"] <- effect[s, "DE"] + beta_i * sum(A_w_col[ind_d == s & ind_o != s, k]) 
      effect[s, "NE"] <- effect[s, "NE"] + beta_i * sum(A_w_col[ind_d != s & ind_o != s, k]) 
      effect[s, "TE"] <- effect[s, "TE"] + beta_i * sum(A_w_col[, k])
    }
    
    # if (delta_i != 0) {
    #   b_mat <- matrix(0, N, length(OD))
    #   temp <- (kronecker(OW, diag(n_d)) + kronecker(diag(n_o), DW))
    #   for(k in 1:length(OD)) {
    #     s <- OD[k] 
    #     b_mat[, k] <- delta_i * temp %*% as.numeric(ind_o == s & ind_d == s) / 2
    #   }
    #   A_w_col <- solve(A, b_mat)
    #   
    #   for(k in 1:length(OD)) {
    #     s <- OD[k]
    #     effect[s, "IE"] <- effect[s, "IE"] + A_w_col[ind_o == s & ind_d == s, k]
    #     effect[s, "OE"] <- effect[s, "OE"] + sum(A_w_col[ind_o == s & ind_d != s, k])
    #     effect[s, "DE"] <- effect[s, "DE"] + sum(A_w_col[ind_d == s & ind_o != s, k])
    #     effect[s, "NE"] <- effect[s, "NE"] + sum(A_w_col[ind_d != s & ind_o != s, k])
    #     effect[s, "TE"] <- effect[s, "TE"] + sum(A_w_col[, k])
    #   }
    # }
  }
  
  
  # computation of the number of local impact
  n_effect <- matrix(0, n_s, 5)
  rownames(n_effect) <- S
  colnames(n_effect) <- c("IE", "OE", "DE", "NE", "TE")
  
  # Total
  n_effect[S, "TE"] <- N
  # intra 
  n_effect[intersect(OD, intra_flows_index), "IE"] <- 1
  # origin
  n_effect[setdiff(S[S %in% O], D), "OE"] <- n_outgoing[setdiff(S[S %in% O], D)] 
  n_effect[intersect(OD, S[S %in% O]), "OE"] <- n_outgoing[intersect(OD, S[S %in% O])] - 1
  # destination
  n_effect[setdiff(S[S %in% D], O), "DE"] <- n_ingoing[setdiff(S[S %in% D], O)] 
  n_effect[intersect(OD, S[S %in% D]), "DE"] <- n_ingoing[intersect(OD, S[S %in% D])] - 1
  # network 
  n_effect[setdiff(S[S %in% O], D), "NE"] <- N - n_outgoing[setdiff(S[S %in% O], D)] 
  n_effect[setdiff(S[S %in% D], O), "NE"] <- N - n_ingoing[setdiff(S[S %in% D], O)] 
  n_effect[OD, "NE"] <- N - n_outgoing[OD] - n_ingoing[OD] + 1
  
  return(list(local_effect = effect, n_effect = n_effect))
}


#########################################################
#####  Compute all effects

compute_all_effect <- function(beta_o, beta_d, beta_i, delta_o, delta_d, 
                               ind_o, ind_d, OW, DW, 
                               rho_o, rho_d, rho_w,
                               type_x = 3, change_x = 1) {
  
  # initialization 
  O <- unique(ind_o)
  n_o <- length(O)
  D <- unique(ind_d)
  n_d <- length(D)
  OD <- intersect(O, D)
  n_od <- length(OD)
  N <- length(ind_o)
  p <- length(beta_o)
  
  # definition of S
  if (type_x == 1) {
    S <- O 
  } else {
    if (type_x == 2) {
      S <- D
    } else {
      if (type_x == 3) {
        # ordering S with respect to O, D and O\cap D
        S <- c(setdiff(O, OD), setdiff(D, OD), OD)  
      }
    }
  }
  n_s <- length(S)
  
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
  lambda <- solve(diag(n_o) - lambda_prime) * change_x
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
  gamma <- solve(diag(n_d) - gamma_prime) * change_x
  rownames(gamma) <- colnames(gamma) <- D
  
  # compute the number of ongoing resp. arriving flows per site
  n_outgoing <- table(ind_o)
  n_ingoing <- table(ind_d)
  
  # we detect the flows that have intra flows 
  intra_flows_index <- ind_o[ind_o == ind_d]
  
  # computation of the local impact
  res_local <- data.frame(site = character(), type = character(), values = numeric())
  for(j in 1:length(S)) {
    s <- S[j]
    if (s %in% O) {
      temp <- beta_o * lambda[ind_o, s]
      for(k in 1:N) {
        temp[k] <- temp[k] + delta_o * sum(OW[, s] * lambda[ind_o[k], ]) 
      }
    } else {
      temp <- rep(0, N)
    }
    if (s %in% D) {
      temp <- temp + beta_d * gamma[ind_d, s]   
      for(k in 1:N) {
        temp[k] <- temp[k] + delta_d * sum(DW[, s] * gamma[ind_d[k], ])
      }
    }
    
    ind_intra <- which(ind_o == s & ind_d == s)
    
    if (length(ind_intra) > 0) {
      res_local <- rbind(res_local, 
              data.frame(site = s, type = "intra", values = temp[ind_intra])
      )
    }
    
    ind_origin <- which(ind_o == s & ind_d != s)
    if (length(ind_origin) > 0) {
      res_local <- rbind(res_local, 
                         data.frame(site = s, type = "origin", values = temp[ind_origin])
      )
    }
    
    ind_dest <- which(ind_d == s & ind_o != s)
    if (length(ind_dest) > 0) {
      res_local <- rbind(res_local, 
                         data.frame(site = s, type = "destination", values = temp[ind_dest])
      )
    }
    res_local <- rbind(res_local, 
                       data.frame(site = s, type = "network", 
                                  values = temp[-c(ind_intra, ind_origin, ind_dest)])
    )
  }
  
  # If s has an intra characteristic, we need to add the effect due to beta_i 
  if (beta_i != 0 & length(OD) > 0) {
    A <- diag(N) - rho_o * kronecker(OW, diag(n_d)) - rho_d * kronecker(diag(n_o), DW) - 
      rho_w * kronecker(OW, DW)
    b_mat <- matrix(0, N, length(OD))
    for(k in 1:length(OD)) {
      b_mat[which(ind_o == OD[k] & ind_d == OD[k]), k] <- 1
    }
    A_w_col <- solve(A, b_mat)
    
    for(k in 1:length(OD)) {
      s <- OD[k]
      temp <- beta_i * A_w_col[, k]

      ind_intra <- which(ind_o == s & ind_d == s)
      res_local[res_local$site == s & res_local$type == "intra", "values"] <- 
        res_local[res_local$site == s & res_local$type == "intra", "values"] + 
        temp[ind_intra]
        
      ind_origin <- which(ind_o == s & ind_d != s)
      res_local[res_local$site == s & res_local$type == "origin", "values"] <- 
        res_local[res_local$site == s & res_local$type == "origin", "values"] + 
        temp[ind_origin]
      
      ind_dest <- which(ind_d == s & ind_o != s)
      res_local[res_local$site == s & res_local$type == "destination", "values"] <- 
        res_local[res_local$site == s & res_local$type == "destination", "values"] + 
        temp[ind_dest]
      
      res_local[res_local$site == s & res_local$type == "network", "values"] <- 
        res_local[res_local$site == s & res_local$type == "network", "values"] + 
        temp[-c(ind_intra, ind_origin, ind_dest)]
    }
  }
  
  # computation of the number of local impact
  n_effect <- matrix(0, n_s, 5)
  rownames(n_effect) <- S
  colnames(n_effect) <- c("IE", "OE", "DE", "NE", "TE")
  
  # Total
  n_effect[S, "TE"] <- N
  # intra 
  n_effect[intersect(OD, intra_flows_index), "IE"] <- 1
  # origin
  n_effect[setdiff(O, D), "OE"] <- n_outgoing[setdiff(O, D)] 
  n_effect[intersect(OD, O), "OE"] <- n_outgoing[intersect(OD, O)] - 1
  # destination
  n_effect[setdiff(D, O), "DE"] <- n_ingoing[setdiff(D, O)] 
  n_effect[intersect(OD, D), "DE"] <- n_ingoing[intersect(OD, D)] - 1
  # network 
  n_effect[setdiff(O, D), "NE"] <- N - n_outgoing[setdiff(O, D)] 
  n_effect[setdiff(D, O), "NE"] <- N - n_ingoing[setdiff(D, O)] 
  n_effect[OD, "NE"] <- N - n_outgoing[OD] - n_ingoing[OD] + 1
  
  return(list(res_local = res_local, n_effect = n_effect))
}

####################################################################

predict_gsim <- function(beta_o, beta_d, beta_i, delta_o, delta_d, 
                         alpha = 0, omega_vec = rep(0, length(ind_o)),
                         ind_o, ind_d,
                         OW, DW, OX, DX,
                         rho_o, rho_d, rho_w) {
  
  # initialization 
  O <- unique(ind_o)
  n_o <- length(O)
  D <- unique(ind_d)
  n_d <- length(D)
  OD <- intersect(O, D)
  n_od <- length(OD)
  N <- length(ind_o)
  p <- length(beta_o)
  is_distance <- !all(omega_vec == 0)
  
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
  term_1 <- beta_o * lambda %*% OX[O]  
  term_2 <- beta_d * gamma %*% DX[D]  
  term_3 <- delta_o * OW %*% lambda %*% OX[O]
  term_4 <- delta_d * DW %*% gamma %*% DX[D]

  # cste effect
  if (alpha != 0) {
    lambda_cste <- apply(lambda, 1, sum)
    temp_alpha <- alpha * lambda_cste[index_o]
  } else {
    temp_alpha <- 0
  }

  # check if we need the full filter matrix 
  if ((beta_i != 0 & length(OD) > 0) | is_distance) {
    A <- diag(N) - rho_o * kronecker(OW, diag(n_d)) - rho_d * kronecker(diag(n_o), DW) - 
      rho_w * kronecker(OW, DW)

  # compute the effect due to the intra
  if (beta_i != 0 & length(OD) > 0) {
    hat_x <- rep(0, N)

    for(k in 1:length(OD)) {
      s <- OD[k]
      hat_x[which(ind_o == OD[k] & ind_d == OD[k])] <- (OX[s] * beta_i)
    }
  }

    # we add the alpha cste term
      hat_x <- hat_x + omega_vec
     # we compute the term
      temp_cste <- solve(A, hat_x)
  } else {
    temp_cste <- rep(0, N)
  }

  
  return(temp_alpha + temp_cste + term_1[ind_o, ] + term_2[ind_d, ] + 
           term_3[ind_o, ] + term_4[ind_d, ])
}


########################################################################

compute_effect_in_log <- function(beta_o, beta_d, beta_i, delta_o, delta_d, 
                                  alpha = 0, omega_vec = rep(0, ind_o),
                                  ind_o, ind_d,
                                  OW, DW, 
                                  rho_o, rho_d, rho_w,
                                  type_x = 3, change_x = 1,
                                  cste = rep(0, length(ind_o)),
                                  OX, DX) {
  
  # initialization 
  O <- unique(ind_o)
  n_o <- length(O)
  D <- unique(ind_d)
  n_d <- length(D)
  OD <- intersect(O, D)
  n_od <- length(OD)
  N <- length(ind_o)
  p <- length(beta_o)
  
  # definition of S
  if (type_x == 1) {
    S <- O 
  } else {
    if (type_x == 2) {
      S <- D
    } else {
      if (type_x == 3) {
        # ordering S with respect to O, D and O\cap D
        S <- c(setdiff(O, OD), setdiff(D, OD), OD)  
      }
    }
  }
  n_s <- length(S)
  
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
  lambda <- solve(diag(n_o) - lambda_prime) * change_x
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
  gamma <- solve(diag(n_d) - gamma_prime) * change_x
  rownames(gamma) <- colnames(gamma) <- D
  
  # computation of the local impact
  effect <- matrix(0, n_s, 5)
  rownames(effect) <- S
  colnames(effect) <- c("IE", "OE", "DE", "NE", "TE")
  effect_o <- effect_d <- effect_delta_o <- effect_delta_d <- effect
  
  # compute the number of ongoing resp. arriving flows per site
  n_outgoing <- table(ind_o)
  n_ingoing <- table(ind_d)
  
  # loop on the intersections
  intra_flows_index <- ind_o[ind_o == ind_d]
  # we detect the flows that have intra flows
  temp <- predict_gsim(beta_o = beta_o, 
                       beta_d = beta_d, 
                       beta_i = beta_i, 
                       delta_o = delta_o, 
                       delta_d = delta_d, 
                       alpha = alpha,
                       omega_vec = omega_vec,
                       ind_o = ind_o, 
                       ind_d = ind_d,
                       OW = OW, 
                       DW = DW, 
                       OX = OX, 
                       DX = DX,
                       rho_o = rho_o, 
                       rho_d = rho_d, 
                       rho_w = rho_w)
  hat_Y <- exp(temp)
  
  # If s has an intra characteristic, we need to add the effect due to beta_i 
  if (beta_i != 0 & length(OD) > 0) {
    A <- diag(N) - rho_o * kronecker(OW, diag(n_d)) - rho_d * kronecker(diag(n_o), DW) - 
      rho_w * kronecker(OW, DW)
    b_mat <- matrix(0, N, length(OD))
    for(k in 1:length(OD)) {
      b_mat[which(ind_o == OD[k] & ind_d == OD[k]), k] <- 1
    }
    A_w_col <- solve(A, b_mat)
    colnames(A_w_col) <- OD 
  }

  # compute intra effect
  for(s in intersect(OD, intra_flows_index)) {
    
    temp <- lambda[s, s] * beta_o + 
      gamma[s, s] * beta_d + 
      sum(lambda[s, ] * OW[, s]) * delta_o + 
      sum(gamma[s, ] * DW[, s]) * delta_d
    
    if (beta_i != 0 & s %in% OD) {
      temp <- temp +  beta_i * A_w_col[ind_o == s & ind_d == s, s]
    }
    
    effect[s, "IE"] <- hat_Y[ind_o == s & ind_d == s] * 
      (exp(temp) - 1)
  }
  
  # loop on the origin sites
  for(s in S[S %in% O]) {

    temp <- lambda[s, s] * beta_o + sum(lambda[s, ] * OW[, s]) *  delta_o
    
    if (s %in% D) {
      temp <- temp + gamma[ind_d[ind_o == s & ind_d != s], s] * beta_d +
        (sapply(D, function(y) sum(gamma[ind_d[ind_o == s], y])) * DW[, s] * delta_d)[setdiff(D, s)]
    }
     
    if (beta_i != 0 & s %in% OD) {
      temp <- temp +  beta_i * A_w_col[ind_o == s & ind_d != s, s]
    }
    
    effect[s, "OE"] <- sum((exp(temp) - 1) *  hat_Y[ind_o == s & ind_d != s] )
    
  }
  
  # loop on the destination sites 
  for(s in S[S %in% D]) {
    
    temp <- gamma[s, s] * beta_d + sum(gamma[s, ] * DW[, s]) * delta_d
    
    if (s %in% O) {
      temp <- temp + lambda[ind_o[ind_d == s & ind_o != s], s] * beta_o +
        (sapply(O, function(y) sum(lambda[ind_o[ind_d == s], y])) * OW[, s] * delta_o)[setdiff(O, s)]
    } 
 
    if (beta_i != 0 & s %in% OD) {
      temp <- temp +  beta_i * A_w_col[ind_d == s & ind_o != s, s]
    }
    
    effect[s, "DE"] <- sum((exp(temp) - 1) *  hat_Y[ind_d == s & ind_o != s] )
    
  }
  
  # loop on the total effects
  for(s in S) {
    temp <- rep(0, N)
    if (s %in% O)
      temp <- temp + lambda[ind_o, s] * beta_o + 
        (lambda[ind_o, ] %*% OW[, s]) * delta_o  
    
    if (s %in% D) {
      temp <- temp + gamma[ind_d, s] * beta_d +
        (gamma[ind_d, ] %*% DW[, s]) * delta_d
    }
    
    if (beta_i != 0 & s %in% OD) {
      temp <- temp +  beta_i * A_w_col[, s]
    }
    
    effect[s, "TE"] <- sum((exp(temp) - 1) *  hat_Y)
    
  }
  
  
  effect[, "NE"] <- effect[, "TE"] - effect[, "IE"] - effect[, "OE"] - effect[, "DE"]
  
  
  # computation of the number of local impact
  n_effect <- matrix(0, n_s, 5)
  rownames(n_effect) <- S
  colnames(n_effect) <- c("IE", "OE", "DE", "NE", "TE")
  
  # Total
  n_effect[S, "TE"] <- N
  # intra 
  n_effect[intersect(OD, intra_flows_index), "IE"] <- 1
  # origin
  n_effect[setdiff(S[S %in% O], D), "OE"] <- n_outgoing[setdiff(S[S %in% O], D)] 
  n_effect[intersect(OD, S[S %in% O]), "OE"] <- n_outgoing[intersect(OD, S[S %in% O])] - 1
  # destination
  n_effect[setdiff(S[S %in% D], O), "DE"] <- n_ingoing[setdiff(S[S %in% D], O)] 
  n_effect[intersect(OD, S[S %in% D]), "DE"] <- n_ingoing[intersect(OD, S[S %in% D])] - 1
  # network 
  n_effect[setdiff(S[S %in% O], D), "NE"] <- N - n_outgoing[setdiff(S[S %in% O], D)] 
  n_effect[setdiff(S[S %in% D], O), "NE"] <- N - n_ingoing[setdiff(S[S %in% D], O)] 
  n_effect[OD, "NE"] <- N - n_outgoing[OD] - n_ingoing[OD] + 1
  
  return(list(local_effect = effect, n_effect = n_effect))
}

#########################################################
# version with a loop
compute_effect_v1 <- function(beta_o, beta_d, beta_i, delta_o, delta_d, 
                              ind_names,
                              ind_o, ind_d,
                              W_o, W_d, W_w,
                              rho_o, rho_d, rho_w,
                              change_x = 1) {
  
  # internal function
  epsilon_when_change_one_unit <- function (which_unit) {
    # add 1 unit 
    ind_o <- (ind_o == which_unit) 
    ind_d <- (ind_d == which_unit) 
    ind_i <- ind_o & ind_d
    
    return(cbind(A_w %*% ind_o, A_w %*% (W_o %*% ind_o), 
                 A_w %*% ind_d, A_w %*% (W_d %*% ind_d),
                 A_w %*% ind_i))
  }
  
  # initialization 
  N <- nrow(W_o)
  p <- length(beta_o)
  # we compute A
  A_w <- solve(diag(N) - rho_o * W_o - rho_d * W_d - rho_w * W_w)
  
  OE <- matrix(0, nrow = p, ncol = 4,
               dimnames = list(names(beta_o), c("delta_i", "delta_o", "delta_d", "delta")))
  DE <- matrix(0, nrow = p, ncol = 4,
               dimnames = list(names(beta_o), c("delta_i", "delta_o", "delta_d", "delta")))
  NE <- matrix(0, nrow = p, ncol = 4,
               dimnames = list(names(beta_o), c("delta_i", "delta_o", "delta_d", "delta")))
  intra <- matrix(0, nrow = p, ncol = 4,
                  dimnames = list(names(beta_o), c("delta_i", "delta_o", "delta_d", "delta")))
  
  # actualization
  beta_o <- beta_o * change_x
  beta_d <- beta_d * change_x  
  beta_i <- beta_i * change_x  
  delta_o <- delta_o * change_x
  delta_d <- delta_d * change_x  
  
  for (k in ind_names) {
    change_Rk <- epsilon_when_change_one_unit(which_unit = as.character(k))
    
    part_o <- change_Rk[, 1]
    part_o_lagg <- change_Rk[, 2]
    part_d <- change_Rk[, 3]
    part_d_lagg <- change_Rk[, 4]    
    part_i <- change_Rk[, 5]
    
    
    # index
    ind_intra <- which(ind_o == as.character(k) & ind_d == as.character(k))
    ind_d_eff <- ind_d == as.character(k) & ind_o != as.character(k)
    ind_o_eff <- ind_o == as.character(k) & ind_d != as.character(k)
    ind_net <- ind_d != as.character(k) & ind_o != as.character(k)
    
    # intra part
    if (length(ind_intra) != 0) {
      intra <- intra + cbind(
        beta_i * part_i[ind_intra], 
        beta_o * part_o[ind_intra] + delta_o * part_o_lagg[ind_intra], 
        beta_d * part_d[ind_intra] + delta_d * part_d_lagg[ind_intra], 
        beta_i * part_i[ind_intra] +
          beta_o * part_o[ind_intra] + delta_o * part_o_lagg[ind_intra] + 
          beta_d * part_d[ind_intra] + delta_d * part_d_lagg[ind_intra])
    }
    # OE part
    oe_sum_part_i <- sum(part_i[ind_o_eff])
    oe_sum_part_o <- sum(part_o[ind_o_eff])
    oe_sum_part_d <- sum(part_d[ind_o_eff])
    oe_sum_part_o_lag <- sum(part_o_lagg[ind_o_eff])
    oe_sum_part_d_lag <- sum(part_d_lagg[ind_o_eff])
    
    OE <- OE + cbind(beta_i * oe_sum_part_i,
                     beta_o * oe_sum_part_o + delta_o * oe_sum_part_o_lag, 
                     beta_d * oe_sum_part_d + delta_d * oe_sum_part_d_lag, 
                     beta_i * oe_sum_part_i +
                       beta_o * oe_sum_part_o + delta_o * oe_sum_part_o_lag + 
                       beta_d * oe_sum_part_d + delta_d * oe_sum_part_d_lag)
    
    # DE part
    de_sum_part_i <- sum(part_i[ind_d_eff])
    de_sum_part_o <- sum(part_o[ind_d_eff])
    de_sum_part_d <- sum(part_d[ind_d_eff])
    de_sum_part_o_lag <- sum(part_o_lagg[ind_d_eff])
    de_sum_part_d_lag <- sum(part_d_lagg[ind_d_eff])
    
    DE <- DE + cbind(
      beta_i * de_sum_part_i, 
      beta_o * de_sum_part_o + delta_o * de_sum_part_o_lag, 
      beta_d * de_sum_part_d + delta_d * de_sum_part_d_lag, 
      beta_i * de_sum_part_i +
        beta_o * de_sum_part_o + delta_o * de_sum_part_o_lag + 
        beta_d * de_sum_part_d + delta_d * de_sum_part_d_lag)
    
    # NE part
    ne_sum_part_i <- sum(part_i[ind_net])
    ne_sum_part_o <- sum(part_o[ind_net])
    ne_sum_part_d <- sum(part_d[ind_net])
    ne_sum_part_o_lag <- sum(part_o_lagg[ind_net])
    ne_sum_part_d_lag <- sum(part_d_lagg[ind_net])
    
    NE <- NE + cbind(beta_i * ne_sum_part_i,
                     beta_o * ne_sum_part_o + delta_o * ne_sum_part_o_lag, 
                     beta_d * ne_sum_part_d + delta_d * ne_sum_part_d_lag, 
                     beta_i * ne_sum_part_i +
                       beta_o * ne_sum_part_o + delta_o * ne_sum_part_o_lag + 
                       beta_d * ne_sum_part_d + delta_d * ne_sum_part_d_lag)
  }
  
  res <- rbind(intra / N, OE / N, DE / N, NE / N, (intra + OE + DE + NE) / N)
  rownames(res) <- c("intra", "OE", "DE", "NE", "total")
  return(res)
}


#########################################################################

# version with matrix computation
compute_effect_v2 <- function(beta_o, beta_d, delta_o, delta_d,
                              ind_names,
                              ind_o, ind_d,
                              W_o, W_d, W_w,
                              rho_o, rho_d, rho_w,
                              change_x = 1) {

  # initialization 
  N <- nrow(W_o)
  p <- length(beta_o)

  # data.frame 
  library(Matrix)
  
  
  o_sparse <- sparse.model.matrix(~.-1, 
        data = data.frame(ind_o = factor(ind_o, levels = ind_names)))
  d_sparse <- sparse.model.matrix(~.-1, 
        data = data.frame(ind_d = factor(ind_d, levels = ind_names)))
  
  # we compute A
  A_w_ind_o <- solve(diag(N) - rho_o * W_o - rho_d * W_d - rho_w * W_w,
                     o_sparse)
  A_w_Wo_ind_o <- solve(diag(N) - rho_o * W_o - rho_d * W_d - rho_w * W_w,
                        W_o %*% o_sparse)
  A_w_ind_d <- solve(diag(N) - rho_o * W_o - rho_d * W_d - rho_w * W_w,
                     d_sparse)
  A_w_Wd_ind_d <- solve(diag(N) - rho_o * W_o - rho_d * W_d - rho_w * W_w,
                        W_d %*% d_sparse)
  
  # we prepare the index for affecting the O, D, N
  intra_sparse <- (o_sparse + d_sparse) / 2 * o_sparse * d_sparse
  o_sparse <- o_sparse - intra_sparse
  d_sparse <- d_sparse - intra_sparse
  n_sparse <- 1 - (intra_sparse + o_sparse + d_sparse)
  
  # actualization
  beta_o <- beta_o * change_x
  beta_d <- beta_d * change_x  
  delta_o <- delta_o * change_x
  delta_d <- delta_d * change_x  
  
  IE_due_to_betao <- sum(A_w_ind_o * intra_sparse) * beta_o + sum(A_w_Wo_ind_o * intra_sparse) * delta_o
  IE_due_to_betad <- sum(A_w_ind_d * intra_sparse) * beta_d + sum(A_w_Wd_ind_d * intra_sparse) * delta_d
  
  OE_due_to_betao <- sum(A_w_ind_o * o_sparse) * beta_o + sum(A_w_Wo_ind_o * o_sparse) * delta_o
  OE_due_to_betad <- sum(A_w_ind_d * o_sparse) * beta_d + sum(A_w_Wd_ind_d * o_sparse) * delta_d
  
  DE_due_to_betao <- sum(A_w_ind_o * d_sparse) * beta_o + sum(A_w_Wo_ind_o * d_sparse) * delta_o
  DE_due_to_betad <- sum(A_w_ind_d * d_sparse) * beta_d + sum(A_w_Wd_ind_d * d_sparse) * delta_d
  
  NE_due_to_betao <- sum(A_w_ind_o * n_sparse) * beta_o + sum(A_w_Wo_ind_o * n_sparse) * delta_o
  NE_due_to_betad <- sum(A_w_ind_d * n_sparse) * beta_d + sum(A_w_Wd_ind_d * n_sparse) * delta_d
  
  res <- list(
    intra = cbind(IE_due_to_betao, IE_due_to_betad, IE_due_to_betao + IE_due_to_betad) / N,
    OE = cbind(OE_due_to_betao, OE_due_to_betad, OE_due_to_betao + OE_due_to_betad) / N,
    DE = cbind(DE_due_to_betao, DE_due_to_betad, DE_due_to_betao + DE_due_to_betad) / N,
    NE = cbind(NE_due_to_betao, NE_due_to_betad, NE_due_to_betao + NE_due_to_betad) / N )
  res$Total <- res$intra + res$OE + res$DE + res$NE
  return(res)
}


# ###################
# ### Test for 
# # We prepare the data :
# # indexes 
# O <- c(paste0("o", 1:4), "od1")
# D <- c(paste0("d", 1:2), "od1")
# # site <- union(o, d)
# n_o <- length(O)
# n_d <- length(D)
# #n <- length(site)
# N <- n_o * n_d
# 
# ind_o <-  rep(O, each = n_d)
# ind_d <- rep(D, times = n_o)
# 
# # spatial weigh matrices
# OW <- matrix(c(0, 1, 0, 0, 0, 
#                1/2, 0, 0, 0, 1/2, 
#                0, 0, 0, 1/2, 1/2, 
#                0, 0, 1, 0, 0, 
#                0, 1/2, 1/2, 0, 0), 
#              n_o, n_o, byrow = T)
# rownames(OW) <- colnames(OW) <- O
# 
# DW <- matrix(c(0, 0, 1, 
#                0, 0, 1, 
#                0.5, 0.5, 0), 
#              n_d, n_d, byrow = T)
# rownames(DW) <- colnames(DW) <- D 
# 
# W_d <- kronecker(diag(n_o), DW)
# W_o <- kronecker(OW, diag(n_d))
# W_w <- kronecker(OW, DW)
# 
# # x matrix
# x <- matrix(c(40, 20, 10, 7, 10, 15, 25), ncol = 1)
# names(x) <- site
# # rho parameters
# rho_d <- 0.4
# rho_o <- 0.4
# rho_w <- -rho_d * rho_o
# 
# # beta parameters
# beta_o <- 0.5
# beta_d <- 1
# 
# A_w <- diag(N) - rho_o * W_o -rho_d * W_d - rho_w * W_w
# # preparation of the grid
# # option 1
# #my_spat <- data.frame(long = rep(1:n_d, each = n_o), 
# #                      lat = rep(1:n_o, times = n_d))
# #my_spat <- st_as_sf(my_spat, coords = c("long", "lat"))
# ## <- my_spat %>% 
# #  st_make_grid(cellsize = c(1, 1), 
# #               n = c(n_d, n_o)) 
# #my_grid <- my_grid %>% 
# #st_sf() 
# 
# 
# 
# compute_effect_v1(beta_o = 0.5, beta_d = 1, 
#                   delta_o = 0, delta_d = 0,
#                   ind_names = S,
#                   ind_o = ind_o, ind_d = ind_d,
#                   W_o = W_o, W_d = W_d, W_w = W_w,
#                   rho_o = 0.4, rho_d = 0.4, rho_w = -0.16)
# 
# compute_effect_v3(beta_o, beta_d, delta_o, delta_d,
#                    ind_o, ind_d, 
#                   OW, DW,
#                   rho_o, rho_d, rho_w,
#                   change_x = 1) 
# ###  Test for computation time 
#  n <- 50
#  N <- n ^ 2
# # # indexes 
#  ind_o <- rep(as.character(1:n), each = n)
#  ind_d <- rep(as.character(1:n), times = n)
# # 
# # # spatial weigh matrices
#  w <- matrix(sample(c(0,1), size = N, prob = c(0.8, 0.2), replace = T), 
#              n, n, byrow = T)
#  diag(w) <- 0
#  # row normalize
#  w <- w / matrix(apply(w, 1, sum), n, n)
#  
# # 
# W_d <- kronecker(diag(n), w)
# W_o <- kronecker(w, diag(n))
# W_w <- kronecker(w, w)
# # 
# # # x matrix
#  x <- matrix(rnorm(n), ncol = 1)
# # 
# # # rho parameters
#  rho_d <- 0.4
#  rho_o <- 0.4
#  rho_w <- -rho_d * rho_o
# # 
# # # beta parameters
#  beta_o <- 0.5
#  beta_d <- 1
#  gamma_od <- -0.5
# # # test
#  x_o <- kronecker(x, rep(1, n))
#  x_d <- kronecker(rep(1, n), x)
#  Y <- x_o * beta_o + x_d * beta_d 
# 
#  # comparison result
#  compute_effect_v1(beta_o = 0.5, beta_d = 1, 
#                    delta_o = 0, delta_d = 0,
#                    ind_names = as.character(1:n),
#                    ind_o = ind_o, ind_d = ind_d,
#                    W_o = W_o, W_d = W_d, W_w = W_w,
#                    rho_o = 0, rho_d = 0, rho_w = 0)
#  
#  compute_effect_v2(beta_o = 0.5, beta_d = 1, 
#                    delta_o = 0, delta_d = 0,
#                    ind_names = as.character(1:n),
#                    ind_o = ind_o, ind_d = ind_d,
#                    W_o = W_o, W_d = W_d, W_w = W_w,
#                    rho_o = 0, rho_d = 0, rho_w = 0)
#  
#  rownames(w) <- colnames(w) <- as.character(1:n)
#    
#  compute_effect_v3(beta_o, beta_d, delta_o, delta_d,
#                    ind_o, ind_d, 
#                    w, w,
#                    rho_o, rho_d, rho_w,
#                    change_x = 1)
#  # comparison time
#  microbenchmark::microbenchmark(compute_effect_v1(beta_o = 0.5, beta_d = 1, 
#                                                   delta_o = 0, delta_d = 0,
#                                                   ind_names = S,
#                                                   ind_o = ind_o, ind_d = ind_d,
#                                                   W_o = W_o, W_d = W_d, W_w = W_w,
#                                                   rho_o = 0.4, rho_d = 0.4, rho_w = -0.16), 
#                                 compute_effect_v2(beta_o = 0.5, beta_d = 1, 
#                                        delta_o = 0, delta_d = 0,
#                                        ind_names = as.character(1:n),
#                                         ind_o = ind_o, ind_d = ind_d,
#                                        W_o = W_o, W_d = W_d, W_w = W_w,
#                                       rho_o = 0, rho_d = 0, rho_w = 0), 
#                                 compute_effect_v3(beta_o, beta_d, delta_o, delta_d,
#                                                   ind_o, ind_d, 
#                                                   w, w,
#                                                   rho_o, rho_d, rho_w,
#                                                   change_x = 1),
#                                 times = 10L)
# 
#  
# 
# compute_effect_v1 <- function(beta_o, beta_d, beta_i, delta_o, delta_d, delta_i,
#                               ind_names,
#                               ind_o, ind_d,
#                               W_o, W_d, W_w,
#                               rho_o, rho_d, rho_w,
#                               change_x = 1) {
#   
#   # internal function
#   epsilon_when_change_one_unit <- function (which_unit) {
#     # add 1 unit 
#     ind_o <- (ind_o == which_unit) 
#     ind_d <- (ind_d == which_unit) 
#     ind_i <- ind_o & ind_d
#     
#     return(cbind(A_w %*% ind_o, A_w %*% (W_o %*% ind_o), 
#                  A_w %*% ind_d, A_w %*% (W_d %*% ind_d),
#                  A_w %*% ind_i, A_w %*% ((W_o + W_d) / 2) %*% ind_i
#     ))
#   }
#   
#   # initialization 
#   N <- nrow(W_o)
#   p <- length(beta_o)
#   # we compute A
#   A_w <- solve(diag(N) - rho_o * W_o - rho_d * W_d - rho_w * W_w)
#   
#   OE <- matrix(0, nrow = p, ncol = 4,
#                dimnames = list(names(beta_o), c("delta_i", "delta_o", "delta_d", "delta")))
#   DE <- matrix(0, nrow = p, ncol = 4,
#                dimnames = list(names(beta_o), c("delta_i", "delta_o", "delta_d", "delta")))
#   NE <- matrix(0, nrow = p, ncol = 4,
#                dimnames = list(names(beta_o), c("delta_i", "delta_o", "delta_d", "delta")))
#   intra <- matrix(0, nrow = p, ncol = 4,
#                   dimnames = list(names(beta_o), c("delta_i", "delta_o", "delta_d", "delta")))
#   
#   # actualization
#   beta_o <- beta_o * change_x
#   beta_d <- beta_d * change_x  
#   beta_i <- beta_i * change_x  
#   delta_o <- delta_o * change_x
#   delta_d <- delta_d * change_x  
#   delta_i <- delta_i * change_x 
#   
#   for (k in ind_names) {
#     change_Rk <- epsilon_when_change_one_unit(which_unit = as.character(k))
#     
#     part_o <- change_Rk[, 1]
#     part_o_lagg <- change_Rk[, 2]
#     part_d <- change_Rk[, 3]
#     part_d_lagg <- change_Rk[, 4]    
#     part_i <- change_Rk[, 5]
#     part_i_lagg <- change_Rk[, 6]    
#     
#     # index
#     ind_intra <- which(ind_o == as.character(k) & ind_d == as.character(k))
#     ind_d_eff <- ind_d == as.character(k) & ind_o != as.character(k)
#     ind_o_eff <- ind_o == as.character(k) & ind_d != as.character(k)
#     ind_net <- ind_d != as.character(k) & ind_o != as.character(k)
#     
#     # intra part
#     if (length(ind_intra) != 0) {
#       intra <- intra + cbind(
#         beta_i * part_i[ind_intra] + delta_i * part_i_lagg[ind_intra], 
#         beta_o * part_o[ind_intra] + delta_o * part_o_lagg[ind_intra], 
#         beta_d * part_d[ind_intra] + delta_d * part_d_lagg[ind_intra], 
#         beta_i * part_i[ind_intra] + delta_i * part_i_lagg[ind_intra] +
#           beta_o * part_o[ind_intra] + delta_o * part_o_lagg[ind_intra] + 
#           beta_d * part_d[ind_intra] + delta_d * part_d_lagg[ind_intra])
#     }
#     # OE part
#     oe_sum_part_i <- sum(part_i[ind_o_eff])
#     oe_sum_part_o <- sum(part_o[ind_o_eff])
#     oe_sum_part_d <- sum(part_d[ind_o_eff])
#     oe_sum_part_i_lag <- sum(part_i_lagg[ind_o_eff])
#     oe_sum_part_o_lag <- sum(part_o_lagg[ind_o_eff])
#     oe_sum_part_d_lag <- sum(part_d_lagg[ind_o_eff])
#     
#     OE <- OE + cbind(beta_i * oe_sum_part_i + delta_i * oe_sum_part_i_lag,
#                      beta_o * oe_sum_part_o + delta_o * oe_sum_part_o_lag, 
#                      beta_d * oe_sum_part_d + delta_d * oe_sum_part_d_lag, 
#                      beta_i * oe_sum_part_i + delta_i * oe_sum_part_i_lag +
#                        beta_o * oe_sum_part_o + delta_o * oe_sum_part_o_lag + 
#                        beta_d * oe_sum_part_d + delta_d * oe_sum_part_d_lag)
#     
#     # DE part
#     de_sum_part_i <- sum(part_i[ind_d_eff])
#     de_sum_part_o <- sum(part_o[ind_d_eff])
#     de_sum_part_d <- sum(part_d[ind_d_eff])
#     de_sum_part_i_lag <- sum(part_i_lagg[ind_d_eff])
#     de_sum_part_o_lag <- sum(part_o_lagg[ind_d_eff])
#     de_sum_part_d_lag <- sum(part_d_lagg[ind_d_eff])
#     
#     DE <- DE + cbind(
#       beta_i * de_sum_part_i + delta_i * de_sum_part_i_lag, 
#       beta_o * de_sum_part_o + delta_o * de_sum_part_o_lag, 
#       beta_d * de_sum_part_d + delta_d * de_sum_part_d_lag, 
#       beta_i * de_sum_part_i + delta_i * de_sum_part_i_lag +
#         beta_o * de_sum_part_o + delta_o * de_sum_part_o_lag + 
#         beta_d * de_sum_part_d + delta_d * de_sum_part_d_lag)
#     
#     # NE part
#     ne_sum_part_i <- sum(part_i[ind_net])
#     ne_sum_part_o <- sum(part_o[ind_net])
#     ne_sum_part_d <- sum(part_d[ind_net])
#     ne_sum_part_i_lag <- sum(part_i_lagg[ind_net])
#     ne_sum_part_o_lag <- sum(part_o_lagg[ind_net])
#     ne_sum_part_d_lag <- sum(part_d_lagg[ind_net])
#     
#     NE <- NE + cbind(beta_i * ne_sum_part_i + delta_i * ne_sum_part_i_lag,
#                      beta_o * ne_sum_part_o + delta_o * ne_sum_part_o_lag, 
#                      beta_d * ne_sum_part_d + delta_d * ne_sum_part_d_lag, 
#                      beta_i * ne_sum_part_i + delta_i * ne_sum_part_i_lag +
#                        beta_o * ne_sum_part_o + delta_o * ne_sum_part_o_lag + 
#                        beta_d * ne_sum_part_d + delta_d * ne_sum_part_d_lag)
#   }
#   
#   res <- rbind(intra / N, OE / N, DE / N, NE / N, (intra + OE + DE + NE) / N)
#   rownames(res) <- c("intra", "OE", "DE", "NE", "total")
#   return(res)
# }
