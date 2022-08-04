####################################################################

predict_gsim <- function(Z, rho, alpha = 0, omega_vec = rep(0, length(ind_o)),
                         OW, DW, OX, DX,
                         ind_o, ind_d) {
  
  # Z is a DF of size p x 5 (names, beta_o, beta_d, delta_o, delta_d, beta_i)
  p <- nrow(Z)
  names_p <- row.names(Z)
  names_o <- names_p[Z[, "beta_o"] != 0]
  names_d <- names_p[Z[, "beta_d"] != 0]  
  names_i <- names_p[Z[, "beta_i"] != 0]
  names_do <- names_p[Z[, "delta_o"] != 0]
  names_dd <- names_p[Z[, "delta_d"] != 0]
  
  # initialization 
  O <- unique(ind_o)
  n_o <- length(O)
  D <- unique(ind_d)
  n_d <- length(D)
  OD <- intersect(O, D)
  n_od <- length(OD)
  N <- length(ind_o)

  is_distance <- !all(omega_vec == 0)
  
  rho_o <- rho[1] 
  rho_d <- rho[2] 
  rho_w <- rho[3]
  
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
  
  # cste effect
  if (alpha != 0) {
    lambda_cste <- apply(lambda, 1, sum)
    temp <- alpha * lambda_cste[ind_o]
  } else {
    temp <- 0
  }
  
  # check if we need the full filter matrix 
  if ((any(Z[, "beta_i"] != 0)  & length(OD) > 0) | is_distance) {
    A <- diag(N) - rho_o * kronecker(OW, diag(n_d)) - rho_d * kronecker(diag(n_o), DW) - 
      rho_w * kronecker(OW, DW)

    # we add the alpha cste term
    temp_2 <- omega_vec
    
    # compute the effect due to the intra
    if (any(Z[, "beta_i"] != 0) & length(OD) > 0) {
      hat_x <- rep(0, N)
      
      for(k in 1:length(OD)) {
        s <- OD[k]
        hat_x[which(ind_o == OD[k] & ind_d == OD[k])] <- sum(as(OX[s, names_i], "matrix") * Z[names_i, "beta_i"])
      }
      temp_2 <- temp_2 + hat_x
    }
    

    # we compute the term
    temp <- temp + solve(A, temp_2)
  } 
  
  if (length(names_o) > 0) {
    term_1 <- as(OX[O, names_o], "matrix") %*% as(Z[names_o, "beta_o"] , "matrix")
    term_1 <- lambda %*% term_1  
    temp <- temp + term_1[ind_o, ]
  } 
  
  if (length(names_d) > 0) {
    term_2 <- as(DX[D, names_d], "matrix") %*% as(Z[names_d, "beta_d"] , "matrix")
    term_2 <- gamma %*% term_2  
    temp <- temp + term_2[ind_d, ] 
  }
  
  if (length(names_do) > 0) {
    term_3 <- as(OX[O, names_do], "matrix") %*% as(Z[names_do, "delta_o"] , "matrix")
    term_3 <- OW %*% lambda %*% term_3
    temp <- temp + term_3[ind_o, ]
  }
  
  if (length(names_dd) > 0) { 
    term_4 <- as(DX[D, names_dd], "matrix") %*% as(Z[names_dd, "delta_d"] , "matrix")
    term_4 <- DW %*% gamma %*% term_4
    temp <- temp + term_4[ind_d, ]
  }
  
  
  return(temp)
}

# #############################
# # compare with
# W_d <- kronecker(diag(n), w)
# W_o <- kronecker(w, diag(n))
# W_w <- kronecker(w, w)
# x_o <- kronecker(as(X, "matrix"), rep(1, n))
# x_d <- kronecker(rep(1, n), as(X, "matrix"))
# x_i <- x_o * as.numeric(index_o == index_d)
# cbind(
# solve(diag(N) - rho_o * W_o - rho_d * W_d - rho_w * W_w, 
#       alpha + x_o * beta_o + x_d * beta_d + x_i * beta_i + 
#         as.numeric(W_o %*% x_o) * delta_o +
#         as.numeric(W_d %*% x_d) * delta_d +
#         gamma_od * g),
# #
# predict_gsim(coeff_x, rho, alpha = alpha, omega_vec = gamma_od * g,
#                          OW = w, DW = w, OX = X, DX = X,
#                          ind_o = index_o, ind_d = index_d)
# )
