
#########################################################
# version with a loop

compute_effect_nc <- function(Z, rho, W_o, W_d, W_w, ind_o, ind_d,
                              change_z = NULL,
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
  N <- length(ind_o)
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
  effect_beta_i <- effect_d <- effect_delta_o <- effect_delta_d <- effect_o
  # computation of the number of local impact
  n_effect <- effect_o
  
  # we detect the flows that have intra flows 
  intra_flows_index <- intersect(OD, intersect(ind_o[ind_o == ind_d], S))
  
  # compute the number of ongoing resp. arriving flows per site
  # index
  for (k in 1:n_s){
    s <- S[k]
    ind_intra <- which(ind_o == s & ind_d == s)
    ind_d_eff <- which(ind_d == s & ind_o != s)
    ind_o_eff <- which(ind_o == s & ind_d != s)
    ind_net <- which(ind_d != s & ind_o != s)
    n_effect[s, ] <- c(length(ind_intra), length(ind_o_eff), length(ind_d_eff), length(ind_net), 
                       length(ind_intra) + length(ind_o_eff) + length(ind_d_eff) + length(ind_net))
  }
  # we compute A
  A <- diag(N) 
  if(!is.null(W_o))
    A <- A - rho_o * W_o
  if(!is.null(W_d))
    A <- A -  rho_d * W_d 
  if(!is.null(W_w))
    A <- A - rho_w * W_w
  
  if(!require("Matrix")) {
    install.packages("Matrix")
    require("Matrix")
  }
  if(class(A)[1] != "Matrix")
    A <- as(A, "dgeMatrix")
 
  if (!Y_in_log) {
    for (k in 1:n_s){
      s <- S[k]
      # add 1 unit 
      change_o <- (ind_o == s) 
      change_d <- (ind_d == s) 
      change_i <- change_o & change_d
    
      # index
      ind_intra <- which(ind_o == s & ind_d == s)
      ind_d_eff <- which(ind_d == s & ind_o != s)
      ind_o_eff <- which(ind_o == s & ind_d != s)
      ind_net <- which(ind_d != s & ind_o != s)
    
      # compute the matricial product
      A_w <- solve(A, cbind(change_o, change_d, (W_o %*% change_o), (W_d %*% change_d), change_i))
      temp_1 <- A_w[, 1] # A_w %*% change_o # A_w[, 1] #
      temp_2 <- A_w[, 2] # A_w %*% change_d # A_w[, 2] # 
      temp_3 <- A_w[, 3] # A_w %*% (W_o %*% change_o) # A_w[, 3] #
      temp_4 <- A_w[, 4] #A_w %*% (W_d %*% change_d) # A_w[, 4] #
      temp_5 <- A_w[, 5] #A_w %*% change_i # A_w[, 5] #
    
      if(length(ind_intra) > 0) {
        effect_o[s, "IE"] <- temp_1[ind_intra]
        effect_d[s, "IE"] <- temp_2[ind_intra]
        effect_delta_o[s, "IE"] <- temp_3[ind_intra]
        effect_delta_d[s, "IE"] <- temp_4[ind_intra]
        effect_beta_i[s, "IE"] <- temp_5[ind_intra]
      }
    
      if(length(ind_o_eff) > 0) {
        effect_o[s, "OE"] <- sum(temp_1[ind_o_eff])
        effect_d[s, "OE"] <- sum(temp_2[ind_o_eff])
        effect_delta_o[s, "OE"] <- sum(temp_3[ind_o_eff])
        effect_delta_d[s, "OE"] <- sum(temp_4[ind_o_eff])
        effect_beta_i[s, "OE"] <- sum(temp_5[ind_o_eff])
      }
    
      if(length(ind_d_eff) > 0) {
        effect_o[s, "DE"] <- sum(temp_1[ind_d_eff])
        effect_d[s, "DE"] <- sum(temp_2[ind_d_eff])
        effect_delta_o[s, "DE"] <- sum(temp_3[ind_d_eff])
        effect_delta_d[s, "DE"] <- sum(temp_4[ind_d_eff])
        effect_beta_i[s, "DE"] <- sum(temp_5[ind_d_eff])
      }
      
      if(length(ind_net) > 0) {
        effect_o[s, "NE"] <- sum(temp_1[ind_net])
        effect_d[s, "NE"] <- sum(temp_2[ind_net])
        effect_delta_o[s, "NE"] <- sum(temp_3[ind_net])
        effect_delta_d[s, "NE"] <- sum(temp_4[ind_net])
        effect_beta_i[s, "NE"] <- sum(temp_5[ind_net])
      }
      
      effect_o[s, "TE"] <- sum(temp_1)
      effect_d[s, "TE"] <- sum(temp_2)
      effect_delta_o[s, "TE"] <- sum(temp_3)
      effect_delta_d[s, "TE"] <- sum(temp_4)
      effect_beta_i[s, "TE"] <- sum(temp_5)
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
      
      n_s <- length(S)
      effect <- matrix(0, n_s, 5)
      rownames(effect) <- S
      colnames(effect) <- c("IE", "OE", "DE", "NE", "TE")
    
      effect <- Z[k, "beta_o"] * effect_o + Z[k, "beta_d"] * effect_d + 
        Z[k, "delta_o"] * effect_delta_o + Z[k, "delta_d"] * effect_delta_d +
        Z[k, "beta_i"] * effect_beta_i
    
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
        names_o <- Z[names_k, "beta_o"] != 0
        names_d <- Z[names_k, "beta_d"] != 0  
        names_i <- Z[names_k, "beta_i"] != 0
        names_do <- Z[names_k, "delta_o"] != 0
        names_dd <- Z[names_k, "delta_d"] != 0
        
        for (j in 1:n_s){
          s <- S[j]
          temp <- rep(0, N)
          if (names_o) {
            temp <- temp +  (ind_o == s) * Z[names_k, "beta_o"]
          }  
          if (names_d) {
            temp <- temp + (ind_d == s) * Z[names_k, "beta_d"]
          } 
          if (names_do) {
            temp <- temp + W_o %*% (ind_o == s) * Z[names_k, "delta_o"]
          }  
          if (names_dd) {
            temp <- temp + W_d %*% (ind_d == s) * Z[names_k, "delta_d"]
          } 
          if (names_i) {
            id_intra <- (ind_o == s) & (ind_o == ind_d)
            temp_intra <- id_intra * Z[names_k, "beta_i"] 
            temp <- temp + temp_intra
          } 
          # check if we need the full filter matrix 
          temp <- solve(A, temp * change_z[k])
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
