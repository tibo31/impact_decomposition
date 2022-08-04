
#########################################################
#####  Compute all effects

compute_all_effect_nc <- function(Z, rho, W_o, W_d, W_w, 
                                  ind_o, ind_d, change_z = NULL,
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

  # computation of the local impact
  res_local <- data.frame(site = character(), type = character(), 
                          var = character(),
                          values = numeric())
  
  for(k in 1:p) {
      names_k <- row.names(Z)[k]
      
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
      names_o <- Z[k, "beta_o"] != 0
      names_d <- Z[k, "beta_d"] != 0  
      names_i <- Z[k, "beta_i"] != 0
      names_do <- Z[k, "delta_o"] != 0
      names_dd <- Z[k, "delta_d"] != 0
      
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
        if (Y_in_log) {
          temp <- hat_y * (exp(temp) - 1)
        }
        # index
        ind_intra <- which(ind_o == s & ind_d == s)
        ind_d_eff <- which(ind_d == s & ind_o != s)
        ind_o_eff <- which(ind_o == s & ind_d != s)
        ind_net <- which(ind_d != s & ind_o != s)
        
        if (length(ind_intra) > 0) {
          res_local <- rbind(res_local, 
                             data.frame(site = s, type = "intra", 
                                        var = names_p[k], values = temp[ind_intra])
          )
        }
        
        if (length(ind_o_eff) > 0) {
          res_local <- rbind(res_local, 
                             data.frame(site = s, type = "origin", 
                                        var = names_p[k], values = temp[ind_o_eff])
          )
        }
        
        if (length(ind_d_eff) > 0) {
          res_local <- rbind(res_local, 
                             data.frame(site = s, type = "destination", 
                                        var = names_p[k], values = temp[ind_d_eff])
          )
        }
        if (length(ind_net) > 0) {
        res_local <- rbind(res_local, 
                           data.frame(site = s, type = "network", 
                                      var = names_p[k], 
                                      values = temp[ind_net])
        )
        }
      }
    }
  
  return(res_local)
}
####################################################################
