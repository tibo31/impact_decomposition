####################################################################

predict_gsim_nc <- function(Z, rho, alpha = 0, omega_vec = rep(0, length(ind_o)),
                         W_o, W_d, W_w, OX, DX,
                         ind_o, ind_d) {
  
  # Z is a DF of size p x 5 (names, beta_o, beta_d, delta_o, delta_d, beta_i)
  p <- nrow(Z)
  names_p <- row.names(Z)
  
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
  
  # we compute X
  temp <- alpha + omega_vec
  names_o <- names_p[Z[, "beta_o"] != 0]
  names_d <- names_p[Z[, "beta_d"] != 0]  
  names_i <- names_p[Z[, "beta_i"] != 0]
  names_do <- names_p[Z[, "delta_o"] != 0]
  names_dd <- names_p[Z[, "delta_d"] != 0]
  if(length(names_o) > 0) {
    temp <- temp + as(OX[ind_o, names_o], "matrix") %*% Z[names_o, "beta_o"]
  }  
  if(length(names_d) > 0) {
    temp <- temp + as(DX[ind_d, names_d], "matrix") %*% Z[names_d, "beta_d"]
  } 
  if(length(names_do) > 0) {
    temp <- temp + W_o %*% as(OX[ind_o, names_o], "matrix") %*% Z[names_o, "delta_o"]
  }  
  if(length(names_dd) > 0) {
    temp <- temp + W_d %*% as(DX[ind_d, names_d], "matrix") %*% Z[names_d, "delta_d"]
  } 
  if(length(names_i) > 0) {
    id_intra <- ind_o == ind_d
    temp_intra <- ifelse(id_intra, as(OX[ind_o, names_o], "matrix") %*% Z[names_o, "beta_i"], 0) 
    temp <- temp + temp_intra
  } 
  # check if we need the full filter matrix 
    res <- solve(A, temp)
  
  return(res)
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
