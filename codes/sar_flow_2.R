sar_flow_2 <- function(x, y, g, W_d = NULL, W_o = NULL, W_w = NULL, 
                       model = "", ndraw = 5500L, nomit = 2500L) {
  
  # x, a data_frame or a matrix with explanatory variable of full size
  # Y, the vector of flows
  # G, the vector of distances
  # model, the model choosen
  
  # initialization
  n <- nrow(x)
  
  # we add the distance variable to x
  x_matrix <- as(x, "matrix")
  x_matrix <- cbind(x_matrix, g)
  
  # we add the constant if not included
  if(! any(apply(x_matrix, 2, function(x) all(x == 1)))) {
    x_matrix <- cbind(1, x_matrix)
    colnames(x_matrix)[1] <- "(intercept)"
  }

  nvars <- ncol(x_matrix)  # number of x + distance
  names_x <- colnames(x_matrix)
  
  # verification
  stopifnot(n == length(y), n == length(g))
  
  # determine the good model 
  if (is.null(W_d) & is.null(W_o) & is.null(W_w)) {
      return(lm(y ~ x + g))
  }
  
  if (!is.null(W_d) & is.null(W_o) & is.null(W_w)) {
    names_rho <- "rho_d"
    y2 <- as.numeric(W_d %*% y)
    zpzty2 <- crossprod(x_matrix, y2) 
    W_2 <- as(W_d, "dgCMatrix")
  }
  
  if (is.null(W_d) & !is.null(W_o) & is.null(W_w)) {
    names_rho <- "rho_o" 
    y2 <- as.numeric(W_o %*% y)
    zpzty2 <- crossprod(x_matrix, y2)
    W_2 <- as(W_o, "dgCMatrix")
  }
  
  if (is.null(W_d) & is.null(W_o) & !is.null(W_w)) {
    names_rho <- "rho_w"  
    y2 <- W_w %*% y 
    zpzty2 <- crossprod(x_matrix, y2)  
    W_2 <- as(W_w, "dgCMatrix")
  }
  
  if (!is.null(W_d) & !is.null(W_o) & is.null(W_w)) {
    if (model != "model_5") {
      names_rho <- c("rho_d", "rho_o")
      y2 <- as.numeric(W_d %*% y)
      zpzty2 <- crossprod(x_matrix, y2)
      W_2 <- as(W_d, "dgCMatrix")
      y3 <- as.numeric(W_o %*% y)
      zpzty3 <- crossprod(x_matrix, y3)
      W_3 <- as(W_o, "dgCMatrix")
      model <- "model_7"
      } else {
        names_rho <- c("rho_od")
        W_2 <- as(0.5 * (W_d + W_o), "dgCMatrix")
        y2 <- as.numeric(W_2 %*% y)
        zpzty2 <- crossprod(x_matrix, y2)
      }
  }
  
  if (!is.null(W_d) & !is.null(W_o) & !is.null(W_w)) {
    if (model == "model_6") {
      names_rho <- c("rho_odw")
      W_2 <- as(1 / 3 * (W_d + W_o + W_w), "dgCMatrix")
      y2 <- as.numeric(W_2 %*% y)
      zpzty2 <- crossprod(x_matrix, y2)
    } else {
      names_rho <- c("rho_d", "rho_o", "rho_w")
      y2 <- W_d %*% y
      zpzty2 <- crossprod(x_matrix, y2) 
      W_2 <- as(W_d, "dgCMatrix")
      y3 <- as.numeric(W_o %*% y)
      zpzty3 <- crossprod(x_matrix, y3)
      W_3 <- as(W_o, "dgCMatrix")
      y4 <- as.numeric(W_w %*% y) 
      zpzty4 <- crossprod(x_matrix, y4) 
      W_4 <- as(W_w, "dgCMatrix")
    }
  }
  
  # number of rho 
  nb_rho <- length(names_rho)
  
  if (nb_rho == 1) {
    m <- 20  # parameter fixed
    traces <- ftrace1(W_2, miter = m)
  }
  
  if (model == "model_8") {
    m <- 20  # parameter fixed
    traces_1 <- ftrace1(W_2, miter = m)
    traces_2 <- ftrace1(W_3, miter = m)    
    rho[3] <- - rho[1] * rho[2]
  }
  
  # initialization of rho
  pvec <- runif(nb_rho)
  rho <- 0.7 * pvec / sum(pvec) 
  
  # to store the results
  bsave <- matrix(0, ndraw, nvars)
  psave <- matrix(0, ndraw, nb_rho)
  ssave <- rep(0, ndraw)
  acc_rate <- matrix(0, ndraw, nb_rho)
  # Parameters related to the prior distribution
  # V <- matrix(1, n, n)
  sige <- 1
  
  cc1 <- 0.2  # tuning parameter, equation (5.28) in Lesage book
  acc1 <- 0
  
  if (nb_rho == 2 | nb_rho == 3) {
    cc2 <- cc1
    acc2 <- 0
  }
  
  if (nb_rho == 3 & model != "model_8") {
    cc3 <- cc1
    acc3 <- 0
  }
  
  rmin <- -1
  rmax <- 1
  
  # x-prime*w
  x_prime_x <- crossprod(x_matrix)
  zpzti <- solve(x_prime_x)
  zpzty1 <- crossprod(x_matrix, y)
  chol_zpzti <- chol(zpzti) 
  
  # Bayesian algorithm
  pb <- txtProgressBar(min = 0, max = ndraw, initial = 0, 
                       char = "=", width = (getOption("width")), 
                       style = 1)
  counter <- 0
  
  for (iter in 1:ndraw) {
    tau <- c(1, -rho)
    h <- sqrt(sige) * chol_zpzti
    
    bdraw1 <- t(h) %*% as.numeric(mvtnorm::rmvnorm(1, mean = rep(0, nvars))) +
      zpzti %*% zpzty1

    bdraw2 <- t(h) %*% as.numeric(mvtnorm::rmvnorm(1, mean = rep(0, nvars))) +
      zpzti %*% zpzty2

    bdraw <- cbind(bdraw1, bdraw2)
    
    if (nb_rho == 2) {
      bdraw3 <- t(h) %*% as.numeric(mvtnorm::rmvnorm(1, mean = rep(0, nvars))) +
        zpzti %*% zpzty3
      
      bdraw <- cbind(bdraw, bdraw3)
    } else {
      if (nb_rho == 3) {
        bdraw3 <- t(h) %*% as.numeric(mvtnorm::rmvnorm(1, mean = rep(0, nvars))) +
          zpzti %*% zpzty3
        
        bdraw4 <- t(h) %*% as.numeric(mvtnorm::rmvnorm(1, mean = rep(0, nvars))) +
          zpzti %*% zpzty4
        
        bdraw <- cbind(bdraw, bdraw3, bdraw4)
      }
    }
    
    beta_draw <- as(bdraw, "matrix") %*% tau
    # update for beta ends here
    
    # update for sige starts here
    y_hat <- x_matrix %*% bdraw
    E1 <- y - y_hat[, 1]
    E2 <- y2 - y_hat[, 2]
    
    Q <- matrix(0, nb_rho + 1, nb_rho + 1)
    Q[1, 1] <- sum(E1^2)
    Q[1, 2] <- Q[2, 1] <- sum(E1 * E2)  
    Q[2, 2] <- sum(E2^2)  
    
    if (nb_rho == 2 | nb_rho == 3) {
      E3 <- y3 - y_hat[, 3]
      Q[1, 3] <- Q[3, 1] <- sum(E1 * E3)
      Q[2, 3] <- Q[3, 2] <- sum(E2 * E3)
      Q[3, 3] <- sum(E3^2)
    }
    
    if (nb_rho == 3) {
      E4 <- y4 - y_hat[, 4]
      Q[1, 4] <- Q[4, 1] <- sum(E1 * E4)
      Q[2, 4] <- Q[4, 2] <- sum(E2 * E4)
      Q[3, 4] <- Q[4, 3] <- sum(E3 * E4)
      Q[4, 4] <- sum(E4^2)
    }
    
    epe <- t(tau) %*% Q %*% tau # Sum of Square residuals (p. 222, in Lesage book) 
    
    nu <- 0 
    d0 <- 0
    nu1 <- n + 2 * nu
    d1 <- 2 * d0 + epe
    chi <- rgamma(1, nu1 * 0.5) * 2
    sige <- as.numeric(d1/chi) 
    # update for sige ends here
    
    # update for rho1, rho2, rho3 starts here
    # update rho1 using metropolis hastings
    #logdet <- as.numeric(determinant(
    #  switch(nb_rho, 
    #         "1" = Diagonal(n) -  rho[1] * W_2, 
    #         "2" = Diagonal(n) -  rho[1] * W_2 - rho[2] * W_3,
    #         "3" = Diagonal(n) -  rho[1] * W_2 - rho[2] * W_3 - rho[3] * W_4), 
    #  logarithm = TRUE)$modulus)
    
    if (nb_rho == 1) {
      logdet <- - sum(traces * rho[1]^(1:m)/(1:m)) 
    } else {
      if (model == "model_8") {
        logdet <- - sum(traces_1 * rho[1]^(1:m)/(1:m)) - sum(traces_2 * rho[2]^(1:m)/(1:m))
      } else {
        if (model == "model_7") {
          logdet <- as.numeric(determinant(Diagonal(n) -  rho[1] * W_2 - rho[2] * W_3, 
                             logarithm = TRUE)$modulus)
        } else {
          logdet <- as.numeric(determinant(Diagonal(n) -  rho[1] * W_2 - rho[2] * W_3 - rho[3] * W_4, 
                                           logarithm = TRUE)$modulus)
        }
      }
    }
    
    rho_x <- logdet - epe / (2 * sige)
      
    accept <- TRUE
    rho1_c <- rho[1] + cc1 * rnorm(1)
    
    while (accept) {
      if ((rho1_c > rmin) & (rho1_c < rmax))
        accept <- FALSE
      else
        rho1_c <- rho[1] + cc1 * rnorm(1)
    }
    
    if (nb_rho == 1)
      rho_temp <- rho1_c
    
    if (nb_rho == 2)
      rho_temp <- c(rho1_c, rho[2])
    
    if (nb_rho == 3)
      rho_temp <- c(rho1_c, rho[2], rho[3])
    
    if (model == "model_8")
      rho_temp[3] <- - rho_temp[1] * rho_temp[2]
    
    if (nb_rho == 1) {
      logdet_y <- - sum(traces * rho_temp[1]^(1:m)/(1:m)) 
    } else {
      if (model == "model_8") {
        logdet_y <- - sum(traces_1 * rho_temp[1]^(1:m)/(1:m)) - sum(traces_2 * rho_temp[2]^(1:m)/(1:m))
      } else {
        if (model == "model_7") {
          logdet_y <- as.numeric(determinant(Diagonal(n) -  rho_temp[1] * W_2 - rho_temp[2] * W_3, 
                                           logarithm = TRUE)$modulus)
        } else {
          logdet_y <- as.numeric(determinant(Diagonal(n) -  rho_temp[1] * W_2 - rho_temp[2] * W_3 - rho_temp[3] * W_4, 
                                           logarithm = TRUE)$modulus)
        }
      }
    }
    
    tau_y <- c(1, -rho_temp)
    epe_y <- t(tau_y) %*% Q %*% tau_y
    rho_y <- logdet_y - epe_y / (2 * sige)
    
    ru <- runif(1)
    
    if ((rho_y - rho_x) > exp(1)) {
      p <- 1
    } else {
      ratio <- exp(rho_y - rho_x)
      p <- min(1, ratio)
    }
    
    if (ru < p) {
      rho[1] <- rho1_c
      acc1 <- acc1 + 1
    }
    
    acc_rate[iter, 1] <- acc1 / iter
    
    if (acc_rate[iter, 1] < 0.4) {
      cc1 <- cc1 / 1.1
    } else {
      if (acc_rate[iter, 1] > 0.6)
        cc1 <- cc1 * 1.1
    }
    
    if (model == "model_8")
      rho[3] <- - rho[1] * rho[2]
    
    # update rho2 using metropolis hastings
    if (nb_rho == 2 | nb_rho == 3) {
      
      if (model == "model_8") {
        logdet_2 <- - sum(traces_1 * rho[1]^(1:m)/(1:m)) - sum(traces_2 * rho[2]^(1:m)/(1:m))
      } else {
        if (model == "model_7") {
          logdet_2 <- as.numeric(determinant(Diagonal(n) -  rho[1] * W_2 - rho[2] * W_3, 
                                             logarithm = TRUE)$modulus)
        } else {
          logdet_2 <- as.numeric(determinant(Diagonal(n) -  rho[1] * W_2 - rho[2] * W_3 - rho[3] * W_4, 
                                             logarithm = TRUE)$modulus)
        }
      }
      
      tau_x <- c(1, -rho)
      epe_x <- t(tau_x) %*% Q %*% tau_x
      rho_x <- logdet_2 - epe_x/(2 * sige)
      accept <- T
      rho2_c <- rho[2] + cc2*rnorm(1)
      
      while (accept) {
        if ((rho2_c > rmin) & (rho2_c < rmax))
          accept <- F
        else
          rho2_c <- rho[2] + cc2*rnorm(1)
      }
      
      if (nb_rho == 2)
        rho_temp <- c(rho[1], rho2_c)
      
      if (nb_rho == 3)
        rho_temp <- c(rho[1], rho2_c, rho[3])
      
      if (model == "model_8") {
        rho_temp[3] <- - rho_temp[1] * rho_temp[2]
        logdet_y <- - sum(traces_1 * rho_temp[1]^(1:m)/(1:m)) - sum(traces_2 * rho_temp[2]^(1:m)/(1:m))
      } else {
        if (model == "model_7") {
          logdet_y <- as.numeric(determinant(Diagonal(n) -  rho_temp[1] * W_2 - rho_temp[2] * W_3, 
                                             logarithm = TRUE)$modulus)
        } else {
          logdet_y <- as.numeric(determinant(Diagonal(n) -  rho_temp[1] * W_2 - rho_temp[2] * W_3 - rho_temp[3] * W_4, 
                                             logarithm = TRUE)$modulus)
        }
      }
      
      tau_y <- c(1, -rho_temp)
      epe_y <- t(tau_y) %*% Q %*% tau_y
      rho_y <- logdet_y - epe_y/(2*sige)

      ru <- runif(1)
      
      if ((rho_y - rho_x) > exp(1)) {
        p <- 1
      } else {
        ratio <- exp(rho_y - rho_x)
        p <- min(1, ratio)
      }
      
      if (ru < p) {
        rho[2] <- rho2_c
        acc2 <- acc2 + 1
      }
      
      acc_rate[iter, 2] <- acc2/iter
      
      if (acc_rate[iter, 2] < 0.4) {
        cc2 <- cc2 / 1.1
      } else {
        if (acc_rate[iter, 2] > 0.6)
          cc2 <- cc2 * 1.1
      }
    }
    
    if (model == "model_8")
      rho[3] <- - rho[1] * rho[2]
    
    # update rho3 using metropolis hastings
      if (nb_rho == 3 & model != "model_8") {
        logdet_2 <- as.numeric(determinant(Diagonal(n) -  rho[1] * W_2 - rho[2] * W_3 - rho[3] * W_4, 
                                                      logarithm = TRUE)$modulus)
        
        tau_x <- c(1, -rho)
        epe_x <- t(tau_x) %*% Q %*% tau_x
        rho_x <- logdet_2 - epe_x/(2 * sige)
        accept <- T
        rho3_c <- rho[3] + cc3*rnorm(1)
        
        while (accept) {
          if ((rho3_c > rmin) & (rho3_c < rmax))
            accept <- F
          else
            rho3_c <- rho[3] + cc3*rnorm(1)
        }
        
        rho_temp <- c(rho[1], rho[2], rho3_c)
        logdet_y <- as.numeric(determinant(Diagonal(n) -  rho_temp[1] * W_2 - rho_temp[2] * W_3 - rho_temp[3] * W_4, 
                                           logarithm = TRUE)$modulus)
        
        tau_y <- c(1, -rho_temp)
        epe_y <- t(tau_y) %*% Q %*% tau_y
        rho_y <- logdet_y - epe_y/(2*sige)
        ru <- runif(1)
        
        if ((rho_y - rho_x) > exp(1)) {
          p <- 1
        } else {
          ratio <- exp(rho_y - rho_x)
          p <- min(1, ratio)
        }
        
        if (ru < p) {
          rho[3] <- rho3_c
          acc3 <- acc3 + 1
        }
        
        acc_rate[iter, 3] <- acc3/iter
        
        if (acc_rate[iter, 3] < 0.4) {
          cc3 <- cc3 / 1.1
        } else {
          if (acc_rate[iter, 3] > 0.6)
            cc3 <- cc3 * 1.1
        }
        
      }
    
    # save the results
    bsave[iter, ] <- as.numeric(beta_draw)
    ssave[iter] <- sige
    psave[iter, ] <- rho
    
    counter <- counter + 1
    setTxtProgressBar(pb, counter)
  }
  
  bsave <- as.matrix(bsave[nomit:ndraw, ])
  psave <- as.matrix(psave[nomit:ndraw, ])
  
  return(list(bsave = bsave,
                psave = psave))
}

summary_bayesian_sar_flow <- function(x) {
  
  res_sar <- cbind(mean = c(apply(x$psave, 2, mean), apply(x$bsave, 2, mean)),
                   median = c(apply(x$psave, 2, median), apply(x$bsave, 2, median)),
                   lower_05 = c(apply(x$psave, 2, quantile, 0.05), 
                                apply(x$bsave, 2, quantile, 0.05)),
                   lower_95 = c(apply(x$psave, 2, quantile, 0.95), 
                                apply(x$bsave, 2, quantile, 0.95)))
  
  res_sar <- cbind(res_sar, 
                   t_stat = res_sar[, "mean"]/c(apply(x$psave, 2, sd), 
                                                apply(x$bsave, 2, sd)))
  
  print(res_sar)
}

#sar_flow_2(x = flows_data[1:60, c("x1_d", "x2_d", "x1_o", "x2_o")], 
#           y = flows_data[1:60, "y_2"], 
#           g = flows_data[1:60, "g"], 
#           W_d = W_d[1:60, 1:60]) 
