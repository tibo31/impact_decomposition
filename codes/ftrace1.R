ftrace1 <- function(w, method = "exact", miter = 10, riter = 50) {
  
  # initialization
  stopifnot(method %in% c("exact", "approx"))
  n <- nrow(w)
  
  if (method == "exact") {
    traces <- numeric(miter)
    traces[1] <- sum(w[0:(n - 1) * n + 1:n])
      w_j <- w
      for (j in 2:miter) {
        w_j <- w %*% w_j 
        traces[j] <- sum(w_j[0:(n - 1) * n + 1:n])
      }
  } else {
    tmat <- matrix(0, miter, riter)
    
    for (j in 1:riter) {
      u <- 2*round(runif(n)) - 1
      wu <- u
      for (i in 1:miter) {
        wu <- w %*% wu
        tmat[i, j] <- sum(wu * u)
      }
    }
    
    traces <- apply(tmat, 1, mean)
    traces[1] <- sum(w[0:(n - 1) * n + 1:n])
    traces[2] <- sum(w * t(w))
  }
  return(traces)
}