library(Matrix)

emsc <- function(x, p = 2, m = NULL, idx = NULL ){
  x <- as.matrix(x)
  if (is.null(m) & is.null(idx))
    m = colMeans(x)
  else{if (is.null(m)) m = x[idx,]
  }
  ns <- dim(x)[1]
  nv <- dim(x)[2]
    emsc_results = matrix(0,ns,nv)
  for (i in 1:ns){
    fit <- multifit(x[i,], m, p)
    emsc_results[i,] <- fit$res/fit$b[1] + m
  }
  return(emsc_results)
}

multifit <- function(y, x = NULL, p = NULL){
  
  if (is.null(p)) 
    x <- t(x) 
  else {m <- length(x)
    background <- matrix(1,m,1)
    if (p > 0){
      for (i in 1:p){
        background <- cbind(background, matrix(seq(0,1,by=1/(m-1))^i))
        }
    }
    x <- cbind(x, background)
  }
  b <- ginv(t(x)%*%x)%*%t(x)%*%y
  fit <- x%*%b
  res <- y - fit
  results <- list(b = b, fit = fit, res = res)
  return(results)
}

snv <- function(x) {
  no_samples <- nrow(x)
  no_vars <- ncol(x)
  meanx <- matrix(rowMeans(x), no_samples, no_vars)
  stdx <- apply(x, 1, sd)
  stdx <- matrix(stdx, no_samples, no_vars)
  results <- x - meanx
  results <- results / stdx
  return(results)
}


mos <- function(data = NULL) {
    ns <- dim(data)[1]
    ms <- matrix(0, ns, 1)
    mf <- matrix(0, ns, 1)
    for (i in 1:ns) {
        diff_data <- diff(data[i,])
        sign_x <- sign(diff_data)
        zcp <- length(which(abs(diff(sign_x)) == 2))
        mf[i] <- norm(data[i,], type = "2") / (norm(as.matrix(diff_data), type = "2") * zcp)
        data_m <- as.matrix(data[i,] - mean(data[i,]))
        ms[i] <- norm(data_m, type = "2") / norm(diff(data_m), type = "2")
    }
    results <- list(mf = mf, ms = ms)
}

baseline_corr <- function(data, lamda = 1e8, p = .001, d = 2) {
    bl <- asysm(t(data), lamda, p, d)
    results <- data - t(bl)
}

asysm <- function(y, lamda = 1e8, p = .001, d = 2) {
    ns <- dim(y)[2]
    ny <- dim(y)[1]
    if (ns == 1) {
        z <- matrix(1, ny, 1)
        m <- length(y)
        w <- matrix(1, m, 1)
        rep <- 1
        iter <- 1
        while (rep & iter <= 50) {
            z <- difsmw(y, lamda, w, d)
            w0 <- w
            w = p * (y > z) + (1 - p) * (y <= z)
            rep = sum(abs(w - w0)) > 0
            iter <- iter+1
        }
    } else {
        z <- matrix(0, ny, ns)
        for (i in 1:ns) {            
            m <- length(y[, i])
            w <- matrix(1, m, 1)
            rep <- 1
            while (rep & iter <= 50) {
                z[, i] <- difsmw(y[, i], lamda, w, d)
                w0 <- w
                w <- p * (y[, i] > z[, i]) + (1 - p) * (y[, i] <= z[, i])
                rep <- sum(abs(w - w0)) > 0
                iter <- iter + 1
            }
        }
    }
    return(z)   
}

difsmw <- function(y, lamda, w, d) {
    m <- length(y)
    W = diag(as.vector(w))
    D <- diff(diag(m), differences = d)
    C <- chol(W + lamda * (t(D) %*% D))
    z <- solve(C) %*% (solve(t(C)) %*% as.matrix((w * y)))
    return(z)
}

gausswin <- function(window = NULL, alpha = NULL) {
    N <- window - 1
    if (is.null(alpha)) alpha = 2.5
    n <- c(0:N) - N / 2
    results <- exp(-(1 / 2) * (2.5 * n / (N / 2)) ^ 2)
    return(results)
}

gaussian_sm <- function(data = NULL, window_width = 7) {
    ns <- dim(data)[1]
    nv <- dim(data)[2]
    half_width <- round(window_width / 2)
    gaussFilter <- gausswin(window_width)
    gaussFilter <- gaussFilter / sum(gaussFilter)
    results <- matrix(0,ns,nv)
    for (i in 1:ns) {
        data_sm <- conv(data[i,], gaussFilter)
        results[i,] <- data_sm[(half_width - 1):(length(data_sm) - half_width)]
    }    
    return(results)
}

normalise_tot <- function(data = NULL) {
    ns <- dim(data)[1]
    nv <- dim(data)[2]
    data_norm <- matrix(0,ns,nv)
    for (i in 1:ns) {
        data_norm[i,] <- data[i,] / sqrt(sum(data[i,]^2))
    }
    return(data_norm)
}

auto_sc <- function(data = NULL) {
  ns <- dim(data)[1]
  nv <- dim(data)[2]
  data_sd <- apply(data, 2, sd)
  scale_factor <- data_sd
  null_var <- NULL
  if (any(data_sd == 0))
    null_var <- which(data_sd == 0)
  data_sd <- t(matrix(as.matrix(data_sd), nv, ns))
  data_centre <- colMeans(data)
  data_mean <- t(matrix(as.matrix(data_centre), nv, ns))
  data_scaled <- (data - data_mean) / data_sd
  if (!is.null(null_var)) 
    data_scaled[, null_var] <- 0
  results <- list(data = data_scaled, scale_factor = scale_factor, data_centre = data_centre)
  return(results)
}

pareto_sc <- function(data = NULL) {
  ns <- dim(data)[1]
  nv <- dim(data)[2]
  data_sd <- apply(data, 2, sd)
  scale_factor <- data_sd
  null_var <- NULL
  if (any(data_sd == 0))
    null_var <- which(data_sd == 0)
  data_sd <- t(matrix(as.matrix(data_sd), nv, ns))
  data_centre <- colMeans(data)
  data_mean <- t(matrix(as.matrix(data_centre), nv, ns))
  data_scaled <- (data - data_mean) / sqrt(data_sd)
  if (!is.null(null_var)) 
    data_scaled[, null_var] <- 0
  results <- list(data = data_scaled, scale_factor = scale_factor, data_centre = data_centre)
  return(results)
}

arpls <- function(y = NULL, lamda = 1e6, ratio = 1e-2){
  if (is.matrix(y) | is.data.frame(y)){
    ns <- dim(y)[1]  
    nv <- dim(y)[2]
    if (ns > 1 & nv > 1) 
      stop("Signal y need to be a vector, not matrix!")  
    else y <- as.vector(y, mode = "numeric")
  }
  N <- length(y)
  D <- diff(sparseMatrix(i = 1:N, j = 1:N, x = 1), difference = 2)
  H <- lamda * (t(D) %*% D)
  w <- as.vector(rep(1, N))
  max_iter <- 100
  k <- 1
  while (k <= max_iter){
    W <- Diagonal(x = w)
    C <- chol(W + H)
    z <- solve(C) %*% (solve(t(C)) %*% as.matrix(w * y))
    d <- y - z
    dn <- d[d < 0]
    m <- mean(dn)
    s <- sd(dn)
    wt <- 1 / (1 + exp(2 * (d - (2 * s - m)) / s))
    if (norm(as.matrix(w - wt), "2") / norm(as.matrix(w), "2") < ratio)
      break
    w <- as.vector(wt)
    k <- k + 1
  }
#  if (k > max_iter)
#     print("Maxium iterations reached, the algorithm probably failed to converge!")
  return(as.vector(z))
}