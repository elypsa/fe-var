library(Matrix)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("src/test.cpp")
source('R/aux.R')

fe_var <- function(Y,
                   X = NULL,
                   fe,
                   lags,
                   irf_hor = 11,
                   conf_levels = c(90, 95),
                   nboot = 2e3,
                   perm = NULL,
                   normalize = T,
                   seed = 2025,
                   tol = 1e-3,
                   maxit = 15,
                   delta_update = .5) {
  # Y: matrix of dependent variables, rows = time, columns = variables
  # X: matrix of lagged dependents, default is NULL and X is created, no const!
  # fe: vector of fixed effects, e.g. country
  # lags: number of lags in VAR
  # irf_hor: number of periods to forecast for IRF
  # nboot: number of bootstrap samples
  # perm: permutation matrix for ordering in cholesky
  # normalize: whether to normalize the IRF
  # seed: random seed for reproducibility
  # tol: tolerance for convergence
  # maxit: maximum number of iterations for convergence

  if (!is.vector(fe) || length(fe) != nrow(Y)) {
    stop("fe must be a vector with the same length as the number of rows in Y.")
  }

  if (!is.numeric(lags) || lags < 1) {
    stop("lags must be a positive integer.")
  }

  if (!is.numeric(irf_hor) || irf_hor < 1) {
    stop("irf_hor must be a positive integer.")
  }

  if (!is.numeric(nboot) || nboot < 1) {
    stop("nboot must be a positive integer.")
  }

  if (is.null(perm)) {
    perm <- diag(ncol(Y))  # Default to identity matrix if no permutation is provided
  } else if (!is.matrix(perm) ||
             nrow(perm) != ncol(Y) || ncol(perm) != ncol(Y)) {
    stop("perm must be a square matrix with dimensions equal to the number of columns in Y.")
  }

  # Create lags if X is not available
  if (is.null(X)) {
    X <- lag_panel(Y, fe, lags)
  }

  if(all(conf_levels < 1)) {
    conf_levels <- conf_levels * 100  # Convert to percentage if given as fraction
    message("conf_levels should be in percentage, e.g. 90 or 95. Converting to percentage.")
  }

  stopifnot(ncol(X) == ncol(Y) * lags)


  # Remove rows with NAs (due to lags)
  valid_rows <- complete.cases(X)
  X <- X[valid_rows, ]
  Y <- Y[valid_rows, ]
  fe <- fe[valid_rows]
  fe_t <- table(fe)

  NT <- nrow(Y)
  Ny <- ncol(Y)
  Nfe <- length(unique(fe))
  Nx <- ncol(X) + Nfe
  Tmax <- max(table(fe))
  ones_vecs <- lapply(table(fe), \(n) rep(1, n))
  D <- bdiag(ones_vecs)
  D <- as.matrix(D)
  XD <- cbind(X, D[, 1:(Nfe - 1)], 1)  # Add fixed effects to X and constant
  X <- cbind(X, 1)

  stopifnot(ncol(XD) == Ny * lags + Nfe)
  stopifnot(nrow(D) == NT)

  bet <- solve(t(XD) %*% XD) %*% t(XD) %*% Y
  res <- Y - XD %*% bet

  S <- t(res) %*% res
  Sigma <- S / (NT - Nx)
  Sigma0 <- Sigma
  bet0 <- bet
  bet_bs_smp <- array(dim = c(Ny * lags, Ny, nboot))
  sigma_bs_smp <- array(dim = c(Ny, Ny, nboot))
  res <- sweep(res, 2, colMeans(res, na.rm = TRUE))  # Detrend residuals

  conv_test <- 10

  bet_bs <- bet0[1:(Ny * lags), ]
  iter <- 0
  while (conv_test > tol && iter < maxit) {
    set.seed(seed)
    print(iter)
    for (i in 1:nboot) {
      rr <- 1 - 2 * (runif(Tmax) > 0.5)
      rr <- unlist(lapply(table(fe), function(x)
        rr[1:x]), use.names = F)
      resb <- rr * res
      boot_var <- sim_var(X, bet_bs, rr, resb, fe_t, Ny, lags)
      Yb <- boot_var$Yb
      Xb <- boot_var$Xb
      X_temp <- cbind(Xb, D[, 1:(Nfe - 1)], 1)
      bet_tmp <- solve(t(X_temp) %*% X_temp) %*% t(X_temp) %*% Yb
      bet_bs_smp[, , i] <- bet_tmp[1:(Ny * lags), ]
      res_temp <- Yb - X_temp %*% bet_tmp
      sigma_bs_smp[, , i] <- t(res_temp) %*% (res_temp) / (NT - Nx)
    }
    betbs_mean <- apply(bet_bs_smp, c(1, 2), mean)
    diff <- bet0[1:(Ny * lags), ] - betbs_mean
    diffvec <- as.vector(diff)
    conv_test <- sqrt(mean(diffvec^2))

    iter <- iter + 1
    if (conv_test > tol) {
      print(conv_test)
      bet_bs <- bet_bs + delta_update * diff
    }
    if (conv_test < tol || iter == maxit) {
      virf_boot <- array(0, dim = c(irf_hor, Ny, Ny, nboot))
      for (j in 1:nboot) {
        virf_boot[, , , j] <- irf(bet_bs_smp[, , j], sigma_bs_smp[, , j] , perm, irf_hor, normalize, lags)
      }
      # calculate CI using conf_levels
      low <- (100 - conf_levels) / 200
      high <- 1 - low
      virf_conf_intervals <- apply(virf_boot, c(1, 2, 3), function(x) {
        quantile(x, probs = c(low, high), na.rm = TRUE)})
      virf_mean <- apply(virf_boot, c(1, 2, 3), mean)
       # first dimension is according to confidence level(s)
      virf <- list(probs = c(low, high), intervals = virf_conf_intervals,
                   virf_mean = virf_mean, virf_boot = virf_boot)
    }
  }
  bc_diff <- bet[1:(Ny * lags), ] - bet_bs
  bet_bc <- rbind(bet_bs, bet[(Ny * lags + 1):NROW(bet), ]) # bias corrected

  # adjust fixed effects
  uhat <- Y - X[, 1:(Ny * lags)] %*% bet_bc[1:(Ny * lags), ]
  D1 <- cbind(D[, 1:(Nfe - 1)], 1)
  fe_bc <- solve(t(D1) %*% D1) %*% t(D1) %*% uhat
  bet_bc[(Ny * lags + 1):NROW(bet_bc), ] <- fe_bc
  res_bc <- Y - XD %*% bet_bc
  sigma_bc <- t(res_bc) %*% res_bc / (NT - Nx)
  res <- res_bc

  # virf$mean <- irf(bet_bc, sigma_bc, perm, irf_hor, normalize, lags)

  return(list(
    bet0 = bet0,
    sigma0 = Sigma0,
    # without bias correction
    bet = bet_bc,
    # bias corrected
    virf = virf,
    sigma_bc = sigma_bc,
    data = list(Y = Y, X = X, fe = fe)
  ))
}

