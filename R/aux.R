get_companion_form <- function(bet, Ny, lags) {
  # bet: vector of coefficients, first Ny*lags are the VAR coefficients
  # Ny: number of variables
  # lags: number of lags in VAR

  rbind(t(bet[1:(Ny*lags),]), cbind(diag(rep(1, Ny*(lags-1))), matrix(0, nrow = Ny*(lags-1), ncol = Ny)))
}

irf <- function(bet, sigma, perm, irf_hor, normalize, lags) {
  Ny <- ncol(sigma)
  Sig_chol <- t(perm) %*% t(chol(perm%*% sigma %*% t(perm)) ) %*% perm
  B_comp <- get_companion_form(bet, Ny, lags)

  cc_dy <- diag(rep(1, Ny*lags))
  msel <- rbind(diag(1, Ny),
                matrix(0, (lags-1)*Ny,Ny))
  if (normalize) {
    diagonal <- diag(diag(Sig_chol))
    Sig_chol <- Sig_chol %*% solve(diagonal)
  }
  virf <- array(0, dim = c(irf_hor, Ny, Ny)) # ordering: horizon, variable, shock

  virf[1,,] <- t(msel)%*%cc_dy%*%msel%*%Sig_chol

  for (i in 2:irf_hor) {
    cc_dy <- B_comp%*%cc_dy
    virf[i,,] <- t(msel)%*%cc_dy%*%msel%*%Sig_chol
  }
  return(virf)
}

lag_panel <- function(Y, fe, n) {
  # Create lagged variables for panel data

  ids <- tapply(1:NROW(Y), fe, function(x)
    x)

  X_ind <- lapply(ids, function(id) {
    if (length(id) < n + 1) {
      stop("Not enough observations for the specified number of lags.")
    }
    Y_temp <- Y[id, ]
    lagged_vars <- embed(Y_temp, n + 1)
    # pad NAs
    lagged_vars <- rbind(matrix(NA, nrow = n, ncol = ncol(lagged_vars)), lagged_vars)
  })
  X <- do.call(rbind, X_ind)
  return(X)
}

plot_irf <- function(virf, variable=1) {
  # irf: list containing the IRF results
  if (!is.list(virf) || !all(c("mean", "intervals", "probs") %in% names(virf))) {
    stop("irf must be a list containing 'mean', 'intervals', and 'probs'.")
  }

  mean_irf <- virf$mean[,,variable]
  intervals <- virf$intervals[,,,variable]
  ci_variables <- list()


  for (i in 1:ncol(mean_irf)) {
    ci_variables[[i]] <- t(intervals[,,i])
  }

  # Plotting the IRF
  par(mfrow = c(1, ncol(mean_irf)), mar = c(4, 4, 2, 1))
  for (i in 1:ncol(mean_irf)) {
    matplot(ci_variables[[i]], type = 'l', lty = 2:(length(virf$probs)/2), col = 'gray')
    lines(mean_irf[, i], type = "l", col = "blue", ylim = range(ci_variables[[i]]),
          xlab = "Horizon", ylab = paste("IRF of variable", i),
          main = paste("Impulse Response Function for Variable", i))
    # lines(intervals[, i, 1], col = "red", lty = 2)
    # lines(intervals[, i, 2], col = "red", lty = 2)
    # legend("topright", legend = paste0(probs * 100, "% CI"), col = "red", lty = 2)
  }
}
