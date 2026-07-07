safe_pmax <- function(X, eps) {
  objectClass <- class(X)[[1]]

  if (objectClass == "gpu.matrix.torch" || objectClass == "gpu.matrix.tensorflow") {
    return(X + (X < eps) * (eps - X))
  }

  pmax(X, eps)
}

snmf <- function(
  V,
  S = diag(ncol(V)),
  k = 10,
  Winit = NULL,
  Hinit = NULL,
  tol = 1e-03,
  Rupdate_iter = 20,
  niter = 100,
  num_initializations = 10,
  dispersion_mode = c("full", "gene", "spot")
) {

  dispersion_mode <- match.arg(dispersion_mode)

  dtype <- "float32"

  objectClass <- class(V)[[1]]
  objectPackage <- attr(class(V), "package")

  if (!is.null(objectPackage) &&
      (objectClass == "gpu.matrix.torch" || objectClass == "gpu.matrix.tensorflow")) {

    if (is.null(Winit)) {
      Winit <- gpu.matrix(
        runif(nrow(V) * k), nrow(V), k,
        dtype = dtype(V),
        type = GPUmatrix:::typeGPUmatrix(V),
        device = GPUmatrix:::device(V)
      )
    }

    if (is.null(Hinit)) {
      Hinit <- gpu.matrix(
        runif(k * ncol(V)), k, ncol(V),
        dtype = dtype(V),
        type = GPUmatrix:::typeGPUmatrix(V),
        device = GPUmatrix:::device(V)
      )
    }

  } else {

    if (is.null(Winit)) {
      Winit <- matrix(runif(nrow(V) * k), nrow(V), k)
    }

    if (is.null(Hinit)) {
      Hinit <- matrix(runif(k * ncol(V)), k, ncol(V))
    }
  }

  controlDimensionNMF(Winit, Hinit, V, k)

  Vold <- V

  if (!is.null(objectPackage) &&
      (objectClass == "gpu.matrix.torch" || objectClass == "gpu.matrix.tensorflow")) {

    R <- list(
      alpha = rep(0, nrow(V)),
      beta = rep(0, ncol(V)),
      phi = gpu.matrix(
        matrix(1e-4, nrow(V), ncol(V)),
        dtype = dtype(V),
        type = GPUmatrix:::typeGPUmatrix(V),
        device = GPUmatrix:::device(V)
      ),
      eps = 1e-10,
      phi.min = 1e-8,
      phi.max = 1e3,
      inner.iter = 10
    )

  } else {

    R <- list(
      alpha = rep(0, nrow(V)),
      beta = rep(0, ncol(V)),
      phi = matrix(1e-4, nrow(V), ncol(V)),
      eps = 1e-10,
      phi.min = 1e-8,
      phi.max = 1e3,
      inner.iter = 10
    )
  }

  initial_iterations <- max(1, floor(niter / 10))

  best_loss <- Inf
  best_W <- NULL
  best_H <- NULL

  for (init_run in seq_len(num_initializations)) {

    if (!is.null(objectPackage) &&
        (objectClass == "gpu.matrix.torch" || objectClass == "gpu.matrix.tensorflow")) {

      W_current <- gpu.matrix(
        runif(nrow(V) * k), nrow(V), k,
        dtype = dtype(V),
        type = GPUmatrix:::typeGPUmatrix(V),
        device = GPUmatrix:::device(V)
      )

      H_current <- gpu.matrix(
        runif(k * ncol(V)), k, ncol(V),
        dtype = dtype(V),
        type = GPUmatrix:::typeGPUmatrix(V),
        device = GPUmatrix:::device(V)
      )

    } else {

      W_current <- matrix(runif(nrow(V) * k), nrow(V), k)
      H_current <- matrix(runif(k * ncol(V)), k, ncol(V))
    }

    for (iter_init in seq_len(initial_iterations)) {

      WHS_current <- W_current %*% H_current %*% S
      WHS_current <- safe_pmax(WHS_current, R$eps)

      W_current <- updateW(V, W_current, H_current, S, WHS_current, R)

      WHS_current <- W_current %*% H_current %*% S
      WHS_current <- safe_pmax(WHS_current, R$eps)

      H_current <- updateH(V, W_current, H_current, S, WHS_current, R)
    }

    V_reconstructed <- W_current %*% H_current %*% S
    current_loss <- mean((V_reconstructed - V)^2)

    if (current_loss < best_loss) {
      best_loss <- current_loss
      best_W <- W_current
      best_H <- H_current
    }
  }

  Winit <- best_W
  Hinit <- best_H

  for (iter in seq_len(niter)) {

    WHS <- Winit %*% Hinit %*% S
    WHS <- safe_pmax(WHS, R$eps)

    Winit <- updateW(V, Winit, Hinit, S, WHS, R)

    WHS <- Winit %*% Hinit %*% S
    WHS <- safe_pmax(WHS, R$eps)

    Hinit <- updateH(V, Winit, Hinit, S, WHS, R)

    if (iter %% Rupdate_iter == 0) {
      WHS <- Winit %*% Hinit %*% S
      WHS <- safe_pmax(WHS, R$eps)

      R <- updateR(V, Winit, Hinit, S, WHS, R, mode = dispersion_mode)
    }

    if (iter %% 100 == 0) {
      cat("Iteration:", iter, "\n")
    }

    if (iter %% 10 == 0) {

      myD <- colSums(Winit)
      myD <- safe_pmax(myD, R$eps)

      Winit <- t(t(Winit) / myD)
      Hinit <- Hinit * myD

      Vnew <- Winit %*% Hinit %*% S

      mse_change <- mean((Vnew - Vold)^2)

      if (is.na(mse_change)) {
        stop("Error in calculating mean squared error. Check dimensions of Vnew and Vold.")
      }

      if (mse_change < tol) {
        message("NMF converged early.")
        return(list(W = Winit, H = Hinit, phi = R$phi, alpha = R$alpha, beta = R$beta, niter=iter + num_initializations*initial_iterations))
      }

      Vold <- Vnew
    }
  }

  warning("Maximum number of iterations reached without convergence. Consider increasing 'niter'.")
  
  return(list(W = Winit, H = Hinit, phi = R$phi, alpha = R$alpha, beta = R$beta, niter = iter + num_initializations*initial_iterations))
}

# Optimized version of updateW function
updateW <- function(V, W, H, S, WHS, R) {

  eps <- if (!is.null(R$eps)) R$eps else 1e-10

  WHS <- safe_pmax(WHS, eps)

  Phi <- R$phi
  Phi <- safe_pmax(Phi, if (!is.null(R$phi.min)) R$phi.min else 1e-8)

  Rpos <- V / WHS

  Q <- (1 + Phi * V) / (1 + Phi * WHS)

  B <- H %*% S

  numerator <- Rpos %*% t(B)
  denominator <- Q %*% t(B)

  W.new <- W * numerator / safe_pmax(denominator, eps)
  W.new <- safe_pmax(W.new, eps)

  return(W.new)
}

# Optimized version of updateH function
updateH <- function(V, W, H, S, WHS, R) {

  eps <- if (!is.null(R$eps)) R$eps else 1e-10

  WHS <- safe_pmax(WHS, eps)

  Phi <- R$phi
  Phi <- safe_pmax(Phi, if (!is.null(R$phi.min)) R$phi.min else 1e-8)

  Rpos <- V / WHS

  Q <- (1 + Phi * V) / (1 + Phi * WHS)

  numerator   <- t(W) %*% Rpos %*% t(S)
  denominator <- t(W) %*% Q %*% t(S)

  H.new <- H * numerator / safe_pmax(denominator, eps)
  H.new <- safe_pmax(H.new, eps)

  return(H.new)
}

updateR <- function(
  V, W, H, S, WHS, R,
  mode = c("full", "gene", "spot")
) {

  mode <- match.arg(mode)

  eps <- if (!is.null(R$eps)) R$eps else 1e-10
  phi.min <- if (!is.null(R$phi.min)) R$phi.min else 1e-8
  phi.max <- if (!is.null(R$phi.max)) R$phi.max else 1e3
  inner.iter <- if (!is.null(R$inner.iter)) R$inner.iter else 10

  m <- nrow(V)  # genes
  n <- ncol(V)  # spots

  WHS <- safe_pmax(WHS, eps)

  # Moment target:
  # Var(V_ij) = mu_ij + phi_ij * mu_ij^2
  D <- (V - WHS)^2 - WHS
  C <- WHS^2

  if (mode == "gene") {

    alpha <- rowSums(D) / safe_pmax(rowSums(C), eps)
    alpha <- pmin(safe_pmax(alpha, phi.min), phi.max)

    Phi <- matrix(alpha, nrow = m, ncol = n, byrow = FALSE)

    R$alpha <- alpha
    R$beta  <- rep(0, n)
    R$phi   <- Phi
    R$theta <- 1 / Phi

    return(R)
  }

  if (mode == "spot") {

    beta <- colSums(D) / safe_pmax(colSums(C), eps)
    beta <- pmin(safe_pmax(beta, phi.min), phi.max)

    Phi <- matrix(beta, nrow = m, ncol = n, byrow = TRUE)

    R$alpha <- rep(0, m)
    R$beta  <- beta
    R$phi   <- Phi
    R$theta <- 1 / Phi

    return(R)
  }

  if (mode == "full") {

    if (is.null(R$alpha)) {
      R$alpha <- rep(0, m)
    }

    if (is.null(R$beta)) {
      R$beta <- rep(0, n)
    }

    alpha <- R$alpha
    beta  <- R$beta

    for (iter in seq_len(inner.iter)) {

      alpha <- (
        rowSums(D) -
        rowSums(C * matrix(beta, nrow = m, ncol = n, byrow = TRUE))
      ) / safe_pmax(rowSums(C), eps)

      beta <- (
        colSums(D) -
        colSums(C * matrix(alpha, nrow = m, ncol = n, byrow = FALSE))
      ) / safe_pmax(colSums(C), eps)

      # Identifiability constraint: mean(beta) = 0
      beta.mean <- mean(beta)

      beta  <- beta - beta.mean
      alpha <- alpha + beta.mean
    }

    Phi <- outer(alpha, beta, "+")

    Phi <- pmin(safe_pmax(Phi, phi.min), phi.max)

    R$alpha <- alpha
    R$beta  <- beta
    R$phi   <- Phi
    R$theta <- 1 / Phi

    return(R)
  }
}

# Function to handle potential negative values very close to zero, setting them to zero
# Note: This function is not currently called in NMFKLMixing, but could be used for numerical stability.
setNegativeZero <- function (x) {
  x[x < 0 & x > -1e-10] <- 0 # Set values between -1e-10 and 0 to 0
  return(x)
}

# Function to control and validate the dimensions of input matrices for NMF
controlDimensionNMF <- function (W, H, V, k) {
  if (ncol(W) != k) {
    stop("W must have 'k' columns. 'k' represents the number of latent components.")
  }
  if (nrow(H) != k) {
    stop("H must have 'k' rows. 'k' represents the number of latent components.")
  }
  if (ncol(V) != ncol(H)) {
    stop("V (original matrix) and H (coefficient matrix) must have the same number of columns.")
  }
  if (nrow(V) != nrow(W)) {
    stop("V (original matrix) and W (basis matrix) must have the same number of rows.")
  }
}

