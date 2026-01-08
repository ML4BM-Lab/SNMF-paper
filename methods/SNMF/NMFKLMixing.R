NMFKLMixing <- function (V, S = diag(ncol(V)), k = 10, Winit = NULL, Hinit = NULL, tol = 1e-03, 
                         niter = 100, num_initializations=10) 
{
  dtype = "float32" # Define the data type for GPU matrices
  
  # Determine if V is a GPU matrix and set initial W and H accordingly
  objectClass <- class(V)[[1]]
  objectPackage <- attr(class(V), "package")
  
  if (!is.null(objectPackage) && (objectClass == "gpu.matrix.torch" || objectClass == "gpu.matrix.tensorflow")) {
    # If V is a GPU matrix, initialize Winit and Hinit as GPU matrices
    if (is.null(Winit)) {
      Winit <- gpu.matrix(runif(nrow(V) * k), nrow(V), k, 
                          dtype = dtype(V), type = GPUmatrix:::typeGPUmatrix(V), 
                          device = GPUmatrix:::device(V))
    }
    if (is.null(Hinit)) {
      Hinit <- gpu.matrix(runif(k * ncol(V)), k, ncol(V), 
                          dtype = dtype(V), type = GPUmatrix:::typeGPUmatrix(V), 
                          device = GPUmatrix:::device(V))
    }
  } else {
    # If V is a regular matrix, initialize Winit and Hinit as regular matrices
    if (is.null(Winit)) {
      Winit <- matrix(runif(nrow(V) * k), nrow(V), k)
    }
    if (is.null(Hinit)) {
      Hinit <- matrix(runif(k * ncol(V)), k, ncol(V))
    }
  }
  
  # Control dimensions of input matrices to ensure compatibility
  controlDimensionNMF(Winit, Hinit, V, k)
  
  Vold <- V # Store the initial V for convergence check
  
  # Initialize S1 and S2 for update rules, as GPU matrices if applicable
  S1 <- gpu.matrix(dtype = dtype, 1, nrow = nrow(Winit), ncol = nrow(S)) %*% S
  S2 <- gpu.matrix(dtype = dtype, 1, nrow = nrow(V), ncol = ncol(V))
  
  # --- Multiple Initializations and Best Selection ---
  initial_iterations <- max(1, floor(niter / 10)) # Number of iterations for each initialization
  
  best_loss <- Inf # Initialize with a very large loss
  best_W <- NULL
  best_H <- NULL
  
  for (init_run in 1:num_initializations) {
    # Re-initialize W and H for each run
    if (!is.null(objectPackage) && (objectClass == "gpu.matrix.torch" || objectClass == "gpu.matrix.tensorflow")) {
      W_current <- gpu.matrix(runif(nrow(V) * k), nrow(V), k, 
                              dtype = dtype(V), type = GPUmatrix:::typeGPUmatrix(V), 
                              device = GPUmatrix:::device(V))
      H_current <- gpu.matrix(runif(k * ncol(V)), k, ncol(V), 
                              dtype = dtype(V), type = GPUmatrix:::typeGPUmatrix(V), 
                              device = GPUmatrix:::device(V))
    } else {
      W_current <- matrix(runif(nrow(V) * k), nrow(V), k)
      H_current <- matrix(runif(k * ncol(V)), k, ncol(V))
    }
    
    # Run a few iterations for the current initialization
    for (iter_init in 1:initial_iterations) {
      R_current <- W_current %*% (H_current %*% S) # Calculate R for update rules
      H_current <- updateH(V, W_current, H_current, S, R_current, S1) # Update H
      W_current <- updateW(V, W_current, H_current, S, R_current, S2) # Update W
    }
    
    # Calculate the loss (e.g., Frobenius norm of the difference) for the current initialization
    V_reconstructed <- W_current %*% (H_current %*% S)
    current_loss <- mean((V_reconstructed - V)^2) # Mean squared error as loss
    
    # Keep track of the best initialization
    if (current_loss < best_loss) {
      best_loss <- current_loss
      best_W <- W_current
      best_H <- H_current
    }
  }
  
  Winit <- best_W # Use the best W from initializations
  Hinit <- best_H # Use the best H from initializations
  
  # --- Main Iteration Loop ---
  for (iter in 1:niter) {
    # Calculate R (reconstructed matrix)
    R <- Winit %*% (Hinit %*% S) 
    
    # Update H and W using the multiplicative update rules
    Hinit <- updateH(V, Winit, Hinit, S, R, S1) # S1 and S2 are precomputed
    Winit <- updateW(V, Winit, Hinit, S, R, S2) 
    
    # Periodically check for convergence and normalize W
    if (iter %% 100 == 0) {
      cat("Iteration:", iter,"\n")
    }
    if (iter %% 10 == 0) {
        
      # Normalize W by column sums and adjust H accordingly
      myD <- colSums(Winit) # Sum of each column in W
      Winit <- t((t(Winit)/myD)) # Normalize W so each column sums to 1
      Hinit <- Hinit * myD # Adjust H to compensate for W normalization
      
      Vnew <- Winit %*% (Hinit %*% S) # Reconstruct V with updated W and H
      
      # Check for convergence based on the change in V
      if(is.na(mean((Vnew - Vold)^2))) {
        browser()
        stop("Error in calculating mean squared error. Check the dimensions of Vnew and Vold.")
      }
      if (mean((Vnew - Vold)^2) < tol) {
        message("NMF converged early.")
        res <- list(W = Winit, H = Hinit) # Return results
        return(res)
      }
      Vold <- Vnew # Update Vold for the next convergence check
    }
  }
  
  # If the loop finishes without converging within `niter` iterations
  if (iter == niter) {
    warning(message = "Maximum number of iterations reached without convergence. Consider increasing 'niter'.")
    res <- list(W = Winit, H = Hinit) # Return the current W and H
    return(res)
  }
}

# Optimized version of updateW function
updateW <- function (V, W, H, S, R, S2) {
  # R is the reconstructed matrix: W %*% H %*% S
  # S2 is a matrix of ones with the same dimensions as V
  
  tHS <- t(H %*% S) # Pre-calculate t(H %*% S) for efficiency
  num <- (V / R) %*% tHS # Numerator of the update rule
  denom <- S2 %*% tHS # Denominator of the update rule
  W <- W * num / denom # Multiplicative update rule for W
  return(W)
}

# Optimized version of updateH function
updateH <- function (V, W, H, S, R, S1) {
  # R is the reconstructed matrix: W %*% H %*% S
  # S1 is (a matrix of ones with dimensions nrow(W) x nrow(S)) %*% S
  
  num <- (t(W) %*% (V / R)) %*% t(S) # Numerator of the update rule
  denom <- t(W) %*% S1 # Denominator of the update rule
  H <- H * num / denom # Multiplicative update rule for H
  return(H)
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



NMFKLMixingFilter <- function (V, S = diag(ncol(V)), k = 10, Winit = NULL, Hinit = NULL, tol = 1e-03, theta = 0,
                         niter = 100, dtype = "float32") {

  # Determine if V is a GPU matrix and set initial W and H accordingly
  objectClass <- class(V)[[1]]
  objectPackage <- attr(class(V), "package")
  
  if (!is.null(objectPackage) && (objectClass == "gpu.matrix.torch" || objectClass == "gpu.matrix.tensorflow")) {
    # If V is a GPU matrix, initialize Winit and Hinit as GPU matrices
    if (is.null(Winit)) {
      Winit <- gpu.matrix(runif(nrow(V) * k), nrow(V), k, 
                          dtype = dtype(V), type = GPUmatrix:::typeGPUmatrix(V), 
                          device = GPUmatrix:::device(V))
    }
    if (is.null(Hinit)) {
      Hinit <- gpu.matrix(runif(k * ncol(V)), k, ncol(V), 
                          dtype = dtype(V), type = GPUmatrix:::typeGPUmatrix(V), 
                          device = GPUmatrix:::device(V))
    }
  } else {
    # If V is a regular matrix, initialize Winit and Hinit as regular matrices
    if (is.null(Winit)) {
      Winit <- matrix(runif(nrow(V) * k), nrow(V), k)
    }
    if (is.null(Hinit)) {
      Hinit <- matrix(runif(k * ncol(V)), k, ncol(V))
    }
  }
  
  # Control dimensions of input matrices to ensure compatibility
  controlDimensionNMF(Winit, Hinit, V, k)
  
  Vold <- V # Store the initial V for convergence check
  
  # Initialize S1 and S2 for update rules, as GPU matrices if applicable
  S1 <- gpu.matrix(dtype = dtype, 1, nrow = nrow(Winit), ncol = nrow(S)) %*% S
  S2 <- gpu.matrix(dtype = dtype, 1, nrow = nrow(V), ncol = ncol(V))
  
  # --- Multiple Initializations and Best Selection ---
  num_initializations <- 10 # Number of random initializations
  initial_iterations <- max(1, floor(niter / 10)) # Number of iterations for each initialization
  
  best_loss <- Inf # Initialize with a very large loss
  best_W <- NULL
  best_H <- NULL
  
  for (init_run in 1:num_initializations) {
    # Re-initialize W and H for each run
    if (!is.null(objectPackage) && (objectClass == "gpu.matrix.torch" || objectClass == "gpu.matrix.tensorflow")) {
      W_current <- gpu.matrix(runif(nrow(V) * k), nrow(V), k, 
                              dtype = dtype(V), type = GPUmatrix:::typeGPUmatrix(V), 
                              device = GPUmatrix:::device(V))
      H_current <- gpu.matrix(runif(k * ncol(V)), k, ncol(V), 
                              dtype = dtype(V), type = GPUmatrix:::typeGPUmatrix(V), 
                              device = GPUmatrix:::device(V))
    } else {
      W_current <- matrix(runif(nrow(V) * k), nrow(V), k)
      H_current <- matrix(runif(k * ncol(V)), k, ncol(V))
    }
    
    # Run a few iterations for the current initialization
    for (iter_init in 1:initial_iterations) {
      R_current <- W_current %*% (H_current %*% S) # Calculate R for update rules
      H_current <- updateH(V, W_current, H_current, S, R_current, S1) # Update H
      dummy <- mean(colSums(H_current))
      H_current <- H_current/dummy
      W_current <- W_current*dummy
      W_current <- updateW(V, W_current, H_current, S, R_current, S2) # Update W
    }
    
    # Calculate the loss (e.g., Frobenius norm of the difference) for the current initialization
    V_reconstructed <- W_current %*% (H_current %*% S)
    current_loss <- mean((V_reconstructed - V)^2) # Mean squared error as loss
    
    # Keep track of the best initialization
    if (current_loss < best_loss) {
      best_loss <- current_loss
      best_W <- W_current
      best_H <- H_current
    }
  }
  
  Winit <- best_W # Use the best W from initializations
  Hinit <- best_H # Use the best H from initializations
  
  # --- Main Iteration Loop ---
  # Here we will use the theta to foster the sparsity
  # For the time being the filtering matrix will not be optimized
  
  
  F <- gpu.matrix(diag(rep(1-theta-theta / (ncol(Winit)-1), ncol(Winit)))) + theta / (ncol(Winit)-1)
  F <- gpu.matrix(F, dtype = dtype)
  
  for (iter in 1:niter) {
    # Calculate R (reconstructed matrix)
    R <- (Winit %*% F) %*% (Hinit %*% S) 
    
    # Update H and W using the multiplicative update rules
    Hinit <- updateH(V, (Winit %*% F), Hinit, S, R, S1) # S1 and S2 are precomputed and are constant
    dummy <- mean(colSums(Hinit))
    Hinit <- Hinit/dummy
    Winit <- Winit*dummy
    Winit <- updateW(V, Winit, F %*% Hinit, S, R, S2) 
    
    # Periodically check for convergence and normalize W
    if (iter %% 100 == 0) {
      cat("Iteration:", iter,"\n")
    }
    if (iter %% 10 == 0) {

      # Normalize W by column sums and adjust H accordingly
      myD <- mean(colSums(Winit)) # Sum of each column in W
      Winit <- Winit/myD # Normalize W so each column sums to 1
      Hinit <- Hinit * myD # Adjust H to compensate for W normalization
      
      Vnew <- (Winit %*% F) %*% (Hinit %*% S) # Reconstruct V with updated W and H
      
      # Check for convergence based on the change in V
      if(is.na(mean((Vnew - Vold)^2))) {
        browser()
        stop("Error in calculating mean squared error. Check the dimensions of Vnew and Vold.")
      }
      if (mean((Vnew - Vold)^2) < tol) {
        message("NMF converged early.")
        res <- list(W = Winit, H = Hinit) # Return results
        return(res)
      }
      Vold <- Vnew # Update Vold for the next convergence check
    }
  }
  
  # If the loop finishes without converging within `niter` iterations
  if (iter == niter) {
    warning(message = "Maximum number of iterations reached without convergence. Consider increasing 'niter'.")
    res <- list(W = Winit, H = Hinit) # Return the current W and H
    return(res)
  }
}

# Optimized version of updateW function
updateW <- function (V, W, H, S, R, S2) {
  # R is the reconstructed matrix: W %*% H %*% S
  # S2 is a matrix of ones with the same dimensions as V
  
  # tHS <- t(H %*% S) # Pre-calculate t(H %*% S) for efficiency
  # num <- (V / R) %*% tHS # Numerator of the update rule
  # denom <- S2 %*% tHS # Denominator of the update rule
  # W <- W * num / denom # Multiplicative update rule for W
  # return(W)
  
  tHS_torch <- torch_t(torch_matmul(H@gm, S@gm))
  num_torch <- torch_matmul(torch_div(V@gm, R@gm), tHS_torch)
  denom_torch <- torch_matmul(S2@gm, tHS_torch)
  W_torch <- torch_mul(W@gm, torch_div(num_torch, denom_torch))
  return(gpu.matrix(W_torch))
}

# Optimized version of updateH function
updateH <- function (V, W, H, S, R, S1) {
  # R is the reconstructed matrix: W %*% H %*% S
  # S1 is (a matrix of ones with dimensions nrow(W) x nrow(S)) %*% S
  
  # num <- (t(W) %*% (V / R)) %*% t(S) # Numerator of the update rule
  # denom <- t(W) %*% S1 # Denominator of the update rule
  # H <- H * num / denom # Multiplicative update rule for H
  # return(H)
  
  denom_torch <- torch_matmul(torch_t(W@gm), S1@gm)
  num_torch <- torch_matmul(torch_matmul(torch_t(W@gm), torch_div(V@gm, R@gm)), torch_t(S@gm))
  H_torch <-torch_mul(H@gm, torch_div(num_torch, denom_torch))
  return(gpu.matrix(H_torch))

}


