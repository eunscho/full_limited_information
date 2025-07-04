generate <- function(conditions, condition_number, rep_set, rep, tau) {
  # Load required libraries
  library(mvtnorm)       # For multivariate normal distributions
  library(covsim)        # For covariance matrix simulations (assumed for takane_simple)
  library(LaplacesDemon) # For non-normal distribution simulations (e.g., skew, kurtosis)
  
  # Set seed for reproducibility based on condition_number, rep_set, and rep
  set.seed(1000 * condition_number + 100 * rep_set + rep)
  
  # Extract values from the conditions matrix
  n <- as.integer(conditions[condition_number, 1]) # Sample size (number of respondents)
  L <- as.integer(conditions[condition_number, 2]) # Number of latent factors
  J <- as.integer(conditions[condition_number, 3]) # Number of items
  fl <- as.double(conditions[condition_number, 4]) # Factor loading value
  xdist <- as.integer(conditions[condition_number, 5]) # Observed score distribution type
  tdist <- as.integer(conditions[condition_number, 6]) # Latent variable distribution (excess kurtosis)
  
  # Determine number of response categories based on observed score distribution
  K <- ifelse(xdist <= 5, 2, 5) # 2 categories for binary, 5 for polytomous responses

  # Calculate items per factor (number of items divided by number of factors)
  ipf <- J / L
  
  # Generate factor loading matrix (assumes takane_simple from covsim package)
  a <- takane_simple(fl)$a # Factor loadings (strength of item-factor relationship)
  
  # Define response categories based on K
  if (K == 2) {
    among <- c(1, 2) # Binary responses (e.g., 1 = No, 2 = Yes)
  } else {
    among <- c(1, 2, 3, 4, 5) # Polytomous responses (e.g., Likert scale 1 to 5)
  }
  
  # Extract skewness and kurtosis for latent variable distribution
  # Assumes get_skewkurt is a user-defined function returning skew and kurtosis
  skew <- get_skewkurt(tdist)$skew
  kurt <- get_skewkurt(tdist)$kurt
  
  # Initialize output matrix (n x J) for simulated responses
  out <- matrix(vector("double", n * J), nrow = n)
  
  # Define covariance matrix for latent variables (assumes 2 factors with correlation 0.6)
  sigma <- matrix(c(1, .6, .6, 1), nrow = 2)
  
  # Generate latent variables (theta) with specified skew and kurtosis
  # rIG: Generates data from inverse Gaussian or similar distribution (LaplacesDemon)
  theta <- rIG(n, sigma, rep(skew, 2), rep(kurt, 2), typeA = "symm")[[1]]
  
  # Generate binary response data (K == 2)
  if (K == 2) {
    # Initialize probability matrix (n x J) for binary responses
    probs <- matrix(vector("double", n * J), nrow = n)
    
    # Calculate probabilities for each respondent (i) and item (j)
    for (i in 1:n) {
      for (j in 1:J) {
        # Single factor case: Use first latent variable (theta[i, 1])
        if (L == 1) {
          probs[i, j] <- pnorm(a * theta[i, 1] - tau) # IRT probability using normal CDF
        } else {
          # Multiple factors: Determine which factor the item loads on
          l <- floor((j - 1)/ipf) + 1 # Factor index for item j
          if (l == 1) {
            probs[i, j] <- pnorm(a * theta[i, 1] - tau) # Use first factor
          } else {
            probs[i, j] <- pnorm(a * theta[i, 2] - tau) # Use second factor
          }
        }
        # Sample binary response (1 or 2) based on calculated probability
        out[i, j] <- sample(among, size = 1, replace = TRUE,
                            prob = c(1 - probs[i, j], probs[i, j]))
      }
    }
    # Generate polytomous response data (K == 5)
  } else {
    # Initialize probability array (n x J x 4) for 4 thresholds
    probs <- array(vector("double", n * J * 4), dim = c(n, J, 4))
    
    # Calculate probabilities for each respondent (i), item (j), and threshold (k)
    for (i in 1:n) {
      for (j in 1:J) {
        for (k in 1:4) {
          if (L == 1) {
            # Single factor: Use first latent variable (theta[i, 1])
            probs[i, j, k] <- pnorm(a * theta[i, 1] - tau[k]) # Probability for threshold k
          } else {
            # Multiple factors: Determine which factor the item loads on
            l <- floor((j - 1)/ipf) + 1 # Factor index for item j
            if (l == 1) {
              probs[i, j, k] <- pnorm(a * theta[i, 1] - tau[k]) # Use first factor
            } else {
              probs[i, j, k] <- pnorm(a * theta[i, 2] - tau[k]) # Use second factor
            }
          }
        }
        # Sample polytomous response (1 to 5) based on threshold probabilities
        out[i, j] <- sample(among, size = 1, replace = TRUE,
                            prob = c(1 - probs[i, j, 1],
                                     probs[i, j, 1] - probs[i, j, 2],
                                     probs[i, j, 2] - probs[i, j, 3],
                                     probs[i, j, 3] - probs[i, j, 4],
                                     probs[i, j, 4]))
      }
    }
  }
  
  # Convert output matrix to data frame and return
  return(data.frame(out))
}
