generate <- function(conditions, condition_number, rep_set, rep, tau) {
  library(mvtnorm)
  library(covsim)
  library(LaplacesDemon)
  set.seed(1000 * condition_number + 100 * rep_set + rep)
  n <- as.integer(conditions[condition_number, 1])
  L <- as.integer(conditions[condition_number, 2]) # number of factors
  J <- as.integer(conditions[condition_number, 3]) # number of items
  fl <- as.double(conditions[condition_number, 4]) # factor loading
  xdist <- as.integer(conditions[condition_number, 5]) # observed score distribution
  tdist <- as.integer(conditions[condition_number, 6]) # latent variable distribution, excess kurtosis
  K <- ifelse(xdist <= 5, 2, 5) # number of categories
  #############################################################################
  # data generation
  #############################################################################
  ipf <- J / L
  a <- takane_simple(fl)$a
  if (K == 2) {
    among <- c(1, 2)
  } else {
    among <- c(1, 2, 3, 4, 5)
  }
  
  skew <- get_skewkurt(tdist)$skew
  kurt <- get_skewkurt(tdist)$kurt

  out <- matrix(vector("double", n * J ), nrow = n)
  sigma <- matrix(c(1, .6, .6, 1), nrow = 2)
  theta <- rIG(n, sigma, rep(skew, 2), rep(kurt, 2), typeA = "symm")[[1]]
  if (K == 2) {
    probs <- matrix(vector("double", n * J), nrow = n)
    for (i in 1:n) {
      for (j in 1:J) {
        if (L == 1) {
          probs[i, j] <- pnorm(a * theta[i, 1] - tau)
        } else {
          l <- floor((j - 1)/ipf) + 1 # What factor does the item load on
          if (l == 1) {
            probs[i, j] <- pnorm(a * theta[i, 1] - tau)
          } else {
            probs[i, j] <- pnorm(a * theta[i, 2] - tau)
          } 
        }
        out[i, j] <- sample(among, size = 1, replace = TRUE,
                            prob = c(1 - probs[i, j], probs[i, j]))
      }
    }
  } else {
    probs <- array(vector("double", n * J * 4), dim = c(n, J, 4))
    for (i in 1:n) {
      for (j in 1:J) {
        for (k in 1:4) {
          if (L == 1) {
            probs[i, j, k] <- pnorm(a * theta[i, 1] - tau[k])
          } else {
            l <- floor((j - 1)/ipf) + 1 # What factor does the item load on
            if (l == 1) {
              probs[i, j, k] <- pnorm(a * theta[i, 1] - tau[k])
            } else {
              probs[i, j, k] <- pnorm(a * theta[i, 2] - tau[k])
            } 
          }

        }
        if (L == 1) {
          probs[i, j, k] <- pnorm(a * theta[i, 1] - tau[k])
        } else {
          l <- floor((j - 1)/ipf) + 1 # What factor does the item load on

        }

        out[i, j] <- sample(among, size = 1, replace = TRUE,
                            prob = c(1 - probs[i, j, 1],
                                     probs[i, j, 1] - probs[i, j, 2],
                                     probs[i, j, 2] - probs[i, j, 3],
                                     probs[i, j, 3] - probs[i, j, 4],
                                     probs[i, j, 4]))
      }
    }
  }
  return(data.frame(out))
}