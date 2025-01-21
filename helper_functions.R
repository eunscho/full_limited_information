takane_simple <- function(loading, tau = NULL, logistic = FALSE) {
  a <- loading / sqrt(1 - loading^2)
  if (is.null(tau)) {
    b <- NULL
  } else {
    b <- tau / sqrt(1 - loading^2)
  }
  if (logistic) {
    D <- 1.702
    a <- D * a
    if (!is.null(tau)) {
      b <- D * b
    }
  }
  return(list(a = a, bs = b))
}

tau2ratio <- function(taus, skew = 0, kurt = 0) {
  library(detectnorm)
  library(tidyverse)
  if (length(taus) == 1) {
    out <- vector("double", 2)
    if (skew == 0 & kurt == 0) {
      out[1] <- pnorm(taus)
      out[2] <- 1- out[1]
    } else {
      n <- 10^6
      x <- rnonnorm(n, skew = skew, kurt = kurt)$dat
      x <- as.data.frame(x)
      colnames(x) <- "x1"
      x_count <- x %>% 
        mutate(count = cut(x1, breaks = c(-Inf, taus, Inf))) %>% 
        count(count)
      out <- x_count[, "n"] / n
    }
  } else {
    out <- vector("double", 5)
    if (skew == 0 & kurt == 0) {
      out[1] <- pnorm(taus[1])
      for (i in 2:4) {
        out[i] <- pnorm(taus[i]) - pnorm(taus[i - 1])
      }
      out[5] <- 1- pnorm(taus[4])
    } else {
      n <- 10^6
      x <- rnonnorm(n, skew = skew, kurt = kurt)$dat
      x <- as.data.frame(x)
      colnames(x) <- "x1"
      x_count <- x %>% 
        mutate(count = cut(x1, breaks = c(-Inf, taus[1], taus[2], taus[3], taus[4], Inf))) %>% 
        count(count)
      out <- x_count[, "n"] / n
    }
  }
  return(out)
}

ratio2tau <- function(ratio, skew = 0, kurt = 0, lambda = .7) {
  if (skew == 0 & kurt == 0) {
    out <- qnorm(cumsum(ratio))[1:(length(ratio) - 1)]
  } else {
    library(covsim)
    n <- 10^7
    rep <- 1:10
    out_rep <- matrix(vector("double", length(rep) * (length(ratio) - 1)), nrow = length(rep))
    for (i in rep) {
      set.seed(i)
      x <- lambda * rIG(n, 
                        sigma.target = diag(2), 
                        skewness = rep(skew,2), 
                        excesskurtosis = rep(kurt, 2))[[1]][, 1] +
        rnorm(n, sd = sqrt(1 - lambda^2))
      out_rep[i, ] <- sort(x)[cumsum(ratio) * n][1:(length(ratio) - 1)]
      out <- apply(out_rep, 2, mean)
    }
  }
  out
}

tau2ratio <- function(taus, skew = 0, kurt = 0) {
  library(detectnorm)
  library(tidyverse)
  if (length(taus) == 1) {
    out <- vector("double", 2)
    if (skew == 0 & kurt == 0) {
      out[1] <- pnorm(taus)
      out[2] <- 1- out[1]
    } else {
      n <- 10^6
      x <- rnonnorm(n, skew = skew, kurt = kurt)$dat
      x <- as.data.frame(x)
      colnames(x) <- "x1"
      x_count <- x %>% 
        mutate(count = cut(x1, breaks = c(-Inf, taus, Inf))) %>% 
        count(count)
      out <- x_count[, "n"] / n
    }
  } else {
    out <- vector("double", 5)
    if (skew == 0 & kurt == 0) {
      out[1] <- pnorm(taus[1])
      for (i in 2:4) {
        out[i] <- pnorm(taus[i]) - pnorm(taus[i - 1])
      }
      out[5] <- 1- pnorm(taus[4])
    } else {
      n <- 10^6
      x <- rnonnorm(n, skew = skew, kurt = kurt)$dat
      x <- as.data.frame(x)
      colnames(x) <- "x1"
      x_count <- x %>% 
        mutate(count = cut(x1, breaks = c(-Inf, taus[1], taus[2], taus[3], taus[4], Inf))) %>% 
        count(count)
      out <- x_count[, "n"] / n
    }
  }
  return(out)
}

mplusfilename <- function(condition_number, rep_set, rep, estimator) {
  syntax <- paste0("c", condition_number, "_s", rep_set, "_r", rep, "_e", estimator)
  return(syntax)
}

dataSection <- function(condition_number, rep_set, rep, estimator) {
  syntax <- paste0(" DATA: \n FILE = \"", 
                   mplusfilename(condition_number, rep_set, rep, estimator),
                   ".dat\";\n")
  return(syntax)
}

savedata <- function(condition_number, rep_set, rep, estimator) {
  syntax <- paste0(" SAVEDATA: \n RESULTS =\"", 
                   mplusfilename(condition_number, rep_set, rep, estimator),
                   ".res\";")
  return(syntax)
}

# estimator 1: DWLS delta, 2: DWLS theta, 3: ULS delta, 4: ULS theta, 5: logistic
estimatorSection <- function(estimator) {
  syntax <- " ANALYSIS: \n"
  if (estimator == 1 | estimator == 2) {
    syntax <- paste0(syntax, " estimator = WLSMV;\n", collapse = "")
  } else if (estimator == 3 | estimator == 4) {
    syntax <- paste0(syntax, " estimator = ULSMV;\n", collapse = "")
  } else {
    syntax <- paste0(syntax, " estimator = mlr; link = probit;\n", collapse = "")
  }
  if (estimator == 1 | estimator == 3) {
    syntax <- paste0(syntax, " parameterization = delta;\n", collapse = "")
  } else if (estimator == 2 | estimator == 4) {
    syntax <- paste0(syntax, " parameterization = theta;\n", collapse = "")
  }
  return(syntax)
}

modelSection <- function(ipf, L, fl, estimator) {
  discrm <- fl / sqrt(1 - fl^2)
  if (estimator == 5 ) {
    starting <- 1.702 * discrm
  } else {
    starting <- discrm
  }
  syntax <- " MODEL: \n"
  if (L == 1) {
    syntax <- paste0(syntax, " f1 by x1-x", ipf, "*", starting, " (p", 1, "-p", ipf, "); \n",
                     "f1@1; \n [f1@0]; \n", 
                     "MODEL CONSTRAINT: \n", 
                     "p1 > 0; \n",
                     collapse = " ")
  } else {
    syntax <- paste0(syntax, " f1 by x1-x", ipf, "*", starting, " (p", 1, "-p", ipf, "); \n", 
                     " f2 by x", ipf + 1, "-x", 2 * ipf, "*", starting, " (p", ipf + 1, "-p", 2 * ipf, "); \n", 
                     # " f1-f3@1; \n [f1-f3@0]; \n",
                     " f1 with f2; \n",
                     " f1-f2@1; [f1-f2@0]; \n",
                     "MODEL CONSTRAINT: \n",
                     "p1 > 0; p", ipf + 1, " > 0;\n",
                     collapse = " ")
  }
  return(syntax)
}

getSyntax <- function(condition_number, rep_set, rep, L, J, fl, estimator) {
  ipf <- J / L
  syntax <- paste0(
    "TITLE:\n",
    " Item per factor = ", ipf, ", number of factor = ", L, ", estimator = ", estimator, "\n",
    dataSection(condition_number, rep_set, rep, estimator),
    " VARIABLE:\n",
    " NAMES ARE ", paste0("x", 1, "-x", (ipf * L), collapse=" "), ";\n",
    " CATEGORICAL ARE ", paste0("x", 1, "-x", (ipf * L), collapse=" "), ";\n",
    estimatorSection(estimator), 
    modelSection(ipf, L, fl, estimator),
    " OUTPUT: cinterval;\n",
    savedata(condition_number, rep_set, rep, estimator)
  )
  return(syntax)
}

fa2irt_est <- function(loading, thresh) {
  discrm <- loading / sqrt(1 - loading^2)
  if (length(dim(thresh)) == 2) {
    J <- nrow(thresh)
    K <- ncol(thresh) + 1
  } else {
    J <- length(thresh)
    thresh <- as.matrix(unlist(thresh))
    K <- 2
  }
  intercept <- matrix(vector("double", J * (K - 1)), nrow = J)
  for (j in 1:J) {
    for (k in 1:(K - 1)) {
      intercept[j, k] <- thresh[j, k] / sqrt(1 - loading[j]^2)
    }
  }
  return(list(discrm = discrm, intercept = intercept))
}

irt2fa_est <- function(discrm, intercept) {
  loading <- discrm / sqrt(1 + discrm^2)
  if (length(dim(intercept)) == 2) {
    J <- nrow(intercept) 
    K <- ncol(intercept) + 1
  } else {
    J <- length(intercept)
    intercept <- as.matrix(unlist(intercept))
    K <- 2
  }
  if (length(discrm) == 1) {
    discrm <- rep(discrm, J)
  }
  thresh <- matrix(vector("double", J * (K - 1)), nrow = J)
  for (j in 1:J) {
    for (k in 1:(K - 1)) {
      thresh[j, k] <- intercept[j, k] / sqrt(1 + discrm[j]^2)
    }
  }
  return(list(loading = loading, thresh = thresh))
}


fa2irt_se <- function(loading, thresh, se_loading, se_thresh) {
  denom <- (1 - loading ^ 2) ^ (3/2)
  se_a <- se_loading / denom
  se_d <- matrix(0, nrow = nrow(se_thresh), ncol = ncol(se_thresh))
  for (i in 1:nrow(se_thresh)) {
    for (j in 1:ncol(se_thresh)) {
      se_d[i, j] <- sqrt((1 - loading[i, 1] ^ 2) ^ 2 * se_thresh[i, j] ^ 2 + thresh[i, j] ^ 2 * loading[i, 1] ^ 2 * se_loading[i, 1] ^ 2) / denom[i, 1]
    }
  }
  return(list(se_a = se_a, se_d = se_d))
}

xdist2param <- function(xdist, fl, skew = 0, kurt = 0) {
  ratio <- list(c(.5, .5),
                c(.6, .4),
                c(.7, .3),
                c(.8, .2),
                c(.9, .1),
                c(.0503, .2101, .4792, .2101, .0503),
                c(.05, .10, .70, .10, .05),
                c(.5565, .2512, .0961, .0481, .0481),
                c(.25, .35, .23, .11, .06),
                c(.20, .20, .20, .20, .20))
  tau <- ratio2tau(ratio[[xdist]], skew = skew, kurt = kurt, lambda = fl)
  out <- takane_simple(fl, tau, logistic = F)
  return(out)
}

get_skewkurt <- function(tdist = 1) {
  mat <- matrix(c(0, 0, # normal
                  0, -1, # platy
                  0, 3.75, # lepto
                  1.25, 0, # skewed
                  1.25, 3.75, # moderate
                  2.5, 7.5), # considerable
                byrow = T, ncol = 2)
  list(skew = mat[tdist, 1], kurt = mat[tdist, 2])
}

calc_tau <- function() {
  xdist <- seq(1, 10)
  tdist <- seq(1, 6)
  lambdas <- c(.5, .8)
  
  out <- list()
  for (i in seq_along(xdist)) {
    for (j in seq_along(tdist)) {
      for (k in seq_along(lambdas)) {
        K <- ifelse(xdist[i] <= 5, 2, 5)
        skew <- get_skewkurt(tdist[j])$skew
        kurt <- get_skewkurt(tdist[j])$kurt
        params <- xdist2param(xdist[i], lambdas[k], skew, kurt)
        a <- params$a
        taus <- params$bs
        
        t1 <- if (K == 2) taus[1] else taus[1]
        t2 <- if (K == 2) NA else taus[2]
        t3 <- if (K == 2) NA else taus[3]
        t4 <- if (K == 2) NA else taus[4]
        
        out <- append(out, list(data.frame(
          xdist = xdist[i],
          tdist = tdist[j],
          lambda = lambdas[k],
          a = a,
          t1 = t1,
          t2 = t2,
          t3 = t3,
          t4 = t4
        )))
      }
    }
  }
  out <- do.call(rbind, out)
  write.csv(out, "taus.csv")
  return(out)
}

get_tau <- function(fl, xdist, tdist) {
  library(tidyverse)
  tau_res <- read_csv("./summary/taus.csv")
  tau <- tau_res %>% filter(lambda == fl & xdist == !!xdist & tdist == !!tdist) %>% 
    dplyr::select(t1:t4)
  if (xdist <= 5) {
    tau <- tau[, 1]
  }
  return(unlist(tau))
}

dIG <- function(x, sigma = 1, skewness = 0, excesskurtosis = 0) {
  # This function is a modification of the code written by Professor Foldnes
  library(PearsonDS)
  out <- vector("double",length(x))
  for (i in 1:length(x)) {
    IG_params <- function(sigma, skewness, excesskurtosis, typeA=c("symm", "triang") ){
      sigma.target <-as.matrix(sigma)
      nvar          <- dim(sigma.target)[2]
      #define functions
      function.skew <- function(IGvalues){
        fval <- numeric(nvar)
        for (i in 1:nvar)
          fval[i] <-   A[i, ]^3 %*% IGvalues/(sum(A[i,]^2)^(3/2))
        
        fval-skewness
      }
      function.kurt <- function(IGvalues){
        fval <- numeric(nvar)
        for (i in 1:nvar)
          fval[i] <-   A[i, ]^4 %*% IGvalues/(sum(A[i,]^2)^(2))
        
        fval-excesskurtosis
      }
      #calculate A
      typeA <- match.arg(typeA)
      if(typeA=="triang")
        A  <- t(chol(sigma.target))
      else
        A  <- lavaan::lav_matrix_symmetric_sqrt(sigma.target)#symmetric
      IGskew        <- nleqslv::nleqslv(x=skewness, function.skew)$x
      IGkurt.excess <- nleqslv::nleqslv(x=excesskurtosis, function.kurt)$x
      parlist       <- list()
      for (i in 1:nvar) parlist[[i]] <- PearsonDS::pearsonFitM(moments=c(mean=0, variance=1, skewness=IGskew[i], 3+IGkurt.excess[i]))
      
      return(list(parlist=parlist, A=A, Ainv=solve(A), Adet = det(A)))
    }
    
    speclist <- IG_params(sigma, skewness, excesskurtosis)
    if (is.vector(x[i]))
      x[i] <- matrix(x[i], nrow = length(x[i]))
    A  <- speclist$A
    Ainv <- speclist$Ainv
    Adet <- speclist$Adet
    parlist <- speclist$parlist
    
    xtransformed <- Ainv %*% x[i]
    
    # the univariate densities
    uni <- lapply(1:length(parlist), function(dim){
      dpearson(xtransformed[dim,], params=parlist[[dim]])
    })
    uni <- do.call(rbind,uni)
    if (is.vector(uni))
      uni <- matrix(uni, nrow = length(uni))
    uniprod <- apply(uni, 2, prod)
    out[i] <- matrix(uniprod/Adet, nrow=1)[[1]]
  }
  return(out)
  
}
