#' Parallel Simulation Function
#'
#' Runs the simulation in parallel across multiple conditions.
#'
#' @param start_con Starting condition number
#' @param end_con Ending condition number
#' @param start_rep Starting replication set
#' @param end_rep Ending replication set
#' @export
sim3 <- function(start_con = 1, end_con = 960) {
  # Define simulation parameters
  n <- c(200, 1000)
  L <- c(1, 2)          # Number of latent variables
  J <- c(6, 20)         # Number of items
  fl <- c(0.5, 0.8)     # Factor loading
  xdist <- seq(1, 10)   # Observed score distribution
  tdist <- seq(1, 6)    # Latent variable distribution, excess kurtosis
  
  # Generate all combinations of conditions
  conditions <- tidyr::crossing(n, L, J, fl, xdist, tdist)
  
  # Define the range of condition numbers to process
  condition_numbers <- start_con:end_con
  
  # Define the simulation name prefix
  name <- "sim3_"
  
  # Define the range of replications
  REPS <- 1000
  ##############################################################################
  # analyze function
  ##############################################################################
  analyze <- function(conditions, condition_number, tau, REPS) {
    # Load necessary libraries
    library(MplusAutomation)
    library(texreg)
    library(tidyverse)
    library(tictoc)
    library(purrr)
    
    # Extract simulation conditions
    n <- as.integer(conditions[condition_number, "n"])
    L <- as.integer(conditions[condition_number, "L"])    # Number of factors
    J <- as.integer(conditions[condition_number, "J"])    # Number of items
    fl <- as.double(conditions[condition_number, "fl"])    # Factor loading
    xdist <- as.integer(conditions[condition_number, "xdist"])  # Observed score distribution
    tdist <- as.integer(conditions[condition_number, "tdist"])  # Latent variable distribution, excess kurtosis
    K <- ifelse(xdist <= 5, 2, 5)                            # Number of categories
    discrm <- takane_simple(fl)$a                             # Discrimination parameter (e.g., Takane method)
    
    #############################################################################
    # Define the analyze_estimator function
    # This function processes Mplus output to extract and calculate necessary parameters
    #############################################################################
    analyze_estimator <- function(L, J, discrm, K, estimator, tau, coef, ci) {
      D <- 1.702
      # original values
      aest_org <- coef %>% filter(stringr::str_detect(paramHeader, "BY")) %>% dplyr::select(est)
      dest_org <- coef %>% filter(paramHeader == "Thresholds") %>% dplyr::select(est)
      ase_org <- coef %>% filter(stringr::str_detect(paramHeader, "BY")) %>% dplyr::select(se)
      dse_org <- coef %>% filter(paramHeader == "Thresholds") %>% dplyr::select(se)
      if (K == 5) {
        dest_org <- matrix(unlist(dest_org), ncol = length(tau), byrow = T)
        dse_org <- matrix(unlist(dse_org), ncol = length(tau), byrow = T)
      }
      if (L == 2) {
        rest_org <- coef %>% filter(stringr::str_detect(paramHeader, "WITH")) %>% dplyr::select(est)
        rse_org <- coef %>% filter(stringr::str_detect(paramHeader, "WITH")) %>% dplyr::select(se)
      } else {
        rest_org <- rse_org <- NA
      }
      
      alow <- ci %>% filter(stringr::str_detect(paramHeader, "BY")) %>% dplyr::select(low2.5)
      aup <- ci %>% filter(stringr::str_detect(paramHeader, "BY")) %>% dplyr::select(up2.5)
      dlow <- ci %>% filter(paramHeader == "Thresholds") %>% dplyr::select(low2.5)
      dup <- ci %>% filter(paramHeader == "Thresholds") %>% dplyr::select(up2.5)
      if (L == 2) {
        rlow <- ci %>% filter(stringr::str_detect(paramHeader, "WITH")) %>% dplyr::select(low2.5)
        rup <- ci %>% filter(stringr::str_detect(paramHeader, "WITH")) %>% dplyr::select(up2.5)
      }
      
      if (K == 5) {
        dlow <- matrix(unlist(dlow), ncol = length(tau), byrow = T)
        dup <- matrix(unlist(dup), ncol = length(tau), byrow = T)
      }
      if (K == 5) {
        taus <- matrix(rep(tau, J), byrow = T, nrow = J)
      } else {
        taus <- tau
      }
      
      # inadmissible solutions
      if (estimator == 1 | estimator == 3) {
        aest_fa <- aest_org
      } else {
        aest_fa <- aest_org / sqrt(1 + aest_org^2)
      } 
      inadmissible <- ifelse(any(aest_fa > .999), TRUE, FALSE)
      
      if (inadmissible) {
        aest <- ase <- acov <- NA
        dest <- dse <- dcov <- NA
        rest <- rse <- rcov <- NA
      } else {
        if (estimator == 1 | estimator == 3) {
          fa_param <- irt2fa_est(discrm, taus)
          a <- loading <- fa_param$loading
          d <- thresh <- as.double(fa_param$thresh)
        } else if (estimator == 2 | estimator == 4) {
          a <- discrm
          d <- taus
        } else if (estimator == 5) {
          a <- discrm
          d <- taus
        } 
        acov <- mean(ifelse(a < alow | a > aup, 0, 1))
        dcov <- mean(ifelse(d < dlow | d > dup, 0, 1))
        if (L == 2) {
          rcov <- mean(ifelse(.6 < rlow | .6 > rup, 0, 1))  
        } else {
          rcov <- NA
        }
        
        
        if (estimator == 1 | estimator == 3) {
          est <- fa2irt_est(unlist(aest_org), unlist(dest_org))
          aest <- data.frame(est$discrm - discrm)
          dest <- data.frame(est$intercept - taus)
          se <- fa2irt_se(aest_org, dest_org, ase_org, dse_org)
          ase <- data.frame((se$se_a))
          dse <- data.frame((se$se_d))
        } else  {
          aest <- aest_org - discrm
          dest <- dest_org - taus
          ase <- ase_org
          dse <- dse_org
        } 
        rest <- rest_org - .6
        rse <- rse_org
      }
      
      out <- list(aest = aest, 
                  ase = ase, acov = acov, 
                  dest = dest, 
                  dse = dse, dcov = dcov,
                  rest = rest, rse = rse, rcov = rcov)
      return(out)
    }
    
    #############################################################################
    # Initialize the results matrix as a matrix of lists
    #############################################################################
    res <- matrix(list(), nrow = REPS, ncol = 5)
    
    #############################################################################
    # Simulation Loop
    #############################################################################
    for (rep in 1:REPS) {
      for (estimator in 1:5) {
        # Generate filenames for Mplus
        mplusfname <- mplusfilename(condition_number, rep_set = 1, rep, estimator)
        datafilename <- paste0(mplusfname, ".dat")
        inputfilename <- paste0(mplusfname, ".inp")
        outputfilename <- paste0(mplusfname, ".out")
        resultfilename <- paste0(mplusfname, ".res")
        
        # Generate and prepare data for Mplus
        data <- generate(conditions, condition_number, rep_set = 1, rep, tau)
        prepareMplusData(data, datafilename)
        
        # Generate Mplus syntax
        syntax <- getSyntax(condition_number, rep_set = 1, rep, L, J, fl, estimator)
        mclm::write_txt(syntax, inputfilename)
        
        # Run Mplus model
        tic()
        runModels(inputfilename)
        toctime <- toc()
        time <- toctime$toc - toctime$tic
        
        # Read Mplus results
        mplusres <- readModels(outputfilename)
        converge <- length(mplusres$warnings) == 0
        
        # Extract coefficients and confidence intervals
        coef <- data.frame(mplusres$parameters$unstandardized)
        ci <- mplusres$parameters$ci.unstandardized
        
        # Safely analyze the estimator
        safe_analyze_estimator <- safely(analyze_estimator)
        res_safe <- safe_analyze_estimator(L, J, discrm, K, estimator, tau, coef, ci)
        
        # Check for errors and convergence
        if (is.null(res_safe$error) & converge) {
          aest <- res_safe$result$aest
          ase <- res_safe$result$ase
          acov <- res_safe$result$acov
          dest <- res_safe$result$dest
          dse <- res_safe$result$dse
          dcov <- res_safe$result$dcov
          rest <- res_safe$result$rest
          rse <- res_safe$result$rse
          rcov <- res_safe$result$rcov
        } else {
          aest <- dest <- ase <- dse <- NA
          acov <- dcov <- rcov <- rest <- rse <- NA
        }
        
        # Determine if the results are proper (not NA)
        proper <- !any(is.na(aest))
        
        # Assign the results to the res matrix using double brackets
        res[[rep, estimator]] <- list(
          estimator = estimator, 
          rep = rep, 
          time = time, 
          proper = proper,
          aest = aest, 
          acov = acov, 
          ase = ase, 
          dest = dest, 
          dcov = dcov, 
          dse = dse,
          rest = rest, 
          rse = rse, 
          rcov = rcov
        )
        
        # Remove temporary files
        file.remove(datafilename)
        file.remove(inputfilename)
        file.remove(outputfilename)
        file.remove(resultfilename)
        
        # Remove unnecessary objects from memory
        rm(data, syntax, mplusres, coef, ci, res_safe)
        
        # Force garbage collection
        gc()
      }
    }
    
    #############################################################################
    # Summarizing Results with Exclusion Strategies
    #############################################################################
    # Initialize lists to store results for partial and full exclusion
    partial_out <- list()
    full_out <- list()
    
    # Determine for each replication if any estimator failed
    for (rep in 1:REPS) {
      # Check if any estimator failed in this replication
      any_failed <- any(sapply(1:5, function(est) {
        if (is.null(res[[rep, est]]$proper)) {
          return(TRUE)
        }
        return(!res[[rep, est]]$proper | any(is.na(res[[rep, est]]$aest)))
      }))
      
      if (any_failed) {
        # Replication has at least one failure
        # Full exclusion applies: exclude all estimators for this replication
        for (estimator in 1:5) {
          partial_out[[length(partial_out) + 1]] <- res[[rep, estimator]]
          full_out[[length(full_out) + 1]] <- list(
            estimator = estimator,
            rep = rep,
            time = res[[rep, estimator]]$time,
            proper = FALSE,  # Mark as excluded
            aest = NA,
            acov = NA,
            ase = NA,
            dest = NA,
            dcov = NA,
            dse = NA,
            rest = NA,
            rse = NA,
            rcov = NA
          )
        }
      } else {
        # All estimators succeeded in this replication
        for (estimator in 1:5) {
          partial_out[[length(partial_out) + 1]] <- res[[rep, estimator]]
          full_out[[length(full_out) + 1]] <- res[[rep, estimator]]
        }
      }
    }
    
    #############################################################################
    # Function to Summarize Results
    #############################################################################
    
    summarize_results <- function(results, exclusion_flag) {
      # Initialize empty tibble
      summarized <- tibble()
      
      for (estimator in 1:5) {
        # Extract results for this estimator
        est_results <- lapply(results, function(x) {
          if (x$estimator == estimator) {
            return(x)
          } else {
            return(NULL)
          }
        })
        est_results <- Filter(Negate(is.null), est_results)
        
        if (length(est_results) == 0) {
          # No data for this estimator under this exclusion strategy
          next
        }
        
        # Extract variables safely
        proper_temp <- sapply(est_results, function(x) if (!is.null(x$proper)) x$proper else FALSE) 
        aest_temp <- sapply(est_results, function(x) if (!is.null(x$aest)) x$aest else NA)
        ase_temp <- sapply(est_results, function(x) if (!is.null(x$ase)) x$ase else NA)
        acov_temp <- sapply(est_results, function(x) if (!is.null(x$acov)) x$acov else NA)
        dest_temp <- sapply(est_results, function(x) if (!is.null(x$dest)) x$dest else NA)
        dse_temp <- sapply(est_results, function(x) if (!is.null(x$dse)) x$dse else NA)
        dcov_temp <- sapply(est_results, function(x) if (!is.null(x$dcov)) x$dcov else NA)
        rest_temp <- sapply(est_results, function(x) if (!is.null(x$rest)) x$rest else NA)
        rse_temp <- sapply(est_results, function(x) if (!is.null(x$rse)) x$rse else NA)
        rcov_temp <- sapply(est_results, function(x) if (!is.null(x$rcov)) x$rcov else NA)
        
        # Calculate proportion of proper solutions (including failures)
        total_cases <- length(proper_temp)
        proper_prop <- sum(proper_temp, na.rm = TRUE) / total_cases
        
        # Convert to numeric, handling cases where multiple values exist
        aest_temp_numeric <- unlist(aest_temp)
        ase_temp_numeric <- unlist(ase_temp)
        dest_temp_numeric <- unlist(dest_temp)
        dse_temp_numeric <- unlist(dse_temp)
        rest_temp_numeric <- unlist(rest_temp)
        rse_temp_numeric <- unlist(rse_temp)
        
        rb <- function(a, b) {
          (a - b) / b
        }
        # Calculate summary statistics, ensuring no errors due to all NAs
        b_aest <- if (all(is.na(aest_temp_numeric))) NA else mean(aest_temp_numeric, na.rm = TRUE)
        rm_aest <- if (all(is.na(aest_temp_numeric))) NA else sqrt(mean(aest_temp_numeric^2, na.rm = TRUE))
        rb_ase <- if (all(is.na(ase_temp_numeric)) | all(is.na(aest_temp_numeric))) NA else rb(mean(ase_temp_numeric, na.rm = TRUE), sd(aest_temp_numeric, na.rm = TRUE))
        acov <- if (all(is.na(acov_temp))) NA else mean(acov_temp, na.rm = TRUE)
        b_dest <- if (all(is.na(dest_temp_numeric))) NA else mean(dest_temp_numeric, na.rm = TRUE)
        rb_dest <- if (all(is.na(dest_temp_numeric))) NA else mean(dest_temp_numeric / tau, na.rm = TRUE)
        rm_dest <- if (all(is.na(dest_temp_numeric))) NA else sqrt(mean(dest_temp_numeric^2, na.rm = TRUE))
        rb_dse <- if (all(is.na(dse_temp_numeric)) | all(is.na(dest_temp_numeric))) NA else rb(mean(dse_temp_numeric, na.rm = TRUE), sd(dest_temp_numeric, na.rm = TRUE))
        dcov <- if (all(is.na(dcov_temp))) NA else mean(dcov_temp, na.rm = TRUE)
        b_rest <- if (all(is.na(rest_temp_numeric))) NA else mean(rest_temp_numeric, na.rm = TRUE)
        rm_rest <- if (all(is.na(rest_temp_numeric))) NA else sqrt(mean(rest_temp_numeric^2, na.rm = TRUE))
        rb_rse <- if (all(is.na(rse_temp_numeric)) | all(is.na(rest_temp_numeric))) NA else rb(mean(rse_temp_numeric, na.rm = TRUE), sd(rest_temp_numeric, na.rm = TRUE))
        rcov <- if (all(is.na(rcov_temp))) NA else mean(rcov_temp, na.rm = TRUE)
        
        # Handle 'tau' as separate columns: tau1, tau2, tau3, tau4
        tau1 <- ifelse(length(tau) >= 1, tau[1], NA)
        tau2 <- ifelse(length(tau) >= 2, tau[2], NA)
        tau3 <- ifelse(length(tau) >= 3, tau[3], NA)
        tau4 <- ifelse(length(tau) >= 4, tau[4], NA)
        
        # Create a tibble for the current estimator and exclusion strategy
        out_temp <- tibble(
          condition_number = condition_number,
          n = n,
          L = L,
          J = J,
          discrm = discrm,
          xdist = xdist,
          tdist = tdist, 
          tau1 = tau1,
          tau2 = tau2,
          tau3 = tau3,
          tau4 = tau4,
          estimator = estimator,
          exclusion = exclusion_flag,
          proper = proper_prop,
          b_aest = b_aest,
          rm_aest = rm_aest,
          rb_ase = rb_ase,
          acov = acov, 
          b_dest = b_dest,
          rb_dest = rb_dest,
          rm_dest = rm_dest,
          rb_dse = rb_dse,
          dcov = dcov, 
          b_rest = b_rest,
          rm_rest = rm_rest,
          rb_rse = rb_rse,
          rcov = rcov
        )
        
        # Combine with summarized tibble
        summarized <- bind_rows(summarized, out_temp)
      }
      
      return(summarized)
      
      rm(aest, ase, acov, dest, dse, dcov, rest, rse, rcov)
      gc()
    }
    
    #############################################################################
    # Summarize Results for Partial and Full Exclusion
    #############################################################################
    
    # Summarize Partial Exclusion Results (exclusion = 0)
    partial_summarized <- summarize_results(partial_out, exclusion_flag = 0)
    
    # Summarize Full Exclusion Results (exclusion = 1)
    full_summarized <- summarize_results(full_out, exclusion_flag = 1)
    
    #############################################################################
    # Combine Both Summaries to Create Final Summarized Results
    #############################################################################
    final_summarized <- bind_rows(partial_summarized, full_summarized)
    
    # Return the summarized results
    return(final_summarized)
  }

  ##############################################################################
  # MAIN SIMULATION
  ##############################################################################
  # Ensure the output directory exists; if not, create it
  output_dir <- "C:/Users/bene/Documents/Simulations/sim3/res"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Set up the parallel processing plan
  # Reserve two core for other tasks by setting workers = availableCores() - 2
  num_workers <- parallel::detectCores() - 2
  future::plan(future::multisession, workers = num_workers)
  
  # Start parallel processing across conditions
  future.apply::future_lapply(condition_numbers, function(condition_number) {
    # Print the current condition being processed
    print(conditions[condition_number, ])
    
    # Start timing for the condition
    tictoc::tic()
    print(paste("Starting condition number", condition_number))
    
    # Define the filename for saving summarized results
    filename <- file.path(output_dir, paste0(name, condition_number, ".csv"))
    
    # Check if the result file already exists to avoid redundant computations
    if (!file.exists(filename)) {
      # Extract parameters for the current condition
      fl_val <- as.double(conditions[condition_number, "fl"])
      xdist_val <- as.integer(conditions[condition_number, "xdist"])
      tdist_val <- as.integer(conditions[condition_number, "tdist"])
      
      # Calculate tau using the user-defined function (ensure this function is defined)
      tau <- get_tau(fl_val, xdist_val, tdist_val)
      
      # Call the analyze function with all replications at once
      summarized_results <- analyze(conditions, condition_number, tau, REPS)
      
      write.csv(summarized_results, file = filename)

      # Print the summarized results for verification
      print(summarized_results)
      tictoc::toc()
    } else {
      print(paste("File", filename, "already exists. Skipping..."))
    }
    
    return(NULL)  # Return NULL as results are saved to disk
  },
  future.chunk.size = 1
  ) # End of parallel processing across conditions
  
  # Optionally, return all_results if needed (currently NULLs)
  return(NULL)
}
