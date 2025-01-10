# The main function and helper functions to obtain oc with proposed method

rm(list=ls())
set.seed(315)

source("./scenarioGenerator.R")

print(paste0("The current w_e is ", w_e))

library(cmdstanr)
model <- cmdstan_model(stan_file = './model.stan')

show_progress <- function(message, newline = TRUE) {
  cat(sprintf("[%s] %s%s", 
              format(Sys.time(), "%H:%M:%S"),
              message,
              if(newline) "\n" else ""))
}

get_scenario_index <- function(scenario) {
  map_abb_to_word <- list(L = "low", M = "moderate", H = "high")
  scenario_chars <- strsplit(scenario, "")[[1]]
  scenario_words <- sapply(scenario_chars, function(x) map_abb_to_word[[x]])
  scenario_str <- paste(scenario_words, collapse = " ")
  which(sapply(scenario_list, function(x) paste(x, collapse = " ")) == scenario_str)
}

getoc <- function(
  ndose = ndose,               # Number of dose levels 
  startdose = startdose,       # Starting dose level
  target_tox = q_T,            # Target toxicity rate
  samplesize = samplesize,     # Total sample size
  cohortsize = csize,          # Cohort size
  we = w_e,                    # Efficacy weight for utility function
  wt = 1 - w_e,                # Toxicity weight for utility function
  
  tox_scenario = "MHM",    # Toxicity scenario (e.g., "MHM" for Moderate-High-Moderate)
  eff_scenario = "HLM",    # Efficacy scenario (e.g., "HLM" for High-Low-Moderate)
  
  warmup = 2000,           # Number of warmup iterations
  iter = 2000,             # Number of sampling iterations
  chains = 2,              # Number of MCMC chains
  thin = 1,                # Thinning factor
  
  interval_width = 0.1,   # Width of probability intervals
  N_star = 6,             # Cutoff for dose exploration control
  c_T = 0.95,             # Probability cutoff for dose elimination
  alpha_u = 1,            # Prior parameter alpha for utility
  beta_u = 1,             # Prior parameter beta for utility
  
  x = log(1:ndose/10),    # Pseudo ascending dose sequence
  a_prior = c(-5, 5),     # Prior range for intercept
  b_prior = c(0, 8)       # Prior range for slope
) {
  
  show_progress("Starting trial simulation...")
  show_progress(sprintf("Using toxicity scenario: %s, efficacy scenario: %s", tox_scenario, eff_scenario))
  
  if(we + wt != 1) stop("Utility weights must sum to 1")
  if(ndose < 2) stop("Number of dose levels must be greater than 1")
  
  show_progress("Generating trial scenarios...")
  trial_results <- generate_trial_scenarios(
    w_e = we, w_t = wt,
    D = ndose, startdose = startdose, 
    samplesize = samplesize, csize = cohortsize,
    alpha_u = alpha_u, beta_u = beta_u,
    q_T = target_tox,
    W = W, E = E,
    p_A_tox = p_A_tox, p_B_tox = p_B_tox, p_C_tox = p_C_tox,
    p_1_eff = p_1_eff, p_2_eff = p_2_eff, p_3_eff = p_3_eff
  )
  
  tox_idx <- get_scenario_index(tox_scenario)
  eff_idx <- get_scenario_index(eff_scenario)
  scenario_idx <- (tox_idx - 1) * 27 + eff_idx
  
  n <- rep(0, ndose)              # Number of patients at each dose
  d <- startdose                  # Current dose
  cohort_no <- 1                  # Current cohort number
  toxdose <- ndose + 1            # Highest non-toxic dose + 1
  st <- 0                         # Trial completion indicator
  st_early <- 0                   # Early stopping indicator
  EU_hat <- rep(0, ndose)         # Expected utilities
  
  tox_response <- vector("list", ndose)
  eff_response <- vector("list", ndose)
  qTP_list <- vector("list", ndose)
  qEP_list <- vector("list", ndose)
  Utility_list <- vector("list", ndose)
  EU_hat_list <- vector("list", ndose)
  
  tox_bounds <- seq(0, 1, interval_width)
  k_star <- which.max(as.numeric(table(cut(target_tox, tox_bounds))))
  
  show_progress("Beginning dose-finding trial...")
  show_progress(sprintf("Starting at dose level %d", startdose))
  
  # Function to handle dose selection when k_U values are tied, following S.1
  select_dose_with_tied_k_U <- function(k_U_candidates, d_current, adm_subset, n_d, n_upper, k_T, is_at_k_star = FALSE) {
    if(length(k_U_candidates) == 2) {
      if(all(k_U_candidates[1] == k_U_candidates[2])) {
        if(d_current == 1 && n_d >= N_star) {
          return(d_current)
        } else if(d_current == ndose) {
          prob_vector <- c(0.7, 0.3)  
          return(sample(adm_subset, 1, prob = prob_vector))
        } else {
          return(sample(adm_subset, 1))
        }
      } else {
        return(adm_subset[which.max(k_U_candidates)])
      }
    } else if(length(k_U_candidates) == 3) {
      if(n_d >= N_star) {
        if(k_U_candidates[1] == k_U_candidates[2]) {
          return(sample(adm_subset[1:2], 1))
        } else {
          return(adm_subset[which.max(k_U_candidates[1:2])])
        }
      } else {
        if(all(k_U_candidates[1] == k_U_candidates[2] && k_U_candidates[2] == k_U_candidates[3])) {
          if(n_upper == 0) {
            return(adm_subset[3])
          } else {
            prob_vector <- if(is_at_k_star) c(0.4, 0.4, 0.2) else c(0.3, 0.3, 0.4)
            return(sample(adm_subset, 1, prob = prob_vector))
          }
        } else if(k_U_candidates[1] == k_U_candidates[2] && k_U_candidates[2] > k_U_candidates[3]) {
          return(d_current)  
        } else if(k_U_candidates[1] < k_U_candidates[2] && k_U_candidates[2] == k_U_candidates[3]) {
          if(n_upper == 0) {
            return(adm_subset[3])
          } else {
            prob_vector <- if(is_at_k_star) c(0.7, 0.3) else c(0.5, 0.5)
            return(sample(adm_subset[2:3], 1, prob = prob_vector))
          }
        } else if(k_U_candidates[1] == k_U_candidates[3] && k_U_candidates[2] < k_U_candidates[3]) {
          if(n_upper == 0) {
            prob_vector <- c(0.4, 0.6)
          } else {
            prob_vector <- if(is_at_k_star) c(0.7, 0.3) else c(0.4, 0.6)
          }
          return(sample(adm_subset[c(1,3)], 1, prob = prob_vector))
        } else {
          return(adm_subset[which.max(k_U_candidates)])
        }
      }
    }
  }
  
  # Main trial loop
  while(st == 0) {
    show_progress(sprintf("Cohort %d: Treating %d patients at dose level %d...", 
                          cohort_no, cohortsize, d), 
                  newline = FALSE)
    
    if(is.null(tox_response[[d]])) {
      tox_response[[d]] <- list()
    } else if(length(tox_response[[d]]) > 0 && 
              floor(length(tox_response[[d]]) / cohortsize) + 1 == cohort_no) {
      tox_response[[d]] <- list()
    }
    
    if(is.null(eff_response[[d]])) {
      eff_response[[d]] <- list()
    } else if(length(eff_response[[d]]) > 0 && 
              floor(length(eff_response[[d]]) / cohortsize) + 1 == cohort_no) {
      eff_response[[d]] <- list()
    }
    
    for(i in 1:cohortsize) {
      tox_i <- matrix(NA, nrow = 3, ncol = 5)
      p_tox_d <- rbind(
        trial_results$scenarios[[scenario_idx]]$toxicity$p_tox_all[d,],
        trial_results$scenarios[[scenario_idx]]$toxicity$p_tox_all[d + ndose,],
        trial_results$scenarios[[scenario_idx]]$toxicity$p_tox_all[d + 2*ndose,]
      )
      tox_i <- t(apply(p_tox_d, 1, rmultinom, n = 1, size = 1))
      tox_response[[d]] <- rlist::list.append(tox_response[[d]], tox_i)
      
      eff_i <- matrix(NA, nrow = 3, ncol = 3)
      p_eff_d <- rbind(
        trial_results$scenarios[[scenario_idx]]$efficacy$p_eff_all[d,],
        trial_results$scenarios[[scenario_idx]]$efficacy$p_eff_all[d + ndose,],
        trial_results$scenarios[[scenario_idx]]$efficacy$p_eff_all[d + 2*ndose,]
      )
      eff_i <- t(apply(p_eff_d, 1, rmultinom, n = 1, size = 1))
      eff_response[[d]] <- rlist::list.append(eff_response[[d]], eff_i)
    }
    
    l <- max(d)
    for(j in 1:l) {
      for(i in 1:length(tox_response[[j]])) {
        qTP_list[[j]][i] <- sum(tox_response[[j]][[i]] * W) / sum(apply(W, 1, max))
      }
      for(i in 1:length(eff_response[[j]])) {
        qEP_list[[j]][i] <- sum(eff_response[[j]][[i]] * E) / sum(apply(E, 1, max))
      }
      for(i in 1:length(qEP_list[[j]])) {
        Utility_list[[j]][i] <- wt * (1 - qTP_list[[j]][i]) + we * qEP_list[[j]][i]
      }
    }
    
    n[d] <- n[d] + cohortsize
    
    qtp_sum <- unlist(sapply(qTP_list, sum))
    stan_data <- list(
      l = l,
      n = as.array(n[1:l]),
      sum_qtp = as.array(qtp_sum[1:l]),
      x = as.array(x[1:l]),
      a = a_prior[1],
      b = a_prior[2],
      c = b_prior[2]
    )
    
    suppressWarnings({
      suppressMessages({
        fit <- model$sample(
          data = stan_data,
          seed = 315,
          chains = chains,
          parallel_chains = chains,
          refresh = 2000,
          iter_sampling = iter,
          iter_warmup = warmup,
          thin = thin,
          show_messages = FALSE
        )
      })
    })
    
    cat(" Completed\n")  
    
    params <- fit$draws(format = "df")
    a_post <- params$alpha
    b_post <- params$beta
    
    show_progress("\nToxicity Safety Check:")
    EqTP <- 1 / (1 + exp(-a_post - b_post * x[d]))
    prob_too_toxic <- mean(EqTP > target_tox)
    
    show_progress(sprintf("Current dose level: %d", d))
    show_progress(sprintf("Target toxicity rate: %.3f", target_tox))
    show_progress(sprintf("Probability cutoff: %.3f", c_T))
    show_progress(sprintf("Current EqTP summary:"))
    show_progress(sprintf("  Mean: %.3f", mean(EqTP)))
    show_progress(sprintf("  Prob(EqTP > target): %.3f", prob_too_toxic))
    
    if(prob_too_toxic > c_T) {
      show_progress(sprintf("WARNING: Probability of excessive toxicity (%.3f) exceeds threshold (%.3f)", 
                            prob_too_toxic, c_T))
      
      if(d == 1) {
        if(length(qTP_list[[1]]) > 0) {
          show_progress(sprintf("Number of patients at dose 1: %d", n[1]))
        }
        
        st_early <- 1
        st <- 1
        show_progress("Trial stopped early due to excessive toxicity at lowest dose")
        break
      } else {
        toxdose <- d
        show_progress(sprintf("Doses above level %d eliminated for excess toxicity", d-1))
      }
    } else {
      show_progress("Toxicity check passed - continuing trial")
    }
    show_progress("")
    
    Mz_post <- as.numeric(table(cut(EqTP, tox_bounds)))
    k_T <- which.max(Mz_post)
    
    d_current <- d
    if(k_T > k_star || toxdose == d_current) {
      if(d_current > 1) d <- d_current - 1
      else d <- 1
    } else {
      u_sum <- unlist(sapply(Utility_list, sum))
      EU_hat[d_current] <- (alpha_u + u_sum[d_current]) / (alpha_u + beta_u + n[d_current])
      
      if(d_current == 1) {
        d_upper <- min(d_current + 1, toxdose - 1)
        adm_subset <- if(d_upper == 1) 1 else c(1, 2)
      } else if(d_current < toxdose - 1) {
        adm_subset <- c(d_current - 1, d_current, min(d_current + 1, toxdose - 1))
      } else {
        adm_subset <- c(d_current - 1, d_current)
      }
      
      k_U <- numeric(length(adm_subset))
      for(j in 1:length(adm_subset)) {
        util_sum <- sum(unlist(Utility_list[[adm_subset[j]]]))
        k_U[j] <- which.max(pbeta(seq(interval_width,1,interval_width), 
                                  alpha_u + util_sum,
                                  beta_u + n[adm_subset[j]] - util_sum) - 
                              pbeta(seq(0,1-interval_width,interval_width),
                                    alpha_u + util_sum,
                                    beta_u + n[adm_subset[j]] - util_sum))
      }
      
      n_upper <- if(d_current < ndose) n[d_current + 1] else 0
      
      d <- select_dose_with_tied_k_U(
        k_U_candidates = k_U,
        d_current = d_current,
        adm_subset = adm_subset,
        n_d = n[d_current],
        n_upper = n_upper,
        k_T = k_T,
        is_at_k_star = (k_T == k_star)
      )
    }

    
    if(st == 0) {  
      show_progress(sprintf("Next cohort will be treated at dose level %d", d))
    }
    
    cohort_no <- cohort_no + 1
    if(sum(n) >= samplesize) {
      show_progress("Maximum sample size reached, ending trial...")
      st <- 1
      break
    }
  }
  
  show_progress("Calculating final results...")
  
  EqTP_final <- numeric(ndose)
  k_EqTP_final <- numeric(ndose)
  
  for(i in 1:ndose) {
    EqTP_final[i] <- mean(1 / (1 + exp(-a_post - b_post * x[i])))
    k_EqTP_final[i] <- which.max(as.numeric(table(cut(EqTP_final[i], tox_bounds))))
  }
  
  MTD <- if(sum(k_EqTP_final > k_star) == ndose) {
    0
  } else {
    which.min(abs(EqTP_final[k_EqTP_final <= k_star] - target_tox))
  }
  
  OBD <- if(st_early == 1) {
    0
  } else {
    which.max(head(EU_hat, MTD))
  }
  
  show_progress(sprintf("Trial completed with %d cohorts", cohort_no - 1))
  show_progress(sprintf("MTD identified as dose level %d", MTD))
  show_progress(sprintf("OBD identified as dose level %d", OBD))
  
  if(st_early == 1) {
    show_progress("Note: Trial terminated early due to excessive toxicity at lowest dose")
  }
  
  show_progress("Final patient allocation by dose level:")
  show_progress(sprintf("Dose levels: %s", paste(1:ndose, collapse = " ")))
  show_progress(sprintf("Patients:    %s", paste(n, collapse = " ")))
  
  invisible(list(
    MTD = MTD,
    OBD = OBD,
    allocation = n,
    early_stop = st_early,
    final_utility = EU_hat,
    EqTP = EqTP_final,
    cohorts = cohort_no - 1,
    scenarios = list(
      toxicity = tox_scenario,
      efficacy = eff_scenario
    )
  ))
}

oc_results <- getoc(
  ndose = ndose,
  startdose = startdose,
  target_tox = q_T,
  samplesize = samplesize,
  cohortsize = csize,
  we = w_e,
  tox_scenario = "HHH",
  eff_scenario = "LLL"
)


