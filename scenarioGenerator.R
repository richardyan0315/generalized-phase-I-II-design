rm(list=ls())

library(rlist)

set.seed(315)

{
  # toxicity-efficacy trade-offs: 
  # when w_e > w_t: prefer high efficacy
  # when w_e < w_t: prefer low toxicity
  
  w_e <- 0.7 # 0.45 / 0.7
  w_t <- 1 - w_e
  
  D <- ndose <- 6       # total dose level number
  x <- log(1:ndose/10)  # pseudo ascending dose sequence
  startdose <- 1        # starting dose level
  samplesize <- 27      # total sample size
  csize <- 3            # cohort size
  
  # prior parameters of EU (Expected Utility) ~ Beta(alpha_u, beta_u)
  alpha_u <- 1
  beta_u <- 1
  
  # Over Toxicity Threshold: Target toxicity probability
  q_T <- 0.3
  
  # severity matrix of toxicity for each type at each grade (grade 0 to 4)
  W <- rbind(
    c(0.1, 0.8, 1.3, 1.8, 2.3),
    c(0.07, 0.5, 1, 1.5, 2),
    c(0.02, 0.05, 0.25, 0.5, 1)
  )
  
  # efficacious weight matrix for each type at each grade (low / middle / high)
  E <- rbind(
    c(0.05, 0.2, 0.47),
    c(0.07, 0.14, 0.41),
    c(0.03, 0.15, 0.45)
  )
  
  # initial probabilities at each grade of each type of toxicity (dose 1 prob)
  p_A_tox <- c(0.891, 0.072, 0.032, 0.004, 0.001)
  p_B_tox <- c(0.868, 0.062, 0.031, 0.029, 0.010)
  p_C_tox <- c(0.817, 0.071, 0.057, 0.050, 0.005)
  
  # initial probabilities at each grade of each type of efficacy (dose 1 prob)
  p_1_eff <- c(0.900, 0.070, 0.030)
  p_2_eff <- c(0.850, 0.130, 0.020)
  p_3_eff <- c(0.870, 0.090, 0.040)
  
  # Define all possible scenario combinations
  scenario_list <- list(
    c("low","low","low"), c("low","low","moderate"), c("low","low","high"),
    c("low","moderate","low"), c("low","moderate","moderate"), c("low","moderate","high"),
    c("low","high","low"), c("low","high","moderate"), c("low","high","high"),
    c("moderate","low","low"), c("moderate","low","moderate"), c("moderate","low","high"),
    c("moderate","moderate","low"), c("moderate","moderate","moderate"), c("moderate","moderate","high"),
    c("moderate","high","low"), c("moderate","high","moderate"), c("moderate","high","high"),
    c("high","low","low"), c("high","low","moderate"), c("high","low","high"),
    c("high","moderate","low"), c("high","moderate","moderate"), c("high","moderate","high"),
    c("high","high","low"), c("high","high","moderate"), c("high","high","high")
  )
}

{
  prefix_scenario <- data.frame()
  MTD_true_current <- numeric()
  
  early_stop_rate <- data.frame()          # early stop rate for each scenario
  poor_allocation_rate <- data.frame()     # percentage of trials with poor dose allocation
  average_patient_number <- data.frame()   # average number of patients per scenario
  average_overdose_number <- data.frame()  # average number of patients receiving overdose
  average_selection_percentage <- data.frame()  # percentage of correct dose selection
  
  rownames_ESR <- vector()  # Early Stop Rate
  rownames_APN <- vector()  # Average Patient Number
  rownames_ASP <- vector()  # Average Selection Percentage
  rownames_csel_MTD_all <- vector()  # Correct Selection of MTD
  rownames_csel_OBD_all <- vector()  # Correct Selection of OBD
  rownames_PA <- vector()  # Poor Allocation
  rownames_OD <- vector()  # Overdose
  
  # DataFrame for storing correct OBD selection rates
  csel_OBD_all <- data.frame()
  
  # Initialize OBD_true matrix (27 x 27 scenarios)
  OBD_true <- matrix(NA, nrow = 27, ncol = 27)
  
  # Set row and column names for OBD_true matrix
  rownames(OBD_true) <- c("Tox_LLL","Tox_LLM","Tox_LLH","Tox_LML","Tox_LMM","Tox_LMH",
                          "Tox_LHL","Tox_LHM","Tox_LHH","Tox_MLL","Tox_MLM","Tox_MLH",
                          "Tox_MML","Tox_MMM","Tox_MMH","Tox_MHL","Tox_MHM","Tox_MHH",
                          "Tox_HLL","Tox_HLM","Tox_HLH","Tox_HML","Tox_HMM","Tox_HMH",
                          "Tox_HHL","Tox_HHM","Tox_HHH")
  colnames(OBD_true) <- c("Eff_LLL","Eff_LLM","Eff_LLH","Eff_LML","Eff_LMM","Eff_LMH",
                          "Eff_LHL","Eff_LHM","Eff_LHH","Eff_MLL","Eff_MLM","Eff_MLH",
                          "Eff_MML","Eff_MMM","Eff_MMH","Eff_MHL","Eff_MHM","Eff_MHH",
                          "Eff_HLL","Eff_HLM","Eff_HLH","Eff_HML","Eff_HMM","Eff_HMH",
                          "Eff_HHL","Eff_HHM","Eff_HHH")
  
  # Lists for storing true probabilities and utilities
  true_tox_prob_list <- vector("list", length = 27)  # True toxicity probability
  true_eff_prob_list <- vector("list", length = 27 * 27)  # True efficacy probability
  EU_true_list <- vector("list", length = 27 * 27)  # True expected utility
  curve_type_eff <- vector("list", length = 27 * 27)  # Efficacy curve types
}

## Define transition matrices and intervals
{
  # Transition matrices for toxicity in various speeds
  # L: low speed, M: moderate speed, H: high speed
  trans_matrix_tox_L <- as.matrix(rbind(
    c(0.9,0.1,0,0,0),
    c(0,0.9,0.1,0,0),
    c(0,0,0.9,0.1,0),
    c(0,0,0,0.9,0.1),
    c(0,0,0,0,1)
  ))
  
  trans_matrix_tox_M <- as.matrix(rbind(
    c(0.75,0.25,0,0,0),
    c(0,0.75,0.25,0,0),
    c(0,0,0.75,0.25,0),
    c(0,0,0,0.75,0.25),
    c(0,0,0,0,1)
  ))
  
  trans_matrix_tox_H <- as.matrix(rbind(
    c(0.45,0.55,0,0,0),
    c(0,0.45,0.55,0,0),
    c(0,0,0.45,0.55,0),
    c(0,0,0,0.45,0.55),
    c(0,0,0,0,1)
  ))
  
  # Transition matrices for efficacy in various speeds
  trans_matrix_eff_L <- as.matrix(rbind(
    c(0.9,0.1,0),
    c(0,0.9,0.1),
    c(0,0,1)
  ))
  
  trans_matrix_eff_M <- as.matrix(rbind(
    c(0.75,0.25,0),
    c(0,0.75,0.25),
    c(0,0,1)
  ))
  
  trans_matrix_eff_H <- as.matrix(rbind(
    c(0.45,0.55,0),
    c(0,0.45,0.55),
    c(0,0,1)
  ))
  
  # uTPI interval settings (including Sensitivity Analysis)
  width_SA <- 0.10  # interval width: 0.05/0.1(default)/0.2
  
  # Equal-width toxicity intervals
  int.bound <- seq(0, 1, 0.1)
  # Location of the highest tolerable toxicity interval
  k_star <- which.max(as.numeric(table(cut(q_T, int.bound))))
  
  # Equal-width efficacy intervals
  u.lower <- seq(0.0, 1 - width_SA, by = width_SA)
  u.upper <- seq(width_SA, 1, by = width_SA)
}

# ---- Helper Functions ----
# Function to calculate probability matrix for all doses
p_all_doses <- function(p_initial, D, trans_matrix) {
  # Initialize probability matrix: rows for doses, columns for grades
  p_all <- matrix(0, nrow = D, ncol = length(p_initial))
  p_all[1, ] <- p_initial
  
  # Calculate probabilities for each dose level
  for(r in 2:D) {
    for(c in 1:length(p_initial)) {
      p_all[r,c] <- p_all[r-1, ] %*% trans_matrix[, c]
    }
  }
  return(p_all)
}

# Function to calculate true MTD
MTD_true <- function(p_A_tox_all, p_B_tox_all, p_C_tox_all, D, W, q_T) {
  W_max <- sum(apply(W, 1, max))  # Maximum possible toxicity weight
  EqTP <- numeric()
  k_EqTP_candidate <- numeric()
  
  # Calculate EqTP for each dose
  for(d in 1:D) {
    dose_d <- rbind(p_A_tox_all[d,],
                    p_B_tox_all[d,],
                    p_C_tox_all[d,])
    EqTP[d] <- sum(dose_d * W) / W_max
    EqTP[d] <- round(EqTP[d], 2)
    
    # Find interval index for current EqTP
    k_EqTP_candidate[d] <- which.max(as.numeric(table(cut(EqTP[d], int.bound))))
  }
  
  # Determine MTD based on toxicity threshold
  if(sum(k_EqTP_candidate > k_star) == D) {
    MTD_true <- 0  # No MTD exists
  } else {
    EqTP_candidate <- EqTP[k_EqTP_candidate <= k_star]
    MTD_true <- which.min(abs(EqTP_candidate - q_T))
  }
  
  return(list(EqTP = EqTP, MTD_true = MTD_true))
}

# Function to calculate equivalent efficacy probability (EqEP)
EqEP <- function(p_1_eff_all, p_2_eff_all, p_3_eff_all, D, E) {
  # Calculate maximum possible efficacy score
  E_max <- sum(apply(E, 1, max))
  EqEP <- numeric()
  
  # Calculate EqEP for each dose level
  for(d in 1:D) {
    dose_d <- rbind(p_1_eff_all[d,],
                    p_2_eff_all[d,],
                    p_3_eff_all[d,])
    EqEP[d] <- sum(dose_d * E) / E_max
    EqEP[d] <- round(EqEP[d], 2)
  }
  
  # Return original EqEP in ascending order
  return(list(EqEP_original = EqEP))
}

# Function to adjust EqEP based on curve types
adjust_efficacy_curves <- function(scenario_idx, original_EqEP) {
  # Define which scenarios belong to which curve type
  peak_scenarios <- c(3, 9, 15, 13, 23, 14, 17)  # Maximum efficacy at different dose levels
  late_activated_scenarios <- c(2, 5, 6, 18, 21, 24, 26)  # Steady then increased
  late_silent_scenarios <- c(1, 8, 12, 25, 27)  # Increase then plateau
  # All other scenarios remain as slope type (original monotonic curve)
  
  EqEP_shaped <- original_EqEP
  curve_type <- NULL
  
  if(scenario_idx %in% peak_scenarios) {
    if(scenario_idx %in% c(3, 9, 15)) {
      # Maximum efficacy at dose 4
      EqEP_shaped <- c(
        original_EqEP[1], 
        original_EqEP[2], 
        original_EqEP[3], 
        original_EqEP[6], 
        original_EqEP[5], 
        original_EqEP[4]
      )
      curve_type <- "peak"
    } else if(scenario_idx %in% c(13, 23)) {
      # Maximum efficacy at dose 3
      EqEP_shaped <- c(
        original_EqEP[1],
        original_EqEP[2],
        original_EqEP[6],
        original_EqEP[5],
        original_EqEP[4],
        original_EqEP[3]
      )
      curve_type <- "peak"
    } else if(scenario_idx %in% c(14, 17)) {
      # Maximum efficacy at dose 2
      EqEP_shaped <- c(
        original_EqEP[2],
        original_EqEP[3],
        original_EqEP[6],
        original_EqEP[5],
        original_EqEP[4],
        original_EqEP[4] - 0.15
      )
      curve_type <- "peak"
    }
  } else if(scenario_idx %in% late_activated_scenarios) {
    # Steady at early doses, then increased
    mean_early <- mean(original_EqEP[1:3])
    EqEP_shaped <- c(
      rep(mean_early, 3),
      original_EqEP[4:6]
    )
    curve_type <- "late-activated"
  } else if(scenario_idx %in% late_silent_scenarios) {
    # Increase then plateau
    plateau_value <- original_EqEP[4]
    EqEP_shaped <- c(
      original_EqEP[1],
      original_EqEP[3],
      original_EqEP[4],
      plateau_value + 0.01,
      plateau_value + 0.02,
      plateau_value + 0.03
    )
    curve_type <- "late-silent"
  } else {
    # Keep original monotonic curve for slope type
    curve_type <- "slope"
  }
  
  EqEP_shaped <- round(EqEP_shaped, 2)
  
  return(list(
    EqEP = EqEP_shaped,
    curve_type = curve_type
  ))
}

# Function to generate toxicity scenarios
Gen_Scn_T <- function(transition_level = c("low", "low", "low"),
                      p_A_tox, p_B_tox, p_C_tox, D, W) {
  # Generate toxicity A probabilities based on transition level
  if(transition_level[1] == "low") {
    p_A_tox_all <- p_all_doses(p_initial = p_A_tox, D = D, 
                               trans_matrix = trans_matrix_tox_L)
  } else if(transition_level[1] == "moderate") {
    p_A_tox_all <- p_all_doses(p_initial = p_A_tox, D = D, 
                               trans_matrix = trans_matrix_tox_M)
  } else if(transition_level[1] == "high") {
    p_A_tox_all <- p_all_doses(p_initial = p_A_tox, D = D, 
                               trans_matrix = trans_matrix_tox_H)
  }
  
  # Generate toxicity B probabilities
  if(transition_level[2] == "low") {
    p_B_tox_all <- p_all_doses(p_initial = p_B_tox, D = D, 
                               trans_matrix = trans_matrix_tox_L)
  } else if(transition_level[2] == "moderate") {
    p_B_tox_all <- p_all_doses(p_initial = p_B_tox, D = D, 
                               trans_matrix = trans_matrix_tox_M)
  } else if(transition_level[2] == "high") {
    p_B_tox_all <- p_all_doses(p_initial = p_B_tox, D = D, 
                               trans_matrix = trans_matrix_tox_H)
  }
  
  # Generate toxicity C probabilities
  if(transition_level[3] == "low") {
    p_C_tox_all <- p_all_doses(p_initial = p_C_tox, D = D, 
                               trans_matrix = trans_matrix_tox_L)
  } else if(transition_level[3] == "moderate") {
    p_C_tox_all <- p_all_doses(p_initial = p_C_tox, D = D, 
                               trans_matrix = trans_matrix_tox_M)
  } else if(transition_level[3] == "high") {
    p_C_tox_all <- p_all_doses(p_initial = p_C_tox, D = D, 
                               trans_matrix = trans_matrix_tox_H)
  }
  
  # Calculate MTD based on EqTP
  MTD_result <- MTD_true(p_A_tox_all, p_B_tox_all, p_C_tox_all, D, W, q_T)
  
  return(list(MTD_true = MTD_result$MTD_true,
              EqTP = MTD_result$EqTP,
              p_tox_all = rbind(p_A_tox_all, p_B_tox_all, p_C_tox_all)))
}

# Function to generate efficacy scenarios
Gen_Scn_E <- function(transition_level = c("low", "low", "low"),
                      p_1_eff, p_2_eff, p_3_eff, D, E) {
  # Generate efficacy 1 probabilities
  if(transition_level[1] == "low") {
    p_1_eff_all <- p_all_doses(p_initial = p_1_eff, D = D, 
                               trans_matrix = trans_matrix_eff_L)
  } else if(transition_level[1] == "moderate") {
    p_1_eff_all <- p_all_doses(p_initial = p_1_eff, D = D, 
                               trans_matrix = trans_matrix_eff_M)
  } else if(transition_level[1] == "high") {
    p_1_eff_all <- p_all_doses(p_initial = p_1_eff, D = D, 
                               trans_matrix = trans_matrix_eff_H)
  }
  
  # Generate efficacy 2 probabilities
  if(transition_level[2] == "low") {
    p_2_eff_all <- p_all_doses(p_initial = p_2_eff, D = D, 
                               trans_matrix = trans_matrix_eff_L)
  } else if(transition_level[2] == "moderate") {
    p_2_eff_all <- p_all_doses(p_initial = p_2_eff, D = D, 
                               trans_matrix = trans_matrix_eff_M)
  } else if(transition_level[2] == "high") {
    p_2_eff_all <- p_all_doses(p_initial = p_2_eff, D = D, 
                               trans_matrix = trans_matrix_eff_H)
  }
  
  # Generate efficacy 3 probabilities
  if(transition_level[3] == "low") {
    p_3_eff_all <- p_all_doses(p_initial = p_3_eff, D = D, 
                               trans_matrix = trans_matrix_eff_L)
  } else if(transition_level[3] == "moderate") {
    p_3_eff_all <- p_all_doses(p_initial = p_3_eff, D = D, 
                               trans_matrix = trans_matrix_eff_M)
  } else if(transition_level[3] == "high") {
    p_3_eff_all <- p_all_doses(p_initial = p_3_eff, D = D, 
                               trans_matrix = trans_matrix_eff_H)
  }
  
  # Calculate original EqEP in ascending order
  EqEP_result <- EqEP(p_1_eff_all, p_2_eff_all, p_3_eff_all, D, E)
  
  return(list(EqEP = EqEP_result$EqEP,
              p_eff_all = rbind(p_1_eff_all, p_2_eff_all, p_3_eff_all)))
}

# Function to record true MTD values for all toxicity scenarios
record_MTD_true_all <- function(p_A_tox, p_B_tox, p_C_tox, W) {
  
  # Initialize storage matrices
  EqTP_all <- matrix(numeric(), nrow = 27, ncol = 6)
  rownames(EqTP_all) <- c("LLL","LLM","LLH","LML","LMM","LMH",
                          "LHL","LHM","LHH","MLL","MLM","MLH",
                          "MML","MMM","MMH","MHL","MHM","MHH",
                          "HLL","HLM","HLH","HML","HMM","HMH",
                          "HHL","HHM","HHH")
  colnames(EqTP_all) <- c("dose1","dose2","dose3","dose4","dose5","dose6")
  MTD_true_all <- numeric(27)
  
  # Calculate MTD and EqTP for each scenario
  for (i in 1:length(scenario_list)) {
    transition_level <- scenario_list[[i]]
    Scn_T <- Gen_Scn_T(transition_level = transition_level, 
                       p_A_tox, p_B_tox, p_C_tox, D = 6, W = W)
    
    EqTP_all[i,] <- Scn_T$EqTP
    MTD_true_all[i] <- Scn_T$MTD_true
  }
  
  return(list(
    EqTP_all = EqTP_all,
    MTD_true_all = MTD_true_all
  ))
}


# ---- Main function for generating trial scenarios ----
generate_trial_scenarios <- function(
  w_e, w_t = 1 - w_e,
  D, startdose, samplesize, csize,
  alpha_u, beta_u,
  q_T,
  W, E,
  p_A_tox, p_B_tox, p_C_tox,
  p_1_eff, p_2_eff, p_3_eff
) {
  # Parameter validation
  if(w_e + w_t != 1) stop("Utility weights must sum to 1")
  if(D < 2) stop("Number of dose levels must be greater than 1")
  
  # Generate scenarios using original functions
  MTD_result <- record_MTD_true_all(p_A_tox, p_B_tox, p_C_tox, W)
  
  # Initialize storage for scenario results
  scenario_results <- list()
  
  # Generate both toxicity and efficacy scenarios
  for(i in 1:27) {
    for(j in 1:27) {
      current_index <- (i-1)*27 + j
      
      # Generate toxicity scenario
      tox_scenario <- Gen_Scn_T(scenario_list[[i]], p_A_tox, p_B_tox, p_C_tox, D, W)
      
      # Generate efficacy scenario
      eff_scenario <- Gen_Scn_E(scenario_list[[j]], p_1_eff, p_2_eff, p_3_eff, D, E)
      
      # Adjust efficacy curves based on curve types
      adjusted_eff <- adjust_efficacy_curves(j, eff_scenario$EqEP)
      eff_scenario$EqEP <- adjusted_eff$EqEP
      eff_scenario$curve_type <- adjusted_eff$curve_type
      
      # Calculate utility for the scenario
      EU_true_candidates <- w_e * eff_scenario$EqEP + w_t * (1 - tox_scenario$EqTP)
      EU_true_candidates <- round(EU_true_candidates, 2)
      
      # Store scenario results
      scenario_results[[current_index]] <- list(
        toxicity = tox_scenario,
        efficacy = eff_scenario,
        EU = EU_true_candidates
      )
    }
  }
  
  # Return results
  return(list(
    scenarios = scenario_results,
    MTD_all = MTD_result$MTD_true_all,
    EqTP_all = MTD_result$EqTP_all,
    parameters = list(
      weights = c(w_e = w_e, w_t = w_t),
      trial_params = list(D = D, 
                          startdose = startdose,
                          samplesize = samplesize, 
                          csize = csize),
      target_tox = q_T
    )
  ))
}

# ---- Function to generate combined results tables ----
generate_summary_tables <- function(trial_results) {
  scenarios <- trial_results$scenarios
  n_scenarios <- length(scenarios)
  
  scenarios_prob_tox <- data.frame(matrix(nrow = 27, ncol = 6))
  colnames(scenarios_prob_tox) <- paste0("pT_", 1:6)
  
  for(i in 1:27) {
    current_tox_scenario <- scenarios[[(i-1)*27 + 1]]$toxicity 
    scenarios_prob_tox[i,] <- current_tox_scenario$EqTP
    
    transition_level_tox <- scenario_list[[i]]
    scenario_tox_abb <- toupper(paste0(substr(transition_level_tox, 1, 1), collapse = ""))
    rownames(scenarios_prob_tox)[i] <- scenario_tox_abb
  }
  
  scenarios_prob_eff <- data.frame(matrix(nrow = 27, ncol = 6))
  colnames(scenarios_prob_eff) <- paste0("pE_", 1:6)
  
  for(j in 1:27) {
    current_eff_scenario <- scenarios[[j]]$efficacy
    scenarios_prob_eff[j,] <- current_eff_scenario$EqEP
    
    transition_level_eff <- scenario_list[[j]]
    scenario_eff_abb <- toupper(paste0(substr(transition_level_eff, 1, 1), collapse = ""))
    rownames(scenarios_prob_eff)[j] <- scenario_eff_abb
  }
  
  scenarios_prob <- data.frame(matrix(nrow = n_scenarios, ncol = 12))
  colnames(scenarios_prob) <- c(paste0("pT_", 1:6), paste0("pE_", 1:6))
  
  for(i in 1:27) {
    for(j in 1:27) {
      current_idx <- (i-1)*27 + j
      current_scenario <- scenarios[[current_idx]]
      
      scenarios_prob[current_idx, 1:6] <- current_scenario$toxicity$EqTP
      scenarios_prob[current_idx, 7:12] <- current_scenario$efficacy$EqEP
    }
  }
  
  w_e <- trial_results$parameters$weights["w_e"]
  w_t <- trial_results$parameters$weights["w_t"]
  
  true_utilities <- data.frame(matrix(nrow = n_scenarios, ncol = 6))
  colnames(true_utilities) <- paste0("dose", 1:6)

  for(i in 1:n_scenarios) {
    current_scenario <- scenarios[[i]]
    true_utilities[i,] <- w_e * current_scenario$efficacy$EqEP + 
      w_t * (1 - current_scenario$toxicity$EqTP)
  }
  
  return(list(
    combined_prob = scenarios_prob,
    toxicity_prob = scenarios_prob_tox,
    efficacy_prob = scenarios_prob_eff,
    utilities = true_utilities
  ))
}

extract_fixed_scenarios <- function(summary_tables, tox_abb, eff_abb) {
  # Validate input lengths
  if(length(tox_abb) != length(eff_abb)) {
    stop("Toxicity and efficacy abbreviation vectors must have the same length")
  }
  
  n_scenarios <- length(tox_abb)
  
  # Extract and format toxicity probabilities
  tox_prob <- matrix(NA, nrow = n_scenarios, ncol = 6)
  for(i in 1:n_scenarios) {
    tox_idx <- which(rownames(summary_tables$toxicity_prob) == tox_abb[i])
    tox_prob[i,] <- as.numeric(summary_tables$toxicity_prob[tox_idx,])
  }
  tox_prob <- as.data.frame(tox_prob)
  colnames(tox_prob) <- paste0("pT_", 1:6)
  rownames(tox_prob) <- paste("Scenario", 1:n_scenarios)
  
  # Extract and format efficacy probabilities
  eff_prob <- matrix(NA, nrow = n_scenarios, ncol = 6)
  for(i in 1:n_scenarios) {
    eff_idx <- which(rownames(summary_tables$efficacy_prob) == eff_abb[i])
    eff_prob[i,] <- as.numeric(summary_tables$efficacy_prob[eff_idx,])
  }
  eff_prob <- as.data.frame(eff_prob)
  colnames(eff_prob) <- paste0("pE_", 1:6)
  rownames(eff_prob) <- paste("Scenario", 1:n_scenarios)
  
  # Calculate utilities
  utilities <- matrix(NA, nrow = n_scenarios, ncol = 6)
  for(i in 1:n_scenarios) {
    tox_idx <- which(rownames(summary_tables$toxicity_prob) == tox_abb[i])
    eff_idx <- which(rownames(summary_tables$efficacy_prob) == eff_abb[i])
    combined_idx <- (tox_idx - 1) * 27 + eff_idx
    utilities[i,] <- as.numeric(summary_tables$utilities[combined_idx,])
  }
  utilities <- as.data.frame(utilities)
  colnames(utilities) <- paste0("dose", 1:6)
  rownames(utilities) <- paste("Scenario", 1:n_scenarios)
  
  return(list(
    toxicity = tox_prob,
    efficacy = eff_prob,
    utilities = utilities,
    scenario_mapping = data.frame(
      scenario = paste("Scenario", 1:n_scenarios),
      tox_type = tox_abb,
      eff_type = eff_abb
    )
  ))
}


# ---- Usage example ----
# Generate trial scenarios
trial_results <- generate_trial_scenarios(    
  w_e = w_e, # 0.45/0.7
  w_t = 1 - w_e,
  
  D, startdose, samplesize, csize,
  alpha_u, beta_u,
  q_T,
  
  W, E,
  p_A_tox, p_B_tox, p_C_tox,
  p_1_eff, p_2_eff, p_3_eff
)

# Generate summary tables
summary_tables <- generate_summary_tables(trial_results)

# random fix scenario combination labels
tox_abb <- c("MHM", "HLH", "HML", "HHM", "LLM", "LLM", "HHH", "HHH", "MMH")
eff_abb <- c("HLM", "LLH", "MLH", "LLH", "LLL", "HMH", "MHM", "LML", "LMM")

fixed_results <- extract_fixed_scenarios(summary_tables, tox_abb, eff_abb)

# Table S.6
print("Toxicity probabilities:")
print(fixed_results$toxicity)

# Table S.7
print("Efficacy probabilities:")
print(fixed_results$efficacy)

# Table S.8 (HE) or Table S.9 (LT)
print("Fixed scenario utilities:")
print(fixed_results$utilities)

# print("Scenario mapping:")
# print(fixed_results$scenario_mapping)


# ----- Table S.5 raw data
tox_A_S5 <- p_all_doses(p_initial = p_A_tox, D = 6, trans_matrix = trans_matrix_tox_H)
tox_B_S5 <- p_all_doses(p_initial = p_B_tox, D = 6, trans_matrix = trans_matrix_tox_H)
tox_C_S5 <- p_all_doses(p_initial = p_C_tox, D = 6, trans_matrix = trans_matrix_tox_H)
eff_1_S5 <- p_all_doses(p_initial = p_1_eff, D = 6, trans_matrix = trans_matrix_eff_H)
eff_2_S5 <- p_all_doses(p_initial = p_2_eff, D = 6, trans_matrix = trans_matrix_eff_H)
eff_3_S5 <- p_all_doses(p_initial = p_3_eff, D = 6, trans_matrix = trans_matrix_eff_H)

create_dose_group <- function(dose_idx, tox_A, tox_B, tox_C, eff_1, eff_2, eff_3) {
  combo_A1 <- c(tox_A[dose_idx,], eff_1[dose_idx,])
  combo_B2 <- c(tox_B[dose_idx,], eff_2[dose_idx,])
  combo_C3 <- c(tox_C[dose_idx,], eff_3[dose_idx,])
  
  dose_matrix <- rbind(combo_A1, combo_B2, combo_C3)
  rownames(dose_matrix) <- c("Type A/1", "Type B/2", "Type C/3")
  colnames(dose_matrix) <- c(paste0("Tox_grade_", 0:4), paste0("Eff_grade_", 1:3))
  
  return(dose_matrix)
}

Table_S5_results <- list()
for(i in 1:6) {
  Table_S5_results[[i]] <- round(create_dose_group(i, 
                                         tox_A_S5, tox_B_S5, tox_C_S5,
                                         eff_1_S5, eff_2_S5, eff_3_S5),3)
}
names(Table_S5_results) <- paste("Dose", 1:6)

# Table S.5 results
Table_S5_results










