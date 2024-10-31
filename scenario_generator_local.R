# Scenario Generation

rm(list=ls())

## path for results saving and packages in use ##
{
  library(rlist)
  # path <- "your path/"
}

set.seed(315)

library(cmdstanr)


## parameters ##
{
  # high efficacy setting (HE)
  # w_e <- 0.7;
  # low toxicity setting (LT)
  w_e <- 0.45;
  # score trade-off
  w_t <- 1 - w_e
  
  ndose <- 6 # total dose level number
  D <- ndose;
  x <- log(1 : ndose/10) # pseudo ascending dose sequence 
  
  startdose <- 1
  samplesize <- 27 
  csize <- 3 # cohort size
  
  # prior of EU ~ Beta(alpha_u, beta_u)
  alpha_u <- 1
  beta_u <- 1
  
  prefix_scenario <- data.frame()
  MTD_true_current <- numeric()
  
  early_stop_rate <- data.frame() # early stop rate for each scenario 
  poor_allocation_rate <- data.frame() 
  average_patient_number <- data.frame()
  average_overdose_number <- data.frame()
  average_selection_percentage <- data.frame()

  rownames_ESR <- vector()
  rownames_APN <- vector()
  rownames_ASP <- vector()
  rownames_csel_MTD_all <- vector()
  rownames_csel_OBD_all <- vector()
  rownames_PA <- vector() # poor allocation
  rownames_OD <- vector() # over dose
  
  csel_OBD_all <- data.frame() # correct selected OBD rate for each scenario
  

  ## OBD_true matrix: it is used to store the true_OBD values for each scenario combination
  
  # ALL scenarios of 27 by 27
  OBD_true <- matrix(NA, nrow = 27, ncol = 27) # length(scenario_list) == 27
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
  Tox_namelist <- rownames(OBD_true) # applied for "multi_test.R" when parallel running
  
  # store the OBD before the simulation getting started
  
  # True toxicity probability - EqTP (marginal toxicity)
  true_tox_prob_list <- vector("list", length = 27)
  
  # True efficacy probability - EqEP (marginal efficacy)
  true_eff_prob_list <- vector("list", length = 27 * 27)
  
  # True expected utility
  EU_true_list <- vector("list", length = 27 * 27)
  
  # EqEP curve types (totally 4 types)
  curve_type_eff <- vector("list", length = 27 * 27)

  # weights for EqTP and EqEP in Utility (pre-set)
  # toxicity-efficacy trade-offs: when w_e > w_t: prefer high efficacy; when w_e < w_t: prefer low toxicity
  
  
  # Over Toxicity Threshold: Target qTP
  q_T <- 0.3
  
  # severity matrix of toxicity for each type at each grade (grade 0 to 4)
  W <- rbind(c(0.1, 0.8, 1.3, 1.8, 2.3), # first dose used to be 0 for all three toxicities
             c(0.07, 0.5, 1, 1.5, 2),
             c(0.02, 0.05, 0.25, 0.5, 1))
  
  # efficacious weight matrix for each type at each grade (low / middle / high)
  E <- rbind(c(0.05, 0.2, 0.47), 
             c(0.07, 0.14, 0.41), 
             c(0.03, 0.15, 0.45))
  
  # initial probabilities at each grade of each type of toxicity (dose 1 prob)
  p_A_tox <- c(0.891, 0.072, 0.032, 0.004, 0.001) 
  p_B_tox <- c(0.868, 0.062, 0.031, 0.029, 0.010)
  p_C_tox <- c(0.817, 0.071, 0.057, 0.050, 0.005)
  
  # initial probabilities at each grade of each type of efficacy (dose 1 prob)
  p_1_eff <- c(0.900, 0.070, 0.030) 
  p_2_eff <- c(0.850, 0.130, 0.020)
  p_3_eff <- c(0.870, 0.090, 0.040)
}

## interval ##
# uTPI settings (including the Sensitivity Analysis)
{
  width_SA <- 0.10 # 0.05 / 0.1 (default) / 0.2
  
  int.bound <- seq(0, 1, 0.1) # equal-width toxicity intervals
  k_star <- which.max(as.numeric(table(cut(q_T, int.bound)))) # location of the highest tolerable toxicity interval
  
  # equal-width efficacy intervals, storing the lower bounds and upper bounds separately
  u.lower <- seq(0.0, 1 - width_SA, by = width_SA)
  u.upper <- seq(width_SA, 1, by = width_SA)
}

## scenario labels ##
{
  # transition rates: scenario combination for toxicity/efficacy 
  scenario_list <- list(c("low","low","low"),
                        c("low","low","moderate"),
                        c("low","low","high"),
                        c("low","moderate","low"),
                        c("low","moderate","moderate"),
                        c("low","moderate","high"),
                        c("low","high","low"),
                        c("low","high","moderate"),
                        c("low","high","high"),
                        c("moderate","low","low"),
                        c("moderate","low","moderate"),
                        c("moderate","low","high"),
                        c("moderate","moderate","low"),
                        c("moderate","moderate","moderate"),
                        c("moderate","moderate","high"),
                        c("moderate","high","low"),
                        c("moderate","high","moderate"),
                        c("moderate","high","high"),
                        c("high","low","low"),
                        c("high","low","moderate"),
                        c("high","low","high"),
                        c("high","moderate","low"),
                        c("high","moderate","moderate"),
                        c("high","moderate","high"),
                        c("high","high","low"),
                        c("high","high","moderate"),
                        c("high","high","high"))
  
  # transition matrix of toxicity in various speed
  trans_matrix_tox_L <- as.matrix(rbind(c(0.9,0.1,0,0,0),
                                        c(0,0.9,0.1,0,0),
                                        c(0,0,0.9,0.1,0),
                                        c(0,0,0,0.9,0.1),
                                        c(0,0,0,0,1)))
  trans_matrix_tox_M <- as.matrix(rbind(c(0.75,0.25,0,0,0), 
                                        c(0,0.75,0.25,0,0),
                                        c(0,0,0.75,0.25,0),
                                        c(0,0,0,0.75,0.25),
                                        c(0,0,0,0,1)))
  trans_matrix_tox_H <- as.matrix(rbind(c(0.45,0.55,0,0,0), 
                                        c(0,0.45,0.55,0,0),
                                        c(0,0,0.45,0.55,0),
                                        c(0,0,0,0.45,0.55),
                                        c(0,0,0,0,1)))
  
  # transition matrix of efficacy in various speed: just for the scenario generation
  trans_matrix_eff_L <- as.matrix(rbind(c(0.9, 0.1, 0),
                                        c(0, 0.9, 0.1),
                                        c(0, 0, 1)))
  trans_matrix_eff_M <- as.matrix(rbind(c(0.75, 0.25, 0),
                                        c(0, 0.75, 0.25),
                                        c(0, 0, 1)))
  trans_matrix_eff_H <- as.matrix(rbind(c(0.45, 0.55, 0), 
                                        c(0, 0.45, 0.55),
                                        c(0, 0, 1)))
}

## FUNCTIONS ##
{
  # probability matrix of observing various grades of toxicity/efficacy at ALL D dose
  p_all_doses <- function(p_initial, D, trans_matrix){
    p_all <- matrix(0, nrow = D, ncol = length(p_initial))
    p_all[1, ] <- p_initial
    
    for(r in 2:D){
      for(c in 1:length(p_initial)){
        p_all[r,c] <- p_all[r-1, ] %*% trans_matrix[, c]
      }
    }
    return(p_all) # A #D(6) by #grade(tox 5 / eff 3) matrix, 
                  # each row for a dose level (1 to D = 6), each column for a grade; 
                  # elements are the probability of observing corresponding dose-grade combination
  }
  
  # generate true MTD
  MTD_true <- function(p_A_tox_all, p_B_tox_all, p_C_tox_all, D, W, q_T){
    W_max <- sum(apply(W, 1, max)) # sum of the maximum of each row of W, W_max is fixed when W is fixed
    EqTP <- numeric()
    k_EqTP_candidate <- numeric() # generate k candidate for each EqTP for various doses
    
    for(d in 1:D){
      dose_d <- rbind(p_A_tox_all[d,],
                      p_B_tox_all[d,],
                      p_C_tox_all[d,])
      EqTP[d] <- sum(dose_d * W) / W_max;
      
      EqTP[d] <- round(EqTP[d], 2);
      
      k_EqTP_candidate[d] <- which.max(as.numeric(table(cut(EqTP[d], int.bound)))) # the k[d] is the interval to which EqTP[d] belongs
    }
    
    if(sum(k_EqTP_candidate > k_star) == D){ # toxicity measures of all doses are greater than the threshold index
      MTD_true <- 0 # no MTD for this scenario
    }
    else{
      EqTP_candidate <- EqTP[k_EqTP_candidate <= k_star] 
      MTD_true <- which.min(abs(EqTP_candidate - q_T)) # choose the true MTD 
    }
    return(list(EqTP = EqTP, 
                MTD_true = MTD_true))
  }
  
  # calculate original EqEP in ascending order
  EqEP <- function(p_1_eff_all, p_2_eff_all, p_3_eff_all, D, E){
    E_max <- sum(apply(E, 1, max)) # sum of the maximum of each row of E, E_max is fixed when E is fixed
    
    EqEP <- numeric()
    #type <- NULL # curve types
    
    for(d in 1:D){
      dose_d <- rbind(p_1_eff_all[d,],
                      p_2_eff_all[d,],
                      p_3_eff_all[d,])
      EqEP[d] <- sum(dose_d * E) / E_max;
      EqEP[d] <- round(EqEP[d], 2);
    }
    
    EqEP_original <- EqEP # the original EqEP in asecending order
    
    # peak
    
    # late-activaed 
    
    # late-slient
    
    # cliff
     
    return(list(EqEP_original = EqEP_original))#, # before specify the curve types
                #EqEP = EqEP, # after specify the curve types
                #type = type))
  }
  
  # generate scenario for toxicity
  Gen_Scn_T <- function(transition_level = c("low", "low", "low"),
                        p_A_tox, p_B_tox, p_C_tox, D, W){
    
    # for toxicity A
    if(transition_level[1] == "low"){
      p_A_tox_all <- p_all_doses(p_initial = p_A_tox, D = D, 
                                 trans_matrix = trans_matrix_tox_L)
    }else if(transition_level[1] == "moderate"){
      p_A_tox_all <- p_all_doses(p_initial = p_A_tox, D = D, 
                                 trans_matrix = trans_matrix_tox_M)
    }else if(transition_level[1] == "high"){
      p_A_tox_all <- p_all_doses(p_initial = p_A_tox, D = D, 
                                 trans_matrix = trans_matrix_tox_H)
    }
    
    # for toxicity B
    if(transition_level[2] == "low"){
      p_B_tox_all <- p_all_doses(p_initial = p_B_tox, D = D, 
                                 trans_matrix = trans_matrix_tox_L)
    }else if(transition_level[2] == "moderate"){
      p_B_tox_all <- p_all_doses(p_initial = p_B_tox, D = D, 
                                 trans_matrix = trans_matrix_tox_M)
    }else if(transition_level[2] == "high"){
      p_B_tox_all <- p_all_doses(p_initial = p_B_tox, D = D, 
                                 trans_matrix = trans_matrix_tox_H)
    }
    
    # for toxicity C
    if(transition_level[3] == "low"){
      p_C_tox_all <- p_all_doses(p_initial = p_C_tox, D = D, 
                                 trans_matrix = trans_matrix_tox_L)
    }else if(transition_level[3] == "moderate"){
      p_C_tox_all <- p_all_doses(p_initial = p_C_tox, D = D, 
                                 trans_matrix = trans_matrix_tox_M)
    }else if(transition_level[3] == "high"){
      p_C_tox_all <- p_all_doses(p_initial = p_C_tox, D = D, 
                                 trans_matrix = trans_matrix_tox_H)
    }
    
    # calculate MTD_true based on EqTP 
    MTD_true <- MTD_true(p_A_tox_all, p_B_tox_all, p_C_tox_all, D, W, q_T)
    
    return(list(MTD_true = MTD_true$MTD_true,
                EqTP = MTD_true$EqTP,
                p_tox_all = rbind(p_A_tox_all, p_B_tox_all, p_C_tox_all)))
    
  }
  
  # generate scenario for efficacy
  Gen_Scn_E <- function(transition_level = c("low", "low", "low"),
                        p_1_eff, p_2_eff, p_3_eff, D, E){
    
    # for efficacy 1
    if(transition_level[1] == "low"){
      p_1_eff_all <- p_all_doses(p_initial = p_1_eff, D = D, 
                                 trans_matrix = trans_matrix_eff_L)
    }else if(transition_level[1] == "moderate"){
      p_1_eff_all <- p_all_doses(p_initial = p_1_eff, D = D, 
                                 trans_matrix = trans_matrix_eff_M)
    }else if(transition_level[1] == "high"){
      p_1_eff_all <- p_all_doses(p_initial = p_1_eff, D = D, 
                                 trans_matrix = trans_matrix_eff_H)
    }
    
    # for efficacy 2
    if(transition_level[2] == "low"){
      p_2_eff_all <- p_all_doses(p_initial = p_2_eff, D = D, 
                                 trans_matrix = trans_matrix_eff_L)
    }else if(transition_level[2] == "moderate"){
      p_2_eff_all <- p_all_doses(p_initial = p_2_eff, D = D, 
                                 trans_matrix = trans_matrix_eff_M)
    }else if(transition_level[2] == "high"){
      p_2_eff_all <- p_all_doses(p_initial = p_2_eff, D = D, 
                                 trans_matrix = trans_matrix_eff_H)
    }
    
    # for efficacy 3
    if(transition_level[3] == "low"){
      p_3_eff_all <- p_all_doses(p_initial = p_3_eff, D = D, 
                                 trans_matrix = trans_matrix_eff_L)
    }else if(transition_level[3] == "moderate"){
      p_3_eff_all <- p_all_doses(p_initial = p_3_eff, D = D, 
                                 trans_matrix = trans_matrix_eff_M)
    }else if(transition_level[3] == "high"){
      p_3_eff_all <- p_all_doses(p_initial = p_3_eff, D = D, 
                                 trans_matrix = trans_matrix_eff_H)
    }
    
    # calculate original EqEP in ascending order
    EqEP <- EqEP(p_1_eff_all, p_2_eff_all, p_3_eff_all, D, E)
    
    
    
    return(list(EqEP = EqEP$EqEP,
                p_eff_all = rbind(p_1_eff_all, p_2_eff_all, p_3_eff_all)))
  }
  
  # record EqTP, true MTD for all toxicity scenarios
  record_MTD_true_all <- function(p_A_tox, p_B_tox, p_C_tox, W){
    scenario_list <- list(c("low","low","low"),
                          c("low","low","moderate"),
                          c("low","low","high"),
                          c("low","moderate","low"),
                          c("low","moderate","moderate"),
                          c("low","moderate","high"),
                          c("low","high","low"),
                          c("low","high","moderate"),
                          c("low","high","high"),
                          c("moderate","low","low"),
                          c("moderate","low","moderate"),
                          c("moderate","low","high"),
                          c("moderate","moderate","low"),
                          c("moderate","moderate","moderate"),
                          c("moderate","moderate","high"),
                          c("moderate","high","low"),
                          c("moderate","high","moderate"),
                          c("moderate","high","high"),
                          c("high","low","low"),
                          c("high","low","moderate"),
                          c("high","low","high"),
                          c("high","moderate","low"),
                          c("high","moderate","moderate"),
                          c("high","moderate","high"),
                          c("high","high","low"),
                          c("high","high","moderate"),
                          c("high","high","high"))
    
    EqTP_all <- matrix(numeric(), nrow = 27, ncol = 6)
    
    rownames(EqTP_all) <- c("LLL","LLM","LLH","LML","LMM","LMH",
                            "LHL","LHM","LHH","MLL","MLM","MLH",
                            "MML","MMM","MMH","MHL","MHM","MHH",
                            "HLL","HLM","HLH","HML","HMM","HMH",
                            "HHL","HHM","HHH")
    colnames(EqTP_all) <- c("dose1","dose2","dose3","dose4","dose5","dose6")
    MTD_true_all <- numeric(27)
    
    for (i in 1: length(scenario_list)){
      transition_level <- scenario_list[[i]]
      
      Scn_T <- Gen_Scn_T(transition_level = transition_level, p_A_tox, p_B_tox, p_C_tox, D = 6, W = W)
      
      EqTP_all[i,] <- Scn_T$EqTP
      MTD_true_all[i] <- Scn_T$MTD_true
    }
    
    return(list(
      EqTP_all = EqTP_all,
      MTD_true_all = MTD_true_all
    ))
    
  }
}

# ----------------- Data Generation ------------------ #

# save p_tox_all and p_eff_all for nine fix scenarios
p_tox_all_fix <- list()
p_eff_all_fix <- list()

map_abb_to_word <- list(M = "moderate", H = "high", L = "low")

tox_abb <- c("MHM", "HLH", "HML", "HHM", "LLM", "LLM", "HHH", "HHH", "MMH")
eff_abb <- c("HLM", "LLH", "MLH", "LLH", "LLL", "HMH", "MHM", "LML", "LMM")

tox_full <- sapply(tox_abb, function(x) {
  words <- strsplit(x, "")[[1]]
  scenario <- sapply(words, function(letter) map_abb_to_word[[letter]])
  paste(scenario, collapse = " ")
})

eff_full <- sapply(eff_abb, function(x) {
  words <- strsplit(x, "")[[1]]
  scenario <- sapply(words, function(letter) map_abb_to_word[[letter]])
  paste(scenario, collapse = " ")
})

scenario_list_str <- sapply(scenario_list, function(x) paste(x, collapse = " "))

indices_tox <- sapply(tox_full, function(x) which(scenario_list_str == x))
indices_eff <- sapply(eff_full, function(x) which(scenario_list_str == x))


# True MTD

record_MTD_true_all <- record_MTD_true_all(p_A_tox, p_B_tox, p_C_tox, W)

# True OBD & True tox/eff probabilities
for (p in 1:27) {
  transition_level_tox <- scenario_list[[p]]
  
  # True MTD
  MTD_true_current[p] <- record_MTD_true_all$MTD_true_all[p]
  
  # TOX scenario generation, with true MTD for current scenario
  Scn_T <- Gen_Scn_T(transition_level = transition_level_tox, 
                     p_A_tox, p_B_tox, p_C_tox, 
                     D = 6, W = W)
  
  true_tox_prob_list[[p]] <- Scn_T$EqTP; # same for every 27 as a set
  
  scenario_tox_abb <- noquote(toupper(paste0(substr(transition_level_tox, 1, 1), collapse = "")));
  
  
  for (q in 1:27){ #length(scenario_list)
    transition_level_eff <- scenario_list[[q]]
    
    # EFF scenario generation
    Scn_E <- Gen_Scn_E(transition_level = transition_level_eff, 
                       p_1_eff, p_2_eff, p_3_eff, 
                       D = 6, E = E)
    
    
    # store the true EqTP / EqEP as the true toxicity / efficacy probability inputs
    current_combination <- (p - 1) * 27 + q; #length(scenario_list)
    
    # use the function to change the order in various shapes
    EqEP_shaped <- NULL
    type <- NULL
    
    # 1. peak
    # maximum efficacy at dose 4
    if(q %in% c(3,9,15)) {
      EqEP_origin <- Scn_E$EqEP; 
      EqEP_shaped <- c(EqEP_origin[1], EqEP_origin[2], EqEP_origin[3], EqEP_origin[6], EqEP_origin[5], EqEP_origin[4]); 
      type <- "peak"}
    # maximum efficacy at dose 3
    else if(q %in% c(13,23)) {
      EqEP_origin <- Scn_E$EqEP; 
      EqEP_shaped <- c(EqEP_origin[1], EqEP_origin[2], EqEP_origin[6], EqEP_origin[5], EqEP_origin[4], EqEP_origin[3]); 
      type <- "peak"}
    # maximum efficacy at dose 2
    else if(q %in% c(14,17)) {
      EqEP_origin <- Scn_E$EqEP; 
      EqEP_shaped <- c(EqEP_origin[2], EqEP_origin[3], EqEP_origin[6], EqEP_origin[5], EqEP_origin[4], EqEP_origin[4] - 0.15); 
      type <- "peak"}
    
    # 2. late-activated
    else if(q %in% c(2,5,6,18,21,24,26)) {
      EqEP_origin <- Scn_E$EqEP; 
      EqEP_shaped <- c(mean(EqEP_origin[1:3]), mean(EqEP_origin[1:3]), mean(EqEP_origin[1:3]), EqEP_origin[4], EqEP_origin[5], EqEP_origin[6]);
      type <- "late-activated"}
    
    # 3. late-silent
    else if(q %in% c(1,8,12,25,27)) {
      EqEP_origin <- Scn_E$EqEP; 
      EqEP_shaped <- c(EqEP_origin[1], EqEP_origin[3], EqEP_origin[4], EqEP_origin[4] + 0.01, EqEP_origin[4] + 0.02, EqEP_origin[4] + 0.03);
      type <- "late-silent"}
    
    # 4. slope
    else {
      EqEP_origin <- Scn_E$EqEP; 
      EqEP_shaped <- EqEP_origin;
      type <- "slope"}
    
    EqEP_shaped <- round(EqEP_shaped, 2)
    
    true_eff_prob_list[[current_combination]] <- EqEP_shaped;
    curve_type_eff[[current_combination]] <- type;
    
    EU_true_candidates <- w_t * (1 - Scn_T$EqTP) + w_e * EqEP_shaped;
    EU_true_candidates <- round(EU_true_candidates, 2);
    
    # True OBD for current scenario
    OBD_true[p, q] <- which.max(head(EU_true_candidates, MTD_true_current[p])) # auto select the lower dose when multiple maximum values exist
    
    # True Utility
    EU_true_list[[current_combination]] <- EU_true_candidates 
    
    # scenario combination list
    scenario_eff_abb <- noquote(toupper(paste0(substr(transition_level_eff, 1, 1), collapse = "")));
    prefix_scenario[current_combination, 1] <- paste0(scenario_tox_abb, "_", scenario_eff_abb)
  }
}



# Save the (pT.true, pE.true)
pT.true <- true_tox_prob_list
pE.true <- true_eff_prob_list
curve_type <- curve_type_eff

MTD_OBD_true <- data.frame("MTD" = MTD_true_current, OBD_true)

# proportion of various types of efficacy curves
sum(curve_type_eff[1:27] == "late-silent") / 27 # 18.5%
sum(curve_type_eff[1:27] == "late-activated") / 27 # 25.9%
sum(curve_type_eff[1:27] == "peak") / 27 # 25.9%
sum(curve_type_eff[1:27] == "slope") / 27 # 29.6%

# ----------------- Table Generation ------------------ #

{
  # print the scenario probabilities (both tox & eff) in table format
  scenarios_prob <- data.frame(matrix(nrow = 729, ncol = 6*2))
  rownames(scenarios_prob) <- unlist(prefix_scenario) #1:729
  original_vector <- rep(c(paste0("pT_", 1:6), paste0("pE_", 1:6)))
  new_vector <- c()
  for (i in 1:6) {
    new_vector <- c(new_vector, original_vector[i], original_vector[i + 6])
  }
  colnames(scenarios_prob) <- new_vector
  for (i in 1:27) {
    for (p in 1:27) {
      scenarios_prob[((i - 1) * 27 + p), seq(1, 11, by = 2)] <- pT.true[[i]]
    }
  }
  for (i in 1:729) {
    scenarios_prob[i, seq(2, 12, by = 2)] <- pE.true[[i]]
  }
  #head(scenarios_prob, 27)
  
  # print the true tox prob only in table format
  scenarios_prob_tox <- data.frame(matrix(nrow = 27, ncol = 6))
  colnames(scenarios_prob_tox) <- paste0("pT_", 1:6)
  rownames_scenarios_prob_tox <- vector()
  for (i in 1:27) {
    scenarios_prob_tox[i, seq(1, 6, by = 1)] <- pT.true[[i]]
    transition_level_tox <- scenario_list[[i]]
    scenario_tox_abb <- noquote(toupper(paste0(substr(transition_level_tox, 1, 1), collapse = "")))
    rownames_scenarios_prob_tox[i] <- scenario_tox_abb
  }
  rownames(scenarios_prob_tox) <- rownames_scenarios_prob_tox
  #head(scenarios_prob_tox, 27)
  
  # print the true eff prob only in table format
  scenarios_prob_eff <- data.frame(matrix(nrow = 27, ncol = 6))
  colnames(scenarios_prob_eff) <- paste0("pE_", 1:6)
  rownames_scenarios_prob_eff <- vector()
  for (i in 1:27) {
    scenarios_prob_eff[i, seq(1, 6, by = 1)] <- pE.true[[i]]
    transition_level_eff <- scenario_list[[i]]
    scenario_eff_abb <- noquote(toupper(paste0(substr(transition_level_eff, 1, 1), collapse = "")))
    rownames_scenarios_prob_eff[i] <- scenario_eff_abb
  }
  rownames(scenarios_prob_eff) <- rownames_scenarios_prob_eff
  #head(scenarios_prob_eff, 27)
  
  # print the true utilities in table format
  true_utilities <- data.frame(matrix(nrow = 729, ncol = 6))
  rownames(true_utilities) <- unlist(prefix_scenario) #1:729
  colnames(true_utilities) <- c("dose1","dose2","dose3","dose4","dose5","dose6")
  for (i in 1:729) {
    true_utilities[i, seq(1, 6, by = 1)] <- EU_true_list[[i]]
  }
  #head(true_utilities, 27)
}



# Here is the end of this file


