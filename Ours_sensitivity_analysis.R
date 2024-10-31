# Sensitivity Analysis: N*, interval width

rm(list=ls())

# Load the true prob and utility
source("./R files/scenario_generator_local.R")

model <- cmdstan_model(stan_file = './R files/qtpi_CC.stan')

# check the w_e (0.7 or 0.45)
print(paste0("The current w_e is ", w_e))

# file name distinguishing
if(w_e > 0.5){INX <- "HE_"} else {INX <- "LT_"}

# Sensitivity Analysis parameter check
print(paste0("The current interval width is ", width_SA)) # change it in the scenario generator.R

N_star <- 6 # 3 / 6 (default) / 9
print(paste0("The current N* is ", N_star))

# path
path_a <- "./Outputs/Ours_SA/"
path_allocation <- "allocation/"
path_MTD_OBD <- "MTD_OBD/"

# fix scenario prefix
tox_abb <-c("MHM", "HLH", "HML", "HHM", "LLM", "LLM", "HHH", "HHH", "MMH")
eff_abb <-c("HLM", "LLH", "MLH", "LLH", "LLL", "HMH", "MHM", "LML", "LMM")

# get the simulation results
fix_tox <- data.frame()
fix_eff <- data.frame()
fix_u <- data.frame()

for (i in 1:length(tox_abb)) {
  fix_tox <- rbind(fix_tox, scenarios_prob_tox[rownames(scenarios_prob_tox) == tox_abb[i],]) # true tox prob (unlist when use)
  fix_eff <- rbind(fix_eff, scenarios_prob_eff[rownames(scenarios_prob_eff) == eff_abb[i],]) # true eff prob
  fix_u <- rbind(fix_u, true_utilities[rownames(true_utilities) == paste0(tox_abb[i], "_", eff_abb[i]),]) # true utility
}

# true MTD & OBD for the fix 9 scenarios
MTD_OBD_true_fix9 <- data.frame()
rownames_fix9 <- vector()
for (i in 1:length(tox_abb)){
  MTD_OBD_true_fix9[i,1] <- MTD_OBD_true[paste0("Tox_", tox_abb[i]), 1]
  MTD_OBD_true_fix9[i,2] <- MTD_OBD_true[paste0("Tox_", tox_abb[i]), paste0("Eff_", eff_abb[i])]
  rownames_fix9[i] <- paste0(tox_abb[i], "_", eff_abb[i])
}
rownames(MTD_OBD_true_fix9) <- rownames_fix9
colnames(MTD_OBD_true_fix9) <- c("MTD", "OBD")


# generate Scn_T / Scn_E for the fix 9 scenarios
transition_level_tox_fix9 <- list()
transition_level_eff_fix9 <- list()

for (i in 1:length(tox_abb)){
  transition_level_tox_fix9[[i]] <- scenario_list[[which(rownames(scenarios_prob_tox) == tox_abb[i])]]
  transition_level_eff_fix9[[i]] <- scenario_list[[which(rownames(scenarios_prob_eff) == eff_abb[i])]]
}

# ----------------- Simulation Starts ------------------ #

set.seed(315)

options(mc.cores = parallel::detectCores()) 

# simulation times for each trial
simN <- 100 # 100

# SA parameter setting
# 1. interval width
#width_SA <- 0.05 # 0.05 / 0.1 (default) / 0.2
#{
#  int.bound <- seq(0, 1, width_SA) # equal-width toxicity intervals
#  k_star <- which.max(as.numeric(table(cut(q_T, int.bound)))) # location of the highest tolerable toxicity interval
  
  # equal-width efficacy intervals, storing the lower bounds and upper bounds separately
#  u.lower <- seq(0.0, 1 - width_SA, by = width_SA)
#  u.upper <- seq(width_SA, 1, by = width_SA)
#}

# 2. N*
#N_star <- 3 # 3 / 6 (default) / 9


{
  for (fix in 1:length(tox_abb)) { # i in 1:9 which is the total number of fix scenarios
    
    current_combination <- fix
    
    Scn_T_fix9 <- Gen_Scn_T(transition_level = transition_level_tox_fix9[[fix]], 
                            p_A_tox, p_B_tox, p_C_tox, 
                            D = 6, W = W)
    Scn_E_fix9 <- Gen_Scn_E(transition_level = transition_level_eff_fix9[[fix]], 
                            p_1_eff, p_2_eff, p_3_eff, 
                            D = 6, E = E)
    
    early_stop <- numeric()
    MTD_OBD <- data.frame()
    patients_allocation <- matrix(NA, nrow = simN, ncol = D)
    PA <- numeric() # poor allocation indicator 
    OD <- numeric() # overdose patient number
    
    for (sim in 1:simN) { 
      set.seed(sim)
      
      n <- rep(0,D)  # to store total patients allocated for each dose in this round
      d <- startdose
      cohort_no. <- 1 # cohort number
      toxdose <- D + 1 # over toxic dose indicator: when toxdose = M the maximum non-toxic dose is M - 1 (used in the dose-finding section)
      st <- 0 # sign of terminated, 0 for not and 1 for so
      st_early <- 0
      
      EU_hat <- rep(0,D)
      
      cohort_sequence <- numeric()
      k_EqTP_candidate <- numeric()
      
      # values lists
      {
        # toxicity observation
        tox_response <- list(d_1 = list(),
                             d_2 = list(),
                             d_3 = list(),
                             d_4 = list(),
                             d_5 = list(),
                             d_6 = list()) 
        # efficacy observation
        eff_response <- list(d_1 = list(),
                             d_2 = list(),
                             d_3 = list(),
                             d_4 = list(),
                             d_5 = list(),
                             d_6 = list()) 
        
        # observed qTP
        qTP_list <- list(d_1 = numeric(), 
                         d_2 = numeric(), 
                         d_3 = numeric(), 
                         d_4 = numeric(), 
                         d_5 = numeric(), 
                         d_6 = numeric()) 
        
        # observed qEP
        qEP_list <- list(d_1 = numeric(), 
                         d_2 = numeric(), 
                         d_3 = numeric(), 
                         d_4 = numeric(), 
                         d_5 = numeric(), 
                         d_6 = numeric()) 
        
        # observed utility
        Utility_list <- list(d_1 = numeric(), 
                             d_2 = numeric(), 
                             d_3 = numeric(), 
                             d_4 = numeric(), 
                             d_5 = numeric(), 
                             d_6 = numeric()) 
        
        # posterior mean of EU
        EU_hat_list <- list(d_1 = numeric(), 
                            d_2 = numeric(), 
                            d_3 = numeric(), 
                            d_4 = numeric(), 
                            d_5 = numeric(), 
                            d_6 = numeric()) 
        
        k_U_vector <- vector("numeric", length = 6)
      }
      
      
      ################## Single Trail ###################
      
      while (st == 0)
      {
        early_stop[sim] <- st_early # early stop status of current scenario
        cohort_sequence[cohort_no.] <- d  # track cohort dose allocation
        
        for (i in 1:csize) { ## loop over each patients' observation in current cohort (csize = 3)
          
          ## TOXICITY ##
          
          # total tox types = 3, grade = 0 to 4, totally 5 (3 by 5)
          tox_i <- matrix(NA, nrow = 3, ncol = 5) 
          # probability of grade 0-4 toxicity for dose d (starting from 1), each kind of toxicity respectively
          # p_tox_all = rbind(p_A_tox_all, p_B_tox_all, p_C_tox_all)))
          p_tox_d <- rbind(Scn_T_fix9$p_tox_all[d,], 
                           Scn_T_fix9$p_tox_all[d + D,],  
                           Scn_T_fix9$p_tox_all[d + 2*D,])
          # MULTILEVEL RESPONSE of i-th patient: ACTUAL grade level for each toxicity of i-th patient, generated based on binomial random
          tox_i <- t(apply(p_tox_d, 1, rmultinom, n = 1, size = 1))
          tox_response[[d]] <- rlist::list.append(tox_response[[d]], tox_i)
          
          
          ## EFFICACY ##
          
          # total eff types = 3, grade = Low / Medium / High, totally 3 (3 by 3)
          eff_i <- matrix(NA, nrow = 3, ncol = 3)
          p_eff_d <- rbind(Scn_E_fix9$p_eff_all[d,],  
                           Scn_E_fix9$p_eff_all[d + D,],
                           Scn_E_fix9$p_eff_all[d + 2*D,])
          eff_i <- t(apply(p_eff_d, 1, rmultinom, n = 1, size = 1))
          eff_response[[d]] <- rlist::list.append(eff_response[[d]], eff_i) 
        }
        
        
        # Observed qTP/ qEP and utility
        {
          ## update observed qTP at each level
          for(j in 1:max(cohort_sequence)){ # to include the possible next level
            for(i in 1:length(tox_response[[j]])){
              qTP_list[[j]][i] <- sum(tox_response[[j]][[i]] * W) / sum(apply(W, 1, max))
            } 
          }
          #qTP_list
          
          ## update observed qEP at each level
          for(j in 1:max(cohort_sequence)){
            for(i in 1:length(eff_response[[j]])){
              qEP_list[[j]][i] <- sum(eff_response[[j]][[i]] * E) / sum(apply(E, 1, max))
            }
          }
          #qEP_list
          
          ## update Utility at each level
          for(j in 1:max(cohort_sequence)){
            for(i in 1:length(qEP_list[[j]])){
              Utility_list[[j]][i] <- w_t * (1 - qTP_list[[j]][i]) + w_e * qEP_list[[j]][i]
            }
          }
          Utility_list 
        }
        
        ## update sample size for current dose d
        n[d] <- n[d] + csize  
        
        ## MCMC for the posterior parameter of a and b in EqTP by stan
        l <- max(cohort_sequence)
        qtp_sum<- unlist(sapply(qTP_list, sum))
        
        data <- list(
          l = l,
          n = as.array(n[1:l]),
          sum_qtp = as.array(qtp_sum[1:l]),
          x = as.array(x[1:l]), # value of dose at dose d
          a = -5,
          b = 5,
          c = 8
        )
        
        # model: the MCMC model result has been generated in "scenario_generator.R"
        
        fit <- model$sample(
          data = data,
          seed = 315,
          chains = 2,
          parallel_chains = 2,
          refresh = 2000, 
          iter_sampling = 2000,
          iter_warmup = 2000, 
          thin = 1
        )
        
        params <- fit$draws(format = "df")
        a_post <- params$alpha 
        b_post <- params$beta
        
        # ------------------ The Dose Finding Procedure -------------------- #
        
        # record the current dose level
        d_current <- d;
        
        # original Admi_set before the whole process getting started
        Admi_set <- c(1:(toxdose - 1)) # now the toxdose == D + 1 as default
        
        cat("The Admi_set now is:", Admi_set, ".")
        
        # Update the toxdose and Admi_set with the toxicity safety rule before each new cohort enrolled
        EqTP_current <- 1 / (1 + exp(-a_post - b_post * x[d_current])) 
        if(mean(EqTP_current > q_T) > 0.95){
          if(d_current == 1){st_early <- 1; st <- 1; early_stop[sim] <- st_early; break;} # Admi_set is empty, early stop the trial
          else {toxdose <- d_current} # update the toxdose, the maximum admirable dose is toxdose-1 
        }
        print(paste0("The mean(EqTP_current > q_T) is: ", mean(EqTP_current > q_T), "; the maximum dose in Admi_set now is: ", toxdose-1))
        
        Mz_post <- as.numeric(table(cut(EqTP_current, int.bound))) # formula (5) 
        k_T <- which.max(Mz_post) # k_T the highest marginal posterior probability
        
        print(paste0("The current k_T is: ", k_T))
        
        # ------------------ DOSE FINDING ------------------- #
        
        if (k_T > k_star | toxdose == d_current){ # get into the over-toxic interval OR the update of Admi_set for the next cohort would be at the beginning when new patients enrolled  
          # and need to decrease the dose to the lower one
          if (d_current > 1) { d_next <- d_current - 1;} 
          else {d_next <- 1;} # even though that now for toxdose == d_current, the d_current would not be 1 (must be early terminated)
          
        } else { # k_T <= k_star, judge the decision accordingly
          
          ## get k_U
          k_U_candidates <- numeric()
          u_sum <- unlist(sapply(Utility_list, sum)) # sum of observed utilities for each group of n_j
          
          ## get EU_hat for current dose
          EU_hat[d_current] <- (alpha_u + u_sum[d_current]) / (alpha_u + beta_u + n[d_current])

          {
            if (d_current == 1){
              d_upper <- min(d_current + 1, toxdose - 1);
              ifelse(d_upper == 1, Admi_subset <- 1, Admi_subset <- c(1, 2));
            }
            else if (d_current < max(Admi_set)){
              d_lower <- d_current - 1; 
              d_upper <- min(d_current + 1, toxdose - 1);
              Admi_subset <- c(d_lower, d_current, d_upper);
            }
            else { # d_current == max(Admi_set), i.e., Admi_set now is (1,2,3,4) and d_current == 4, the Admi_subset should be (3,4)
              Admi_subset <- c(d_current-1, max(Admi_subset));
            }
          }
          
          cat("Admi Subset now is: ", Admi_subset, " | d_current is: ", d_current, " | ")
          
          # Calculate the k_U's for the doses among the Admi_subset
          for (j in Admi_subset[1]:Admi_subset[length(Admi_subset)]){ 
            k_U_vector[j] <- which.max(pbeta(u.upper, alpha_u + u_sum[j], beta_u + n[j] - u_sum[j]) - pbeta(u.lower, alpha_u + u_sum[j], beta_u + n[j] - u_sum[j]))
          }
          k_U_candidates <- k_U_vector[Admi_subset]
          
          cat("The k_U_candidates for Admi Subset are: ", k_U_candidates, " | ")
          
          # ----------- with the k_U's, execute the decision procedure --------- #
          
          if (k_T == k_star) {
            
            # judge the possible {d-1, d} OR {d-1, d, d+1} allowance based on the n[d_current]
            print(paste0("Now k_T == k_star. The patient number of current dose is: ", n[d_current], "."))
            
            if(length(k_U_candidates) == 2) { 
              
              if(k_U_candidates[1] == k_U_candidates[2]) { # When equal: Decisions cannot be made directly based on the size relationship of k_U's.
                if(d_current == 1 & n[d_current] >= N_star) {d_next <- d_current} # when k_T == k_star AND n[d_current] >= N_star, subset is {d-1, d}, the maximum dose is d_current
                else if(d_current == 6) { # d_current == 6, d_lower == 5 (which should be the most common scenario)
                  indicator <- sample(0:1, 1, prob = c(0.3,0.7)); # double check the lower dose to avoid getting stuck at the end of the dose series
                  d_next <- ifelse(indicator == 1, Admi_subset[1], Admi_subset[2]);
                } 
                else {indicator <- sample(0:1, 1); # Make an indiscriminate equal-probability selection because equal k_U values represent approximate equality in probability.
                d_next <- ifelse(indicator == 1, Admi_subset[1], Admi_subset[2]);
                }
              } 
              else { # k_U_candidates[1] != k_U_candidates[2]
                if(d_current == 1 & n[d_current] >= N_star) {d_next <- d_current}  # when k_T == k_star AND n[d_current] >= N_star, subset is {d-1, d}, the maximum dose is d_current
                else {d_next <- Admi_subset[which.max(k_U_candidates)] } 
              }
            } 
            
            else { # length(k_U_candidates) == 3, which means Admi subset are 3 distinct doses: 2 3 4 // 3 4 5 etc.
              
              if (n[d_current] >= N_star) { #subset must be (d-1, d)
                if (k_U_candidates[1] == k_U_candidates[2]) {
                  indicator <- sample(0:1, 1); 
                  d_next <- ifelse(indicator == 1, Admi_subset[1], Admi_subset[2]) # Equally likely to choose either d_lower OR d_current.
                } 
                else {d_next <- Admi_subset[which.max(k_U_candidates[1:2])]}
              }
              else { # n[d_current] < N_star, consider 3 subset doses (d-1, d, d+1) based on k_T == k_star
                
                print(paste0("The patient number of upper dose is: ", n[d_upper]))
                # when there are 3 doses in the Admi subset there are four special cases where a unique maximum k_U is not produced, in addition to the normal situation where a unique maximum k_U exists:
                # 1. same-same-same: 
                if (k_U_candidates[1] == k_U_candidates[2] & k_U_candidates[2] == k_U_candidates[3]) {
                  d_next <- ifelse(n[d_upper] == 0, 
                                   Admi_subset[3], # n[d_current] < N_star & n[d_upper] == 0, encourage aggressive dose exploration
                                   sample(Admi_subset, 1, prob = c(0.4,0.4,0.2))) # otherwise randomly select, give lower prob for upper dose (k_T == k_star) 
                }
                # 2. same-same-low: When d-1 and d are approximately equivalent and the number of n[d_current] is low, maintain d_current in the next cohort to reconfirm by increasing the number of patients.
                else if (k_U_candidates[1] == k_U_candidates[2] & k_U_candidates[2] > k_U_candidates[3]) {d_next <- d_current}
                # 3. low-same-same: 
                else if (k_U_candidates[1] < k_U_candidates[2] & k_U_candidates[2] == k_U_candidates[3]) {
                  d_next <- ifelse(n[d_upper] == 0, 
                                   Admi_subset[3], # when n[d_current] < N_star & n[d_upper] == 0, encourage aggressive dose exploration
                                   sample(Admi_subset[2:3], 1, prob = c(0.7,0.3))) # otherwise randomly select, give lower prob for upper dose (k_T == k_star)
                }
                # 4. same-low-same: d-1 and d+1 are approximately equivalent
                else if (k_U_candidates[1] == k_U_candidates[3] & k_U_candidates[2] < k_U_candidates[3]) {
                  d_next <- ifelse(n[d_upper] == 0, 
                                   sample(Admi_subset[c(1,3)], 1, prob = c(0.4,0.6)), # 当n[d_current] < N_star & n[d_upper] == 0, encourage aggressive dose exploration
                                   sample(Admi_subset[c(1,3)], 1, prob = c(0.7,0.3))) # otherwise randomly select, give lower prob for upper dose (k_T == k_star)
                }
                # *. unique maximum exists
                else {d_next <- Admi_subset[which.max(k_U_candidates)]}
                
              }
            }
          }
          
          else { # k_T < k_star: the current dose is at the safe interval
            
            if(length(k_U_candidates) == 2) {
              if(k_U_candidates[1] == k_U_candidates[2]){
                if(d_current == 6) { 
                  indicator <- sample(0:1, 1, prob = c(0.3,0.7)); # double check the lower dose to avoid getting stuck at the end of the dose series
                  d_next <- ifelse(indicator == 1, Admi_subset[1], Admi_subset[2]);
                }
                else {
                  indicator <- sample(0:1, 1); 
                  d_next <- ifelse(indicator == 1, Admi_subset[1], Admi_subset[2])
                }
              }
              else { # k_U_candidates[1] != k_U_candidates[2]
                d_next <- Admi_subset[which.max(k_U_candidates)]
              }
            }
            else { # length(k_U_candidates) == 3, which means Admi subset are 3 distinct doses: 2 3 4 // 3 4 5 etc.
              
              print(paste0("The patient number of upper dose is: ", n[d_upper]))
              
              if (k_U_candidates[1] == k_U_candidates[2] & k_U_candidates[2] == k_U_candidates[3]){
                d_next <- ifelse(n[d_upper] == 0, 
                                 Admi_subset[3], # when n[d_upper] == 0, encourage aggressive dose exploration
                                 sample(Admi_subset, 1, prob = c(0.3, 0.3, 0.4))) # otherwise randomly select, give same prob for each dose (k_T < k_star) 
              }
              else if (k_U_candidates[1] == k_U_candidates[2] & k_U_candidates[2] > k_U_candidates[3]){d_next <- sample(Admi_subset[1:2], 1, prob = c(0.4,0.6))}
              else if (k_U_candidates[1] < k_U_candidates[2] & k_U_candidates[2] == k_U_candidates[3]){
                d_next <- ifelse(n[d_upper] == 0, 
                                 Admi_subset[3], # when n[d_upper] == 0, encourage aggressive dose exploration
                                 sample(Admi_subset[2:3], 1)) # otherwise randomly select, give same prob for d and d+1 (k_T < k_star) 
              }
              else if (k_U_candidates[1] == k_U_candidates[3] & k_U_candidates[2] < k_U_candidates[3]){
                d_next <- ifelse(n[d_upper] == 0, 
                                 sample(Admi_subset[c(1,3)], 1, prob = c(0.2,0.8)), # when n[d_upper] == 0, encourage aggressive dose exploration
                                 sample(Admi_subset[c(1,3)], 1, prob = c(0.4,0.6))) # otherwise randomly select, give same prob for each dose (k_T < k_star) 
              }
              else {d_next <- Admi_subset[which.max(k_U_candidates)]}
            }
          }
        }
        
        d <- d_next
        
        ## enroll the next cohort of patients
        cohort_no. <- cohort_no. + 1
        
        ## stop if the maximum number of patients reached
        if (sum(n) >= samplesize){
          st <- 1;
          break
        } else {
          print(paste0("next dose is: ", d))
        }
        
        
      }
      
      ######################## Single Trial Summary #######################
      patients_allocation[sim, ] <- n
      cat("dose allocation #", sim, ": ", n)
      
      # poor allocation (PA) indicator for each round of simulation
      PA[sim] <- ifelse(n[MTD_OBD_true_fix9[fix, 2]] < samplesize*0.2, 1, 0)
      
      # Over dose patient number for each round of simulation
      OD[sim] <- sum(n[which(fix_tox[fix, ] > q_T)])
      
      ## MTD & OBD
      EqTP <- apply(sapply(x, function(x) {1 / (1 + exp(-a_post - b_post * x))}), 2, mean) # Expected (mean of) qTP for each dose level
      
      for(z in 1:length(EqTP)){
        k_EqTP_candidate[z] <- which.max(as.numeric(table(cut(EqTP[z], int.bound))))
      }
      
      #print(k_EqTP_candidate)
      
      if(sum(k_EqTP_candidate > k_star) == length(EqTP)){ 
        MTD_OBD[sim, 1] <- 0 # no MTD 
      } else {
        EqTP_candidate <- EqTP[k_EqTP_candidate <= k_star] 
        MTD_OBD[sim, 1] <- which.min(abs(EqTP_candidate - q_T)) # MTD #should give it a k_EqTP_candidate similar to the true MTD definition
      }
      
      # test
      #print(EU_hat)
      
      if(early_stop[sim] == 1){
        MTD_OBD[sim, 2] <- 0; 
      } else {
        MTD_OBD[sim, 2] <- which.max(head(EU_hat, MTD_OBD[sim, 1])) # OBD
      }
    }
    
    
    #################### All Trials End, Output Results ####################    
    # patient allocation results for all 100 trials of current scenario
    colnames(patients_allocation) <- c("dose_1","dose_2","dose_3","dose_4","dose_5","dose_6")
    write.table(patients_allocation, paste0(path_a, path_allocation, INX, "_patients_allocation_fix_", fix, ".txt"), append = F)
    
    # average patient number for each dose
    average_patient_number[current_combination, 1:6] <- colMeans(patients_allocation)
    colnames(average_patient_number) <- c("dose_1","dose_2","dose_3","dose_4","dose_5","dose_6")
    rownames_APN[current_combination] <- paste0(tox_abb[fix], "_", eff_abb[fix]);
    
    # average dose selection percentage for each dose
    for(s in 0:6){
      average_selection_percentage[current_combination, s+1] <- (sum(MTD_OBD[,2] == s) / simN) * 100
    }
    colnames(average_selection_percentage) <- c("early_stop","dose_1","dose_2","dose_3","dose_4","dose_5","dose_6")
    rownames_ASP[current_combination] <- paste0(tox_abb[fix], "_", eff_abb[fix]);
    
    # MTD & OBD selection results for all trials of current scenario
    colnames(MTD_OBD) <- c("MTD", "OBD")
    write.table(MTD_OBD, paste0(path_a, path_MTD_OBD, INX, "_MTD_OBD_selection_fix_", fix,".txt"), append = F)
    
    # early stop rate for current scenario
    early_stop_rate[current_combination, 1] <- round(mean(early_stop)*100, 2);
    rownames_ESR[current_combination] <- paste0(tox_abb[fix], "_", eff_abb[fix]);
    
    # correct OBD selection percentage for current scenario
    correct_selection_OBD_percentage <- sum(MTD_OBD[,2] == MTD_OBD_true_fix9[fix, 2]) / nrow(MTD_OBD)

    csel_OBD_all[current_combination, 1] <- correct_selection_OBD_percentage * 100;
    csel_OBD_all[current_combination, 2] <- MTD_OBD_true_fix9[fix, 2]; #OBD_true[, q];
    rownames_csel_OBD_all[current_combination] <- paste0(tox_abb[fix], "_", eff_abb[fix]);
    
    # poor allocation result
    # when less than 20% of the total patients allocated to the OBD, this round of simulation is poor allocated
    poor_allocation_rate[current_combination, 1] <- round(mean(PA)*100, 2)
    rownames_PA[current_combination] <- paste0(tox_abb[fix], "_", eff_abb[fix]);
    
    # overdose patient number
    # over-toxic dose: true qTP > q_T == 0.3, those who allocated to these doses are treated as overdose
    average_overdose_number[current_combination, 1] <- round(mean(OD), 2)
    rownames_OD[current_combination] <- paste0(tox_abb[fix], "_", eff_abb[fix]);
    
  }
  
  # average patient number
  rownames(average_patient_number) <- rownames_APN
  colnames(average_patient_number) <- c("dose_1","dose_2","dose_3","dose_4","dose_5","dose_6")
  #print(average_patient_number)
  write.table(average_patient_number, paste0(path_a, N_star, "_", width_SA, "/", INX, "average_patient_number_fix.txt"), append = F, row.names = TRUE)
  
  # average selection percentage
  rownames(average_selection_percentage) <- rownames_ASP
  colnames(average_selection_percentage) <- c("early_stop","dose_1","dose_2","dose_3","dose_4","dose_5","dose_6")
  #print(average_selection_percentage)
  write.table(average_selection_percentage, paste0(path_a, N_star, "_", width_SA, "/", INX, "average_selection_percentage_fix.txt"), append = F, row.names = TRUE)
  
  # early stop rate output
  rownames(early_stop_rate) <- rownames_ESR
  colnames(early_stop_rate) <- "early_stop_rate(%)"
  #print(early_stop_rate)
  write.table(early_stop_rate, paste0(path_a, N_star, "_", width_SA, "/", INX, "early_stop_rate_fix.txt"), append = F, row.names = TRUE)
  
  # poor allocation rate output
  rownames(poor_allocation_rate) <- rownames_PA
  colnames(poor_allocation_rate) <- "poor_allocation_rate(%)"
  write.table(poor_allocation_rate, paste0(path_a, N_star, "_", width_SA, "/", INX, "poor_allocation_rate_fix.txt"), append = F, row.names = TRUE)
  
  # overdose patient numbers output
  rownames(average_overdose_number) <- rownames_OD
  colnames(average_overdose_number) <- "average_overdose_numbers"
  write.table(average_overdose_number, paste0(path_a, N_star, "_", width_SA, "/", INX, "average_overdose_numbers_fix.txt"), append = F, row.names = TRUE)
  
  # correct selected OBD rate output
  rownames(csel_OBD_all) <- rownames_csel_OBD_all
  colnames(csel_OBD_all) <- c("correct_selected_OBD_rate(%)", "OBD_true")
  #print(csel_OBD_all)
  write.table(csel_OBD_all, paste0(path_a, N_star, "_", width_SA, "/", INX, "correct_selected_OBD_rate_fix.txt"), append = F, row.names = TRUE)
}


# ------------------- Dose Finding END ---------------------- #





# ------------ Visualization Sensitivity Analysis ------------ #

# A.efficacy width = 0.10, N* = 3, 6, 9
# A.1 OBD correct selection rate
OBD_rate_N3 <- read.table("./Outputs/Ours_SA/3_0.1/HE_correct_selected_OBD_rate_fix.txt")
OBD_rate_N6 <- read.table("./Outputs/Ours_SA/6_0.1_default/HE_correct_selected_OBD_rate_fix.txt")
OBD_rate_N9 <- read.table("./Outputs/Ours_SA/9_0.1/HE_correct_selected_OBD_rate_fix.txt")

OBD_rate_N_all <- data.frame(matrix(nrow = 9, ncol = 4))
OBD_rate_N_all[,1:3] <- cbind(OBD_rate_N3[,1], OBD_rate_N6[,1], OBD_rate_N9[,1])
OBD_rate_N_all[,4] <- c("1","2","3","4","5","6","7","8","9")
colnames(OBD_rate_N_all) <- c("N_star=3", "N_star=6", "N_star=9", "scenario")
rownames(OBD_rate_N_all) <- seq(1,9,1)

library(tidyr)
OBD_rate_N_all_long <- OBD_rate_N_all %>% pivot_longer(cols = -scenario, names_to = "N_star", values_to = "OBD_rate")

library(ggplot2)
AA1 <- ggplot(OBD_rate_N_all_long, aes(x = scenario, y = OBD_rate, color = N_star, shape = N_star)) +
  geom_point(size = 5) +
  scale_color_manual(values = c("N_star=3" = "#45639e", "N_star=6" = "#abcc87", "N_star=9" = "purple")) +
  scale_shape_manual(values = c("N_star=3" = 15, "N_star=6" = 16, "N_star=9" = 4)) +
  labs(title = "Sensitivity Analysis: N_star",
       x = "Scenario",
       y = "OBD correct selection rate(%)",
       color = "N_star",
       shape = "N_star") + 
  scale_y_continuous(limits = c(0,80), breaks = seq(0, 80, 10), expand = c(0,0)) + 
  #theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + 
  theme_minimal() + 
  theme(axis.text = element_text(family = "Helvetica", size = 20, colour = "black"),
        axis.title = element_text(family = "Helvetica", size = 20, colour = "black"),
        legend.text = element_text(family = "Helvetica", size = 20, colour = "black"),
        legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_blank(),
        #panel.grid = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.4))

# ---------------------------------------------------------------- #

# A.2 Number of patients allocated to the OBD
pts_N3 <- read.table("./Outputs/Ours_SA/3_0.1/HE_average_patient_number_fix.txt")
pts_N6 <- read.table("./Outputs/Ours_SA/6_0.1_default/HE_average_patient_number_fix.txt")
pts_N9 <- read.table("./Outputs/Ours_SA/9_0.1/HE_average_patient_number_fix.txt")

pts_N_all <- data.frame(matrix(nrow = 9, ncol = 5))
pts_N_all[,4] <- OBD_rate_N3[,2]
pts_N_all[,5] <- c("1","2","3","4","5","6","7","8","9")
colnames(pts_N_all) <- c("N_star=3", "N_star=6", "N_star=9", "true_OBD", "scenario")
rownames(pts_N_all) <- seq(1,9,1)

for(i in 1:9){
  pts_N_all[i,1] <- pts_N3[i, pts_N_all[,4][i]]
  pts_N_all[i,2] <- pts_N6[i, pts_N_all[,4][i]]
  pts_N_all[i,3] <- pts_N9[i, pts_N_all[,4][i]]
}

pts_N_all <- pts_N_all[,-4]

library(tidyr)
pts_N_all_long <- pts_N_all %>% pivot_longer(cols = -scenario, names_to = "N_star", values_to = "No.OBD.pts")

library(ggplot2)
AA2 <- ggplot(pts_N_all_long, aes(x = scenario, y = No.OBD.pts, color = N_star, shape = N_star)) +
  geom_point(size = 5) +
  scale_color_manual(values = c("N_star=3" = "#45639e", "N_star=6" = "#abcc87", "N_star=9" = "purple")) +
  scale_shape_manual(values = c("N_star=3" = 15, "N_star=6" = 16, "N_star=9" = 4)) +
  labs(title = "Sensitivity Analysis: N_star",
       x = "Scenario",
       y = "No. of patients allocated to OBD",
       color = "N_star",
       shape = "N_star") + 
  scale_y_continuous(limits = c(0,15), breaks = seq(0, 15, 3), expand = c(0,0)) + 
  #theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + 
  theme_minimal() + 
  theme(axis.text = element_text(family = "Helvetica", size = 20, colour = "black"),
        axis.title = element_text(family = "Helvetica", size = 20, colour = "black"),
        legend.text = element_text(family = "Helvetica", size = 20, colour = "black"),
        legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_blank(),
        #panel.grid = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.4))

# ---------------------------------------------------------------- #

# B. efficacy width = 0.05, 0.10, 0.20  N* = 6
# B.1 OBD correct selection rate

OBD_width_0.05 <- read.table("./Outputs/Ours_SA/6_0.05/HE_correct_selected_OBD_rate_fix.txt")
OBD_width_0.10 <- read.table("./Outputs/Ours_SA/6_0.1_default/HE_correct_selected_OBD_rate_fix.txt")
OBD_width_0.20 <- read.table("./Outputs/Ours_SA/6_0.2/HE_correct_selected_OBD_rate_fix.txt")

OBD_width_all <- data.frame(matrix(nrow = 9, ncol = 4))
OBD_width_all[,1:3] <- cbind(OBD_width_0.05[,1], OBD_width_0.10[,1], OBD_width_0.20[,1])
OBD_width_all[,4] <- c("1","2","3","4","5","6","7","8","9")
colnames(OBD_width_all) <- c("width=0.05", "width=0.10", "width=0.20", "scenario")
rownames(OBD_width_all) <- seq(1,9,1)

OBD_width_all_long <- OBD_width_all %>% pivot_longer(cols = -scenario, names_to = "width", values_to = "OBD_rate")

BB1 <- ggplot(OBD_width_all_long, aes(x = scenario, y = OBD_rate, color = width, shape = width)) +
  geom_point(size = 5) +
  scale_color_manual(values = c("width=0.05" = "#45639e", "width=0.10" = "#abcc87", "width=0.20" = "purple")) +
  scale_shape_manual(values = c("width=0.05" = 15, "width=0.10" = 16, "width=0.20" = 4)) +
  labs(title = "Sensitivity Analysis: efficacy interval width",
       x = "Scenario",
       y = "OBD correct selection rate(%)",
       color = "width",
       shape = "width") + 
  scale_y_continuous(limits = c(0,80), breaks = seq(0, 80, 10), expand = c(0,0)) + 
  #theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + 
  theme_minimal() + 
  theme(axis.text = element_text(family = "Helvetica", size = 20, colour = "black"),
        axis.title = element_text(family = "Helvetica", size = 20, colour = "black"),
        legend.text = element_text(family = "Helvetica", size = 20, colour = "black"),
        legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_blank(),
        #panel.grid = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.4))

# --------------------------------------------------------------------- #
# B.2 Number of patients allocated to the OBD
pts_width_0.05 <- read.table("./Outputs/Ours_SA/6_0.05/HE_average_patient_number_fix.txt")
pts_width_0.10 <- read.table("./Outputs/Ours_SA/6_0.1_default/HE_average_patient_number_fix.txt")
pts_width_0.20 <- read.table("./Outputs/Ours_SA/6_0.2/HE_average_patient_number_fix.txt")

pts_width_all <- data.frame(matrix(nrow = 9, ncol = 5))
pts_width_all[,4] <- OBD_rate_N3[,2]
pts_width_all[,5] <- c("1","2","3","4","5","6","7","8","9")
colnames(pts_width_all) <- c("width=0.05", "width=0.10", "width=0.20", "true_OBD", "scenario")
rownames(pts_width_all) <- seq(1,9,1)

for(i in 1:9){
  pts_width_all[i,1] <- pts_width_0.05[i, pts_width_all[,4][i]]
  pts_width_all[i,2] <- pts_width_0.10[i, pts_width_all[,4][i]]
  pts_width_all[i,3] <- pts_width_0.20[i, pts_width_all[,4][i]]
}

pts_width_all <- pts_width_all[,-4]

pts_width_all_long <- pts_width_all %>% pivot_longer(cols = -scenario, names_to = "width", values_to = "No.OBD.pts")

BB2 <- ggplot(pts_width_all_long, aes(x = scenario, y = No.OBD.pts, color = width, shape = width)) +
  geom_point(size = 5) +
  scale_color_manual(values = c("width=0.05" = "#45639e", "width=0.10" = "#abcc87", "width=0.20" = "purple")) +
  scale_shape_manual(values = c("width=0.05" = 15, "width=0.10" = 16, "width=0.20" = 4)) +
  labs(title = "Sensitivity Analysis: efficacy interval width",
       x = "Scenario",
       y = "No. of patients allocated to OBD",
       color = "width",
       shape = "width") + 
  scale_y_continuous(limits = c(0,15), breaks = seq(0, 15, 3), expand = c(0,0)) + 
  #theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + 
  theme_minimal() + 
  theme(axis.text = element_text(family = "Helvetica", size = 20, colour = "black"),
        axis.title = element_text(family = "Helvetica", size = 20, colour = "black"),
        legend.text = element_text(family = "Helvetica", size = 20, colour = "black"),
        legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_blank(),
        #panel.grid = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.4))

# merge as one plot together
library(ggpubr)
ggarrange(AA1, AA2,
          #labels = c("A1", "A2"), 
          ncol = 2, nrow = 1, 
          font.label = list(size = 20, family = "Helvetica", color = "black"),
          common.legend = TRUE, legend = "top")

ggarrange(BB1, BB2,
          #labels = c("B1", "B2"),
          ncol = 2, nrow = 1, 
          font.label = list(size = 20, family = "Helvetica", color = "black"),
          common.legend = TRUE, legend = "top")



