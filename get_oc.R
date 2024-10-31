# get simulation oc results

rm(list=ls())

# Load the true prob and utility
source("./scenario_generator_local.R")

# check the w_e
print(paste0("The current w_e is ", w_e))

# file name distinguishing
if(w_e > 0.5){INX <- "HE_"} else {INX <- "LT_"}


set.seed(315)

#library(rlist)
#library(cmdstanr)
options(mc.cores = parallel::detectCores()) 

# simulation times for each trial
simN <- 100 #  1000


# -------- path for saving simulation results -------- #
p <- commandArgs(trailingOnly = TRUE)
p <- as.numeric(p)
path_a <- paste0("./your.path_", p, "/")  

# check the value of p
print(paste0("The current subset is: ", p))
print(paste0("The current path is: ", path_a))

path_allocation <- "allocation/"
path_MTD_OBD <- "MTD_OBD/"


# data in source file
MTD_OBD_true_p <- MTD_OBD_true[p,]
transition_level_tox <- scenario_list[[p]]
# TOX scenario generation for current p
Scn_T <- Gen_Scn_T(transition_level = transition_level_tox, 
                   p_A_tox, p_B_tox, p_C_tox, 
                   D = 6, W = W)

# -------------- Simulation --------------- # 

for (q in 1:27) { #length(scenario_list)
  
  current_combination <- q  # (p - 1) * length(scenario_list) + q
  transition_level_eff <- scenario_list[[q]] # transition level of efficacy
  
  # EFF scenario generation
  Scn_E <- Gen_Scn_E(transition_level = transition_level_eff, 
                     p_1_eff, p_2_eff, p_3_eff, 
                     D = 6, E = E)
  
  early_stop <- numeric()
  MTD_OBD <- data.frame()
  patients_allocation <- matrix(NA, nrow = simN, ncol = D) # total patients allocated for each dose in each simulation round
  PA <- numeric() # poor allocation indicator 
  OD <- numeric() # overdose patient number
  
  ### START SUMULATION ###
  
  for (sim in 1:simN) { 
    set.seed(sim)
    
    # print(paste0("Running the No.", sim, " simulation for the current scenario combination...")) # which would be shown during the running to check the progress
    
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
      
      # k_U
      #k_U_list <- list(d_1 = numeric(), 
      #                 d_2 = numeric(), 
      #                 d_3 = numeric(), 
      #                 d_4 = numeric(), 
      #                 d_5 = numeric(), 
      #                 d_6 = numeric()) 
      
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
        p_tox_d <- rbind(Scn_T$p_tox_all[d,], 
                         Scn_T$p_tox_all[d + D,],  
                         Scn_T$p_tox_all[d + 2*D,])
        # MULTILEVEL RESPONSE of i-th patient: ACTUAL grade level for each toxicity of i-th patient, generated based on binomial random
        tox_i <- t(apply(p_tox_d, 1, rmultinom, n = 1, size = 1))
        tox_response[[d]] <- rlist::list.append(tox_response[[d]], tox_i)
        
        
        ## EFFICACY ##
        
        # total eff types = 3, grade = Low / Medium / High, totally 3 (3 by 3)
        eff_i <- matrix(NA, nrow = 3, ncol = 3)
        p_eff_d <- rbind(Scn_E$p_eff_all[d,],  
                         Scn_E$p_eff_all[d + D,],
                         Scn_E$p_eff_all[d + 2*D,])
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
      
      # model: the MCMC model result has been generated in "scenario_generator_local.R"
      
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
      
      # ------------------ Dose Finding Procedure -------------------- #
      
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
          else { 
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
            
            if(k_U_candidates[1] == k_U_candidates[2]) { 
              if(d_current == 1 & n[d_current] >= 6) {d_next <- d_current} 
              else if(d_current == 6) { 
                indicator <- sample(0:1, 1, prob = c(0.3,0.7));
                d_next <- ifelse(indicator == 1, Admi_subset[1], Admi_subset[2]);
              } 
              else {indicator <- sample(0:1, 1); 
              d_next <- ifelse(indicator == 1, Admi_subset[1], Admi_subset[2]);
              }
            } 
            else { 
              if(d_current == 1 & n[d_current] >= 6) {d_next <- d_current}  
              else {d_next <- Admi_subset[which.max(k_U_candidates)] } 
            }
          } 
          
          else { 
            
            if (n[d_current] >= 6) { 
              if (k_U_candidates[1] == k_U_candidates[2]) {
                indicator <- sample(0:1, 1); 
                d_next <- ifelse(indicator == 1, Admi_subset[1], Admi_subset[2]) 
              } 
              else {d_next <- Admi_subset[which.max(k_U_candidates[1:2])]}
            }
            else { 
              
              print(paste0("The patient number of upper dose is: ", n[d_upper]))
            
              if (k_U_candidates[1] == k_U_candidates[2] & k_U_candidates[2] == k_U_candidates[3]) {
                d_next <- ifelse(n[d_upper] == 0, 
                                 Admi_subset[3], 
                                 sample(Admi_subset, 1, prob = c(0.4,0.4,0.2))) 
              }
            
              else if (k_U_candidates[1] == k_U_candidates[2] & k_U_candidates[2] > k_U_candidates[3]) {d_next <- d_current}

              else if (k_U_candidates[1] < k_U_candidates[2] & k_U_candidates[2] == k_U_candidates[3]) {
                d_next <- ifelse(n[d_upper] == 0, 
                                 Admi_subset[3], 
                                 sample(Admi_subset[2:3], 1, prob = c(0.7,0.3))) 
              }

              else if (k_U_candidates[1] == k_U_candidates[3] & k_U_candidates[2] < k_U_candidates[3]) {
                d_next <- ifelse(n[d_upper] == 0, 
                                 sample(Admi_subset[c(1,3)], 1, prob = c(0.4,0.6)), 
                                 sample(Admi_subset[c(1,3)], 1, prob = c(0.7,0.3))) 
              }

              else {d_next <- Admi_subset[which.max(k_U_candidates)]}
              
            }
          }
        }
        
        else { 
          
          if(length(k_U_candidates) == 2) {
            if(k_U_candidates[1] == k_U_candidates[2]){
              if(d_current == 6) { 
                indicator <- sample(0:1, 1, prob = c(0.3,0.7)); 
                d_next <- ifelse(indicator == 1, Admi_subset[1], Admi_subset[2]);
              }
              else {
                indicator <- sample(0:1, 1); 
                d_next <- ifelse(indicator == 1, Admi_subset[1], Admi_subset[2])
              }
            }
            else { 
              d_next <- Admi_subset[which.max(k_U_candidates)]
            }
          }
          else { 
            
            print(paste0("The patient number of upper dose is: ", n[d_upper]))
            
            if (k_U_candidates[1] == k_U_candidates[2] & k_U_candidates[2] == k_U_candidates[3]){
              d_next <- ifelse(n[d_upper] == 0, 
                               Admi_subset[3], 
                               sample(Admi_subset, 1, prob = c(0.3, 0.3, 0.4))) 
            }
            else if (k_U_candidates[1] == k_U_candidates[2] & k_U_candidates[2] > k_U_candidates[3]){d_next <- sample(Admi_subset[1:2], 1, prob = c(0.4,0.6))}
            else if (k_U_candidates[1] < k_U_candidates[2] & k_U_candidates[2] == k_U_candidates[3]){
              d_next <- ifelse(n[d_upper] == 0, 
                               Admi_subset[3],
                               sample(Admi_subset[2:3], 1)) 
            }
            else if (k_U_candidates[1] == k_U_candidates[3] & k_U_candidates[2] < k_U_candidates[3]){
              d_next <- ifelse(n[d_upper] == 0, 
                               sample(Admi_subset[c(1,3)], 1, prob = c(0.2,0.8)), 
                               sample(Admi_subset[c(1,3)], 1, prob = c(0.4,0.6))) 
            }
            else {d_next <- Admi_subset[which.max(k_U_candidates)]}
          }
        }
      }
      
      d <- d_next
      print(paste0("next dose is: ", d))
      
      cohort_no. <- cohort_no. + 1
      
      if (sum(n) >= samplesize){
        st <- 1
        break
      }
      
      
    }
    
    ######################## Single Trial Summary #######################
    patients_allocation[sim, ] <- n
    cat("dose allocation #", sim, ": ", n)
    
    # poor allocation (PA) indicator for each round of simulation
    PA[sim] <- ifelse(n[MTD_OBD_true_p[, q+1]] < samplesize*0.2, 1, 0)
    
    # Over dose patient number for each round of simulation
    OD[sim] <- sum(n[which(true_tox_prob_list[[p]] > q_T)])
    
    ## MTD & OBD
    
    EqTP <- apply(sapply(x, function(x) {1 / (1 + exp(-a_post - b_post * x))}), 2, mean) # Expected (mean of) qTP for each dose level
    
    for(i in 1:length(EqTP)){
      k_EqTP_candidate[i] <- which.max(as.numeric(table(cut(EqTP[i], int.bound))))
    }
    

    if(sum(k_EqTP_candidate > k_star) == length(EqTP)){ 
      MTD_OBD[sim, 1] <- 0 # no MTD 
    } else {
      EqTP_candidate <- EqTP[k_EqTP_candidate <= k_star] 
      MTD_OBD[sim, 1] <- which.min(abs(EqTP_candidate - q_T)) # MTD 
    }
    
    
    if(early_stop[sim] == 1){
      MTD_OBD[sim, 2] <- 0; 
    } else {
      MTD_OBD[sim, 2] <- which.max(head(EU_hat, MTD_OBD[sim, 1])) # OBD
    }
  }
  
  
  #################### All Trials End, Output Results ####################
  scenario_tox_abb <- noquote(toupper(paste0(substr(transition_level_tox, 1, 1), collapse = "")));
  scenario_eff_abb <- noquote(toupper(paste0(substr(transition_level_eff, 1, 1), collapse = "")));
  prefix_scenario[current_combination, 1] <- paste0(scenario_tox_abb, "_", scenario_eff_abb)
  
  # patient allocation results for all 100 trials of current scenario
  colnames(patients_allocation) <- c("dose_1","dose_2","dose_3","dose_4","dose_5","dose_6")
  write.table(patients_allocation, paste0(path_a, path_allocation, INX, scenario_tox_abb, "_", scenario_eff_abb, "_patients_allocation.txt"), append = F)
  
  # average patient number for each dose
  average_patient_number[current_combination, 1:6] <- colMeans(patients_allocation)
  colnames(average_patient_number) <- c("dose_1","dose_2","dose_3","dose_4","dose_5","dose_6")
  rownames_APN[current_combination] <- paste0(scenario_tox_abb, "_", scenario_eff_abb);
  
  # average dose selection percentage for each dose
  for(i in 0:6){
    average_selection_percentage[current_combination, i+1] <- (sum(MTD_OBD[,2] == i) / simN) * 100
  }
  colnames(average_selection_percentage) <- c("early_stop","dose_1","dose_2","dose_3","dose_4","dose_5","dose_6")
  rownames_ASP[current_combination] <- paste0(scenario_tox_abb, "_", scenario_eff_abb);
  
  # MTD & OBD selection results for all trials of current scenario
  colnames(MTD_OBD) <- c("MTD", "OBD")
  write.table(MTD_OBD, paste0(path_a, path_MTD_OBD, INX, scenario_tox_abb, "_", scenario_eff_abb, "_MTD_OBD_selection.txt"), append = F)
  
  # early stop rate for current scenario
  # print(paste0("The early stop rate for current scenario is: ", mean(early_stop)*100, "%"))
  early_stop_rate[current_combination, 1] <- round(mean(early_stop)*100, 2);
  rownames_ESR[current_combination] <- paste0(scenario_tox_abb, "_", scenario_eff_abb);
  
  # correct OBD selection percentage for current scenario
  correct_selection_OBD_percentage <- sum(MTD_OBD[,2] == MTD_OBD_true_p[, q + 1]) / nrow(MTD_OBD)
  #correct_selection_OBD_percentage <- sum(MTD_OBD[,2] == OBD_true[, q]) / nrow(MTD_OBD)
  #print(paste0("The correct selection rate of OBD for current scenario is: ", correct_selection_OBD_percentage*100, "%. ",
  #             "The true OBD is: dose ", OBD_true[p, q],"."))
  
  csel_OBD_all[current_combination, 1] <- correct_selection_OBD_percentage * 100;
  csel_OBD_all[current_combination, 2] <- MTD_OBD_true_p[, q + 1]; #OBD_true[, q];
  rownames_csel_OBD_all[current_combination] <- paste0(scenario_tox_abb, "_", scenario_eff_abb);
  
  # poor allocation result
  # when less than 20% of the total patients allocated to the OBD, this round of simulation is poor allocated
  poor_allocation_rate[current_combination, 1] <- round(mean(PA)*100, 2)
  rownames_PA[current_combination] <- paste0(scenario_tox_abb, "_", scenario_eff_abb);
  
  # overdose patient number
  # over-toxic dose: true qTP > q_T == 0.3, those who allocated to these doses are treated as overdose
  average_overdose_number[current_combination, 1] <- round(mean(OD), 2)
  rownames_OD[current_combination] <- paste0(scenario_tox_abb, "_", scenario_eff_abb);
  
}



# prefix record output
# colnames(prefix_scenario) <- "prefix_of_scenario"
# write.table(prefix_scenario, paste0(path, "prefix_scenario.txt"), append = F, row.names = TRUE)

# average patient number
rownames(average_patient_number) <- rownames_APN
colnames(average_patient_number) <- c("dose_1","dose_2","dose_3","dose_4","dose_5","dose_6")
#print(average_patient_number)
write.table(average_patient_number, paste0(path_a, INX, "average_patient_number_", p ,".txt"), append = F, row.names = TRUE)

# average selection percentage
rownames(average_selection_percentage) <- rownames_ASP
colnames(average_selection_percentage) <- c("early_stop","dose_1","dose_2","dose_3","dose_4","dose_5","dose_6")
#print(average_selection_percentage)
write.table(average_selection_percentage, paste0(path_a, INX, "average_selection_percentage_", p ,".txt"), append = F, row.names = TRUE)

# early stop rate output
rownames(early_stop_rate) <- rownames_ESR
colnames(early_stop_rate) <- "early_stop_rate(%)"
#print(early_stop_rate)
write.table(early_stop_rate, paste0(path_a, INX, "early_stop_rate_", p ,".txt"), append = F, row.names = TRUE)

# poor allocation rate output
rownames(poor_allocation_rate) <- rownames_PA
colnames(poor_allocation_rate) <- "poor_allocation_rate(%)"
write.table(poor_allocation_rate, paste0(path_a, INX, "poor_allocation_rate_", p ,".txt"), append = F, row.names = TRUE)

# overdose patient numbers output
rownames(average_overdose_number) <- rownames_OD
colnames(average_overdose_number) <- "average_overdose_numbers"
write.table(average_overdose_number, paste0(path_a, INX, "average_overdose_numbers_", p ,".txt"), append = F, row.names = TRUE)

# correct selected OBD rate output
rownames(csel_OBD_all) <- rownames_csel_OBD_all
colnames(csel_OBD_all) <- c("correct_selected_OBD_rate(%)","OBD_true")
#print(csel_OBD_all)
write.table(csel_OBD_all, paste0(path_a, INX, "correct_selected_OBD_rate_", p ,".txt"), append = F, row.names = TRUE)


# The whole trial ends here




