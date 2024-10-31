# Updated 2024-10-29 Vancouver

# Now is for "HE_" scenarios for all 5 methods

# Organize the 9 selected fix scenario for "fix scenario simulation"
# and result visualization

# *** The Result Figures are in this file

rm(list=ls())
# load the scenario data
source("./R files/scenario_generator_local.R")
# load the oc simulation results of all 5 methods
source("./R files/OBD_selection_visualization_local.R")


{
  library(ggplot2)
  library(dplyr)
  library(tidyverse)
  library(gridExtra)
  library(cowplot)
  library(ggpubr)
}


#MTD_OBD_true
#scenarios_prob_tox
#scenarios_prob_eff

# check the w_e
print(paste0("The current w_e is ", w_e))


# ---------------- Selecting Fix Scenario ---------------- #

# Aiming: top 2 leading OBD, last 2 PA and OD

# ------- OBD correct rate (%) ------- #
{
  Selecting_Dataset_A <- as.data.frame(matrix(nrow = 729, ncol = 6))
  rownames(Selecting_Dataset_A) <- prefix_scenario[,1]
  colnames(Selecting_Dataset_A) <- c("OBD_true", "Iso", "TEPI", "uTPI", "BOIN12", "Proposed")
  
  Selecting_Dataset_A[,1] <- Iso_correct_selected_OBD_rate[, 2]
  Selecting_Dataset_A[,2] <- Iso_correct_selected_OBD_rate[, 1]
  Selecting_Dataset_A[,3] <- TEPI_correct_selected_OBD_rate[, 1]
  Selecting_Dataset_A[,4] <- uTPI_correct_selected_OBD_rate[, 4]
  Selecting_Dataset_A[,5] <- BOIN12_correct_selected_OBD_rate[, 4]
  Selecting_Dataset_A[,6] <- Ours_correct_selected_OBD_rate[, 1]
  
  head(Selecting_Dataset_A)
}

# ------------------------------- #

# -------- OBD patients --------- #
{
  Selecting_Dataset_B <- as.data.frame(matrix(nrow = 729, ncol = 5))
  rownames(Selecting_Dataset_B) <- prefix_scenario[,1]
  colnames(Selecting_Dataset_B) <- c("Iso", "TEPI", "uTPI", "BOIN12", "Proposed")
  
  for (i in 1:5) {
    Selecting_Dataset_B[,i] <- combined_data_OBD_patient_check[, i]
  }
  
  head(Selecting_Dataset_B) 
}
# ------------------------------- #

# ----- Poor Allocation Rate(%) ----- #
{
  Selecting_Dataset_C <- as.data.frame(matrix(nrow = 729, ncol = 5))
  rownames(Selecting_Dataset_C) <- prefix_scenario[,1]
  colnames(Selecting_Dataset_C) <- c("Iso", "TEPI", "uTPI", "BOIN12", "Proposed")
  
  Selecting_Dataset_C[,1] <- Iso_poor_allocation_rate[, 1]
  Selecting_Dataset_C[,2] <- TEPI_poor_allocation_rate[, 1]
  Selecting_Dataset_C[,3] <- uTPI_poor_allocation_rate[, 7]
  Selecting_Dataset_C[,4] <- BOIN12_poor_allocation_rate[, 7]
  Selecting_Dataset_C[,5] <- Ours_poor_allocation_rate[, 1]
  
  head(Selecting_Dataset_C)
}
# ------------------------------- #

# ------ Overdose patients ------- #
{
  Selecting_Dataset_D <- as.data.frame(matrix(nrow = 729, ncol = 5))
  rownames(Selecting_Dataset_D) <- prefix_scenario[,1]
  colnames(Selecting_Dataset_D) <- c("Iso", "TEPI", "uTPI", "BOIN12", "Proposed")
  
  Selecting_Dataset_D[,1] <- Iso_overdose_rate[, 1]
  Selecting_Dataset_D[,2] <- TEPI_overdose_rate[, 1]
  Selecting_Dataset_D[,3] <- uTPI_overdose_rate[, 6]
  Selecting_Dataset_D[,4] <- BOIN12_overdose_rate[, 6]
  Selecting_Dataset_D[,5] <- Ours_overdose_rate[, 1]
  
  tail(Selecting_Dataset_D, n = 30) 
}
# ------------------------------- #

# threshold (to select qualified scenarios)
subset_A <- rownames(Selecting_Dataset_A[Selecting_Dataset_A$Ours == apply(Selecting_Dataset_A[ ,2:6], 1, max), ])
subset_B <- rownames(Selecting_Dataset_B[Selecting_Dataset_B$Ours == apply(Selecting_Dataset_B[ ,1:5], 1, max), ])
#subset_C <- rownames(Selecting_Dataset_C[Selecting_Dataset_C$Ours == apply(Selecting_Dataset_C[ ,1:5], 1, min), ])
#subset_D <- rownames(Selecting_Dataset_D[Selecting_Dataset_D$Ours == apply(Selecting_Dataset_D[ ,1:5], 1, min), ])

# In C and D, identify the two lowest values for each scenario that meet specific criteria (this step cannot simultaneously retrieve the corresponding scenario)."
lowest_2_C <- t(apply(Selecting_Dataset_C[, 1:5], 1, function(x) sort(x)[1:2]))
lowest_2_D <- t(apply(Selecting_Dataset_D[which(Selecting_Dataset_D$Ours != 0 & Selecting_Dataset_D$BOIN12 >= 4), 1:5], 1, function(x) sort(x)[1:2]))

TF_indicator_C <- vector() # Identify the row names (scenarios) where the 'Ours' results are among the two lowest values.
for (i in 1:nrow(lowest_2_C)) {
  TF_indicator_C[i] <- Selecting_Dataset_C$Ours[i] %in% lowest_2_C[i,]
}
lowest_2_CC <- rownames(lowest_2_C[which(TF_indicator_C == TRUE), ])

TF_indicator_D <- vector() 
for (i in 1:nrow(lowest_2_D)) {
  TF_indicator_D[i] <- Selecting_Dataset_D[which(Selecting_Dataset_D$Ours != 0 & Selecting_Dataset_D$BOIN12 >= 4), "Proposed"][i] %in% lowest_2_D[i,]
}
lowest_2_DD <- rownames(lowest_2_D[which(TF_indicator_D == TRUE), ])

subset_C <- lowest_2_CC
subset_D <- lowest_2_DD

subset_ALL <- Reduce(intersect, list(subset_A, subset_B, subset_C, subset_D))
length(subset_ALL)

# The final qualified scenario subsets for all 4 performance 
Selecting_Dataset_A[subset_ALL, ]
Selecting_Dataset_B[subset_ALL, ]
Selecting_Dataset_C[subset_ALL, ]
Selecting_Dataset_D[subset_ALL, ]

# qualified tox and eff prob
scenarios_prob_subset <- scenarios_prob[rownames(scenarios_prob) %in% subset_ALL,]
scenarios_prob_subset_tox <- scenarios_prob_subset[,c(1,3,5,7,9,11)]
scenarios_prob_subset_eff <- scenarios_prob_subset[,c(2,4,6,8,10,12)]
# qualified true utilities
true_utilities_subset <- true_utilities[rownames(true_utilities) %in% subset_ALL,]


# 11.07: record the final qualified scenario subset prefix
subset_ALL[c(1,7,13,39)] # best performance for Ours in all criteria (MTD are all 3)
# higher MTD
# lower MTD

# tox: LLM / HHH
# eff: HMH / LLL / MHM


# ----------------- Selecting Fix Scenario END ------------------ #



# select qualified fix scenario, regarding the OBD selection rate
threshold_OBD_rate <- 20 # : with at least xxx % correct OBD select rate

{
  vector1 <- which(Iso_correct_selected_OBD_rate[,1] >= threshold_OBD_rate)
  vector2 <- which(TEPI_correct_selected_OBD_rate[,1] >= threshold_OBD_rate)
  vector3 <- which(uTPI_correct_selected_OBD_rate$best_dose_avg >= threshold_OBD_rate)
  vector4 <- which(BOIN12_correct_selected_OBD_rate$best_dose_avg >= threshold_OBD_rate)
  vector5 <- which(Ours_correct_selected_OBD_rate[, 1] >= threshold_OBD_rate)
  # get the qualified intersection of all 5 methods: row numbers
  intersection <- Reduce(intersect, list(vector1, vector2, vector3, vector4, vector5))
  # get the qualified intersection of all 5 methods: scenario names tox_abb eff_abb
  qualified_fix_scenario_list <- rownames(Iso_correct_selected_OBD_rate[intersection,])
  
  Iso_correct_selected_OBD_rate[qualified_fix_scenario_list, ]
  TEPI_correct_selected_OBD_rate[qualified_fix_scenario_list, ]
  uTPI_correct_selected_OBD_rate[qualified_fix_scenario_list, ]
  BOIN12_correct_selected_OBD_rate[qualified_fix_scenario_list, ]
  Ours_correct_selected_OBD_rate[qualified_fix_scenario_list, ]
  
  qualified_OBD_data <- cbind(Iso_correct_selected_OBD_rate[qualified_fix_scenario_list, ], TEPI_correct_selected_OBD_rate[qualified_fix_scenario_list, ],
                              uTPI_correct_selected_OBD_rate[qualified_fix_scenario_list, c(4,5)], BOIN12_correct_selected_OBD_rate[qualified_fix_scenario_list, c(4,5)],
                              Ours_correct_selected_OBD_rate[qualified_fix_scenario_list, ])
  qualified_OBD_data <- qualified_OBD_data[,-c(2,4,6,8)]
  colnames(qualified_OBD_data) <- c("Iso", "TEPI", "uTPI", "BOIN12", "Proposed", "OBD_true") 
}

dim(qualified_OBD_data) # total qualified scenarios

# true prob and true utility of qualified scenario
#scenarios_prob[which(rownames(scenarios_prob) %in% rownames(qualified_OBD_data)),] # all true prob

#scenarios_prob[which(rownames(scenarios_prob) %in% rownames(qualified_OBD_data)), c(1,3,5,7,9,11)] # true prob tox
#scenarios_prob[which(rownames(scenarios_prob) %in% rownames(qualified_OBD_data)), c(2,4,6,8,10,12)] # true prob eff

#true_utilities[which(rownames(true_utilities) %in% rownames(qualified_OBD_data)),]

#View(true_utilities[which(rownames(true_utilities) %in% rownames(qualified_OBD_data)),])

#scenarios_prob_eff
#scenarios_prob_tox






# ---------------- fix scenario selection (HE, w_e = 0.7, w_t = 0.3 / LT, w_e = 0.45, w_t = 0.55) --------------------- #
{
  
  colors_barplots <- rev(c("orange", "#45639e", "#57b3ab", "#abcc87", "#d1c76e", "#7e91b1"))
  
  # set the chosen fix scenario list
  #tox_abb <- c("LLL", "HMH", "MLM", "MHL", "MMH", "MHH", "LLM", "HHM", "LML")
  #eff_abb <- c("LML", "LLL", "LLL", "MHM", "LLL", "MLH", "HML", "HHH", "MLL")
  
  #tox_abb <- c("LLL", "LHL", "LHH", "MHH", "HMM", "LLM", "MMH", "HLH", "HMH") # MML the last, #6 possible: LML / LLM
  #eff_abb <- c("LLL", "MMM", "LLL", "LML", "LMH", "MHM", "HHL", "MLL", "HMM") # HLM the last
  
  #tox_abb <- c("LLL", "LHL", "LHH", "MHH", "HMM", "LLM", "MMH", "LLM", "LLL") # MML the last, #6 possible: LML / LLM
  #eff_abb <- c("LLL", "MMM", "LLL", "LML", "LMH", "MHM", "HHL", "MLL", "HHM") # HLM the last
  
  tox_abb <-c("MHM", "HLH", "HML", "HHM", "LLM", "LLM", "HHH", "HHH", "MMH")
  eff_abb <-c("HLM", "LLH", "MLH", "LLH", "LLL", "HMH", "MHM", "LML", "LMM")
  
  # get the simulation results
  fix_tox <- data.frame()
  fix_eff <- data.frame()
  fix_u <- data.frame()
  
  # Method #1
  for (i in 1:length(tox_abb)) {
    fix_tox <- rbind(fix_tox, scenarios_prob_tox[rownames(scenarios_prob_tox) == tox_abb[i],]) # true tox prob (unlist when use)
    fix_eff <- rbind(fix_eff, scenarios_prob_eff[rownames(scenarios_prob_eff) == eff_abb[i],]) # true eff prob
    fix_u <- rbind(fix_u, true_utilities[rownames(true_utilities) == paste0(tox_abb[i], "_", eff_abb[i]),]) # true utility
    #MTD_OBD_true[rownames(MTD_OBD_true) == paste0("Tox_", tox_abb[i]), 1] # true MTD
    #MTD_OBD_true[rownames(MTD_OBD_true) == paste0("Tox_", tox_abb[i]), colnames(MTD_OBD_true) == paste0("Eff_", tox_abb[i])] # true OBD 
  }
  
  # Method #2 (11.07 Updated)
  #fix_tox <- scenarios_prob_subset_tox # true tox prob (unlist when use)
  #fix_eff <- scenarios_prob_subset_eff # true eff prob
  #fix_u <-  true_utilities_subset# true utility
  
  
  colnames(fix_tox) <- c("dose1","dose2","dose3","dose4","dose5","dose6")
  colnames(fix_eff) <- c("dose1","dose2","dose3","dose4","dose5","dose6")
  
  Dose <- seq(1:6)
  fix_tox <- t(fix_tox)
  fix_tox <- cbind(Dose, fix_tox)
  
  fix_eff <- t(fix_eff)
  fix_eff <- cbind(Dose, fix_eff)
  
  fix_u <- t(fix_u)
  fix_u <- cbind(Dose, fix_u)
  
  all_plot_list_fix <- list()
}


# generate the probability and utility curve plots of chosen scenarios 

for (i in 1:length(tox_abb)) { #length(tox_abb) / length(subset_ALL)
  local_i <- i
  
  y_fix_tox <- fix_tox[,local_i+1]
  y_fix_eff <- fix_eff[,local_i+1]
  y_fix_u <- fix_u[,local_i+1]
  
  p <- ggplot() +
    geom_point(data = as.data.frame(fix_tox), aes_(x = Dose, y = y_fix_tox, color = "Toxicity"), size = 3) +
    geom_point(data = as.data.frame(fix_eff), aes_(x = Dose, y = y_fix_eff, color = "Efficacy"), size = 3) +
    geom_point(data = as.data.frame(fix_u), aes_(x = Dose, y = y_fix_u, color = "Utility"), size = 3) + 
    
    geom_line(data = as.data.frame(fix_tox), aes_(x = Dose, y = y_fix_tox, color = "Toxicity"), linewidth = 1, linetype = "dotdash") +
    geom_line(data = as.data.frame(fix_eff), aes_(x = Dose, y = y_fix_eff, color = "Efficacy"), linewidth = 1, linetype = "longdash") +
    geom_line(data = as.data.frame(fix_u), aes_(x = Dose, y = y_fix_u, color = "Utility"), linewidth = 1) +
    
    geom_hline(yintercept = 0.3, linetype = "dashed", color = "black", linewidth = 0.6) + 
    
    scale_color_manual(values = c("Toxicity" = "red", "Efficacy" = "lightblue", "Utility" = "purple")) +     
    labs(x = "Dose Level", y = "Probability or Utility", title = paste0("Scenario ",i)) +
    scale_x_continuous(breaks = c(1,2,3,4,5,6)) +
    scale_y_continuous(breaks = seq(0.0, 1.0, 0.1)) +
    theme_minimal() + 
    theme(axis.text = element_text(family = "Helvetica", size = 20, colour = "black"),
          axis.title = element_text(family = "Helvetica", size = 20, colour = "black"),
          legend.text = element_text(family = "Helvetica", size = 20, colour = "black"),
          legend.title = element_blank(),
          legend.position = "none",
          plot.title = element_text(family = "Helvetica", size = 24, colour = "black"),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black", linewidth = 0.4))
  
  all_plot_list_fix[[i]] <- p
}

# combine all scenarios
#combined_plots_fix <- grid.arrange(grobs = all_plot_list_fix, ncol = 3)

# Method #2 
APLF_1 <- ggarrange(all_plot_list_fix[[1]], all_plot_list_fix[[2]], all_plot_list_fix[[3]], 
                    all_plot_list_fix[[4]], all_plot_list_fix[[5]], all_plot_list_fix[[6]], 
                    all_plot_list_fix[[7]], all_plot_list_fix[[8]], all_plot_list_fix[[9]], 
                    common.legend = TRUE, legend = "bottom", nrow = 3, ncol = 3)#, 
#font.label = list(size = 20, family = "Helvetica", color = "black")) # use for the label but not the legend
# APLF_2 <- ggarrange(all_plot_list_fix[[10]], all_plot_list_fix[[11]], all_plot_list_fix[[12]], 
#               all_plot_list_fix[[13]], all_plot_list_fix[[14]], all_plot_list_fix[[15]], 
#               all_plot_list_fix[[16]], all_plot_list_fix[[17]], all_plot_list_fix[[18]], 
#               common.legend = TRUE, legend = "bottom", nrow = 3, ncol = 3)
# 
# APLF_3 <- ggarrange(all_plot_list_fix[[19]], all_plot_list_fix[[20]], all_plot_list_fix[[21]], 
#               all_plot_list_fix[[22]], all_plot_list_fix[[23]], all_plot_list_fix[[24]], 
#               all_plot_list_fix[[25]], all_plot_list_fix[[26]], all_plot_list_fix[[27]], 
#               common.legend = TRUE, legend = "bottom", nrow = 3, ncol = 3)
# 
# APLF_4 <- ggarrange(all_plot_list_fix[[28]], all_plot_list_fix[[29]], all_plot_list_fix[[30]], 
#               all_plot_list_fix[[31]], all_plot_list_fix[[32]], all_plot_list_fix[[33]], 
#               all_plot_list_fix[[34]], all_plot_list_fix[[35]], all_plot_list_fix[[36]], 
#               common.legend = TRUE, legend = "bottom", nrow = 3, ncol = 3)
# 
# APLF_5 <- ggarrange(all_plot_list_fix[[37]], all_plot_list_fix[[38]], all_plot_list_fix[[39]], 
#               all_plot_list_fix[[40]], all_plot_list_fix[[41]], all_plot_list_fix[[42]],
#               common.legend = TRUE, legend = "bottom", nrow = 3, ncol = 2)
# 
# ggarrange(APLF_1, APLF_2, #APLF_3, APLF_4, #APLF_5,
#           common.legend = TRUE, legend = "bottom", nrow = 1, ncol = 2)
# ggarrange(APLF_3, APLF_4, 
#           common.legend = TRUE, legend = "bottom", nrow = 1, ncol = 2)

# Method #1

ggarrange(all_plot_list_fix[[1]], all_plot_list_fix[[2]], all_plot_list_fix[[3]], 
          all_plot_list_fix[[4]], all_plot_list_fix[[5]], all_plot_list_fix[[6]], 
          all_plot_list_fix[[7]], all_plot_list_fix[[8]], all_plot_list_fix[[9]], 
          common.legend = TRUE, legend = "bottom", nrow = 3, ncol = 3)


# -------------------------------------------------------------------------------- #

# The following updates with gBOIN-ET results (10.29)
gBOIN_ET_results_fix_HE <- readRDS("./Outputs/gBOIN-ET_results_fix_HE.rds") # scenarios 1 to 9
gBOIN_ET_results_fix_LT <- readRDS("./Outputs/gBOIN-ET_results_fix_LT.rds") # scenarios 1 to 9

gBOIN_ET_results_all_HE <- readRDS("./Outputs/gBOIN-ET_results_all_HE.rds") # all 729 random scenarios
gBOIN_ET_results_all_LT <- readRDS("./Outputs/gBOIN-ET_results_all_LT.rds") # all 729 random scenarios

# Poor Allocation by manually updated gboinet function
# access by: gBOIN_ET_poorAllocation_fix_HE[[i]]$poor_allocation_rate
gBOIN_ET_poorAllocation_fix_HE <- readRDS("./Outputs/gBOIN-ET_poorAllocation_fix_HE.rds")
gBOIN_ET_poorAllocation_fix_LT <- readRDS("./Outputs/gBOIN-ET_poorAllocation_fix_LT.rds")

gBOIN_ET_poorAllocation_all_HE <- readRDS("./Outputs/gBOIN-ET_poorAllocation_all_HE.rds")
gBOIN_ET_poorAllocation_all_LT <- readRDS("./Outputs/gBOIN-ET_poorAllocation_all_LT.rds")

# True fix MTD
MTD_true_fix_HE <- c(3, 3, 3, 2, 6, 6, 2, 2, 4)
MTD_true_fix_LT <- c(3, 3, 3, 2, 6, 6, 2, 2, 4)
  
# True fix OBD
OBD_true_fix_HE <- c(3, 3, 3, 2, 3, 6, 2, 1, 4)
OBD_true_fix_LT <- c(1, 1, 2, 1, 2, 6, 1, 1, 1)


# ---- Data Cleaning
# A. OBD correct selection rate (HE and LT)
gBOIN_ET_correct_selected_OBD_rate_HE <- gBOIN_ET_correct_selected_OBD_rate_LT <- data.frame()
for (i in 1:9) {
  gBOIN_ET_correct_selected_OBD_rate_HE <- rbind(gBOIN_ET_correct_selected_OBD_rate_HE, gBOIN_ET_results_fix_HE[[i]]$prop.select)
  gBOIN_ET_correct_selected_OBD_rate_LT <- rbind(gBOIN_ET_correct_selected_OBD_rate_LT, gBOIN_ET_results_fix_LT[[i]]$prop.select)
}
colnames(gBOIN_ET_correct_selected_OBD_rate_HE) <- colnames(gBOIN_ET_correct_selected_OBD_rate_LT) <- 
  c("Dose 1", "Dose 2", "Dose 3", "Dose 4", "Dose 5", "Dose 6")
# correct OBD selection percentage results:
# sapply(seq_along(OBD_true_fix_HE), function(i) gBOIN_ET_correct_selected_OBD_rate_HE[i, OBD_true_fix_HE[i]])
# sapply(seq_along(OBD_true_fix_LT), function(i) gBOIN_ET_correct_selected_OBD_rate_LT[i, OBD_true_fix_LT[i]])



# B. Number of patients allocated to OBD (HE and LT)
gBOIN_ET_correct_OBD_pts_HE <- gBOIN_ET_correct_OBD_pts_LT <- data.frame()
for (i in 1:9) {
  gBOIN_ET_correct_OBD_pts_HE <- rbind(gBOIN_ET_correct_OBD_pts_HE, gBOIN_ET_results_fix_HE[[i]]$n.patient)
  gBOIN_ET_correct_OBD_pts_LT <- rbind(gBOIN_ET_correct_OBD_pts_LT, gBOIN_ET_results_fix_LT[[i]]$n.patient)
}
colnames(gBOIN_ET_correct_OBD_pts_HE) <- colnames(gBOIN_ET_correct_OBD_pts_LT) <- 
  c("Dose 1", "Dose 2", "Dose 3", "Dose 4", "Dose 5", "Dose 6")
# Number of patients allocated to OBD results:
# sapply(seq_along(OBD_true_fix_HE), function(i) gBOIN_ET_correct_OBD_pts_HE[i, OBD_true_fix_HE[i]])
# sapply(seq_along(OBD_true_fix_LT), function(i) gBOIN_ET_correct_OBD_pts_LT[i, OBD_true_fix_LT[i]])



# C. Poor Allocation Rate (HE and LT)
# access by: gBOIN_ET_poorAllocation_fix_HE[[i]]$poor_allocation_rate



# D. Average number of overdose patients (HE and LT)
# gBOIN_ET_correct_OBD_pts_HE
# mapply(function(row, threshold) { if (threshold < 6) { sum(gBOIN_ET_correct_OBD_pts_HE[row, (threshold+1):6]) } else { 0.0 }}, row = 1:9, threshold = MTD_true_fix_HE)






# --------------------- OBD correct selection rate visualization ----------------------- #

# generate the data.frame to store the OBD correct selection rate of chosen fix scenarios
OBD_rate_all_fix <- as.data.frame(matrix(nrow = length(tox_abb), ncol = 6)) # row as scenario 1-9, column as methods A to F

rownames(OBD_rate_all_fix) <- seq(1:9)
colnames(OBD_rate_all_fix) <- c("Iso", "TEPI", "uTPI", "BOIN12", "gBOIN-ET", "Proposed")

# load the correct OBD selection percentage results
for (i in 1:length(tox_abb)) {
  OBD_rate_all_fix[i,1] <- Iso_correct_selected_OBD_rate[which(rownames(Iso_correct_selected_OBD_rate) == paste0(tox_abb[i], "_", eff_abb[i])), 1]
  OBD_rate_all_fix[i,2] <- TEPI_correct_selected_OBD_rate[which(rownames(TEPI_correct_selected_OBD_rate) == paste0(tox_abb[i], "_", eff_abb[i])), 1]
  OBD_rate_all_fix[i,3] <- uTPI_correct_selected_OBD_rate[which(rownames(uTPI_correct_selected_OBD_rate) == paste0(tox_abb[i], "_", eff_abb[i])), 4]
  OBD_rate_all_fix[i,4] <- BOIN12_correct_selected_OBD_rate[which(rownames(BOIN12_correct_selected_OBD_rate) == paste0(tox_abb[i], "_", eff_abb[i])), 4]
  
  OBD_rate_all_fix[i,6] <- Ours_correct_selected_OBD_rate[which(rownames(Ours_correct_selected_OBD_rate) == paste0(tox_abb[i], "_", eff_abb[i])), 1]
}

# add gBOIN-ET results into the data.frame as a new column (notice the HE or LT case)
OBD_rate_all_fix[,5] <- sapply(seq_along(OBD_true_fix_HE), function(i) gBOIN_ET_correct_selected_OBD_rate_HE[i, OBD_true_fix_HE[i]])
# OBD_rate_all_fix[,5] <- sapply(seq_along(OBD_true_fix_LT), function(i) gBOIN_ET_correct_selected_OBD_rate_LT[i, OBD_true_fix_LT[i]])


# generate the barplots to visualize the OBD correct selection rate

# modify the data to visualize
OBD_rate_all_fix$scenario <- rownames(OBD_rate_all_fix)
OBD_rate_all_fix <- OBD_rate_all_fix %>%
  pivot_longer(cols = -scenario, names_to = "Methods", values_to = "OBD_rate") %>%
  mutate(Methods = factor(Methods, levels = c("Iso", "TEPI", "uTPI", "BOIN12", "gBOIN-ET", "Proposed")))
# generate the color vector
#library(RColorBrewer)
#colors_OBD_rate_barplots <- brewer.pal(n = 5, name = "Set1")
#colors_OBD_rate_barplots <- c("#13a983","#00bdcd","#cc340c","#006b7b","#f88421")
# generate the barplot
A <- ggplot(OBD_rate_all_fix, aes(x = scenario, y = OBD_rate, fill = Methods)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) + # , color = "black", linewidth = 0.4
  labs(title = "Selection of OBD (%)", x = "Scenario", y = "Selection of OBD (%)") +
  scale_fill_manual(values = colors_barplots) +
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 20), expand = c(0,0)) + 
  geom_text(aes(label = OBD_rate), position = position_dodge(width = 0.8), vjust = -0.5, size = 6) + 
  theme_minimal() + 
  theme(axis.text = element_text(family = "Helvetica", size = 30, colour = "black"),
        axis.title = element_text(family = "Helvetica", size = 30, colour = "black"),
        legend.text = element_text(family = "Helvetica", size = 30, colour = "black"),
        legend.title = element_blank(),
        #legend.position = "none",
        plot.title = element_text(family = "Helvetica", size = 30, colour = "black"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.4))
# -------------------------------------------------------------------------------- #


# --------------------- Number of patients allocated to OBD visualization --------------------------------- #

# generate the data.frame to store the number of patient of OBD of chosen fix scenarios
OBD_patients_all_fix <- as.data.frame(matrix(nrow = 9, ncol = 6))
rownames(OBD_patients_all_fix) <- seq(1:9)
colnames(OBD_patients_all_fix) <- c("Iso", "TEPI", "uTPI", "BOIN12", "gBOIN-ET", "Proposed")

# get the same qualified subset for OBD patient numbers data
# qualified_OBD_patients_data <- combined_data_OBD_patient_check[which(rownames(combined_data_OBD_patient_check) %in% rownames(qualified_OBD_data)),]


# load the number of patient of OBD results
for (i in 1:length(tox_abb)) {   
  OBD_patients_all_fix[i, 1:4] <- combined_data_OBD_patient_check[which(rownames(combined_data_OBD_patient_check) == paste0(tox_abb[i], "_", eff_abb[i])), 1:4]
  OBD_patients_all_fix[i, 6] <- combined_data_OBD_patient_check[which(rownames(combined_data_OBD_patient_check) == paste0(tox_abb[i], "_", eff_abb[i])), 5]
  
  OBD_patients_all_fix[i, c(1:4, 6)] <- round(OBD_patients_all_fix[i, c(1:4, 6)], 2)
}
# add gBOIN-ET results into the data.frame as a new column (notice the HE or LT case)
OBD_patients_all_fix[, 5] <- sapply(seq_along(OBD_true_fix_HE), function(i) gBOIN_ET_correct_OBD_pts_HE[i, OBD_true_fix_HE[i]])
# OBD_patients_all_fix[, 5] <- sapply(seq_along(OBD_true_fix_LT), function(i) gBOIN_ET_correct_OBD_pts_LT[i, OBD_true_fix_LT[i]])

# generate the barplots to visualize the average patients allocated to OBD 

# modify the data to visualize
OBD_patients_all_fix$scenario <- rownames(OBD_patients_all_fix)
OBD_patients_all_fix <- OBD_patients_all_fix %>%
  pivot_longer(cols = -scenario, names_to = "Methods", values_to = "OBD_patients") %>%
  mutate(Methods = factor(Methods, levels = c("Iso", "TEPI", "uTPI", "BOIN12", "gBOIN-ET", "Proposed")))
# generate the color vector
#library(RColorBrewer)
#colors_OBD_rate_barplots <- brewer.pal(n = 5, name = "Set1")
#colors_OBD_patients_barplots <- c("#13a983","#00bdcd","#cc340c","#006b7b","#f88421")
#colors_OBD_patients_barplots <- rev(c("#26355d", "#45639e", "#57b3ab", "#abcc87", "#d1c76e"))
# generate the barplot
B <- ggplot(OBD_patients_all_fix, aes(x = scenario, y = OBD_patients, fill = Methods)) +
  
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) + # for HE
  # geom_bar(stat = "identity", position = position_dodge(width = 0.8)) + # , color = "black", linewidth = 0.4
  
  labs(title = "Average Number of Patients at OBD", x = "Scenario", y = "No. of Patients at OBD") +
  scale_fill_manual(values = colors_barplots) +
  
  scale_y_continuous(limits = c(0,15), breaks = seq(0, 15, 3), expand = c(0,0)) + # for HE
  # scale_y_continuous(limits = c(0,20), breaks = seq(0, 20, 5), expand = c(0,0)) +  # for LT
  
  geom_text(aes(label = OBD_patients), position = position_dodge(width = 0.9), vjust = -0.5, size = 6) + # for HE
  # geom_text(aes(label = OBD_patients), position = position_dodge(width = 0.8), vjust = -0.5, size = 6) + 

  theme_minimal() + 
  theme(axis.text = element_text(family = "Helvetica", size = 30, colour = "black"),
        axis.title = element_text(family = "Helvetica", size = 30, colour = "black"),
        legend.text = element_text(family = "Helvetica", size = 30, colour = "black"),
        legend.title = element_blank(),
        #legend.position = "none",
        plot.title = element_text(family = "Helvetica", size = 30, colour = "black"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.4))

# -------------------------------------------------------------------------------- #






# --------------------- Poor Allocation rate visualization ----------------------- #

# generate the data.frame to store the  poor allocation rate of chosen fix scenarios
PA_rate_all_fix <- as.data.frame(matrix(nrow = length(tox_abb), ncol = 6)) # row as scenario 1-9, column as methods A to F
rownames(PA_rate_all_fix) <- seq(1:9)
colnames(PA_rate_all_fix) <- c("Iso", "TEPI", "uTPI", "BOIN12", "gBOIN-ET", "Proposed")

# load the poor allocation rate results
for (i in 1:length(tox_abb)) {
  PA_rate_all_fix[i,1] <- Iso_poor_allocation_rate[which(rownames(Iso_poor_allocation_rate) == paste0(tox_abb[i], "_", eff_abb[i])), 1]
  PA_rate_all_fix[i,2] <- TEPI_poor_allocation_rate[which(rownames(TEPI_poor_allocation_rate) == paste0(tox_abb[i], "_", eff_abb[i])), 1]
  PA_rate_all_fix[i,3] <- uTPI_poor_allocation_rate[which(rownames(uTPI_poor_allocation_rate) == paste0(tox_abb[i], "_", eff_abb[i])), 7]
  PA_rate_all_fix[i,4] <- BOIN12_poor_allocation_rate[which(rownames(BOIN12_poor_allocation_rate) == paste0(tox_abb[i], "_", eff_abb[i])), 7]
  
  PA_rate_all_fix[i,5] <- (gBOIN_ET_poorAllocation_fix_HE[[i]]$poor_allocation_rate) * 100
  # PA_rate_all_fix[i,5] <- (gBOIN_ET_poorAllocation_fix_LT[[i]]$poor_allocation_rate) * 100
  
  PA_rate_all_fix[i,6] <- Ours_poor_allocation_rate[which(rownames(Ours_poor_allocation_rate) == paste0(tox_abb[i], "_", eff_abb[i])), 1]
}


# generate the barplots to visualize the poor allocation rate

# modify the data to visualize
PA_rate_all_fix$scenario <- rownames(PA_rate_all_fix)
PA_rate_all_fix <- PA_rate_all_fix %>%
  pivot_longer(cols = -scenario, names_to = "Methods", values_to = "PA_rate") %>%
  mutate(Methods = factor(Methods, levels = c("Iso", "TEPI", "uTPI", "BOIN12", "gBOIN-ET", "Proposed")))
# generate the color vector
#library(RColorBrewer)
#colors_OBD_rate_barplots <- brewer.pal(n = 5, name = "Set1")
#colors_OBD_rate_barplots <- c("#13a983","#00bdcd","#cc340c","#006b7b","#f88421")
#colors_PA_rate_barplots <- rev(c("#26355d", "#45639e", "#57b3ab", "#abcc87", "#d1c76e"))
# generate the barplot
C <- ggplot(PA_rate_all_fix, aes(x = scenario, y = PA_rate, fill = Methods)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) + # , color = "black", linewidth = 0.4
  labs(title = "Poor Allocation Rate (%)", x = "Scenario", y = "Poor Allocation Rate (%)") +
  scale_fill_manual(values = colors_barplots) +
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 20), expand = c(0,0)) + 
  geom_text(aes(label = PA_rate), position = position_dodge(width = 0.8), vjust = -0.5, size = 6) + 
  theme_minimal() + 
  theme(axis.text = element_text(family = "Helvetica", size = 30, colour = "black"),
        axis.title = element_text(family = "Helvetica", size = 30, colour = "black"),
        legend.text = element_text(family = "Helvetica", size = 30, colour = "black"),
        legend.title = element_blank(),
        #legend.position = "none",
        plot.title = element_text(family = "Helvetica", size = 30, colour = "black"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.4))

# -------------------------------------------------------------------------------- #




# -------------------------- Overdose rate visualization ------------------------- #

# generate the data.frame to store the Overdose rate of chosen fix scenarios
OD_rate_all_fix <- as.data.frame(matrix(nrow = length(tox_abb), ncol = 6)) # row as scenario 1-9, column as methods A to F
rownames(OD_rate_all_fix) <- seq(1:9)
colnames(OD_rate_all_fix) <- c("Iso", "TEPI", "uTPI", "BOIN12", "gBOIN-ET", "Proposed")

# load the overdose rate results
for (i in 1:length(tox_abb)) {
  OD_rate_all_fix[i,1] <- Iso_overdose_rate[which(rownames(Iso_overdose_rate) == paste0(tox_abb[i], "_", eff_abb[i])), 1]
  OD_rate_all_fix[i,2] <- TEPI_overdose_rate[which(rownames(TEPI_overdose_rate) == paste0(tox_abb[i], "_", eff_abb[i])), 1]
  OD_rate_all_fix[i,3] <- uTPI_overdose_rate[which(rownames(uTPI_overdose_rate) == paste0(tox_abb[i], "_", eff_abb[i])), 6]
  OD_rate_all_fix[i,4] <- BOIN12_overdose_rate[which(rownames(BOIN12_overdose_rate) == paste0(tox_abb[i], "_", eff_abb[i])), 6]
  
  OD_rate_all_fix[i,6] <- Ours_overdose_rate[which(rownames(Ours_overdose_rate) == paste0(tox_abb[i], "_", eff_abb[i])), 1]
}
# add gBOIN-ET results (notice HE or LT)
OD_rate_all_fix[,5] <- mapply(function(row, threshold) { if (threshold < 6) { sum(gBOIN_ET_correct_OBD_pts_HE[row, (threshold+1):6]) } else { 0.0 }}, row = 1:9, threshold = MTD_true_fix_HE)
# OD_rate_all_fix[,5] <- mapply(function(row, threshold) { if (threshold < 6) { sum(gBOIN_ET_correct_OBD_pts_LT[row, (threshold+1):6]) } else { 0.0 }}, row = 1:9, threshold = MTD_true_fix_LT)


# generate the barplots to visualize the OBD correct selection rate

# modify the data to visualize
OD_rate_all_fix$scenario <- rownames(OD_rate_all_fix)
OD_rate_all_fix <- OD_rate_all_fix %>%
  pivot_longer(cols = -scenario, names_to = "Methods", values_to = "OD_rate") %>%
  mutate(Methods = factor(Methods, levels = c("Iso", "TEPI", "uTPI", "BOIN12", "gBOIN-ET", "Proposed")))
# generate the color vector
#library(RColorBrewer)
#colors_OBD_rate_barplots <- brewer.pal(n = 5, name = "Set1")
#colors_OBD_rate_barplots <- c("#13a983","#00bdcd","#cc340c","#006b7b","#f88421")
#colors_OD_rate_barplots <- rev(c("#26355d", "#45639e", "#57b3ab", "#abcc87", "#d1c76e"))
# generate the barplot
D <- ggplot(OD_rate_all_fix, aes(x = scenario, y = OD_rate, fill = Methods)) +
  
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +   # for HE
  # geom_bar(stat = "identity", position = position_dodge(width = 0.8)) + # , color = "black", linewidth = 0.4
  
  labs(title = "Average Number of Overdose Patients", x = "Scenario", y = "No. of Overdose Patients") +
  scale_fill_manual(values = colors_barplots) +
  
  scale_y_continuous(limits = c(0,30), breaks = seq(0, 30, 5), expand = c(0,0)) + # for HE
  # scale_y_continuous(limits = c(0,20), breaks = seq(0, 20, 5), expand = c(0,0)) + # for LT
  
  geom_text(aes(label = OD_rate), position = position_dodge(width = 0.9), vjust = -0.5, size = 6) + # for HE
  # geom_text(aes(label = OD_rate), position = position_dodge(width = 0.8), vjust = -0.5, size = 6) + 
  
  theme_minimal() + 
  theme(axis.text = element_text(family = "Helvetica", size = 30, colour = "black"),
        axis.title = element_text(family = "Helvetica", size = 30, colour = "black"),
        legend.text = element_text(family = "Helvetica", size = 30, colour = "black"),
        legend.title = element_blank(),
        #legend.position = "none",
        plot.title = element_text(family = "Helvetica", size = 30, colour = "black"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.4))

# --------------------------------- END --------------------------------- #





# -------------------- merge all the barplots together ------------------------ #
# ggarrange(A, B, 
#           labels = c("A", "B"), ncol = 1, nrow = 2, 
#           font.label = list(size = 24, family = "Helvetica", color = "black"), 
#           common.legend = TRUE, legend = "top")
# 
# ggarrange(C, D, 
#           labels = c("C", "D"), ncol = 1, nrow = 2, 
#           font.label = list(size = 24, family = "Helvetica", color = "black"), 
#           common.legend = TRUE, legend = "top")
{
  A <- A + guides(fill = guide_legend(nrow = 1))
  B <- B + guides(fill = guide_legend(nrow = 1))
  C <- C + guides(fill = guide_legend(nrow = 1))
  D <- D + guides(fill = guide_legend(nrow = 1))
}



ggarrange(A, B, C, D, 
          labels = c("A", "B", "C", "D"), ncol = 1, nrow = 4, 
          font.label = list(size = 30, family = "Helvetica", color = "black"), 
          common.legend = TRUE, legend = "top")

# end here










