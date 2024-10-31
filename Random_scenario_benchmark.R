rm(list=ls())

# ------------- load the packages ------------ #
{
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(tidyverse)
  library(gridExtra)
  library(cowplot)
  library(ggpubr)
}


# ---------- File names for all 5 methods are "HE_" (based on w_e) ----------- #

# to get the true OBD of each scenario from "OBD_true" data.frame
source("./R files/scenario_generator_local.R")
source("./R files/OBD_selection_visualization_local.R")


# check the w_e & width_SA
print(paste0("The current w_e is ", w_e))
print(paste0("The current interval width is ", width_SA))

# file name distinguishing based on w_e
if(w_e > 0.5){INX <- "HE_"} else {INX <- "LT_"}

MTD_true_vector <- MTD_true_current 
OBD_true_vector <- as.vector(t(OBD_true))


# sub-group definition
prefix_scenario_3 <- rownames(scenarios_prob_tox)

{
  # Three types of tox and eff, each with three different transition speeds. 
  # For tox and eff, there are a total of 27 combinations for each. 
  # These 27 different transition rate combinations are categorized by 'dominant rate' into three major groups: L, M, and H, containing 7, 10, and 10 combinations respectively. 
  # Thus, 729 specific random scenarios are divided among combinations of nine major categories (LL, LM, ...).
  
  L_count <- str_count(prefix_scenario_3, "L")
  M_count <- str_count(prefix_scenario_3, "M")
  H_count <- str_count(prefix_scenario_3, "H")
  
  H_group <- prefix_scenario_3[(H_count >= 1 & L_count == 0) | H_count >= 2] # 10
  L_group <- prefix_scenario_3[L_count >= 2] # 7
  M_group <- setdiff(setdiff(prefix_scenario_3, H_group), L_group) # 10
  
  L_group_tox <- prefix_scenario[substr(prefix_scenario[,1], 1, 3) %in% L_group, 1]
  M_group_tox <- prefix_scenario[substr(prefix_scenario[,1], 1, 3) %in% M_group, 1]
  H_group_tox <- prefix_scenario[substr(prefix_scenario[,1], 1, 3) %in% H_group, 1]
  
  L_group_eff <- prefix_scenario[substr(prefix_scenario[,1], nchar(prefix_scenario[,1]) - 2, nchar(prefix_scenario[,1])) %in% L_group, 1]
  M_group_eff <- prefix_scenario[substr(prefix_scenario[,1], nchar(prefix_scenario[,1]) - 2, nchar(prefix_scenario[,1])) %in% M_group, 1]
  H_group_eff <- prefix_scenario[substr(prefix_scenario[,1], nchar(prefix_scenario[,1]) - 2, nchar(prefix_scenario[,1])) %in% H_group, 1]
  
  LL_group <- intersect(L_group_tox, L_group_eff) # 49
  LM_group <- intersect(L_group_tox, M_group_eff) # 70
  LH_group <- intersect(L_group_tox, H_group_eff) # 70
  ML_group <- intersect(M_group_tox, L_group_eff) # 70
  MM_group <- intersect(M_group_tox, M_group_eff) # 100
  MH_group <- intersect(M_group_tox, H_group_eff) # 100
  HL_group <- intersect(H_group_tox, L_group_eff) # 70
  HM_group <- intersect(H_group_tox, M_group_eff) # 100
  HH_group <- intersect(H_group_tox, H_group_eff) # 100
}



# --------- Load the data --------- #

# The other method results
{
  RS_OBD_rate <- combined_data_OBD_rate
  RS_OBD_pts <- combined_data_OBD_patient
  RS_PA_rate <- combined_data_PA_rate
  RS_OD_pts <- combined_data_OD_rate
  
  colnames(RS_OBD_rate) <- c("scenarios", "OBD_rate", "design")
  colnames(RS_OD_pts) <- c("scenarios", "OD_pts", "design")
  
  RS_OBD_rate[,1] <- RS_OBD_pts[,1] <- RS_PA_rate[,1] <- RS_OD_pts[,1] <- prefix_scenario[,1]
  
  prefix_index <- c("LL", "LM", "LH", "ML", "MM", "MH", "HL", "HM", "HH")
  
  design_index <- c("Iso", "TEPI", "uTPI", "BOIN12", "gBOIN-ET", "Proposed")
  # design_index <- c("Iso", "TEPI", "uTPI", "BOIN12", "Proposed")
  
  # Generate variable names with changing suffixes and assign values in a loop
  base_name_1 <- "RS_OBD_rate_"
  base_name_2 <- "RS_OBD_pts_"
  base_name_3 <- "RS_PA_rate_"
  base_name_4 <- "RS_OD_pts_"
  
  RS_OBD_rate_mean9 <- data.frame(matrix(nrow = length(prefix_index) + 1, ncol = 6))
  RS_OBD_pts_mean9 <- data.frame(matrix(nrow = length(prefix_index) + 1, ncol = 6))
  RS_PA_rate_mean9 <- data.frame(matrix(nrow = length(prefix_index) + 1, ncol = 6))
  RS_OD_pts_mean9 <- data.frame(matrix(nrow = length(prefix_index) + 1, ncol = 6))
  
  rownames(RS_OBD_rate_mean9) <- rownames(RS_OBD_pts_mean9) <- rownames(RS_PA_rate_mean9) <- rownames(RS_OD_pts_mean9) <- c(prefix_index, "Overall")
  colnames(RS_OBD_rate_mean9) <- colnames(RS_OBD_pts_mean9) <- colnames(RS_PA_rate_mean9) <- colnames(RS_OD_pts_mean9) <- design_index
  
  
  for (i in 1:length(prefix_index)) {
    suffix <- prefix_index[i]
    
    dynamic_variable_1 <- paste0(base_name_1, suffix)
    dynamic_variable_2 <- paste0(base_name_2, suffix)
    dynamic_variable_3 <- paste0(base_name_3, suffix)
    dynamic_variable_4 <- paste0(base_name_4, suffix)
    
    assign(dynamic_variable_1, RS_OBD_rate[RS_OBD_rate$scenarios %in% get(paste0(suffix, "_group")),])
    assign(dynamic_variable_2, RS_OBD_pts[RS_OBD_pts$scenarios %in% get(paste0(suffix, "_group")),])
    assign(dynamic_variable_3, RS_PA_rate[RS_PA_rate$scenarios %in% get(paste0(suffix, "_group")),])
    assign(dynamic_variable_4, RS_OD_pts[RS_OD_pts$scenarios %in% get(paste0(suffix, "_group")),])
    
    for (j in 1:6) {
      RS_OBD_rate_mean9[i,j] <- round(mean(subset(get(dynamic_variable_1), get(dynamic_variable_1)[,3] == design_index[j])[ ,2]), 2)
      RS_OBD_pts_mean9[i,j] <- round(mean(subset(get(dynamic_variable_2), get(dynamic_variable_2)[,3] == design_index[j])[ ,2]), 2)
      RS_PA_rate_mean9[i,j] <- round(mean(subset(get(dynamic_variable_3), get(dynamic_variable_3)[,3] == design_index[j])[ ,2]), 2)
      RS_OD_pts_mean9[i,j] <- round(mean(subset(get(dynamic_variable_4), get(dynamic_variable_4)[,3] == design_index[j])[ ,2]), 2)
    }
    
  }
  
  RS_OBD_rate_mean9[10, ] <- round(colMeans(RS_OBD_rate_mean9[-10,]), 2)
  RS_OBD_pts_mean9[10, ] <- round(colMeans(RS_OBD_pts_mean9[-10,]), 2)
  RS_PA_rate_mean9[10, ] <- round(colMeans(RS_PA_rate_mean9[-10,]), 2)
  RS_OD_pts_mean9[10, ] <- round(colMeans(RS_OD_pts_mean9[-10,]), 2)
  
}

  
# ---------------- Visualization --------------- #
{
  colors_barplots <- rev(c("orange", "#45639e", "#57b3ab", "#abcc87", "#d1c76e", "#7e91b1"))
  
  
  RS_OBD_rate_mean9$scenario <- rownames(RS_OBD_rate_mean9)
  RS_OBD_rate_mean9 <- RS_OBD_rate_mean9 %>%
    pivot_longer(cols = -scenario, names_to = "Methods", values_to = "OBD_rate") %>%
    mutate(Methods = factor(Methods, levels = c("Iso", "TEPI", "uTPI", "BOIN12", "gBOIN-ET", "Proposed"))) %>%
    mutate(scenario = factor(scenario, levels = c("LL", "LM", "LH", "ML", "MM", "MH", "HL", "HM", "HH", "Overall")))
  
  RS_1 <- ggplot(RS_OBD_rate_mean9, aes(x = scenario, y = OBD_rate, fill = Methods)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) + # , color = "black", linewidth = 0.4
    labs(title = "Selection of OBD (%)", x = "Random Scenario Subset", y = "Selection of OBD (%)") +
    scale_fill_manual(values = colors_barplots) +
    
    # scale_y_continuous(limits = c(0,60), breaks = seq(0, 60, 10), expand = c(0,0)) + # for HE
    scale_y_continuous(limits = c(0,80), breaks = seq(0, 80, 20), expand = c(0,0)) + # for LT
    
    geom_text(aes(label = OBD_rate), position = position_dodge(width = 0.9), vjust = -0.5, size = 6) + # for HE
    # geom_text(aes(label = OBD_rate), position = position_dodge(width = 0.8), vjust = -0.5, size = 6) + 
    
    theme_minimal() + 
    theme(axis.text = element_text(family = "Helvetica", size = 30, colour = "black"),
          axis.title = element_text(family = "Helvetica", size = 30, colour = "black"),
          legend.text = element_text(family = "Helvetica", size = 30, colour = "black"),
          legend.title = element_blank(),
          #legend.position = "none",
          plot.title = element_text(family = "Helvetica", size = 30, colour = "black"),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black", linewidth = 0.4))
  
  
  
  RS_OBD_pts_mean9$scenario <- rownames(RS_OBD_pts_mean9)
  RS_OBD_pts_mean9 <- RS_OBD_pts_mean9 %>%
    pivot_longer(cols = -scenario, names_to = "Methods", values_to = "OBD_pts") %>%
    mutate(Methods = factor(Methods, levels = c("Iso", "TEPI", "uTPI", "BOIN12", "gBOIN-ET", "Proposed"))) %>%
    mutate(scenario = factor(scenario, levels = c("LL", "LM", "LH", "ML", "MM", "MH", "HL", "HM", "HH", "Overall")))
  
  RS_2 <- ggplot(RS_OBD_pts_mean9, aes(x = scenario, y = OBD_pts, fill = Methods)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) + # , color = "black", linewidth = 0.4
    labs(title = "Average Number of Patients at OBD", x = "Random Scenario Subset", y = "No. of Patients at OBD") +
    scale_fill_manual(values = colors_barplots) +
    
    # scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, 2), expand = c(0,0)) + # for HE
    scale_y_continuous(limits = c(0,15), breaks = seq(0, 15, 3), expand = c(0,0)) + # for LT
    
    geom_text(aes(label = OBD_pts), position = position_dodge(width = 0.9), vjust = -0.5, size = 6) + # for HE
    # geom_text(aes(label = OBD_pts), position = position_dodge(width = 0.8), vjust = -0.5, size = 6) + 
    
    theme_minimal() + 
    theme(axis.text = element_text(family = "Helvetica", size = 30, colour = "black"),
          axis.title = element_text(family = "Helvetica", size = 30, colour = "black"),
          legend.text = element_text(family = "Helvetica", size = 30, colour = "black"),
          legend.title = element_blank(),
          #legend.position = "none",
          plot.title = element_text(family = "Helvetica", size = 30, colour = "black"),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black", linewidth = 0.4))
  
  
  
  RS_PA_rate_mean9$scenario <- rownames(RS_PA_rate_mean9)
  RS_PA_rate_mean9 <- RS_PA_rate_mean9 %>%
    pivot_longer(cols = -scenario, names_to = "Methods", values_to = "PA_rate") %>%
    mutate(Methods = factor(Methods, levels = c("Iso", "TEPI", "uTPI", "BOIN12", "gBOIN-ET", "Proposed"))) %>%
    mutate(scenario = factor(scenario, levels = c("LL", "LM", "LH", "ML", "MM", "MH", "HL", "HM", "HH", "Overall")))
  
  RS_3 <- ggplot(RS_PA_rate_mean9, aes(x = scenario, y = PA_rate, fill = Methods)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) + # , color = "black", linewidth = 0.4
    labs(title = "Poor Allocation Rate (%)", x = "Random Scenario Subset", y = "Poor Allocation Rate (%)") +
    scale_fill_manual(values = colors_barplots) +
    
    scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 20), expand = c(0,0)) + # for HE
    # scale_y_continuous(limits = c(0,80), breaks = seq(0, 80, 20), expand = c(0,0)) + 
    
    geom_text(aes(label = PA_rate), position = position_dodge(width = 0.9), vjust = -0.5, size = 6) + 
    theme_minimal() + 
    theme(axis.text = element_text(family = "Helvetica", size = 30, colour = "black"),
          axis.title = element_text(family = "Helvetica", size = 30, colour = "black"),
          legend.text = element_text(family = "Helvetica", size = 30, colour = "black"),
          legend.title = element_blank(),
          #legend.position = "none",
          plot.title = element_text(family = "Helvetica", size = 30, colour = "black"),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black", linewidth = 0.4))
  
  
  RS_OD_pts_mean9$scenario <- rownames(RS_OD_pts_mean9)
  RS_OD_pts_mean9 <- RS_OD_pts_mean9 %>%
    pivot_longer(cols = -scenario, names_to = "Methods", values_to = "OD_pts") %>%
    mutate(Methods = factor(Methods, levels = c("Iso", "TEPI", "uTPI", "BOIN12", "gBOIN-ET", "Proposed"))) %>%
    mutate(scenario = factor(scenario, levels = c("LL", "LM", "LH", "ML", "MM", "MH", "HL", "HM", "HH", "Overall")))
  
  RS_4 <- ggplot(RS_OD_pts_mean9, aes(x = scenario, y = OD_pts, fill = Methods)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) + # , color = "black", linewidth = 0.4
    labs(title = "Average Number of Overdose Patients", x = "Random Scenario Subset", y = "No. of Overdose Patients") +
    scale_fill_manual(values = colors_barplots) +
    
    scale_y_continuous(limits = c(0,15), breaks = seq(0, 15, 3), expand = c(0,0)) + # for HE
    # scale_y_continuous(limits = c(0,12), breaks = seq(0, 12, 3), expand = c(0,0)) + 
    
    geom_text(aes(label = OD_pts), position = position_dodge(width = 0.9), vjust = -0.5, size = 6) + 
    theme_minimal() + 
    theme(axis.text = element_text(family = "Helvetica", size = 30, colour = "black"),
          axis.title = element_text(family = "Helvetica", size = 30, colour = "black"),
          legend.text = element_text(family = "Helvetica", size = 30, colour = "black"),
          legend.title = element_blank(),
          #legend.position = "none",
          plot.title = element_text(family = "Helvetica", size = 30, colour = "black"),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black", linewidth = 0.4))
  
  
  RS_1 <- RS_1 + guides(fill = guide_legend(nrow = 1))
  RS_2 <- RS_2 + guides(fill = guide_legend(nrow = 1))
  RS_3 <- RS_3 + guides(fill = guide_legend(nrow = 1))
  RS_4 <- RS_4 + guides(fill = guide_legend(nrow = 1))
  
}


# ----------------- merging into one plot --------------------- #


ggarrange(RS_1, RS_2, RS_3, RS_4,  
          labels = c("A", "B", "C", "D"), ncol = 1, nrow = 4, 
          font.label = list(size = 30, family = "Helvetica", color = "black"), 
          common.legend = TRUE, legend = "top")


# ------------------------- END --------------------------- #





