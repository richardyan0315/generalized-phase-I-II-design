# The Result Figures are generated in this file:

# Fix Scenarios: Figures 5 to 8
# Random Scenarios: Figures 9 and S.2
# Sensitivity Analysis: Figures 10 and 11

rm(list=ls())
source("./simulationResultsCleaning.R")

{
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(tidyverse)
  library(gridExtra)
  library(cowplot)
  library(ggpubr)
  library(tidyr)
}


colors_barplots <- rev(c("orange", "#45639e", "#57b3ab", "#abcc87", "#d1c76e", "#7e91b1"))


# ---------------- Fix Scenarios (HE, w_e = 0.7, w_t = 0.3 / LT, w_e = 0.45, w_t = 0.55) --------------------- #
{
  fix_tox <- data.frame()
  fix_eff <- data.frame()
  fix_u <- data.frame()
  
  scenarios_prob_tox <- summary_tables$toxicity_prob
  scenarios_prob_eff <- summary_tables$efficacy_prob
  true_utilities <- round(summary_tables$utilities, 2)
  
  rownames(true_utilities) <- paste0(
    expand.grid(rownames(scenarios_prob_tox), rownames(scenarios_prob_eff))[,2],
    "_",
    expand.grid(rownames(scenarios_prob_tox), rownames(scenarios_prob_eff))[,1]
  )
  
  for (i in 1:length(tox_abb)) {
    fix_tox <- rbind(fix_tox, scenarios_prob_tox[rownames(scenarios_prob_tox) == tox_abb[i],]) # true tox prob (unlist when use)
    fix_eff <- rbind(fix_eff, scenarios_prob_eff[rownames(scenarios_prob_eff) == eff_abb[i],]) # true eff prob
    fix_u <- rbind(fix_u, true_utilities[rownames(true_utilities) == paste0(tox_abb[i], "_", eff_abb[i]),]) # true utility
  }
  
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


# Figure 5 and 7: generate the probability and utility curve sub-plots of fix scenarios 
for (i in 1:length(tox_abb)) {
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

# Figure 5 and 7
ggarrange(all_plot_list_fix[[1]], all_plot_list_fix[[2]], all_plot_list_fix[[3]], 
          all_plot_list_fix[[4]], all_plot_list_fix[[5]], all_plot_list_fix[[6]], 
          all_plot_list_fix[[7]], all_plot_list_fix[[8]], all_plot_list_fix[[9]], 
          common.legend = TRUE, legend = "bottom", nrow = 3, ncol = 3)


# -------------------------------------------------------------------------------- #


# Figure 6 and 8, A to D
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

B <- ggplot(OBD_patients_all_fix, aes(x = scenario, y = OBD_patients, fill = Methods)) +
  
  geom_bar(stat = "identity", position = position_dodge(width = ifelse(INX == "HE_", 0.9, 0.8))) + 

  labs(title = "Average Number of Patients at OBD", x = "Scenario", y = "No. of Patients at OBD") +
  scale_fill_manual(values = colors_barplots) +

  {if(INX == "HE_") 
    scale_y_continuous(limits = c(0,15),
                       breaks = seq(0, 15, 3),
                       expand = c(0,0))
    else 
      scale_y_continuous(limits = c(0,20),
                         breaks = seq(0, 20, 5),
                         expand = c(0,0))
  } + 
  
  {if(INX == "HE_")
    geom_text(aes(label = OBD_patients), 
              position = position_dodge(width = 0.9), 
              vjust = -0.5, 
              size = 6)
    else 
      geom_text(aes(label = OBD_patients), 
                position = position_dodge(width = 0.8), 
                vjust = -0.5, 
                size = 6)
  } + 
  
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

D <- ggplot(OD_rate_all_fix, aes(x = scenario, y = OD_rate, fill = Methods)) +
  
  {if(INX == "HE_")
    geom_bar(stat = "identity", position = position_dodge(width = 0.9))
    else 
      geom_bar(stat = "identity", position = position_dodge(width = 0.8))
  } + 

  labs(title = "Average Number of Overdose Patients", x = "Scenario", y = "No. of Overdose Patients") +
  scale_fill_manual(values = colors_barplots) +
  
  {if(INX == "HE_")
    scale_y_continuous(limits = c(0,30), breaks = seq(0, 30, 5), expand = c(0,0))
    else 
      scale_y_continuous(limits = c(0,20), breaks = seq(0, 20, 5), expand = c(0,0))
  } + 
  
  {if(INX == "HE_")
    geom_text(aes(label = OD_rate), position = position_dodge(width = 0.9), vjust = -0.5, size = 6) 
    else 
      geom_text(aes(label = OD_rate), position = position_dodge(width = 0.8), vjust = -0.5, size = 6)
  } + 
  
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

# Figure 6 and 8
ggarrange(A, B, C, D, 
          labels = c("A", "B", "C", "D"), ncol = 1, nrow = 4, 
          font.label = list(size = 30, family = "Helvetica", color = "black"), 
          common.legend = TRUE, legend = "top")

# ----------------------------- Fix Scenarios END -------------------------------- #



# ---------------- Random Scenarios (HE, w_e = 0.7, w_t = 0.3 / LT, w_e = 0.45, w_t = 0.55) --------------------- #

RS_OBD_rate_mean9$scenario <- rownames(RS_OBD_rate_mean9)
RS_OBD_rate_mean9 <- RS_OBD_rate_mean9 %>%
  pivot_longer(cols = -scenario, names_to = "Methods", values_to = "OBD_rate") %>%
  mutate(Methods = factor(Methods, levels = c("Iso", "TEPI", "uTPI", "BOIN12", "gBOIN-ET", "Proposed"))) %>%
  mutate(scenario = factor(scenario, levels = c("LL", "LM", "LH", "ML", "MM", "MH", "HL", "HM", "HH", "Overall")))

RS_1 <- ggplot(RS_OBD_rate_mean9, aes(x = scenario, y = OBD_rate, fill = Methods)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) + 
  labs(title = "Selection of OBD (%)", x = "Random Scenario Subset", y = "Selection of OBD (%)") +
  scale_fill_manual(values = colors_barplots) +
  
  {if(INX == "HE_")
    scale_y_continuous(limits = c(0,60), breaks = seq(0, 60, 10), expand = c(0,0))
    else 
      scale_y_continuous(limits = c(0,80), breaks = seq(0, 80, 20), expand = c(0,0))
  } + 
  
  {if(INX == "HE_")
    geom_text(aes(label = OBD_rate), position = position_dodge(width = 0.9), vjust = -0.5, size = 6) 
    else 
      geom_text(aes(label = OBD_rate), position = position_dodge(width = 0.8), vjust = -0.5, size = 6) 
  } + 

  theme_minimal() + 
  theme(axis.text = element_text(family = "Helvetica", size = 30, colour = "black"),
        axis.title = element_text(family = "Helvetica", size = 30, colour = "black"),
        legend.text = element_text(family = "Helvetica", size = 30, colour = "black"),
        legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_text(family = "Helvetica", size = 30, colour = "black"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.4))



RS_OBD_pts_mean9$scenario <- rownames(RS_OBD_pts_mean9)
RS_OBD_pts_mean9 <- RS_OBD_pts_mean9 %>%
  pivot_longer(cols = -scenario, names_to = "Methods", values_to = "OBD_pts") %>%
  mutate(Methods = factor(Methods, levels = c("Iso", "TEPI", "uTPI", "BOIN12", "gBOIN-ET", "Proposed"))) %>%
  mutate(scenario = factor(scenario, levels = c("LL", "LM", "LH", "ML", "MM", "MH", "HL", "HM", "HH", "Overall")))

RS_2 <- ggplot(RS_OBD_pts_mean9, aes(x = scenario, y = OBD_pts, fill = Methods)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) + 
  labs(title = "Average Number of Patients at OBD", x = "Random Scenario Subset", y = "No. of Patients at OBD") +
  scale_fill_manual(values = colors_barplots) +
  
  {if(INX == "HE_")
    scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, 2), expand = c(0,0))  
    else 
      scale_y_continuous(limits = c(0,15), breaks = seq(0, 15, 3), expand = c(0,0)) 
  } + 
  
  {if(INX == "HE_")
    geom_text(aes(label = OBD_pts), position = position_dodge(width = 0.9), vjust = -0.5, size = 6)
    else 
      geom_text(aes(label = OBD_pts), position = position_dodge(width = 0.8), vjust = -0.5, size = 6)
  } + 
  
  theme_minimal() + 
  theme(axis.text = element_text(family = "Helvetica", size = 30, colour = "black"),
        axis.title = element_text(family = "Helvetica", size = 30, colour = "black"),
        legend.text = element_text(family = "Helvetica", size = 30, colour = "black"),
        legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_text(family = "Helvetica", size = 30, colour = "black"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.4))



RS_PA_rate_mean9$scenario <- rownames(RS_PA_rate_mean9)
RS_PA_rate_mean9 <- RS_PA_rate_mean9 %>%
  pivot_longer(cols = -scenario, names_to = "Methods", values_to = "PA_rate") %>%
  mutate(Methods = factor(Methods, levels = c("Iso", "TEPI", "uTPI", "BOIN12", "gBOIN-ET", "Proposed"))) %>%
  mutate(scenario = factor(scenario, levels = c("LL", "LM", "LH", "ML", "MM", "MH", "HL", "HM", "HH", "Overall")))

RS_3 <- ggplot(RS_PA_rate_mean9, aes(x = scenario, y = PA_rate, fill = Methods)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) + 
  labs(title = "Poor Allocation Rate (%)", x = "Random Scenario Subset", y = "Poor Allocation Rate (%)") +
  scale_fill_manual(values = colors_barplots) +
  
  {if(INX == "HE_")
    scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 20), expand = c(0,0))
    else 
      scale_y_continuous(limits = c(0,80), breaks = seq(0, 80, 20), expand = c(0,0))
  } + 
  
  geom_text(aes(label = PA_rate), position = position_dodge(width = 0.9), vjust = -0.5, size = 6) + 
  theme_minimal() + 
  theme(axis.text = element_text(family = "Helvetica", size = 30, colour = "black"),
        axis.title = element_text(family = "Helvetica", size = 30, colour = "black"),
        legend.text = element_text(family = "Helvetica", size = 30, colour = "black"),
        legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_text(family = "Helvetica", size = 30, colour = "black"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.4))


RS_OD_pts_mean9$scenario <- rownames(RS_OD_pts_mean9)
RS_OD_pts_mean9 <- RS_OD_pts_mean9 %>%
  pivot_longer(cols = -scenario, names_to = "Methods", values_to = "OD_pts") %>%
  mutate(Methods = factor(Methods, levels = c("Iso", "TEPI", "uTPI", "BOIN12", "gBOIN-ET", "Proposed"))) %>%
  mutate(scenario = factor(scenario, levels = c("LL", "LM", "LH", "ML", "MM", "MH", "HL", "HM", "HH", "Overall")))

RS_4 <- ggplot(RS_OD_pts_mean9, aes(x = scenario, y = OD_pts, fill = Methods)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) + 
  labs(title = "Average Number of Overdose Patients", x = "Random Scenario Subset", y = "No. of Overdose Patients") +
  scale_fill_manual(values = colors_barplots) +
  
  {if(INX == "HE_")
    scale_y_continuous(limits = c(0,15), breaks = seq(0, 15, 3), expand = c(0,0)) 
    else 
      scale_y_continuous(limits = c(0,12), breaks = seq(0, 12, 3), expand = c(0,0))
  } + 
  
  geom_text(aes(label = OD_pts), position = position_dodge(width = 0.9), vjust = -0.5, size = 6) + 
  theme_minimal() + 
  theme(axis.text = element_text(family = "Helvetica", size = 30, colour = "black"),
        axis.title = element_text(family = "Helvetica", size = 30, colour = "black"),
        legend.text = element_text(family = "Helvetica", size = 30, colour = "black"),
        legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_text(family = "Helvetica", size = 30, colour = "black"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.4))


RS_1 <- RS_1 + guides(fill = guide_legend(nrow = 1))
RS_2 <- RS_2 + guides(fill = guide_legend(nrow = 1))
RS_3 <- RS_3 + guides(fill = guide_legend(nrow = 1))
RS_4 <- RS_4 + guides(fill = guide_legend(nrow = 1))

# Figure 9 and S.2
ggarrange(RS_1, RS_2, RS_3, RS_4,  
          labels = c("A", "B", "C", "D"), ncol = 1, nrow = 4, 
          font.label = list(size = 30, family = "Helvetica", color = "black"), 
          common.legend = TRUE, legend = "top")

# --------------------------- Random Scenarios END -------------------------------- #


# --------------------------- Sensitivity Analysis -------------------------------- #
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

# --------------------------------------------------------------------- #
# Figures 10 and 11
ggarrange(AA1, AA2,
          ncol = 2, nrow = 1, 
          font.label = list(size = 20, family = "Helvetica", color = "black"),
          common.legend = TRUE, legend = "top")

ggarrange(BB1, BB2,
          ncol = 2, nrow = 1, 
          font.label = list(size = 20, family = "Helvetica", color = "black"),
          common.legend = TRUE, legend = "top")














