rm(list=ls())
source("./scenarioGenerator.R")

# check the w_e
print(paste0("The current w_e is ", w_e))
# file name distinguishing based on w_e
if(w_e > 0.5){INX <- "HE_"} else {INX <- "LT_"}

#---- simulation results loading
# 0.a true OBD
ifelse(INX == "HE_", 
       OBD_true <- readRDS("./Outputs/OBD_true_random_HE.rds"),
       OBD_true <- readRDS("./Outputs/OBD_true_random_LT.rds"))

MTD_true_vector <- MTD_true_current <- trial_results$MTD_all
OBD_true_vector <- as.vector(t(OBD_true)) # transfer from matrix to vector

# 0.b True fix scenario MTD
MTD_true_fix_HE <- c(3, 3, 3, 2, 6, 6, 2, 2, 4)
MTD_true_fix_LT <- c(3, 3, 3, 2, 6, 6, 2, 2, 4)

# 0.c True fix scenario OBD
OBD_true_fix_HE <- c(3, 3, 3, 2, 3, 6, 2, 1, 4)
OBD_true_fix_LT <- c(1, 1, 2, 1, 2, 6, 1, 1, 1)

# 0.d gBOIN-ET results loading and wrangling
gBOIN_ET_results_fix_HE <- readRDS("./Outputs/Benchmark oc/gBOIN-ET/gBOIN-ET_results_fix_HE.rds") # scenarios 1 to 9
gBOIN_ET_results_fix_LT <- readRDS("./Outputs/Benchmark oc/gBOIN-ET/gBOIN-ET_results_fix_LT.rds") # scenarios 1 to 9
gBOIN_ET_results_all_HE <- readRDS("./Outputs/Benchmark oc/gBOIN-ET/gBOIN-ET_results_all_HE.rds") # all 729 random scenarios
gBOIN_ET_results_all_LT <- readRDS("./Outputs/Benchmark oc/gBOIN-ET/gBOIN-ET_results_all_LT.rds") # all 729 random scenarios

gBOIN_ET_poorAllocation_fix_HE <- readRDS("./Outputs/Benchmark oc/gBOIN-ET/gBOIN-ET_poorAllocation_fix_HE.rds") # scenarios 1 to 9
gBOIN_ET_poorAllocation_fix_LT <- readRDS("./Outputs/Benchmark oc/gBOIN-ET/gBOIN-ET_poorAllocation_fix_LT.rds") # scenarios 1 to 9
gBOIN_ET_poorAllocation_all_HE <- readRDS("./Outputs/Benchmark oc/gBOIN-ET/gBOIN-ET_poorAllocation_all_HE.rds") # all 729 random scenarios
gBOIN_ET_poorAllocation_all_LT <- readRDS("./Outputs/Benchmark oc/gBOIN-ET/gBOIN-ET_poorAllocation_all_LT.rds") # all 729 random scenarios

gBOIN_ET_correct_selected_OBD_rate_HE <- gBOIN_ET_correct_selected_OBD_rate_LT <- data.frame()
for (i in 1:9) {
  gBOIN_ET_correct_selected_OBD_rate_HE <- rbind(gBOIN_ET_correct_selected_OBD_rate_HE, gBOIN_ET_results_fix_HE[[i]]$prop.select)
  gBOIN_ET_correct_selected_OBD_rate_LT <- rbind(gBOIN_ET_correct_selected_OBD_rate_LT, gBOIN_ET_results_fix_LT[[i]]$prop.select)
}
colnames(gBOIN_ET_correct_selected_OBD_rate_HE) <- colnames(gBOIN_ET_correct_selected_OBD_rate_LT) <- 
  c("Dose 1", "Dose 2", "Dose 3", "Dose 4", "Dose 5", "Dose 6")

gBOIN_ET_correct_OBD_pts_HE <- gBOIN_ET_correct_OBD_pts_LT <- data.frame()
for (i in 1:9) {
  gBOIN_ET_correct_OBD_pts_HE <- rbind(gBOIN_ET_correct_OBD_pts_HE, gBOIN_ET_results_fix_HE[[i]]$n.patient)
  gBOIN_ET_correct_OBD_pts_LT <- rbind(gBOIN_ET_correct_OBD_pts_LT, gBOIN_ET_results_fix_LT[[i]]$n.patient)
}
colnames(gBOIN_ET_correct_OBD_pts_HE) <- colnames(gBOIN_ET_correct_OBD_pts_LT) <- 
  c("Dose 1", "Dose 2", "Dose 3", "Dose 4", "Dose 5", "Dose 6")

gBOIN_ET_npts_all_HE <- data.frame()
gBOIN_ET_prop_csel_all_HE <- data.frame()
gBOIN_ET_poorAlloc_all_HE <- c()

gBOIN_ET_npts_all_LT <- data.frame()
gBOIN_ET_prop_csel_all_LT <- data.frame()
gBOIN_ET_poorAlloc_all_LT <- c()

for (i in 1:729) {
  gBOIN_ET_npts_all_HE <- rbind(gBOIN_ET_npts_all_HE, gBOIN_ET_results_all_HE[[i]]$n.patient)
  gBOIN_ET_prop_csel_all_HE <- rbind(gBOIN_ET_prop_csel_all_HE, gBOIN_ET_results_all_HE[[i]]$prop.select)
  gBOIN_ET_poorAlloc_all_HE <- c(gBOIN_ET_poorAlloc_all_HE, gBOIN_ET_poorAllocation_all_HE[[i]]$poor_allocation_rate)
  
  gBOIN_ET_npts_all_LT <- rbind(gBOIN_ET_npts_all_LT, gBOIN_ET_results_all_LT[[i]]$n.patient)
  gBOIN_ET_prop_csel_all_LT <- rbind(gBOIN_ET_prop_csel_all_LT, gBOIN_ET_results_all_LT[[i]]$prop.select)
  gBOIN_ET_poorAlloc_all_LT <- c(gBOIN_ET_poorAlloc_all_LT, gBOIN_ET_poorAllocation_all_LT[[i]]$poor_allocation_rate)
  
}

colnames(gBOIN_ET_npts_all_HE) <- colnames(gBOIN_ET_prop_csel_all_HE) <- colnames(gBOIN_ET_npts_all_LT) <- colnames(gBOIN_ET_prop_csel_all_LT) <- 
  c("Dose 1", "Dose 2", "Dose 3", "Dose 4", "Dose 5", "Dose 6")

# ------------------------- end --------------------------- #


### Fix Scenarios

# 1. OBD correct selection rate
# benchmark methods
Iso_correct_selected_OBD_rate <- read.table(paste0("./Outputs/Benchmark oc/Iso/", INX, "correct_selected_OBD_rate.txt"))
TEPI_correct_selected_OBD_rate <- read.table(paste0("./Outputs/Benchmark oc/TEPI/", INX, "correct_selected_OBD_rate.txt"))
uTPI_correct_selected_OBD_rate <- read.table(paste0("./Outputs/Benchmark oc/uTPI/", INX, "other_oc_uTPI.txt"))
BOIN12_correct_selected_OBD_rate <- read.table(paste0("./Outputs/Benchmark oc/BOIN12/", INX, "other_oc_BOIN12.txt"))

# proposed method
file_list_OBD_rate <- list()
for (i in 1:27) {
  file_path_OBD_rate <- paste0("./Outputs/Ours_Results/OBD/", INX, "correct_selected_OBD_rate_", i, ".txt")
  data_OBD_rate <- read.table(file_path_OBD_rate)
  file_list_OBD_rate[[i]] <- data_OBD_rate
}
Ours_correct_selected_OBD_rate <- do.call(rbind, file_list_OBD_rate)

# ------------------------- end --------------------------- #

# 2. Number of patient allocated to correct OBD
# benchmark methods
Iso_OBD_allocated_patients <- read.table(paste0("./Outputs/Benchmark oc/Iso/", INX, "average_patient_number.txt"))
TEPI_OBD_allocated_patients <- read.table(paste0("./Outputs/Benchmark oc/TEPI/", INX, "average_patient_number.txt"))
uTPI_OBD_allocated_patients <- read.table(paste0("./Outputs/Benchmark oc/uTPI/", INX, "patient_allocation_average_uTPI.txt"))
BOIN12_OBD_allocated_patients <- read.table(paste0("./Outputs/Benchmark oc/BOIN12/", INX, "patient_allocation_average_BOIN12.txt"))

# proposed methods
file_list_OBD_allocated_patients <- list()
for (i in 1:27) {
  file_path_OBD_allocated_patients <- paste0("./Outputs/Ours_Results/Number_of_Patient_at_OBD/", INX, "average_patient_number_", i,".txt")
  data_OBD_allocated_patients <- read.table(file_path_OBD_allocated_patients)
  file_list_OBD_allocated_patients[[i]] <- data_OBD_allocated_patients
}
Ours_OBD_allocated_patients <- do.call(rbind, file_list_OBD_allocated_patients)

# data wrangling
Iso_OBD_allocated_patients <- cbind(Iso_OBD_allocated_patients, OBD_true_vector)
TEPI_OBD_allocated_patients <- cbind(TEPI_OBD_allocated_patients, OBD_true_vector)
uTPI_OBD_allocated_patients <- cbind(uTPI_OBD_allocated_patients, OBD_true_vector)
BOIN12_OBD_allocated_patients <- cbind(BOIN12_OBD_allocated_patients, OBD_true_vector)
Ours_OBD_allocated_patients <- cbind(Ours_OBD_allocated_patients, OBD_true_vector)

Iso_OBD_patient <- numeric(nrow(Iso_OBD_allocated_patients))
TEPI_OBD_patient <- numeric(nrow(TEPI_OBD_allocated_patients))
uTPI_OBD_patient <- numeric(nrow(uTPI_OBD_allocated_patients))
BOIN12_OBD_patient <- numeric(nrow(BOIN12_OBD_allocated_patients))
Ours_OBD_patient <- numeric(nrow(Ours_OBD_allocated_patients))
gBOIN_ET_OBD_patient <- numeric(nrow(Ours_OBD_allocated_patients))

for (i in 1:nrow(Iso_OBD_allocated_patients)) {
  col_index <- OBD_true_vector[i]  
  
  Iso_OBD_patient[i] <- Iso_OBD_allocated_patients[i, col_index]  
  TEPI_OBD_patient[i] <- TEPI_OBD_allocated_patients[i, col_index]  
  uTPI_OBD_patient[i] <- uTPI_OBD_allocated_patients[i, col_index]  
  BOIN12_OBD_patient[i] <- BOIN12_OBD_allocated_patients[i, col_index]  
  Ours_OBD_patient[i] <- Ours_OBD_allocated_patients[i, col_index]  
  
}

# gBOIN-ET
ifelse(INX == "HE_",
       gBOIN_ET_OBD_patient <- sapply(seq_along(OBD_true_vector), function(i) gBOIN_ET_npts_all_HE[i, OBD_true_vector[i]]),
       gBOIN_ET_OBD_patient <- sapply(seq_along(OBD_true_vector), function(i) gBOIN_ET_npts_all_LT[i, OBD_true_vector[i]])
)

combined_data_OBD_patient_check <- cbind(Iso_OBD_patient, TEPI_OBD_patient, uTPI_OBD_patient, BOIN12_OBD_patient, gBOIN_ET_OBD_patient, Ours_OBD_patient, OBD_true_vector)
rownames(combined_data_OBD_patient_check) <- rownames(Iso_OBD_allocated_patients)
colnames(combined_data_OBD_patient_check) <- c("Iso", "TEPI", "uTPI", "BOIN12", "gBOIN-ET", "Proposed", "OBD_true")

# ------------------------- end --------------------------- #

# 3. Poor Allocation Rate
# benchmark methods
Iso_poor_allocation_rate <- read.table(paste0("./Outputs/Benchmark oc/Iso/", INX, "poor_allocation_rate.txt"))
TEPI_poor_allocation_rate <- read.table(paste0("./Outputs/Benchmark oc/TEPI/", INX, "poor_allocation_rate.txt"))
uTPI_poor_allocation_rate <- read.table(paste0("./Outputs/Benchmark oc/uTPI/", INX, "other_oc_uTPI.txt"))
BOIN12_poor_allocation_rate <- read.table(paste0("./Outputs/Benchmark oc/BOIN12/", INX, "other_oc_BOIN12.txt"))

# proposed method
file_list_PA_rate <- list()
for (i in 1:27) {
  file_path_PA_rate <- paste0("./Outputs/Ours_Results/Poor_Allocation/", INX, "poor_allocation_rate_", i, ".txt")
  data_PA_rate <- read.table(file_path_PA_rate)
  file_list_PA_rate[[i]] <- data_PA_rate
}
Ours_poor_allocation_rate <- do.call(rbind, file_list_PA_rate)

Iso_PA_rate <- Iso_poor_allocation_rate[, 1]
TEPI_PA_rate <- TEPI_poor_allocation_rate[, 1]
uTPI_PA_rate <- uTPI_poor_allocation_rate$poor_allocation_rate
BOIN12_PA_rate <- BOIN12_poor_allocation_rate$poor_allocation_rate
Ours_PA_rate <- Ours_poor_allocation_rate[, 1]

# gBOIN-ET
ifelse(INX == "HE_",
       gBOIN_ET_PA_rate <- gBOIN_ET_poorAlloc_all_HE * 100,
       gBOIN_ET_PA_rate <- gBOIN_ET_poorAlloc_all_LT * 100 
)

# ------------------------- end --------------------------- #

# 4. Overdose Patient in Average
# benchmark methods
Iso_overdose_rate <- read.table(paste0("./Outputs/Benchmark oc/Iso/", INX, "average_overdose_numbers.txt"))
TEPI_overdose_rate <- read.table(paste0("./Outputs/Benchmark oc/TEPI/", INX, "average_overdose_numbers.txt"))
uTPI_overdose_rate <- read.table(paste0("./Outputs/Benchmark oc/uTPI/", INX, "other_oc_uTPI.txt")) 
BOIN12_overdose_rate <- read.table(paste0("./Outputs/Benchmark oc/BOIN12/", INX, "other_oc_BOIN12.txt")) 

# proposed method
file_list_OD_rate <- list()
for (i in 1:27) {
  file_path_OD_rate <- paste0("./Outputs/Ours_Results/Overdose/", INX, "average_overdose_numbers_", i, ".txt")
  data_OD_rate <- read.table(file_path_OD_rate)
  file_list_OD_rate[[i]] <- data_OD_rate
}
Ours_overdose_rate <- do.call(rbind, file_list_OD_rate)

Iso_OD_rate <- Iso_overdose_rate[, 1]
TEPI_OD_rate <- TEPI_overdose_rate[, 1]
uTPI_OD_rate <- (uTPI_overdose_rate$overdose)
BOIN12_OD_rate <- (BOIN12_overdose_rate$overdose)
Ours_OD_rate <- Ours_overdose_rate[, 1]

# gBOIN-ET
ifelse(INX == "HE_",
       gBOIN_ET_OD_rate <- mapply(function(row, threshold) { if (threshold < 6) { sum(gBOIN_ET_npts_all_HE[row, (threshold+1):6]) } else { 0.0 }}, row = 1:729, threshold = MTD_true_current),
       gBOIN_ET_OD_rate <- mapply(function(row, threshold) { if (threshold < 6) { sum(gBOIN_ET_npts_all_LT[row, (threshold+1):6]) } else { 0.0 }}, row = 1:729, threshold = MTD_true_current)
)

# ------------------------- end --------------------------- #


#---- Visualization: A
# generate the data.frame to store the OBD correct selection rate of chosen fix scenarios
OBD_rate_all_fix <- as.data.frame(matrix(nrow = length(tox_abb), ncol = 6)) 
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
ifelse(INX == "HE_",
       OBD_rate_all_fix[,5] <- sapply(seq_along(OBD_true_fix_HE), function(i) gBOIN_ET_correct_selected_OBD_rate_HE[i, OBD_true_fix_HE[i]]),
       OBD_rate_all_fix[,5] <- sapply(seq_along(OBD_true_fix_LT), function(i) gBOIN_ET_correct_selected_OBD_rate_LT[i, OBD_true_fix_LT[i]])
)

# modify the data to visualize
OBD_rate_all_fix$scenario <- rownames(OBD_rate_all_fix)
OBD_rate_all_fix <- OBD_rate_all_fix %>%
  pivot_longer(cols = -scenario, names_to = "Methods", values_to = "OBD_rate") %>%
  mutate(Methods = factor(Methods, levels = c("Iso", "TEPI", "uTPI", "BOIN12", "gBOIN-ET", "Proposed")))


#---- Visualization: B
OBD_patients_all_fix <- as.data.frame(matrix(nrow = 9, ncol = 6))
rownames(OBD_patients_all_fix) <- seq(1:9)
colnames(OBD_patients_all_fix) <- c("Iso", "TEPI", "uTPI", "BOIN12", "gBOIN-ET", "Proposed")

# load the number of patient of OBD results
for (i in 1:length(tox_abb)) {   
  OBD_patients_all_fix[i, 1:4] <- combined_data_OBD_patient_check[which(rownames(combined_data_OBD_patient_check) == paste0(tox_abb[i], "_", eff_abb[i])), 1:4]
  OBD_patients_all_fix[i, 6] <- combined_data_OBD_patient_check[which(rownames(combined_data_OBD_patient_check) == paste0(tox_abb[i], "_", eff_abb[i])), 6]
  
  OBD_patients_all_fix[i, c(1:4, 6)] <- round(OBD_patients_all_fix[i, c(1:4, 6)], 2)
}

# add gBOIN-ET results into the data.frame as a new column (notice the HE or LT case)
ifelse(INX == "HE_",
       OBD_patients_all_fix[, 5] <- sapply(seq_along(OBD_true_fix_HE), function(i) gBOIN_ET_correct_OBD_pts_HE[i, OBD_true_fix_HE[i]]),
       OBD_patients_all_fix[, 5] <- sapply(seq_along(OBD_true_fix_LT), function(i) gBOIN_ET_correct_OBD_pts_LT[i, OBD_true_fix_LT[i]])
)

# modify the data to visualize
OBD_patients_all_fix$scenario <- rownames(OBD_patients_all_fix)
OBD_patients_all_fix <- OBD_patients_all_fix %>%
  pivot_longer(cols = -scenario, names_to = "Methods", values_to = "OBD_patients") %>%
  mutate(Methods = factor(Methods, levels = c("Iso", "TEPI", "uTPI", "BOIN12", "gBOIN-ET", "Proposed")))


#---- Visualization: C
# generate the data.frame to store the  poor allocation rate of chosen fix scenarios
PA_rate_all_fix <- as.data.frame(matrix(nrow = length(tox_abb), ncol = 6)) 
rownames(PA_rate_all_fix) <- seq(1:9)
colnames(PA_rate_all_fix) <- c("Iso", "TEPI", "uTPI", "BOIN12", "gBOIN-ET", "Proposed")

# load the poor allocation rate results
for (i in 1:length(tox_abb)) {
  PA_rate_all_fix[i,1] <- Iso_poor_allocation_rate[which(rownames(Iso_poor_allocation_rate) == paste0(tox_abb[i], "_", eff_abb[i])), 1]
  PA_rate_all_fix[i,2] <- TEPI_poor_allocation_rate[which(rownames(TEPI_poor_allocation_rate) == paste0(tox_abb[i], "_", eff_abb[i])), 1]
  PA_rate_all_fix[i,3] <- uTPI_poor_allocation_rate[which(rownames(uTPI_poor_allocation_rate) == paste0(tox_abb[i], "_", eff_abb[i])), 7]
  PA_rate_all_fix[i,4] <- BOIN12_poor_allocation_rate[which(rownames(BOIN12_poor_allocation_rate) == paste0(tox_abb[i], "_", eff_abb[i])), 7]
  
  ifelse(INX == "HE_",
         PA_rate_all_fix[i,5] <- (gBOIN_ET_poorAllocation_fix_HE[[i]]$poor_allocation_rate) * 100,
         PA_rate_all_fix[i,5] <- (gBOIN_ET_poorAllocation_fix_LT[[i]]$poor_allocation_rate) * 100
  )
  
  PA_rate_all_fix[i,6] <- Ours_poor_allocation_rate[which(rownames(Ours_poor_allocation_rate) == paste0(tox_abb[i], "_", eff_abb[i])), 1]
}

# modify the data to visualize
PA_rate_all_fix$scenario <- rownames(PA_rate_all_fix)
PA_rate_all_fix <- PA_rate_all_fix %>%
  pivot_longer(cols = -scenario, names_to = "Methods", values_to = "PA_rate") %>%
  mutate(Methods = factor(Methods, levels = c("Iso", "TEPI", "uTPI", "BOIN12", "gBOIN-ET", "Proposed")))


#---- Visualization: D
# generate the data.frame to store the Overdose rate of chosen fix scenarios
OD_rate_all_fix <- as.data.frame(matrix(nrow = length(tox_abb), ncol = 6)) 
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
ifelse(INX == "HE_",
       OD_rate_all_fix[,5] <- mapply(function(row, threshold) { if (threshold < 6) { sum(gBOIN_ET_correct_OBD_pts_HE[row, (threshold+1):6]) } else { 0.0 }}, row = 1:9, threshold = MTD_true_fix_HE),
       OD_rate_all_fix[,5] <- mapply(function(row, threshold) { if (threshold < 6) { sum(gBOIN_ET_correct_OBD_pts_LT[row, (threshold+1):6]) } else { 0.0 }}, row = 1:9, threshold = MTD_true_fix_LT)  
)

# modify the data to visualize
OD_rate_all_fix$scenario <- rownames(OD_rate_all_fix)
OD_rate_all_fix <- OD_rate_all_fix %>%
  pivot_longer(cols = -scenario, names_to = "Methods", values_to = "OD_rate") %>%
  mutate(Methods = factor(Methods, levels = c("Iso", "TEPI", "uTPI", "BOIN12", "gBOIN-ET", "Proposed")))


# --------------------------------- END --------------------------------- #


### Random Scenarios

# sub-group definition
prefix_scenario_3 <- rownames(summary_tables$toxicity_prob)

prefix_scenario <- paste0(
  expand.grid(rownames(summary_tables$toxicity_prob), rownames(summary_tables$efficacy_prob))[,2],
  "_",
  expand.grid(rownames(summary_tables$toxicity_prob), rownames(summary_tables$efficacy_prob))[,1]
)

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
  
  L_group_tox <- prefix_scenario[substr(prefix_scenario, 1, 3) %in% L_group]
  M_group_tox <- prefix_scenario[substr(prefix_scenario, 1, 3) %in% M_group]
  H_group_tox <- prefix_scenario[substr(prefix_scenario, 1, 3) %in% H_group]
  
  L_group_eff <- prefix_scenario[substr(prefix_scenario, nchar(prefix_scenario) - 2, nchar(prefix_scenario)) %in% L_group]
  M_group_eff <- prefix_scenario[substr(prefix_scenario, nchar(prefix_scenario) - 2, nchar(prefix_scenario)) %in% M_group]
  H_group_eff <- prefix_scenario[substr(prefix_scenario, nchar(prefix_scenario) - 2, nchar(prefix_scenario)) %in% H_group]
  
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
  # create the data frames
  Iso_OBD_rate <- Iso_correct_selected_OBD_rate[, 1]
  TEPI_OBD_rate <- TEPI_correct_selected_OBD_rate[, 1]
  uTPI_OBD_rate <- uTPI_correct_selected_OBD_rate$best_dose_avg
  BOIN12_OBD_rate <- BOIN12_correct_selected_OBD_rate$best_dose_avg
  Ours_OBD_rate <- Ours_correct_selected_OBD_rate[, 1]
  
  ifelse(INX == "HE_",
         gBOIN_ET_OBD_rate <- sapply(seq_along(OBD_true_vector), function(i) gBOIN_ET_prop_csel_all_HE[i, OBD_true_vector[i]]),
         gBOIN_ET_OBD_rate <- sapply(seq_along(OBD_true_vector), function(i) gBOIN_ET_prop_csel_all_LT[i, OBD_true_vector[i]]))

  data_OBD_rate_iso <- data.frame(scenarios = 1:length(Iso_OBD_rate), correct_selected_OBD_rate = Iso_OBD_rate, group = "Iso")
  data_OBD_rate_tepi <- data.frame(scenarios = 1:length(TEPI_OBD_rate), correct_selected_OBD_rate = TEPI_OBD_rate, group = "TEPI")
  data_OBD_rate_utpi <- data.frame(scenarios = 1:length(uTPI_OBD_rate), correct_selected_OBD_rate = uTPI_OBD_rate, group = "uTPI")
  data_OBD_rate_boin12 <- data.frame(scenarios = 1:length(BOIN12_OBD_rate), correct_selected_OBD_rate = BOIN12_OBD_rate, group = "BOIN12")
  data_OBD_rate_ours <- data.frame(scenarios = 1:length(Ours_OBD_rate), correct_selected_OBD_rate = Ours_OBD_rate, group = "Proposed")
  data_OBD_rate_gboinet <- data.frame(scenarios = 1:length(gBOIN_ET_OBD_rate), correct_selected_OBD_rate = gBOIN_ET_OBD_rate, group = "gBOIN-ET")
  
  # merge all data
  combined_data_OBD_rate <- rbind(data_OBD_rate_iso, data_OBD_rate_tepi, data_OBD_rate_utpi, data_OBD_rate_boin12, data_OBD_rate_gboinet,
                                  data_OBD_rate_ours)
  
  # create the data frames
  data_OBD_patient_iso <- data.frame(scenarios = 1:length(Iso_OBD_patient), OBD_patient = Iso_OBD_patient, group = "Iso")
  data_OBD_patient_tepi <- data.frame(scenarios = 1:length(TEPI_OBD_patient), OBD_patient = TEPI_OBD_patient, group = "TEPI")
  data_OBD_patient_utpi <- data.frame(scenarios = 1:length(uTPI_OBD_patient), OBD_patient = uTPI_OBD_patient, group = "uTPI")
  data_OBD_patient_boin12 <- data.frame(scenarios = 1:length(BOIN12_OBD_patient), OBD_patient = BOIN12_OBD_patient, group = "BOIN12")
  data_OBD_patient_gboinet <- data.frame(scenarios = 1:length(gBOIN_ET_OBD_patient), OBD_patient = gBOIN_ET_OBD_patient, group = "gBOIN-ET")
  data_OBD_patient_ours <- data.frame(scenarios = 1:length(Ours_OBD_patient), OBD_patient = Ours_OBD_patient, group = "Proposed")
  
  # merge all data
  combined_data_OBD_patient <- rbind(data_OBD_patient_iso, data_OBD_patient_tepi, data_OBD_patient_utpi, data_OBD_patient_boin12, data_OBD_patient_gboinet, data_OBD_patient_ours)
  
  # create the data frames
  data_PA_rate_iso <- data.frame(scenarios = 1:length(Iso_PA_rate), PA_rate = Iso_PA_rate, group = "Iso")
  data_PA_rate_tepi <- data.frame(scenarios = 1:length(TEPI_PA_rate), PA_rate = TEPI_PA_rate, group = "TEPI")
  data_PA_rate_utpi <- data.frame(scenarios = 1:length(uTPI_PA_rate), PA_rate = uTPI_PA_rate, group = "uTPI")
  data_PA_rate_boin12 <- data.frame(scenarios = 1:length(BOIN12_PA_rate), PA_rate = BOIN12_PA_rate, group = "BOIN12")
  data_PA_rate_gboinet <- data.frame(scenarios = 1:length(gBOIN_ET_PA_rate), PA_rate = gBOIN_ET_PA_rate, group = "gBOIN-ET")
  data_PA_rate_ours <- data.frame(scenarios = 1:length(Ours_PA_rate), PA_rate = Ours_PA_rate, group = "Proposed")

  # merge all data
  combined_data_PA_rate <- rbind(data_PA_rate_iso, data_PA_rate_tepi, data_PA_rate_utpi, data_PA_rate_boin12, data_PA_rate_gboinet, data_PA_rate_ours)
  
  # create the data frames
  data_OD_rate_iso <- data.frame(scenarios = 1:length(Iso_OD_rate), OD_rate = Iso_OD_rate, group = "Iso")
  data_OD_rate_tepi <- data.frame(scenarios = 1:length(TEPI_OD_rate), OD_rate = TEPI_OD_rate, group = "TEPI")
  data_OD_rate_utpi <- data.frame(scenarios = 1:length(uTPI_OD_rate), OD_rate = uTPI_OD_rate, group = "uTPI")
  data_OD_rate_boin12 <- data.frame(scenarios = 1:length(BOIN12_OD_rate), OD_rate = BOIN12_OD_rate, group = "BOIN12")
  data_OD_rate_gboinet <- data.frame(scenarios = 1:length(gBOIN_ET_OD_rate), OD_rate = gBOIN_ET_OD_rate, group = "gBOIN-ET")
  data_OD_rate_ours <- data.frame(scenarios = 1:length(Ours_OD_rate), OD_rate = Ours_OD_rate, group = "Proposed")
  
  # merge all data
  combined_data_OD_rate <- rbind(data_OD_rate_iso, data_OD_rate_tepi, data_OD_rate_utpi, data_OD_rate_boin12, data_OD_rate_gboinet, data_OD_rate_ours)
  
  
  
  RS_OBD_rate <- combined_data_OBD_rate
  RS_OBD_pts <- combined_data_OBD_patient
  RS_PA_rate <- combined_data_PA_rate
  RS_OD_pts <- combined_data_OD_rate
  
  colnames(RS_OBD_rate) <- c("scenarios", "OBD_rate", "design")
  colnames(RS_OD_pts) <- c("scenarios", "OD_pts", "design")
  
  RS_OBD_rate[,1] <- RS_OBD_pts[,1] <- RS_PA_rate[,1] <- RS_OD_pts[,1] <- prefix_scenario
  
  prefix_index <- c("LL", "LM", "LH", "ML", "MM", "MH", "HL", "HM", "HH")
  
  design_index <- c("Iso", "TEPI", "uTPI", "BOIN12", "gBOIN-ET", "Proposed")

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

# --------------------------- Sensitivity Analysis Data -------------------------------- #

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

OBD_rate_N_all_long <- OBD_rate_N_all %>% pivot_longer(cols = -scenario, names_to = "N_star", values_to = "OBD_rate")

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

pts_N_all_long <- pts_N_all %>% pivot_longer(cols = -scenario, names_to = "N_star", values_to = "No.OBD.pts")

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

