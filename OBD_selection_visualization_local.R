rm(list=ls())
# to get the true OBD of each scenario from "OBD_true" data.frame
source("./R files/scenario_generator_local.R")

library(ggplot2)

# The data and visualization in this file is for AVERAGE analysis

# *** The DATA in this file will be used in "fix scenario selection.R" as the source file 
# *** The PLOTS A1 to D1 are testing results

# check the w_e
print(paste0("The current w_e is ", w_e))

# file name distinguishing based on w_e
if(w_e > 0.5){INX <- "HE_"} else {INX <- "LT_"}

OBD_true_vector <- as.vector(t(OBD_true)) # transfer from matrix to vector



# --------- OBD correct selection rate : 729 scenarios------------ #
# benchmark methods
Iso_correct_selected_OBD_rate <- read.table(paste0("./Outputs/Benchmark oc/Iso/", INX, "correct_selected_OBD_rate.txt"))
TEPI_correct_selected_OBD_rate <- read.table(paste0("./Outputs/Benchmark oc/TEPI/", INX, "correct_selected_OBD_rate.txt"))
uTPI_correct_selected_OBD_rate <- read.table(paste0("./Outputs/Benchmark oc/uTPI/", INX, "other_oc_uTPI.txt"))
BOIN12_correct_selected_OBD_rate <- read.table(paste0("./Outputs/Benchmark oc/BOIN12/", INX, "other_oc_BOIN12.txt"))

gBOIN_ET_results_all_HE <- readRDS("./Outputs/gBOIN-ET_results_all_HE.rds")
gBOIN_ET_results_all_LT <- readRDS("./Outputs/gBOIN-ET_results_all_LT.rds")
gBOIN_ET_poorAllocation_all_HE <- readRDS("./Outputs/gBOIN-ET_poorAllocation_all_HE.rds")
gBOIN_ET_poorAllocation_all_LT <- readRDS("./Outputs/gBOIN-ET_poorAllocation_all_LT.rds")


# proposed method
file_list_OBD_rate <- list()
for (i in 1:27) {
  file_path_OBD_rate <- paste0("./Outputs/Ours_Results/OBD/", INX, "correct_selected_OBD_rate_", i, ".txt")
  data_OBD_rate <- read.table(file_path_OBD_rate)
  file_list_OBD_rate[[i]] <- data_OBD_rate
}
Ours_correct_selected_OBD_rate <- do.call(rbind, file_list_OBD_rate)


# ----- gBOIN-ET data cleaning
# --- HE
gBOIN_ET_npts_all_HE <- data.frame()
gBOIN_ET_prop_csel_all_HE <- data.frame()
gBOIN_ET_poorAlloc_all_HE <- c()

# ---LT
gBOIN_ET_npts_all_LT <- data.frame()
gBOIN_ET_prop_csel_all_LT <- data.frame()
gBOIN_ET_poorAlloc_all_LT <- c()

for (i in 1:729) {
  # --- HE
  gBOIN_ET_npts_all_HE <- rbind(gBOIN_ET_npts_all_HE, gBOIN_ET_results_all_HE[[i]]$n.patient)
  gBOIN_ET_prop_csel_all_HE <- rbind(gBOIN_ET_prop_csel_all_HE, gBOIN_ET_results_all_HE[[i]]$prop.select)
  gBOIN_ET_poorAlloc_all_HE <- c(gBOIN_ET_poorAlloc_all_HE, gBOIN_ET_poorAllocation_all_HE[[i]]$poor_allocation_rate)

  # --- LT
  gBOIN_ET_npts_all_LT <- rbind(gBOIN_ET_npts_all_LT, gBOIN_ET_results_all_LT[[i]]$n.patient)
  gBOIN_ET_prop_csel_all_LT <- rbind(gBOIN_ET_prop_csel_all_LT, gBOIN_ET_results_all_LT[[i]]$prop.select)
  gBOIN_ET_poorAlloc_all_LT <- c(gBOIN_ET_poorAlloc_all_LT, gBOIN_ET_poorAllocation_all_LT[[i]]$poor_allocation_rate)

}

colnames(gBOIN_ET_npts_all_HE) <- colnames(gBOIN_ET_prop_csel_all_HE) <- colnames(gBOIN_ET_npts_all_LT) <- colnames(gBOIN_ET_prop_csel_all_LT) <- 
  c("Dose 1", "Dose 2", "Dose 3", "Dose 4", "Dose 5", "Dose 6")



# ----- Data Cleaning by oc

Iso_OBD_rate <- Iso_correct_selected_OBD_rate[, 1]
TEPI_OBD_rate <- TEPI_correct_selected_OBD_rate[, 1]
uTPI_OBD_rate <- uTPI_correct_selected_OBD_rate$best_dose_avg
BOIN12_OBD_rate <- BOIN12_correct_selected_OBD_rate$best_dose_avg
Ours_OBD_rate <- Ours_correct_selected_OBD_rate[, 1]

# gBOIN_ET_OBD_rate <- sapply(seq_along(OBD_true_vector), function(i) gBOIN_ET_prop_csel_all_HE[i, OBD_true_vector[i]]) # for HE
gBOIN_ET_OBD_rate <- sapply(seq_along(OBD_true_vector), function(i) gBOIN_ET_prop_csel_all_LT[i, OBD_true_vector[i]]) # for LT



# create the data frames
data_OBD_rate_iso <- data.frame(scenarios = 1:length(Iso_OBD_rate), correct_selected_OBD_rate = Iso_OBD_rate, group = "Iso")
data_OBD_rate_tepi <- data.frame(scenarios = 1:length(TEPI_OBD_rate), correct_selected_OBD_rate = TEPI_OBD_rate, group = "TEPI")
data_OBD_rate_utpi <- data.frame(scenarios = 1:length(uTPI_OBD_rate), correct_selected_OBD_rate = uTPI_OBD_rate, group = "uTPI")
data_OBD_rate_boin12 <- data.frame(scenarios = 1:length(BOIN12_OBD_rate), correct_selected_OBD_rate = BOIN12_OBD_rate, group = "BOIN12")
data_OBD_rate_ours <- data.frame(scenarios = 1:length(Ours_OBD_rate), correct_selected_OBD_rate = Ours_OBD_rate, group = "Proposed")

data_OBD_rate_gboinet <- data.frame(scenarios = 1:length(gBOIN_ET_OBD_rate), correct_selected_OBD_rate = gBOIN_ET_OBD_rate, group = "gBOIN-ET")

# merge all data
combined_data_OBD_rate <- rbind(data_OBD_rate_iso, data_OBD_rate_tepi, data_OBD_rate_utpi, data_OBD_rate_boin12, data_OBD_rate_gboinet,
                                data_OBD_rate_ours)
#combined_data_OBD_rate <- rbind(data_OBD_rate_iso, data_OBD_rate_tepi, data_OBD_rate_utpi, data_OBD_rate_boin12, data_OBD_rate_gboinet)


# ------------------------- end --------------------------- #




# ------------- Number of patient allocated to correct OBD : 729 scenarios -------------- #

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

for (i in 1:nrow(Iso_OBD_allocated_patients)) { # 729 fixed, could use any one method nrow() to get
  col_index <- OBD_true_vector[i]  # get OBD dose in true OBD column
  
  Iso_OBD_patient[i] <- Iso_OBD_allocated_patients[i, col_index]  
  TEPI_OBD_patient[i] <- TEPI_OBD_allocated_patients[i, col_index]  
  uTPI_OBD_patient[i] <- uTPI_OBD_allocated_patients[i, col_index]  
  BOIN12_OBD_patient[i] <- BOIN12_OBD_allocated_patients[i, col_index]  
  Ours_OBD_patient[i] <- Ours_OBD_allocated_patients[i, col_index]  
  
}
#gBOIN-ET
# gBOIN_ET_OBD_patient <- sapply(seq_along(OBD_true_vector), function(i) gBOIN_ET_npts_all_HE[i, OBD_true_vector[i]]) # for HE
gBOIN_ET_OBD_patient <- sapply(seq_along(OBD_true_vector), function(i) gBOIN_ET_npts_all_LT[i, OBD_true_vector[i]]) # for LT


# check 
# with Ours
combined_data_OBD_patient_check <- cbind(Iso_OBD_patient, TEPI_OBD_patient, uTPI_OBD_patient, BOIN12_OBD_patient, gBOIN_ET_OBD_patient, Ours_OBD_patient, OBD_true_vector)
# without Ours
#combined_data_OBD_patient_check <- cbind(Iso_OBD_patient, TEPI_OBD_patient, uTPI_OBD_patient, BOIN12_OBD_patient, gBOIN_ET_OBD_patient, OBD_true_vector)

rownames(combined_data_OBD_patient_check) <- rownames(Iso_OBD_allocated_patients)
# with Ours
colnames(combined_data_OBD_patient_check) <- c("Iso", "TEPI", "uTPI", "BOIN12", "gBOIN-ET", "Proposed", "OBD_true")
# without Ours
#colnames(combined_data_OBD_patient_check) <- c("Iso", "TEPI", "uTPI", "BOIN12", "gBOIN-ET", "OBD_true")


# create the data frames
data_OBD_patient_iso <- data.frame(scenarios = 1:length(Iso_OBD_patient), OBD_patient = Iso_OBD_patient, group = "Iso")
data_OBD_patient_tepi <- data.frame(scenarios = 1:length(TEPI_OBD_patient), OBD_patient = TEPI_OBD_patient, group = "TEPI")
data_OBD_patient_utpi <- data.frame(scenarios = 1:length(uTPI_OBD_patient), OBD_patient = uTPI_OBD_patient, group = "uTPI")
data_OBD_patient_boin12 <- data.frame(scenarios = 1:length(BOIN12_OBD_patient), OBD_patient = BOIN12_OBD_patient, group = "BOIN12")
data_OBD_patient_gboinet <- data.frame(scenarios = 1:length(gBOIN_ET_OBD_patient), OBD_patient = gBOIN_ET_OBD_patient, group = "gBOIN-ET")
data_OBD_patient_ours <- data.frame(scenarios = 1:length(Ours_OBD_patient), OBD_patient = Ours_OBD_patient, group = "Proposed")


# merge all data
# with Ours
combined_data_OBD_patient <- rbind(data_OBD_patient_iso, data_OBD_patient_tepi, data_OBD_patient_utpi, data_OBD_patient_boin12, data_OBD_patient_gboinet, data_OBD_patient_ours)
# without Ours
#combined_data_OBD_patient <- rbind(data_OBD_patient_iso, data_OBD_patient_tepi, data_OBD_patient_utpi, data_OBD_patient_boin12)

# ------------------------- end --------------------------- #


# --------- Poor Allocation Rate : 729 scenarios------------ #
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

# gBOIN_ET_PA_rate <- gBOIN_ET_poorAlloc_all_HE * 100 # for HE
gBOIN_ET_PA_rate <- gBOIN_ET_poorAlloc_all_LT * 100 # for LT


# create the data frames
data_PA_rate_iso <- data.frame(scenarios = 1:length(Iso_PA_rate), PA_rate = Iso_PA_rate, group = "Iso")
data_PA_rate_tepi <- data.frame(scenarios = 1:length(TEPI_PA_rate), PA_rate = TEPI_PA_rate, group = "TEPI")
data_PA_rate_utpi <- data.frame(scenarios = 1:length(uTPI_PA_rate), PA_rate = uTPI_PA_rate, group = "uTPI")
data_PA_rate_boin12 <- data.frame(scenarios = 1:length(BOIN12_PA_rate), PA_rate = BOIN12_PA_rate, group = "BOIN12")
data_PA_rate_gboinet <- data.frame(scenarios = 1:length(gBOIN_ET_PA_rate), PA_rate = gBOIN_ET_PA_rate, group = "gBOIN-ET")
data_PA_rate_ours <- data.frame(scenarios = 1:length(Ours_PA_rate), PA_rate = Ours_PA_rate, group = "Proposed")



# merge all data
combined_data_PA_rate <- rbind(data_PA_rate_iso, data_PA_rate_tepi, data_PA_rate_utpi, data_PA_rate_boin12, data_PA_rate_gboinet, data_PA_rate_ours)
#combined_data_PA_rate <- rbind(data_PA_rate_iso, data_PA_rate_tepi, data_PA_rate_utpi, data_PA_rate_boin12, data_PA_rate_gboinet)

# ------------------------- end --------------------------- #


# --------- Overdose Patient in Average : 729 scenarios------------ #
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

# gBOIN_ET_OD_rate <- mapply(function(row, threshold) { if (threshold < 6) { sum(gBOIN_ET_npts_all_HE[row, (threshold+1):6]) } else { 0.0 }}, row = 1:729, threshold = MTD_true_current)
gBOIN_ET_OD_rate <- mapply(function(row, threshold) { if (threshold < 6) { sum(gBOIN_ET_npts_all_LT[row, (threshold+1):6]) } else { 0.0 }}, row = 1:729, threshold = MTD_true_current)


# create the data frames
data_OD_rate_iso <- data.frame(scenarios = 1:length(Iso_OD_rate), OD_rate = Iso_OD_rate, group = "Iso")
data_OD_rate_tepi <- data.frame(scenarios = 1:length(TEPI_OD_rate), OD_rate = TEPI_OD_rate, group = "TEPI")
data_OD_rate_utpi <- data.frame(scenarios = 1:length(uTPI_OD_rate), OD_rate = uTPI_OD_rate, group = "uTPI")
data_OD_rate_boin12 <- data.frame(scenarios = 1:length(BOIN12_OD_rate), OD_rate = BOIN12_OD_rate, group = "BOIN12")
data_OD_rate_gboinet <- data.frame(scenarios = 1:length(gBOIN_ET_OD_rate), OD_rate = gBOIN_ET_OD_rate, group = "gBOIN-ET")
data_OD_rate_ours <- data.frame(scenarios = 1:length(Ours_OD_rate), OD_rate = Ours_OD_rate, group = "Proposed")

# merge all data
combined_data_OD_rate <- rbind(data_OD_rate_iso, data_OD_rate_tepi, data_OD_rate_utpi, data_OD_rate_boin12, data_OD_rate_gboinet, data_OD_rate_ours)
#combined_data_OD_rate <- rbind(data_OD_rate_iso, data_OD_rate_tepi, data_OD_rate_utpi, data_OD_rate_boin12, data_OD_rate_gboinet)


# ------------------------- end --------------------------- #







