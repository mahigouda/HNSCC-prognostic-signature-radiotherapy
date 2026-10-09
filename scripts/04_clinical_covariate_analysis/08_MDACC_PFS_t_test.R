## =============================================================================
## Project : Prognostic gene-expression signature for radiotherapy outcome in HPV-negative HNSCC
## Module  : 04 | Clinical covariates: multivariable Cox, risk groups vs. clinical variables
## Script  : 08_MDACC_PFS_t_test.R
## Purpose : MDACC cohort, PFS signature: multivariable Cox models (signature genes + clinical
##           covariates), forest plots, median-split high/low risk groups, Kaplan-Meier curves
##           and t-tests of risk groups/risk scores against clinical variables.
## Input   : data/clinical_covariates/merged_df_MDACC_PFS_for_analysis.txt
## Output  : results/04_clinical_covariate_analysis/MDACC_PFS_*.txt and plots
## =============================================================================

## ---- Packages ---------------------------------------------------------------

library(Biobase)
library(purrr)
library(scales)
library(readr)
library(survival)
library(caret)
library(dplyr)
library(ggplot2)
library(SummarizedExperiment)
library(limma)
library(edgeR)
library(tibble)
library(stringr)
library(Matrix)
library(glmnet)
library(ggpubr)
library(survminer)
library(devtools)
library(AnnotationDbi)
library(MASS)
library(tidyr)
library(cowplot)
library(pheatmap)
library(grid)
library(leaps)
library(GGally)
library(corrplot)
library(broom)
library(biomaRt)
library(dbplyr)
library(rlang)
library(DBI)
library(httr)
library(metadat)
library(numDeriv)
library(metafor)
library(plyr)
library(forestmodel)
library(car)
library(forestplot)

## ---- Paths ------------------------------------------------------------------
## Run from the repository root. Inputs are read from data/clinical_covariates/
## (see data/README.md); outputs are written to results/04_clinical_covariate_analysis/.
data_dir <- file.path(getwd(), "data", "clinical_covariates")
out_dir  <- file.path(getwd(), "results", "04_clinical_covariate_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)

## Sub data
file_name <- file.path(data_dir, "merged_df_MDACC_PFS_for_analysis.txt")
signature_df <- read.table(file_name, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = "\t")

# Map the numeric stages to their descriptive names
signature_df$Stages <- mapvalues(signature_df$Stages,
                                 from = c('1', '2', '3', '4'),
                                 to = c('STAGE I', 'STAGE II', 'STAGE III', 'STAGE IV'))

# Map Vital status to numeric
signature_df$Vital_status <- mapvalues(signature_df$Vital_status,
                                       from = c('Alive', 'Dead'),
                                       to = c('0', '1'))

# Remove rows with NA values in any column
signature_df <- na.omit(signature_df)

# specify the range of columns more safely as follows:
num_cols <- ncol(signature_df)
gene_expression_cols <- signature_df[, 1:(num_cols - 7)]

# Now calculate the median for each gene
signature_df_median <- apply(gene_expression_cols, 2, median)

# categorize age or other continuous variables if needed
# For example, categorizing age into groups
signature_df$age_group <- cut(signature_df$Age,
                              breaks=c(0, 50, 65, 80, Inf), labels=c("0-50", "51-65", "66-80", "80+"))

# Convert OS_Months from months to years
signature_df$OS_Months <- signature_df$OS_Months / 12

# Convert categorical variables to factors
signature_df$Sex <- as.factor(signature_df$Sex)
signature_df$Smoking_status <- as.factor(signature_df$Smoking_status)
signature_df$Stages <- as.factor(signature_df$Stages)
signature_df$Treatment <- as.factor(signature_df$Treatment)
signature_df$age_group <- as.factor(signature_df$age_group)

# Check for NA, NaN, or Inf in your dataset and handle them
signature_df <- na.omit(signature_df) # This removes rows with NA values

# Inspect the Vital_status column
table(signature_df$Vital_status) # Check for any non-numeric values

# Convert Vital_status to a character type
signature_df$Vital_status <- as.character(signature_df$Vital_status)

# Filter out rows where Vital_status is "."
signature_df <- signature_df[signature_df$Vital_status != ".", ]

# Now, convert Vital_status back to numeric
signature_df$Vital_status <- as.numeric(signature_df$Vital_status)

# Check the table again to ensure the periods are removed
table(signature_df$Vital_status)

# Save signature_df as a text file in the working directory
write.table(signature_df, file = "MDACC_PFS_signature_df_r_output.txt", sep = "\t", row.names = FALSE)

# List of gene columns to convert to numeric
genes_to_convert = c("ADAM8", "RAC2", "CSF2", "PPP1R18", "HBEGF", "ITGA5",
                     "LRRC59", "PLAU", "FOSL1", "PLAUR")

# Replace all periods with nothing and convert to numeric
for (gene in genes_to_convert) {
  signature_df[[gene]] <- as.numeric(gsub("\\.", "", signature_df[[gene]]))
}


# Convert Vital_status to a Surv object
surv_obj <- Surv(time = signature_df$OS_Months, event = signature_df$Vital_status)

# Fit a Cox proportional hazards model
cox_model <- coxph(surv_obj ~ ADAM8 + RAC2 +
                     CSF2 + PPP1R18 + HBEGF + ITGA5 + LRRC59 + PLAU + FOSL1 + PLAUR + Age +
                     Sex + Smoking_status + Stages, data = signature_df)
summary(cox_model)

# Example: Fit a Cox model with fewer variables
cox_model_simple <- coxph(surv_obj ~ Age + Sex + Smoking_status, data = signature_df)
summary(cox_model_simple)

# Example: Expanding the model to include gene variables
extended_cox_model <- coxph(surv_obj ~ Age + Sex + Smoking_status + ADAM8 + RAC2 +
                              CSF2 + PPP1R18 + HBEGF + ITGA5 + LRRC59 + PLAU + FOSL1 + PLAUR, data = signature_df)
summary(extended_cox_model)

# Extract model coefficients, confidence intervals
coef_df <- summary(extended_cox_model)$coefficients

# Display the column names of coef_df
colnames(coef_df)

# Prepare the data for the forest plot
# Calculate the 95% confidence intervals
forest_data <- data.frame(
  Variable = rownames(coef_df),
  HR = exp(coef_df[, "coef"]), # Hazard Ratios
  lower = exp(coef_df[, "coef"] - 1.96 * coef_df[, "se(coef)"]), # Lower bound
  upper = exp(coef_df[, "coef"] + 1.96 * coef_df[, "se(coef)"])  # Upper bound
)

# Create the forest plot
forestplot(labeltext = forest_data$Variable,
           mean = forest_data$HR,
           lower = forest_data$lower,
           upper = forest_data$upper,
           xlog = TRUE,  # Use log scale for Hazard Ratios
           title = "Forest Plot of Cox Proportional Hazards Model",
           xlab = "Hazard Ratios (log scale)")


# Vital status continuous outcome variable in the signature_df dataset
linear_model <- lm(Vital_status ~ Age + Sex + Smoking_status +
                     ADAM8 + RAC2 + CSF2 + PPP1R18 + HBEGF + ITGA5 + LRRC59
                   + PLAU + FOSL1 + PLAUR, data = signature_df)

# View the summary of the linear model
summary(linear_model)

forest_model(linear_model)

## Median High Vs low risk
# Assuming 'extended_cox_model' is your final Cox model
# Calculate risk scores
risk_scores <- predict(extended_cox_model, type = "risk")

# Check the number of rows in the dataset and the length of risk_scores
nrow(signature_df)
length(risk_scores)

# Assuming the first 251 rows of signature_df were used in the Cox model
signature_df <- signature_df[1:74, ]

# Exclude observations from signature_df that were not used in the Cox model
signature_df <- signature_df[!is.na(risk_scores), ]

# Add risk scores to the dataset
signature_df$risk_score <- risk_scores

# Median split for high vs low risk
signature_df$risk_group <- ifelse(signature_df$risk_score > median(signature_df$risk_score), "High", "Low")

surv_obj <- Surv(signature_df$OS_Months, signature_df$Vital_status)
km_fit <- survfit(surv_obj ~ risk_group, data = signature_df)

signature_df$risk_group <- as.factor(signature_df$risk_group)


# Plot survival curves
plot(km_fit, main = "Survival Curves for High vs. Low Risk Groups", xlab = "Time", ylab = "Survival Probability", col = 1:2)
legend("topright", legend = levels(signature_df$risk_group), col = 1:2, lty = 1)

# Log-rank test
survdiff(surv_obj ~ risk_group, data = signature_df)

# Log-rank test
survdiff(surv_obj ~ risk_group, data = signature_df)


# Making a copy of signature_df
signature_df_1 <- signature_df

# Assuming the first 13 columns are gene columns
gene_columns <- names(signature_df_1)[1:10]

# Create new columns for each gene to indicate risk category
for (gene in gene_columns) {
  new_col_name <- paste(gene, "Risk", sep = "_")
  signature_df_1 <- signature_df_1 %>%
    mutate(!!new_col_name := ifelse(get(gene) >= some_threshold, "High", "Low"))
}


for (col_name in colnames(signature_df_1)) {
  # Skip for non-continuous variables or already existing risk categories
  if (col_name %in% c("Age", "Sex", "Smoking_status", "Stages", "Treatment", "age_group", "risk_group")) next

  risk_col_name <- paste0(col_name, "_Risk")
  signature_df_1[[risk_col_name]] <- categorize_risk(signature_df_1[[col_name]])
}


signature_df_1 <- signature_df_1 %>%
  mutate(
    Stages_Risk = ifelse(Stages %in% c("STAGE III", "STAGE IV"), "High", "Low"),
    Treatment_Risk = ifelse(Treatment %in% c("cR_S", "IC_S", "IC_S_cRT", "IC_S_RT", "S", "S_cRT", "S_IC_S_cRT", "S_RT"), "High", "Low"),
    age_group_Risk = ifelse(age_group %in% c("0-50", "51-65", "66-80", "80+"), "High", "Low")
  )

# Summary for Stages
table(signature_df_1$Stages, signature_df_1$Stages_Risk)

# Summary for Radiation Therapy
table(signature_df_1$Treatment, signature_df_1$Treatment_Risk)

# Summary for Age Group
table(signature_df_1$age_group, signature_df_1$age_group_Risk)

signature_df_1$Stages_Risk <- as.factor(signature_df_1$Stages_Risk)
signature_df_1$Treatment_Risk <- as.factor(signature_df_1$Treatment_Risk)
signature_df_1$age_group_Risk <- as.factor(signature_df_1$age_group_Risk)


levels(signature_df_1$Stages_Risk)
levels(signature_df_1$Treatment_Risk)
levels(signature_df_1$age_group_Risk)

## Linear model
linear_model_1 <- lm(Vital_status ~ risk_group + age_group + Sex + Smoking_status + Stages + Treatment +
                       ADAM8_Risk + RAC2_Risk + CSF2_Risk + PPP1R18_Risk + HBEGF_Risk + ITGA5_Risk +
                       LRRC59_Risk + PLAU_Risk + FOSL1_Risk +
                       PLAUR_Risk,
                     data = signature_df_1)
summary(linear_model_1)
forest_model(linear_model_1)


# make a copy of signatur_df_1
signature_df_2 <- signature_df_1

# Assuming columns 1 to 10 are your gene expression values
gene_expression_cols <- 1:10

# Checking for any negative or zero values in gene expression columns
if(any(signature_df_2[gene_expression_cols] <= 0)) {
  print("There are non-positive values in gene expression columns")
}


# Apply log2 transformation after adding a small constant
signature_df_2[gene_expression_cols] <- log2(signature_df_2[gene_expression_cols] + 1)

# Scale the log-transformed data
signature_df_2[gene_expression_cols] <- scale(signature_df_2[gene_expression_cols])

# Calculate the mean for each row
signature_df_2$Gene_Risk_Score <- rowMeans(signature_df_2[gene_expression_cols], na.rm = TRUE)

# Relabel the column
signature_df_2 <- dplyr::rename(signature_df_2, PFS_Signature_risk_score = Gene_Risk_Score)


# make a copy of signatur_df_1
signature_df_3 <- signature_df_2

# Classify age_groups into <65 and >65
signature_df_3 <- signature_df_3 %>%
  mutate(new_age_group = ifelse(age_group %in% c("0-50", "51-65"), "<65",
                                ifelse(age_group %in% c("66-80", "80+"), ">65", NA)))

# Check the first few rows of the modified dataframe
head(signature_df_3)

# Update the encoding for age_group to reflect the new classification
signature_df_3 <- signature_df_3 %>%
  mutate(age_group_encoded = ifelse(new_age_group == "<65", 0, 1)) # 0 for "<65", 1 for ">65"

# Encode Sex and Stages as numeric
signature_df_3$Sex_encoded <- ifelse(signature_df_3$Sex == "F", 0, 1)  # Example: 0 for Female, 1 for Male
signature_df_3$Stages_encoded <- as.numeric(as.factor(signature_df_3$Stages))

# Update your categories to reflect the new age group classification
age_groups <- c("<65", ">65")
stages <- c("STAGE I", "STAGE II", "STAGE III", "STAGE IV")
sexes <- c("F", "M")
risk_groups <- c("High", "Low")
treatments <- c("cR_S", "IC_S", "IC_S_cRT", "IC_S_RT", "S", "S_cRT", "S_IC_S_cRT", "S_RT")
smoking_status <- c("Current", "Former", "Never_Smoker")


# Example for one group: Age < 65 and High Risk
subset_data <- subset(signature_df_3, new_age_group == "<65" & risk_group == "High")


# Check the data type of PFS_Signature_risk_score
if ("PFS_Signature_risk_score" %in% names(signature_df_3)) {
  data_type <- class(signature_df_3$PFS_Signature_risk_score)
  print(paste("Data type of PFS_Signature_risk_score:", data_type))
} else {
  print("PFS_Signature_risk_score does not exist in the dataset.")
}


# Categorizing PFS_Signature_risk_score
median_score <- median(signature_df_3$PFS_Signature_risk_score, na.rm = TRUE)
signature_df_3$Categorized_Risk_Score <- ifelse(signature_df_3$PFS_Signature_risk_score >= median_score, "High", "Low")

# Convert to factor
signature_df_3$Categorized_Risk_Score <- factor(signature_df_3$Categorized_Risk_Score)

# Create a contingency table
contingency_table <- table(signature_df_3$new_age_group, signature_df_3$Categorized_Risk_Score)

# Print the contingency table for inspection
print("Contingency Table for Categorized OS_Signature_risk_score and new_age_group")
print(contingency_table)


# Initialize a data frame to store the t-test results for age groups
t_test_results_age <- data.frame(
  Age_Group = character(),
  Mean_Difference = numeric(),
  P_Value = numeric(),
  Test_Valid = logical(),
  stringsAsFactors = FALSE
)

# Loop through each age group
for (age in age_groups) {  # Using the age_groups vector defined earlier
  # Filter data for the specific age group and High risk
  high_risk_scores <- signature_df_3$PFS_Signature_risk_score[signature_df_3$new_age_group == age & signature_df_3$Categorized_Risk_Score == "High"]

  # Filter data for the specific age group and Low risk
  low_risk_scores <- signature_df_3$PFS_Signature_risk_score[signature_df_3$new_age_group == age & signature_df_3$Categorized_Risk_Score == "Low"]

  # Check if there are enough data points in both groups
  if (length(high_risk_scores) > 1 && length(low_risk_scores) > 1) {
    # Perform t-test
    t_test_result <- t.test(high_risk_scores, low_risk_scores, var.equal = TRUE) # assuming equal variances

    # Add the results to the data frame
    t_test_results_age <- rbind(t_test_results_age, data.frame(
      Age_Group = age,
      Mean_Difference = mean(high_risk_scores) - mean(low_risk_scores),
      P_Value = t_test_result$p.value,
      Test_Valid = TRUE
    ))
  } else {
    # Not enough data to perform t-test
    t_test_results_age <- rbind(t_test_results_age, data.frame(
      Age_Group = age,
      Mean_Difference = NA,
      P_Value = NA,
      Test_Valid = FALSE
    ))
  }
}

# Print the t-test results table for age groups
print(t_test_results_age)


# Save results_table as a text file
write.table(t_test_results_age, file = "MDACC_PFS_signature_genes_risk_scores_age_group_t_test.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# make a copy of signatur_df_2
signature_df_4 <- signature_df_2

# Classify age_groups into <65 and >65
signature_df_4 <- signature_df_4 %>%
  mutate(new_age_group = ifelse(age_group %in% c("0-50", "51-65"), "<65",
                                ifelse(age_group %in% c("66-80", "80+"), ">65", NA)))

# Classify age_groups into <65 and >65
signature_df_4 <- signature_df_4 %>%
  mutate(new_stages = ifelse(Stages %in% c("STAGE I", "STAGE II"), "Stage I/II",
                             ifelse(Stages %in% c("STAGE III", "Stage IV"), "Stage III/IV", NA)))


# Update the encoding for age_group to reflect the new classification
signature_df_4 <- signature_df_4 %>%
  mutate(age_group_encoded = ifelse(new_age_group == "<65", 0, 1)) # 0 for "<65", 1 for ">65"

# Encode Sex and Stages as numeric
signature_df_4$Sex_encoded <- ifelse(signature_df_4$Sex == "F", 0, 1)  # Example: 0 for Female, 1 for Male
signature_df_4$new_stages_encoded <- as.numeric(as.factor(signature_df_4$new_stages))


# Update your categories to reflect the new age group classification
age_groups <- c("<65", ">65")
stages <- c("Stage I/II", "Stage III/IV")
sexes <- c("F", "M")
risk_groups <- c("High", "Low")
treatments <- c("cR_S", "IC_S", "IC_S_cRT", "IC_S_RT", "S", "S_cRT", "S_IC_S_cRT", "S_RT")

print(table(signature_df_4$new_stages))
print(table(signature_df_4$risk_group))

# Categorizing PFS_Signature_risk_score
median_score <- median(signature_df_4$PFS_Signature_risk_score, na.rm = TRUE)
signature_df_4$Categorized_Risk_Score <- ifelse(signature_df_4$PFS_Signature_risk_score >= median_score, "High", "Low")

# Convert to factor
signature_df_4$Categorized_Risk_Score <- factor(signature_df_4$Categorized_Risk_Score)

# Create a contingency table
contingency_table <- table(signature_df_4$new_stages, signature_df_4$Categorized_Risk_Score)


# Initialize a data frame to store the t-test results for stages
t_test_results_stages <- data.frame(
  Stage_Group = character(),
  Mean_Difference = numeric(),
  P_Value = numeric(),
  Test_Valid = logical(),
  stringsAsFactors = FALSE
)

# Loop through each stage group
for (stage in c("Stage I/II", "Stage III/IV")) {
  # Filter data for the specific stage group and High risk
  high_risk_scores <- signature_df_4$PFS_Signature_risk_score[signature_df_4$new_stages == stage & signature_df_4$Categorized_Risk_Score == "High"]

  # Filter data for the specific stage group and Low risk
  low_risk_scores <- signature_df_4$PFS_Signature_risk_score[signature_df_4$new_stages == stage & signature_df_4$Categorized_Risk_Score == "Low"]

  # Check if there are enough data points in each group
  if (length(high_risk_scores) >= 2 && length(low_risk_scores) >= 2) {
    # Perform t-test
    t_test_result <- t.test(high_risk_scores, low_risk_scores)

    # Add the results to the data frame
    t_test_results_stages <- rbind(t_test_results_stages, data.frame(
      Stage_Group = stage,
      Mean_Difference = t_test_result$estimate[1] - t_test_result$estimate[2],
      P_Value = t_test_result$p.value,
      Test_Valid = TRUE
    ))
  } else {
    # Not enough data to perform t-test
    t_test_results_stages <- rbind(t_test_results_stages, data.frame(
      Stage_Group = stage,
      Mean_Difference = NA,
      P_Value = NA,
      Test_Valid = FALSE
    ))
  }
}

# Print the t-test results table for stages
print(t_test_results_stages)


# Save results_table as a text file
write.table(t_test_results_stages, file = "MDACC_PFS_signature_genes_risk_scores_stages_t_test.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# make a copy of signatur_df_4
signature_df_5 <- signature_df_4


# Initialize a data frame to store the t-test results for sex
t_test_results_sex <- data.frame(
  Sex_Group = character(),
  Mean_Difference = numeric(),
  P_Value = numeric(),
  Test_Valid = logical(),
  stringsAsFactors = FALSE
)

# Loop through each sex group
for (sex_code in unique(signature_df_5$Sex_encoded)) {
  sex_label <- ifelse(sex_code == 0, "Female", "Male")

  # Filter data for the specific sex group and High risk
  high_risk_scores <- signature_df_5$PFS_Signature_risk_score[signature_df_5$Sex_encoded == sex_code & signature_df_5$Categorized_Risk_Score == "High"]

  # Filter data for the specific sex group and Low risk
  low_risk_scores <- signature_df_5$PFS_Signature_risk_score[signature_df_5$Sex_encoded == sex_code & signature_df_5$Categorized_Risk_Score == "Low"]

  # Perform t-test if there are enough data points in both groups
  if (length(high_risk_scores) > 1 && length(low_risk_scores) > 1) {
    t_test_result <- t.test(high_risk_scores, low_risk_scores, var.equal = TRUE)

    # Add the results to the data frame
    t_test_results_sex <- rbind(t_test_results_sex, data.frame(
      Sex_Group = sex_label,
      Mean_Difference = mean(high_risk_scores) - mean(low_risk_scores),
      P_Value = t_test_result$p.value,
      Test_Valid = TRUE
    ))
  } else {
    # Not enough data to perform t-test
    t_test_results_sex <- rbind(t_test_results_sex, data.frame(
      Sex_Group = sex_label,
      Mean_Difference = NA,
      P_Value = NA,
      Test_Valid = FALSE
    ))
  }
}

# Print the t-test results table for sex
print(t_test_results_sex)


# Save results_table as a text file
write.table(t_test_results_sex, file = "MDACC_PFS_signature_genes_risk_scores_Sex_t_test.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# make a copy of signatur_df_4
signature_df_6 <- signature_df_4

signature_df_6$Smoking_Status_encoded <- as.numeric(as.factor(signature_df_6$Smoking_status))

smoking_status <- c("Current", "Former", "Never_Smoker")


# Initialize a data frame to store the t-test results for smoking status
t_test_results_smoking <- data.frame(
  Smoking_Status = character(),
  Mean_Difference = numeric(),
  P_Value = numeric(),
  Test_Valid = logical(),
  stringsAsFactors = FALSE
)

# Loop through each smoking status
for (status in smoking_status) {
  # Filter data for the specific smoking status group and High risk
  high_risk_scores <- signature_df_6$PFS_Signature_risk_score[signature_df_6$Smoking_status == status & signature_df_6$Categorized_Risk_Score == "High"]

  # Filter data for the specific smoking status group and Low risk
  low_risk_scores <- signature_df_6$PFS_Signature_risk_score[signature_df_6$Smoking_status == status & signature_df_6$Categorized_Risk_Score == "Low"]

  # Check if there are enough data points in both groups
  if (length(high_risk_scores) > 1 && length(low_risk_scores) > 1) {
    # Perform t-test
    t_test_result <- t.test(high_risk_scores, low_risk_scores, var.equal = TRUE) # assuming equal variances

    # Add the results to the data frame
    t_test_results_smoking <- rbind(t_test_results_smoking, data.frame(
      Smoking_Status = status,
      Mean_Difference = mean(high_risk_scores) - mean(low_risk_scores),
      P_Value = t_test_result$p.value,
      Test_Valid = TRUE
    ))
  } else {
    # Not enough data to perform t-test
    t_test_results_smoking <- rbind(t_test_results_smoking, data.frame(
      Smoking_Status = status,
      Mean_Difference = NA,
      P_Value = NA,
      Test_Valid = FALSE
    ))
  }
}

# Print the t-test results table for smoking status
print(t_test_results_smoking)


# Save results_table as a text file
write.table(t_test_results_smoking, file = "MDACC_PFS_signature_genes_risk_scores_Smoking_status_t_test.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# make a copy of signatur_df_6
signature_df_7 <- signature_df_6

# Classify treatment into RT and nonRT
signature_df_7 <- signature_df_7 %>%
  mutate(new_treatment = ifelse(Treatment %in% c("cR_S", "IC_S", "S"), "nonRT",
                                ifelse(Treatment %in% c("IC_S_cRT", "IC_S_RT", "S_cRT",
                                                        "S_IC_S_cRT", "S_RT"), "RT", NA)))

# Correctly update the encoding for treatment to reflect the new classification
signature_df_7$new_treatment_encoded <- ifelse(signature_df_7$new_treatment == "nonRT", 0, 1)

# Convert Treatment to a numeric variable in the dataframe

# Convert Treatment to a numeric variable in the dataframe

table(signature_df_7$new_treatment_encoded)


# Assuming the Categorized_Risk_Score column has been created
# Categorize the risk score based on the median score
median_score <- median(signature_df_7$PFS_Signature_risk_score, na.rm = TRUE)
signature_df_7$Categorized_Risk_Score <- ifelse(signature_df_7$PFS_Signature_risk_score >= median_score, "High", "Low")

# Initialize a data frame to store the t-test results for new treatment groups
t_test_results_treatment <- data.frame(
  Treatment_Group = character(),
  Mean_Difference = numeric(),
  P_Value = numeric(),
  Test_Valid = logical(),
  stringsAsFactors = FALSE
)

# Loop through each treatment group
for (treatment in unique(signature_df_7$new_treatment)) {
  # Skip if treatment is NA
  if (is.na(treatment)) next

  # Filter data for the specific treatment group and High risk
  high_risk_scores <- signature_df_7$PFS_Signature_risk_score[signature_df_7$new_treatment == treatment & signature_df_7$Categorized_Risk_Score == "High"]

  # Filter data for the specific treatment group and Low risk
  low_risk_scores <- signature_df_7$PFS_Signature_risk_score[signature_df_7$new_treatment == treatment & signature_df_7$Categorized_Risk_Score == "Low"]

  # Check if there are enough data points in both groups
  if (length(high_risk_scores) > 1 && length(low_risk_scores) > 1) {
    # Perform t-test
    t_test_result <- t.test(high_risk_scores, low_risk_scores, var.equal = TRUE) # assuming equal variances

    # Add the results to the data frame
    t_test_results_treatment <- rbind(t_test_results_treatment, data.frame(
      Treatment_Group = treatment,
      Mean_Difference = mean(high_risk_scores, na.rm = TRUE) - mean(low_risk_scores, na.rm = TRUE),
      P_Value = t_test_result$p.value,
      Test_Valid = TRUE
    ))
  } else {
    # Not enough data to perform t-test
    t_test_results_treatment <- rbind(t_test_results_treatment, data.frame(
      Treatment_Group = treatment,
      Mean_Difference = NA,
      P_Value = NA,
      Test_Valid = FALSE
    ))
  }
}

# Print the t-test results table for new treatment groups
print(t_test_results_treatment)

# Save results_table as a text file
write.table(t_test_results_treatment, file = "MDACC_PFS_signature_genes_risk_scores_Treatment_t_test.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
