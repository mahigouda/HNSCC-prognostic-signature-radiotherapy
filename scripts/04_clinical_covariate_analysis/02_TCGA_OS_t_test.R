## =============================================================================
## Project : Prognostic gene-expression signature for radiotherapy outcome in HPV-negative HNSCC
## Module  : 04 | Clinical covariates: multivariable Cox, risk groups vs. clinical variables
## Script  : 02_TCGA_OS_t_test.R
## Purpose : TCGA cohort, OS signature: multivariable Cox models (signature genes + clinical
##           covariates), forest plots, median-split high/low risk groups, Kaplan-Meier curves
##           and t-tests of risk groups/risk scores against clinical variables.
## Input   : data/clinical_covariates/TCGA_OS_signature_clinical_df.txt
## Output  : results/04_clinical_covariate_analysis/TCGA_OS_*.txt and plots
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
library(knitr)

## ---- Paths ------------------------------------------------------------------
## Run from the repository root. Inputs are read from data/clinical_covariates/
## (see data/README.md); outputs are written to results/04_clinical_covariate_analysis/.
data_dir <- file.path(getwd(), "data", "clinical_covariates")
out_dir  <- file.path(getwd(), "results", "04_clinical_covariate_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)

## Sub data
file_name <- file.path(data_dir, "TCGA_OS_signature_clinical_df.txt")
signature_df <- read.table(file_name, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = "\t")

# Map the numeric stages to their descriptive names
signature_df$Stages <- mapvalues(signature_df$Stages,
                                 from = c('1', '2', '3', '4', '5'),
                                 to = c('STAGE I', 'STAGE II', 'STAGE III', 'STAGE IVA', 'STAGE IVB'))


# specify the range of columns more safely as follows:
num_cols <- ncol(signature_df)
gene_expression_cols <- signature_df[, 1:(num_cols - 7)]

# Now calculate the median for each gene
signature_df_median <- apply(gene_expression_cols, 2, median)

# categorize age or other continuous variables if needed
# For example, categorizing age into groups
signature_df$age_group <- cut(signature_df$Age,
                              breaks=c(0, 50, 65, 80, Inf), labels=c("0-50", "51-65", "66-80", "80+"))

# Fit the Cox model
# Assuming 'OS_Months' is your survival time and 'Vital_status' is the event indicator
cox_model <- coxph(Surv(OS_Months, Vital_status) ~ ., data = signature_df)

# Get the summary
summary_cox <- summary(cox_model)


# Check the table again to ensure the periods are removed
table(signature_df$Vital_status)

# Convert Vital_status to a Surv object
surv_obj <- Surv(time = signature_df$OS_Months, event = signature_df$Vital_status)

# Fit a Cox proportional hazards model

# Example: Fit a Cox model with fewer variables
cox_model_simple <- coxph(surv_obj ~ Age + Stages + Radiation_therapy, data = signature_df)
summary(cox_model_simple)

# Example: Expanding the model to include gene variables
extended_cox_model <- coxph(surv_obj ~ Age + Stages + Radiation_therapy + ODC1 + TUBB + ETV4 +
                              + NUDT15 + AK4 + PPIB + PLAU + LRRC59 + CTSV + HBEGF + PPIA + MAP4K4 + AXL, data = signature_df)
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
linear_model <- lm(Vital_status ~ age_group + Stages + Radiation_therapy + ODC1 + TUBB + ETV4 + NUDT15 + AK4 + PPIB + PLAU + LRRC59 + CTSV + HBEGF + PPIA + MAP4K4 + AXL, data = signature_df)

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
signature_df <- signature_df[1:291, ]

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


# Making a copy of signature_df
signature_df_2 <- signature_df

signature_df_2$Radiation_therapy <- mapvalues(signature_df_2$Radiation_therapy,
                                              from = c('1', '0'),
                                              to = c('Yes', 'No'))
## Categorize the genes according to risk scores
# Assuming the first 8 columns are gene columns
gene_columns <- names(signature_df_2)[1:13]

# Initialize a list to store results
gene_stats <- list()

for (gene in gene_columns) {
  # Subset data for High and Low risk groups
  high_risk_data <- subset(signature_df_2, risk_group == "High", select = c(gene, "risk_score"))
  low_risk_data <- subset(signature_df_2, risk_group == "Low", select = c(gene, "risk_score"))

  # Calculate statistics - for example, mean risk score
  high_risk_mean <- mean(high_risk_data$risk_score, na.rm = TRUE)
  low_risk_mean <- mean(low_risk_data$risk_score, na.rm = TRUE)

  # Store the results
  gene_stats[[gene]] <- list(High_Risk_Mean = high_risk_mean, Low_Risk_Mean = low_risk_mean)
}

# To view results for a specific gene
gene_stats[["ODC1"]]


# Making a copy of signature_df
signature_df_3 <- signature_df

signature_df_3$Radiation_therapy <- mapvalues(signature_df_3$Radiation_therapy,
                                              from = c('1', '0'),
                                              to = c('Yes', 'No'))

# Assuming the first 13 columns are gene columns
gene_columns <- names(signature_df_3)[1:13]

# Loop through each row (patient) in the dataset
for (i in 1:nrow(signature_df_3)) {
  # Determine the risk group of the patient
  risk_group <- signature_df_3$risk_group[i]

  # Replace gene expression values with risk group label
  signature_df_3[i, gene_columns] <- risk_group
}

# Replace gene expression values based on risk group
signature_df_3 <- signature_df_3 %>%
  mutate(across(all_of(gene_columns), ~ risk_group))


## Linear model for signature_df_3
# Vital status continuous outcome variable in the signature_df_3 dataset
linear_model_signature_df_3 <- lm(Vital_status ~ age_group + Stages + Radiation_therapy + ODC1 + TUBB +
                                    ETV4 + NUDT15 + AK4 + PPIB + PLAU + LRRC59 + CTSV + HBEGF + PPIA +
                                    MAP4K4 + AXL, data = signature_df_3)
# View the summary of the linear model
summary(linear_model_signature_df_3)
forest_model(linear_model_signature_df_3)


# Calculating Variance Inflation Factor
# Check for constant variables
sapply(signature_df_3[, c("ODC1", "TUBB", "ETV4", "NUDT15", "AK4", "PPIB", "PLAU", "LRRC59", "CTSV", "HBEGF", "PPIA", "MAP4K4", "AXL")], function(x) length(unique(x)))

# Check the factor levels
str(signature_df_3$ODC1)

# Simplify the model by removing potentially problematic variables
# Example: removing some gene variables
simplified_model <- lm(Vital_status ~ risk_score + age_group + Stages + Radiation_therapy + ODC1 + TUBB + ETV4, data = signature_df_3)
summary(simplified_model)


# Making a copy of signature_df
signature_df_4 <- signature_df

# Check the unique values and their type in Radiation_therapy column
unique(signature_df_4$Radiation_therapy)


# Assuming the first 13 columns are gene columns
gene_columns <- names(signature_df_4)[1:13]

# Create new columns for each gene to indicate risk category
for (gene in gene_columns) {
  new_col_name <- paste(gene, "Risk", sep = "_")
  signature_df_3 <- signature_df_3 %>%
    mutate(!!new_col_name := ifelse(get(gene) >= some_threshold, "High", "Low"))
}


# Making a copy of signature_df
signature_df_6 <- signature_df_4

# Define a function to categorize continuous variables
categorize_risk <- function(x) {
  median_val <- median(x, na.rm = TRUE)
  ifelse(x >= median_val, "High", "Low")
}

# Apply the function to continuous variables
for (col_name in colnames(signature_df_6)) {
  if (col_name %in% c("Stages", "Radiation_therapy", "age_group", "risk_group", "Stages_Risk", "Radiation_therapy_Risk", "age_group_Risk")) next

  risk_col_name <- paste0(col_name, "_Risk")
  signature_df_6[[risk_col_name]] <- categorize_risk(signature_df_6[[col_name]])
}

signature_df_6 <- signature_df_6 %>%
  mutate(Stages_Risk = ifelse(Stages %in% c("STAGE IVA", "STAGE IVB"), "High", "Low"),
         Radiation_therapy_Risk = ifelse(Radiation_therapy %in% c("Yes", "No"), "High", "Low"),
         age_group_Risk = ifelse(age_group %in% c("0-50", "51-65", "66-80", "80+"), "High", "Low"))

# Summary for Stages
table(signature_df_6$Stages, signature_df_6$Stages_Risk)

# Summary for Radiation Therapy
table(signature_df_6$Radiation_therapy, signature_df_6$Radiation_therapy_Risk)

# Summary for Age Group
table(signature_df_6$age_group, signature_df_6$age_group_Risk)

signature_df_6$Stages_Risk <- as.factor(signature_df_6$Stages_Risk)
signature_df_6$Radiation_therapy_Risk <- as.factor(signature_df_6$Radiation_therapy_Risk)
signature_df_6$age_group_Risk <- as.factor(signature_df_6$age_group_Risk)
signature_df_6$ODC1_Risk <- as.factor(signature_df_6$ODC1_Risk) # Assuming ODC1_Risk is a binary variable

levels(signature_df_6$ODC1_Risk)
levels(signature_df_6$Stages_Risk)
levels(signature_df_6$Radiation_therapy_Risk)
levels(signature_df_6$age_group_Risk)


# Assuming the first 13 columns are gene columns
gene_columns <- names(signature_df_6)[1:13]

# Create new columns for each gene to indicate risk category
for (gene in gene_columns) {
  new_col_name <- paste(gene, "Risk", sep = "_")
  signature_df_6 <- signature_df_6 %>%
    mutate(!!new_col_name := ifelse(get(gene) >= some_threshold, "High", "Low"))
}

# Define a function to categorize continuous variables
categorize_risk <- function(x) {
  median_val <- median(x, na.rm = TRUE)
  ifelse(x >= median_val, "High", "Low")
}

# Loop through columns in signature_df_6
for (col_name in colnames(signature_df_6)) {
  # Skip non-numeric columns or already existing risk categories
  if (col_name %in% c("Age", "Stages", "age_group", "Radiation_therapy", "risk_group") ||
      !is.numeric(signature_df_6[[col_name]])) {
    next
  }

  # Continue with your operations for other columns
  # For example, creating a risk category for each gene
  risk_col_name <- paste0(col_name, "_Risk")
  # Apply categorize_risk to numeric columns
  signature_df_6[[risk_col_name]] <- categorize_risk(signature_df_6[[col_name]])
}


# columns ODC1 to AXL are gene expression values
gene_expression_cols <- 1:13

# Checking for any negative or zero values in gene expression columns
if (any(signature_df_6[, gene_expression_cols] <= 0, na.rm = TRUE)) {
  print("There are non-positive values in gene expression columns")
}

# Apply log2 transformation after adding a small constant
signature_df_6[, gene_expression_cols] <- log2(signature_df_6[, gene_expression_cols] + 1)

# Scale the log-transformed data
signature_df_6[, gene_expression_cols] <- scale(signature_df_6[, gene_expression_cols])

# Calculate the mean for each row
signature_df_6$Gene_Risk_Score <- rowMeans(signature_df_6[, gene_expression_cols], na.rm = TRUE)

# Relabel the column
signature_df_6 <- dplyr::rename(signature_df_6, OS_Signature_risk_score = Gene_Risk_Score)


# Making a copy of signature_df_6
signature_df_7 <- signature_df_6

# Classify age_groups into <65 and >65
signature_df_7 <- signature_df_7 %>%
  mutate(new_age_group = ifelse(age_group %in% c("0-50", "51-65"), "<65",
                                ifelse(age_group %in% c("66-80", "80+"), ">65", NA)))

# Check the first few rows of the modified dataframe
head(signature_df_7)

# Update the encoding for age_group to reflect the new classification
signature_df_7 <- signature_df_7 %>%
  mutate(age_group_encoded = ifelse(new_age_group == "<65", 0, 1)) # 0 for "<65", 1 for ">65"

# Encode Sex and Stages as numeric

# Update your categories to reflect the new age group classification
sexes <- c("F", "M")
# Define your categories
age_groups <- c("<65", ">65")
stages <- c("STAGE I", "STAGE II", "STAGE III", "STAGE IVA", "STAGE IVB")
Radiation_therapies <- c("Yes", "No")
risk_groups <- c("High", "Low")


# Encode Sex and Stages as numeric
signature_df_7$Sex_encoded <- ifelse(signature_df_7$Sex == "F", 0, 1)  # Example: 0 for Female, 1 for Male
signature_df_7$Stages_encoded <- as.numeric(as.factor(signature_df_7$Stages))


# Example for one group: Age < 65 and High Risk
subset_data <- subset(signature_df_7, new_age_group == "<65" & risk_group == "High")


# Categorizing OS_Signature_risk_score
median_score <- median(signature_df_7$OS_Signature_risk_score, na.rm = TRUE)
signature_df_7$Categorized_Risk_Score <- ifelse(signature_df_7$OS_Signature_risk_score >= median_score, "High", "Low")


# Initialize a data frame to store the results
t_test_results <- data.frame(
  Age_Group = character(),
  Mean_Difference = numeric(),
  P_Value = numeric(),
  Test_Valid = logical(),
  Count = integer(),
  stringsAsFactors = FALSE
)

# Loop through each age group
for (age in c("<65", ">65")) {
  # Filter data for the specific age group and High risk
  high_risk_scores <- signature_df_7$OS_Signature_risk_score[signature_df_7$new_age_group == age & signature_df_7$Categorized_Risk_Score == "High"]

  # Filter data for the specific age group and Low risk
  low_risk_scores <- signature_df_7$OS_Signature_risk_score[signature_df_7$new_age_group == age & signature_df_7$Categorized_Risk_Score == "Low"]

  # Check if there are enough data points in each group
  if (length(high_risk_scores) >= 2 && length(low_risk_scores) >= 2) {
    # Perform t-test
    t_test_result <- t.test(high_risk_scores, low_risk_scores)

    # Add the results to the data frame
    t_test_results <- rbind(t_test_results, data.frame(
      Age_Group = age,
      Mean_Difference = t_test_result$estimate[1] - t_test_result$estimate[2],
      P_Value = t_test_result$p.value,
      Test_Valid = TRUE
    ))
  } else {
    # Not enough data to perform t-test
    t_test_results <- rbind(t_test_results, data.frame(
      Age_Group = age,
      Mean_Difference = NA,
      P_Value = NA,
      Count = nrow(subset_data),
      Test_Valid = FALSE
    ))
  }
}

# Print the t-test results table
print(t_test_results)


# Save results_table as a text file
write.table(t_test_results, file = "TCGA_OS_signature_genes_risk_scores_age_group_t_test.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# make a copy of signatur_df_2
signature_df_8 <- signature_df_7

# Classify age_groups into <65 and >65
signature_df_8 <- signature_df_8 %>%
  mutate(new_age_group = ifelse(age_group %in% c("0-50", "51-65"), "<65",
                                ifelse(age_group %in% c("66-80", "80+"), ">65", NA)))

# Classify stages
signature_df_8 <- signature_df_8 %>%
  mutate(new_stages = ifelse(Stages %in% c("STAGE I", "STAGE II"), "Stage I/II",
                             ifelse(Stages %in% c("STAGE III", "STAGE IVA", "STAGE IVB"), "Stage III/IV", NA)))


# Update the encoding for age_group to reflect the new classification
signature_df_8 <- signature_df_8 %>%
  mutate(age_group_encoded = ifelse(new_age_group == "<65", 0, 1)) # 0 for "<65", 1 for ">65"

# Encode Sex and Stages as numeric
signature_df_8$Sex_encoded <- ifelse(signature_df_8$Sex == "F", 0, 1)  # Example: 0 for Female, 1 for Male
signature_df_8$new_stages_encoded <- as.numeric(as.factor(signature_df_8$new_stages))


# Update your categories to reflect the new age group classification
age_groups <- c("<65", ">65")
stages <- c("Stage I/II", "Stage III/IV")
sexes <- c("F", "M")
risk_groups <- c("High", "Low")
Radiation_therapies <- c("Yes", "No")


print(table(signature_df_8$new_stages))
print(table(signature_df_8$risk_group))

# Categorizing OS_Signature_risk_score
median_score <- median(signature_df_8$OS_Signature_risk_score, na.rm = TRUE)
signature_df_8$Categorized_Risk_Score <- ifelse(signature_df_8$OS_Signature_risk_score >= median_score, "High", "Low")


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
  high_risk_scores <- signature_df_8$OS_Signature_risk_score[signature_df_8$new_stages == stage & signature_df_8$Categorized_Risk_Score == "High"]

  # Filter data for the specific stage group and Low risk
  low_risk_scores <- signature_df_8$OS_Signature_risk_score[signature_df_8$new_stages == stage & signature_df_8$Categorized_Risk_Score == "Low"]

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
write.table(t_test_results_stages, file = "TCGA_OS_signature_genes_risk_scores_stages_t_test.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# make a copy of signatur_df_8
signature_df_9 <- signature_df_8


# Initialize a data frame to store the t-test results for sex
t_test_results_sex <- data.frame(
  Sex_Group = character(),
  Mean_Difference = numeric(),
  P_Value = numeric(),
  Test_Valid = logical(),
  stringsAsFactors = FALSE
)

# Loop through each sex group
for (sex_code in unique(signature_df_9$Sex_encoded)) {
  sex_label <- ifelse(sex_code == 0, "Female", "Male")

  # Filter data for the specific sex group and High risk
  high_risk_scores <- signature_df_9$OS_Signature_risk_score[signature_df_9$Sex_encoded == sex_code & signature_df_9$Categorized_Risk_Score == "High"]

  # Filter data for the specific sex group and Low risk
  low_risk_scores <- signature_df_9$OS_Signature_risk_score[signature_df_9$Sex_encoded == sex_code & signature_df_9$Categorized_Risk_Score == "Low"]

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
write.table(t_test_results_sex, file = "TCGA_OS_signature_genes_risk_scores_Sex_t_test.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# make a copy of signatur_df_9
signature_df_10 <- signature_df_9

# Classify treatment into RT and nonRT
signature_df_10 <- signature_df_10 %>%
  mutate(new_treatment = ifelse(Radiation_therapy %in% c("0"), "nonRT",
                                ifelse(Radiation_therapy %in% c("1"), "RT", NA)))

# Correctly update the encoding for treatment to reflect the new classification
signature_df_10$new_treatment_encoded <- ifelse(signature_df_10$new_treatment == "nonRT", 0, 1)

table(signature_df_10$new_treatment_encoded)

# Initialize a data frame to store the t-test results for new_treatment_encoded
t_test_results_treatment <- data.frame(
  Treatment_Group = character(),
  Mean_Difference = numeric(),
  P_Value = numeric(),
  Test_Valid = logical(),
  stringsAsFactors = FALSE
)

# Loop through each treatment group
for (treatment_code in unique(signature_df_10$new_treatment_encoded)) {
  treatment_label <- ifelse(treatment_code == 0, "nonRT", "RT")

  # Filter data for the specific treatment group and High risk
  high_risk_scores <- signature_df_10$OS_Signature_risk_score[signature_df_10$new_treatment_encoded == treatment_code & signature_df_10$risk_group == "High"]

  # Filter data for the specific treatment group and Low risk
  low_risk_scores <- signature_df_10$OS_Signature_risk_score[signature_df_10$new_treatment_encoded == treatment_code & signature_df_10$risk_group == "Low"]

  # Perform t-test if there are enough data points in both groups
  if (length(high_risk_scores) > 1 && length(low_risk_scores) > 1) {
    t_test_result <- t.test(high_risk_scores, low_risk_scores, var.equal = TRUE)

    # Add the results to the data frame
    t_test_results_treatment <- rbind(t_test_results_treatment, data.frame(
      Treatment_Group = treatment_label,
      Mean_Difference = mean(high_risk_scores) - mean(low_risk_scores),
      P_Value = t_test_result$p.value,
      Test_Valid = TRUE
    ))
  } else {
    # Not enough data to perform t-test
    t_test_results_treatment <- rbind(t_test_results_treatment, data.frame(
      Treatment_Group = treatment_label,
      Mean_Difference = NA,
      P_Value = NA,
      Test_Valid = FALSE
    ))
  }
}

# Print the t-test results table for new_treatment_encoded
print(t_test_results_treatment)


# Save results_table as a text file
write.table(t_test_results_treatment, file = "TCGA_OS_signature_genes_risk_scores_Treatment_t_test.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# make a copy of signatur_df_10
signature_df_11 <- signature_df_10

unique(signature_df_11$anatomic_site)

# Convert anatomic_site to a factor
signature_df_11$anatomic_site_encoded <- as.numeric(factor(signature_df_11$anatomic_site))


# Initialize a data frame to store the t-test results for anatomic_site
t_test_results_anatomic_site <- data.frame(
  Anatomic_Site = character(),
  Mean_Difference = numeric(),
  P_Value = numeric(),
  Test_Valid = logical(),
  stringsAsFactors = FALSE
)

# Get the levels of the encoded anatomic_site
anatomic_sites <- levels(factor(signature_df_11$anatomic_site))

# Loop through each anatomic site
for (site in anatomic_sites) {
  # Filter data for the specific anatomic site and High risk
  high_risk_scores <- signature_df_11$OS_Signature_risk_score[signature_df_11$anatomic_site == site & signature_df_11$risk_group == "High"]

  # Filter data for the specific anatomic site and Low risk
  low_risk_scores <- signature_df_11$OS_Signature_risk_score[signature_df_11$anatomic_site == site & signature_df_11$risk_group == "Low"]

  # Perform t-test if there are enough data points in both groups
  if (length(high_risk_scores) > 1 && length(low_risk_scores) > 1) {
    t_test_result <- t.test(high_risk_scores, low_risk_scores, var.equal = TRUE)

    # Add the results to the data frame
    t_test_results_anatomic_site <- rbind(t_test_results_anatomic_site, data.frame(
      Anatomic_Site = site,
      Mean_Difference = mean(high_risk_scores) - mean(low_risk_scores),
      P_Value = t_test_result$p.value,
      Test_Valid = TRUE
    ))
  } else {
    # Not enough data to perform t-test
    t_test_results_anatomic_site <- rbind(t_test_results_anatomic_site, data.frame(
      Anatomic_Site = site,
      Mean_Difference = NA,
      P_Value = NA,
      Test_Valid = FALSE
    ))
  }
}

# Print the t-test results table for anatomic_site
print(t_test_results_anatomic_site)


# Save results_table as a text file
write.table(t_test_results_anatomic_site, file = "TCGA_OS_signature_genes_risk_scores_anatomic_site_t_test.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# make a copy of signatur_df_11
signature_df_12 <- signature_df_11

unique(signature_df_12$new_tumor_event)

# Initialize a data frame to store the t-test results for new_tumor_event
t_test_results_new_tumor_event <- data.frame(
  New_Tumor_Event_Group = character(),
  Mean_Difference = numeric(),
  P_Value = numeric(),
  Test_Valid = logical(),
  stringsAsFactors = FALSE
)

# Loop through each category of new_tumor_event
for (event in unique(signature_df_12$new_tumor_event)) {
  event_label <- ifelse(event == 1, "Event Occurred", "No Event")

  # Filter data for the specific new_tumor_event and High risk
  high_risk_scores <- signature_df_12$OS_Signature_risk_score[signature_df_12$new_tumor_event == event & signature_df_12$risk_group == "High"]

  # Filter data for the specific new_tumor_event and Low risk
  low_risk_scores <- signature_df_12$OS_Signature_risk_score[signature_df_12$new_tumor_event == event & signature_df_12$risk_group == "Low"]

  # Perform t-test if there are enough data points in both groups
  if (length(high_risk_scores) > 1 && length(low_risk_scores) > 1) {
    t_test_result <- t.test(high_risk_scores, low_risk_scores, var.equal = TRUE)

    # Add the results to the data frame
    t_test_results_new_tumor_event <- rbind(t_test_results_new_tumor_event, data.frame(
      New_Tumor_Event_Group = event_label,
      Mean_Difference = mean(high_risk_scores) - mean(low_risk_scores),
      P_Value = t_test_result$p.value,
      Test_Valid = TRUE
    ))
  } else {
    # Not enough data to perform t-test
    t_test_results_new_tumor_event <- rbind(t_test_results_new_tumor_event, data.frame(
      New_Tumor_Event_Group = event_label,
      Mean_Difference = NA,
      P_Value = NA,
      Test_Valid = FALSE
    ))
  }
}

# Print the t-test results table for new_tumor_event
print(t_test_results_new_tumor_event)


# Save results_table as a text file
write.table(t_test_results_new_tumor_event, file = "TCGA_OS_signature_genes_risk_scores_primary_vs_relapse_t_test.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
