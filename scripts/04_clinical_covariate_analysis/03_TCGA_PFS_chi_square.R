## =============================================================================
## Project : Prognostic gene-expression signature for radiotherapy outcome in HPV-negative HNSCC
## Module  : 04 | Clinical covariates: multivariable Cox, risk groups vs. clinical variables
## Script  : 03_TCGA_PFS_chi_square.R
## Purpose : TCGA cohort, PFS signature: multivariable Cox models (signature genes + clinical
##           covariates), forest plots, median-split high/low risk groups, Kaplan-Meier curves
##           and chi-square tests of risk groups/risk scores against clinical variables.
## Input   : data/clinical_covariates/TCGA_PFS_signature_clinical_df.txt
## Output  : results/04_clinical_covariate_analysis/TCGA_PFS_*.txt and plots
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
file_name <- file.path(data_dir, "TCGA_PFS_signature_clinical_df.txt")
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
# Progress_free_survival_Months is survival time and 'Vital_status' is the event indicator

# Get the summary
summary_cox <- summary(cox_model)


# Check the table again to ensure the periods are removed
table(signature_df$Vital_status)

# Convert Vital_status to a Surv object
surv_obj <- Surv(time = signature_df$Progress_free_survival_Months, event = signature_df$Vital_status)


# Example: Fit a Cox model with fewer variables
cox_model_simple <- coxph(surv_obj ~ Age + Sex + Stages + Radiation_therapy, data = signature_df)
summary(cox_model_simple)

# Example: Expanding the model to include gene variables
extended_cox_model <- coxph(surv_obj ~ Age + Sex + Stages + Radiation_therapy + ADAM8 + RAC2 +
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


## Median High Vs low risk
# Assuming 'extended_cox_model' is your final Cox model
# Calculate risk scores
risk_scores <- predict(extended_cox_model, type = "risk")

# Check the number of rows in the dataset and the length of risk_scores
nrow(signature_df)
length(risk_scores)

# Assuming the first 251 rows of signature_df were used in the Cox model
signature_df <- signature_df[1:289, ]

# Exclude observations from signature_df that were not used in the Cox model
signature_df <- signature_df[!is.na(risk_scores), ]

# Add risk scores to the dataset
signature_df$risk_score <- risk_scores

# Median split for high vs low risk
signature_df$risk_group <- ifelse(signature_df$risk_score > median(signature_df$risk_score), "High", "Low")

surv_obj <- Surv(signature_df$Progress_free_survival_Months, signature_df$Vital_status)
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
gene_columns <- names(signature_df_2)[1:10]

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
gene_stats[["ADAM8"]]


# Making a copy of signature_df
signature_df_3 <- signature_df

signature_df_3$Radiation_therapy <- mapvalues(signature_df_3$Radiation_therapy,
                                              from = c('1', '0'),
                                              to = c('Yes', 'No'))

# Assuming the first 13 columns are gene columns
gene_columns <- names(signature_df_3)[1:10]

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


# Calculating Variance Inflation Factor
# Check for constant variables
sapply(signature_df_3[, c("ADAM8", "RAC2", "CSF2", "PPP1R18", "HBEGF", "ITGA5", "LRRC59", "PLAU", "FOSL1", "PLAUR")], function(x) length(unique(x)))

# Check the factor levels
str(signature_df_3$ADAM8)

# Simplify the model by removing potentially problematic variables
# Example: removing some gene variables
simplified_model <- lm(Vital_status ~ risk_score + age_group + Sex + Stages + Radiation_therapy + ADAM8 + RAC2 +
                         CSF2, data = signature_df_3)
summary(simplified_model)


# Making a copy of signature_df
signature_df_4 <- signature_df

signature_df_4$Radiation_therapy <- mapvalues(signature_df_3$Radiation_therapy,
                                              from = c('1', '0'),
                                              to = c('Yes', 'No'))

# Assuming the first 13 columns are gene columns
gene_columns <- names(signature_df_4)[1:10]

# Create new columns for each gene to indicate risk category
for (gene in gene_columns) {
  new_col_name <- paste(gene, "Risk", sep = "_")
  signature_df_4 <- signature_df_4 %>%
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
  if (col_name %in% c("Stages", "Sex", "Radiation_therapy", "age_group", "risk_group", "Stages_Risk", "Radiation_therapy_Risk", "age_group_Risk")) next

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
signature_df_6$ADAM8_Risk <- as.factor(signature_df_6$ADAM8_Risk)

levels(signature_df_6$ADAM8_Risk)
levels(signature_df_6$Stages_Risk)
levels(signature_df_6$Radiation_therapy_Risk)
levels(signature_df_6$age_group_Risk)


# Assuming the first 10 columns are gene columns
gene_columns <- names(signature_df_6)[1:10]

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


# gene expression values
gene_expression_cols <- 1:10

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
signature_df_6 <- dplyr::rename(signature_df_6, PFS_Signature_risk_score = Gene_Risk_Score)


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


# Check the data type of OS_Signature_risk_score
if ("PFS_Signature_risk_score" %in% names(signature_df_7)) {
  data_type <- class(signature_df_7$PFS_Signature_risk_score)
  print(paste("Data type of PFS_Signature_risk_score:", data_type))
} else {
  print("PFS_Signature_risk_score does not exist in the dataset.")
}


# Categorizing PFS_Signature_risk_score
median_score <- median(signature_df_7$PFS_Signature_risk_score, na.rm = TRUE)
signature_df_7$Categorized_Risk_Score <- ifelse(signature_df_7$PFS_Signature_risk_score >= median_score, "High", "Low")

# Convert to factor
signature_df_7$Categorized_Risk_Score <- factor(signature_df_7$Categorized_Risk_Score)

# Create a contingency table
contingency_table <- table(signature_df_7$new_age_group, signature_df_7$Categorized_Risk_Score)

# Print the contingency table for inspection
print("Contingency Table for Categorized PFS_Signature_risk_score and new_age_group")
print(contingency_table)

# Perform Chi-square test
if (all(contingency_table >= 5)) {
  test_result <- chisq.test(contingency_table)
  print(test_result)
} else {
  print("Chi-square test not valid due to low expected counts in some cells.")
}


# Loop through each combination of age group and risk group
for (age in c("<65", ">65")) {
  for (risk in c("High", "Low")) {
    # Filter data for the specific combination
    subset_data <- subset(signature_df_7, new_age_group == age & risk_group == risk)

    # Create a contingency table with Categorized_Risk_Score
    contingency_table <- table(subset_data$Categorized_Risk_Score)

    # Print the contingency table for inspection
    print(paste("Contingency Table for", age, "Age Group and", risk, "Risk Group"))
    print(contingency_table)

    # Perform Chi-square test if valid
    if (all(contingency_table >= 5)) {
      test_result <- chisq.test(contingency_table)
      print(test_result)
    } else {
      print("Chi-square test not valid due to low expected counts in some cells.")
    }
  }
}

# Initialize a data frame to store the results
chi_square_results <- data.frame(
  Age_Group = character(),
  Risk_Group = character(),
  Chi_Square_Statistic = numeric(),
  P_Value = numeric(),
  Count = integer(),
  Test_Valid = logical(),
  stringsAsFactors = FALSE
)

# Loop through each combination of age group and risk group
for (age in c("<65", ">65")) {
  for (risk in c("High", "Low")) {
    # Filter data for the specific combination
    subset_data <- subset(signature_df_7, new_age_group == age & risk_group == risk)

    # Create a contingency table with Categorized_Risk_Score
    contingency_table <- table(subset_data$Categorized_Risk_Score)

    if (all(contingency_table >= 5)) {
      # Perform Chi-square test
      test_result <- chisq.test(contingency_table)
      chi_sq_statistic <- test_result$statistic
      p_value <- test_result$p.value
      test_valid <- TRUE
    } else {
      # Set NA values if the test is not valid
      chi_sq_statistic <- NA
      p_value <- NA
      test_valid <- FALSE
    }

    # Add the results to the data frame
    chi_square_results <- rbind(chi_square_results, data.frame(
      Age_Group = age,
      Risk_Group = risk,
      Chi_Square_Statistic = chi_sq_statistic,
      P_Value = p_value,
      Count = nrow(subset_data),
      Test_Valid = test_valid
    ))
  }
}

# Print the results table
print(chi_square_results)

# Save results_table as a text file
write.table(chi_square_results, file = "TCGA_PFS_signature_genes_risk_scores_age_group.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# Create a combined contingency table for new_age_group and risk_group
combined_contingency_table_age <- table(signature_df_7$new_age_group, signature_df_7$risk_group)

# Perform Chi-square test on the combined contingency table
combined_test_result_age <- chisq.test(combined_contingency_table_age)

# Store the results in a data frame
combined_chi_square_results_age <- data.frame(
  Chi_Square_Statistic = combined_test_result_age$statistic,
  P_Value = combined_test_result_age$p.value,
  Test_Valid = !any(combined_contingency_table_age < 5) # Test is valid if all expected frequencies are >= 5
)

# Print the combined results
print(combined_chi_square_results_age)

# Save results_table as a text file
write.table(combined_chi_square_results_age, file = "TCGA_PFS_signature_genes_risk_scores_age_group_chi_square.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# Create a combined contingency table for new_age_group and risk_group
combined_contingency_table_age <- table(signature_df_7$new_age_group, signature_df_7$risk_group)

# Perform Chi-square test on the combined contingency table
combined_test_result_age <- chisq.test(combined_contingency_table_age)

# Store the results in a data frame
combined_chi_square_results_age <- data.frame(
  Chi_Square_Statistic = combined_test_result_age$statistic,
  P_Value = combined_test_result_age$p.value,
  Test_Valid = !any(combined_contingency_table_age < 5) # Test is valid if all expected frequencies are >= 5
)

# Print the combined results
print(combined_chi_square_results_age)


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
median_score <- median(signature_df_8$PFS_Signature_risk_score, na.rm = TRUE)
signature_df_8$Categorized_Risk_Score <- ifelse(signature_df_8$PFS_Signature_risk_score >= median_score, "High", "Low")

# Convert to factor
signature_df_8$Categorized_Risk_Score <- factor(signature_df_8$Categorized_Risk_Score)

# Create a contingency table
contingency_table <- table(signature_df_8$new_stages, signature_df_8$Categorized_Risk_Score)


# Initialize a data frame to store the results
chi_square_results_stages <- data.frame(
  Stage_Group = character(),
  Risk_Group = character(),
  Chi_Square_Statistic = numeric(),
  P_Value = numeric(),
  Test_Valid = logical(),
  Count = integer(),
  stringsAsFactors = FALSE
)

# Loop through each combination of stage group and risk group
for (stage in c("Stage I/II", "Stage III/IV")) {
  for (risk in c("High", "Low")) {
    # Filter data for the specific combination
    subset_data <- subset(signature_df_8, new_stages == stage & risk_group == risk)

    # Create a contingency table with Categorized_Risk_Score
    contingency_table <- table(subset_data$Categorized_Risk_Score)

    if (all(contingency_table >= 5)) {
      # Perform Chi-square test
      test_result <- chisq.test(contingency_table)
      chi_sq_statistic <- test_result$statistic
      p_value <- test_result$p.value
      test_valid <- TRUE
    } else {
      # Set NA values if the test is not valid
      chi_sq_statistic <- NA
      p_value <- NA
      test_valid <- FALSE
    }

    # Add the results to the data frame
    chi_square_results_stages <- rbind(chi_square_results_stages, data.frame(
      Stage_Group = stage,
      Risk_Group = risk,
      Chi_Square_Statistic = chi_sq_statistic,
      P_Value = p_value,
      Count = nrow(subset_data),
      Test_Valid = test_valid
    ))
  }
}

# Print the results table
print(chi_square_results_stages)


# Save results_table as a text file
write.table(chi_square_results_stages, file = "TCGA_PFS_signature_genes_risk_scores_stages.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# Create a combined contingency table for new_stages and risk_group
combined_contingency_table_stages <- table(signature_df_8$new_stages, signature_df_8$risk_group)

# Perform Chi-square test on the combined contingency table
combined_test_result_stages <- chisq.test(combined_contingency_table_stages)

# Store the results in a data frame
combined_chi_square_results_stages <- data.frame(
  Chi_Square_Statistic = combined_test_result_stages$statistic,
  P_Value = combined_test_result_stages$p.value,
  Test_Valid = !any(combined_contingency_table_stages < 5) # Test is valid if all expected frequencies are >= 5
)

# Print the combined results
print(combined_chi_square_results_stages)

# Save results_table as a text file
write.table(combined_chi_square_results_stages, file = "TCGA_PFS_signature_genes_risk_scores_stages_chi_square.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# Create a combined contingency table for new_stages and risk_group
combined_contingency_table_stages <- table(signature_df_8$new_stages, signature_df_8$risk_group)

# Perform Chi-square test on the combined contingency table
combined_test_result_stages <- chisq.test(combined_contingency_table_stages)

# Store the results in a data frame
combined_chi_square_results_stages <- data.frame(
  Chi_Square_Statistic = combined_test_result_stages$statistic,
  P_Value = combined_test_result_stages$p.value,
  Test_Valid = !any(combined_contingency_table_stages < 5) # Test is valid if all expected frequencies are >= 5
)

# Print the combined results
print(combined_chi_square_results_stages)

# Save results_table as a text file
write.table(combined_chi_square_results_stages, file = "TCGA_PFS_signature_genes_risk_scores_stages_chi_square.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# make a copy of signatur_df_8
signature_df_9 <- signature_df_8


# Initialize a data frame to store the results
chi_square_results_sex <- data.frame(
  Sex_Group = character(),
  Risk_Group = character(),
  Chi_Square_Statistic = numeric(),
  P_Value = numeric(),
  Test_Valid = logical(),
  Count = integer(),
  stringsAsFactors = FALSE
)

# Loop through each combination of sex group and risk group
for (sex in unique(signature_df_9$Sex_encoded)) { # Use unique to get all sex groups
  for (risk in c("High", "Low")) {
    # Filter data for the specific combination
    subset_data <- subset(signature_df_9, Sex_encoded == sex & risk_group == risk)

    # Create a contingency table with Categorized_Risk_Score
    contingency_table <- table(subset_data$Categorized_Risk_Score)

    if (all(contingency_table >= 5)) {
      # Perform Chi-square test
      test_result <- chisq.test(contingency_table)
      chi_sq_statistic <- test_result$statistic
      p_value <- test_result$p.value
      test_valid <- TRUE
    } else {
      # Set NA values if the test is not valid
      chi_sq_statistic <- NA
      p_value <- NA
      test_valid <- FALSE
    }

    # Add the results to the data frame
    chi_square_results_sex <- rbind(chi_square_results_sex, data.frame(
      Sex_Group = ifelse(sex == 0, "Female", "Male"),
      Risk_Group = risk,
      Chi_Square_Statistic = chi_sq_statistic,
      P_Value = p_value,
      Count = nrow(subset_data),
      Test_Valid = test_valid
    ))
  }
}

# Print the results table
print(chi_square_results_sex)


# Save results_table as a text file
write.table(chi_square_results_sex, file = "TCGA_PFS_signature_genes_risk_scores_Sex.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# Create a combined contingency table for Sex_encoded and risk_group
combined_contingency_table_sex <- table(signature_df_9$Sex_encoded, signature_df_9$risk_group)

# Perform Chi-square test on the combined contingency table
combined_test_result_sex <- chisq.test(combined_contingency_table_sex)

# Store the results in a data frame
combined_chi_square_results_sex <- data.frame(
  Chi_Square_Statistic = combined_test_result_sex$statistic,
  P_Value = combined_test_result_sex$p.value,
  Test_Valid = !any(combined_contingency_table_sex < 5) # Test is valid if all expected frequencies are >= 5
)

# Print the combined results
print(combined_chi_square_results_sex)

# Save results_table as a text file
write.table(combined_chi_square_results_sex, file = "TCGA_PFS_signature_genes_risk_scores_Sex_chi_square.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# make a copy of signatur_df_9
signature_df_10 <- signature_df_9

# Classify treatment into RT and nonRT
signature_df_10 <- signature_df_10 %>%
  mutate(new_treatment = ifelse(Radiation_therapy %in% c("No"), "nonRT",
                                ifelse(Radiation_therapy %in% c("Yes"), "RT", NA)))

# Correctly update the encoding for treatment to reflect the new classification
signature_df_10$new_treatment_encoded <- ifelse(signature_df_10$new_treatment == "nonRT", 0, 1)

# Convert Treatment to a numeric variable in the dataframe

# Convert Treatment to a numeric variable in the dataframe

table(signature_df_10$new_treatment_encoded)

# Initialize a data frame to store the results
chi_square_results_treatment <- data.frame(
  Treatment_Group = character(),
  Risk_Group = character(),
  Chi_Square_Statistic = numeric(),
  P_Value = numeric(),
  Test_Valid = logical(),
  Count = integer(),
  stringsAsFactors = FALSE
)

# Loop through each combination of treatment group and risk group
for (treatment in unique(signature_df_10$new_treatment_encoded)) {
  for (risk in c("High", "Low")) {
    # Filter data for the specific combination
    subset_data <- subset(signature_df_10, new_treatment_encoded == treatment & risk_group == risk)

    # Create a contingency table with Categorized_Risk_Score
    contingency_table <- table(subset_data$Categorized_Risk_Score)

    if (all(contingency_table >= 5)) {
      # Perform Chi-square test
      test_result <- chisq.test(contingency_table)
      chi_sq_statistic <- test_result$statistic
      p_value <- test_result$p.value
      test_valid <- TRUE
    } else {
      # Set NA values if the test is not valid
      chi_sq_statistic <- NA
      p_value <- NA
      test_valid <- FALSE
    }

    # Add the results to the data frame
    chi_square_results_treatment <- rbind(chi_square_results_treatment, data.frame(
      Treatment_Group = ifelse(treatment == 0, "nonRT", "RT"),
      Risk_Group = risk,
      Chi_Square_Statistic = chi_sq_statistic,
      P_Value = p_value,
      Count = nrow(subset_data),
      Test_Valid = test_valid
    ))
  }
}

# Print the results table
print(chi_square_results_treatment)


# Save results_table as a text file
write.table(chi_square_results_treatment, file = "TCGA_PFS_signature_genes_risk_scores_Treatment.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# Create a combined contingency table for new_treatment_encoded and risk_group
combined_contingency_table <- table(signature_df_10$new_treatment_encoded, signature_df_10$risk_group)

# Perform Chi-square test on the combined contingency table
combined_test_result <- chisq.test(combined_contingency_table)

# Store the results in a data frame
combined_chi_square_results <- data.frame(
  Chi_Square_Statistic = combined_test_result$statistic,
  P_Value = combined_test_result$p.value,
  Test_Valid = !any(combined_contingency_table < 5) # Test is valid if all expected frequencies are >= 5
)

# Print the combined results
print(combined_chi_square_results)

# Save results_table as a text file
write.table(combined_chi_square_results, file = "TCGA_PFS_signature_genes_risk_scores_Treatment_chi_square.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# Create a combined contingency table for new_treatment_encoded and risk_group
combined_contingency_table <- table(signature_df_10$new_treatment_encoded, signature_df_10$risk_group)

# Perform Chi-square test on the combined contingency table
combined_test_result <- chisq.test(combined_contingency_table)

# Store the results in a data frame
combined_chi_square_results <- data.frame(
  Chi_Square_Statistic = combined_test_result$statistic,
  P_Value = combined_test_result$p.value,
  Test_Valid = !any(combined_contingency_table < 5) # Test is valid if all expected frequencies are >= 5
)

# Print the combined results
print(combined_chi_square_results)


# make a copy of signatur_df_10
signature_df_11 <- signature_df_10

unique(signature_df_11$anatomic_site)

# Convert anatomic_site to a factor
signature_df_11$anatomic_site_encoded <- as.numeric(factor(signature_df_11$anatomic_site))

# Initialize a data frame to store the results for anatomic_site
chi_square_results_anatomic_site <- data.frame(
  Anatomic_Site = character(),
  Risk_Group = character(),
  Chi_Square_Statistic = numeric(),
  P_Value = numeric(),
  Test_Valid = logical(),
  Count = integer(),
  stringsAsFactors = FALSE
)

# Get the levels of the encoded anatomic_site
anatomic_sites <- levels(factor(signature_df_11$anatomic_site))

# Loop through each combination of anatomic site and risk group
for (site in anatomic_sites) {
  for (risk in c("High", "Low")) {
    # Filter data for the specific combination
    subset_data <- subset(signature_df_11, anatomic_site == site & risk_group == risk)

    # Create a contingency table with Categorized_Risk_Score
    contingency_table <- table(subset_data$Categorized_Risk_Score)

    if (all(contingency_table >= 5)) {
      # Perform Chi-square test
      test_result <- chisq.test(contingency_table)
      chi_sq_statistic <- test_result$statistic
      p_value <- test_result$p.value
      test_valid <- TRUE
    } else {
      # Set NA values if the test is not valid
      chi_sq_statistic <- NA
      p_value <- NA
      test_valid <- FALSE
    }

    # Add the results to the data frame
    chi_square_results_anatomic_site <- rbind(chi_square_results_anatomic_site, data.frame(
      Anatomic_Site = site,
      Risk_Group = risk,
      Chi_Square_Statistic = chi_sq_statistic,
      P_Value = p_value,
      Count = nrow(subset_data),
      Test_Valid = test_valid
    ))
  }
}

# Print the results table
print(chi_square_results_anatomic_site)

# Save results_table as a text file
write.table(chi_square_results_anatomic_site, file = "TCGA_PFS_signature_genes_risk_scores_anatomic_site.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# Create a combined contingency table for new_treatment_encoded and risk_group
combined_contingency_table <- table(signature_df_11$anatomic_site_encoded, signature_df_11$risk_group)

# Perform Chi-square test on the combined contingency table
combined_test_result <- chisq.test(combined_contingency_table)

# Check if all expected frequencies are >= 5
test_valid <- all(combined_test_result$expected >= 5)

# Store the results in a data frame
combined_chi_square_results <- data.frame(
  Chi_Square_Statistic = combined_test_result$statistic,
  P_Value = combined_test_result$p.value,
  Test_Valid = test_valid
)

# Print the combined results
print(combined_chi_square_results)

# Save results_table as a text file
write.table(combined_chi_square_results, file = "TCGA_PFS_signature_genes_risk_scores_anatomic_site_chi_square.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# make a copy of signatur_df_11
signature_df_12 <- signature_df_11

unique(signature_df_12$new_tumor_event)


# Initialize a data frame to store the results for new_tumor_event
chi_square_results_new_tumor_event <- data.frame(
  New_Tumor_Event = character(),
  Risk_Group = character(),
  Chi_Square_Statistic = numeric(),
  P_Value = numeric(),
  Test_Valid = logical(),
  Count = integer(),
  stringsAsFactors = FALSE
)

# Loop through each risk group
for (risk in c("High", "Low")) {
  for (event in unique(signature_df_12$new_tumor_event)) {
    # Filter data for the specific combination
    subset_data <- subset(signature_df_12, new_tumor_event == event & risk_group == risk)

    # Create a contingency table with Categorized_Risk_Score
    contingency_table <- table(subset_data$Categorized_Risk_Score)

    if (all(contingency_table >= 5)) {
      # Perform Chi-square test
      test_result <- chisq.test(contingency_table)
      chi_sq_statistic <- test_result$statistic
      p_value <- test_result$p.value
      test_valid <- TRUE
    } else {
      # Set NA values if the test is not valid
      chi_sq_statistic <- NA
      p_value <- NA
      test_valid <- FALSE
    }

    # Add the results to the data frame
    chi_square_results_new_tumor_event <- rbind(chi_square_results_new_tumor_event, data.frame(
      New_Tumor_Event = ifelse(event == 1, "YES", "NO"),
      Risk_Group = risk,
      Chi_Square_Statistic = chi_sq_statistic,
      P_Value = p_value,
      Count = nrow(subset_data),
      Test_Valid = test_valid
    ))
  }
}

# Print the results table
print(chi_square_results_new_tumor_event)


# Save results_table as a text file
write.table(chi_square_results_new_tumor_event, file = "TCGA_PFS_signature_genes_risk_scores_primary_vs_relapse.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# Create a combined contingency table for new_tumor_event and risk_group
combined_contingency_table <- table(signature_df_12$new_tumor_event, signature_df_12$risk_group)

# Perform Chi-square test on the combined contingency table
combined_test_result <- chisq.test(combined_contingency_table)

# Check if all expected frequencies are >= 5
test_valid <- all(combined_test_result$expected >= 5)

# Store the results in a data frame
combined_chi_square_results <- data.frame(
  Chi_Square_Statistic = combined_test_result$statistic,
  P_Value = combined_test_result$p.value,
  Test_Valid = test_valid
)

# Print the combined results
print(combined_chi_square_results)

# Save results_table as a text file
write.table(combined_chi_square_results, file = "TCGA_PFS_signature_genes_risk_scores_primary_vs_relapse_chi_square.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# Save results_table as a text file
write.table(signature_df_12, file = "TCGA_PFS_signature_genes_data_for_KM_analysis.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# make a copy of signatur_df_12
signature_df_13 <- signature_df_12

# Renaming columns
colnames(signature_df_13)[colnames(signature_df_13) == "new_age_group"] <- "age"
colnames(signature_df_13)[colnames(signature_df_13) == "new_stages"] <- "stages"
colnames(signature_df_13)[colnames(signature_df_13) == "new_treatment"] <- "treatment"
colnames(signature_df_13)[colnames(signature_df_13) == "new_tumor_event"] <- "recurrence"

signature_df_13$recurrence <- mapvalues(signature_df_13$recurrence,
                                        from = c('1', '0'),
                                        to = c('Yes', 'No'))

## Linear model
linear_model_1 <- lm(Vital_status ~ risk_group + age + Sex + stages + treatment + anatomic_site + recurrence +
                       ADAM8_Risk + RAC2_Risk + CSF2_Risk + PPP1R18_Risk + HBEGF_Risk + ITGA5_Risk +
                       LRRC59_Risk + PLAU_Risk + FOSL1_Risk +
                       PLAUR_Risk,
                     data = signature_df_13)
summary(linear_model_1)
forest_model(linear_model_1)
