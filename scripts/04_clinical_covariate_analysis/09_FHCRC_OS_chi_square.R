## =============================================================================
## Project : Prognostic gene-expression signature for radiotherapy outcome in HPV-negative HNSCC
## Module  : 04 | Clinical covariates: multivariable Cox, risk groups vs. clinical variables
## Script  : 09_FHCRC_OS_chi_square.R
## Purpose : FHCRC cohort, OS signature: multivariable Cox models (signature genes + clinical
##           covariates), forest plots, median-split high/low risk groups, Kaplan-Meier curves
##           and chi-square tests of risk groups/risk scores against clinical variables.
## Input   : data/clinical_covariates/merged_df_FHCRC_for_analysis.txt
## Output  : results/04_clinical_covariate_analysis/FHCRC_OS_*.txt and plots
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
library(ggfortify)
library(ranger)
library(knitr)
library(gt)
library(zoo)
library(lmtest)

## ---- Paths ------------------------------------------------------------------
## Run from the repository root. Inputs are read from data/clinical_covariates/
## (see data/README.md); outputs are written to results/04_clinical_covariate_analysis/.
data_dir <- file.path(getwd(), "data", "clinical_covariates")
out_dir  <- file.path(getwd(), "results", "04_clinical_covariate_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)


## Sub data
file_name <- file.path(data_dir, "merged_df_FHCRC_for_analysis.txt")
signature_df <- read.table(file_name, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = "\t")

# Map Vital status to numeric
signature_df$Vital_status <- mapvalues(signature_df$Vital_status,
                                       from = c('Alive', 'Dead'),
                                       to = c('0', '1'))

# Remove "fu time: " prefix and convert to numeric
signature_df$OS_Months <- as.numeric(gsub("fu time: ", "", signature_df$OS_Months))

# Check the result
head(signature_df$OS_Months)

# Remove rows with NA values in any column
signature_df <- na.omit(signature_df)

# specify the range of columns more safely as follows:
num_cols <- ncol(signature_df)
gene_expression_cols <- signature_df[, 1:(num_cols - 6)]

# Now calculate the median for each gene
signature_df_median <- apply(gene_expression_cols, 2, median)

head(signature_df$Age)

# Function to calculate midpoint of age ranges
calculate_midpoint <- function(age_range) {
  age_parts <- as.numeric(unlist(strsplit(age_range, "-")))
  return(mean(age_parts))
}

# Apply this function to the Age column
signature_df$NumericAge <- sapply(signature_df$Age, calculate_midpoint)

## Categorize numeric ages
signature_df$age_group <- cut(signature_df$NumericAge,
                              breaks = c(0, 50, 65, 80, Inf),
                              labels = c("0-50", "51-65", "66-80", "80+"))

# Check the first few values of the new age_group column
head(signature_df$age_group)


# Remove the NumericAge column
signature_df <- subset(signature_df, select = -NumericAge)

# Set NumericAge column to NULL
signature_df$NumericAge <- NULL

# Current column names
current_cols <- colnames(signature_df)

# Create a vector of column names without OS_Months and Vital_status
cols_without_outcome <- current_cols[!current_cols %in% c("OS_Months", "Vital_status")]

# Append OS_Months and Vital_status to the end
new_col_order <- c(cols_without_outcome, "OS_Months", "Vital_status")

# Reorder the columns in signature_df
signature_df <- signature_df[, new_col_order]

# Verify the new column order
colnames(signature_df)

# Adjusted analysis columns
analysis_cols <- colnames(signature_df)[1:17]
analysis_cols <- analysis_cols[!analysis_cols %in% c("OS_Months", "Vital_status")]

# Confirm the adjusted columns
print(analysis_cols)

# Initialize a dataframe to store results
results_df <- data.frame(Variable = character(), HR = numeric(),
                         CI = character(), pvalue = numeric(), stringsAsFactors = FALSE)

unique(signature_df$Vital_status)

# Convert to numeric
signature_df$Vital_status <- as.numeric(signature_df$Vital_status)

table(signature_df$Vital_status)

# List of gene columns to convert to numeric
genes_to_convert = c("ODC1", "TUBB", "ETV4", "NUDT15", "AK4", "PPIB", "PLAU", "LRRC59", "CTSV", "HBEGF", "PPIA", "MAP4K4", "AXL")

# Replace all periods with nothing and convert to numeric
for (gene in genes_to_convert) {
  signature_df[[gene]] <- as.numeric(gsub("\\.", "", signature_df[[gene]]))
}


# Convert Vital_status to a Surv object
surv_obj <- Surv(time = signature_df$OS_Months, event = signature_df$Vital_status)

# Fit a Cox proportional hazards model

# Example: Fit a Cox model with fewer variables
# Example: Fit a Cox model with fewer variables
cox_model_simple <- coxph(surv_obj ~ Age + Sex + Stages, data = signature_df)
summary(cox_model_simple)


# Example: Expanding the model to include gene variables
extended_cox_model <- coxph(surv_obj ~ Age + Sex + Stages + ODC1 + TUBB + ETV4 +
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


# Step 1: Clean the numerical data
numerical_cols <- c("ODC1", "TUBB", "ETV4", "NUDT15", "AK4", "PPIB", "PLAU", "LRRC59",
                    "CTSV", "HBEGF", "PPIA", "MAP4K4", "AXL")

# Replace periods with commas and convert to numeric
signature_df[numerical_cols] <- lapply(signature_df[numerical_cols], function(x) {
  as.numeric(gsub("\\.", "", x))
})

# Step 2: Convert 'Age' to a numeric variable (using midpoints for example)
signature_df$NumericAge <- sapply(signature_df$Age, function(age_range) {
  ages <- as.numeric(unlist(strsplit(age_range, "-")))
  mean(ages)
})

# Step 3: Convert categorical variables to factors
signature_df$Sex <- as.factor(signature_df$Sex)
signature_df$Stages <- as.factor(signature_df$Stages)
signature_df$Disease_specific <- as.factor(signature_df$Disease_specific)
signature_df$age_group <- as.factor(signature_df$age_group)

# Step 4: Convert 'Vital_status' to numeric
signature_df$Vital_status <- as.numeric(signature_df$Vital_status)


# Adjust the formula as needed to include the variables you're interested in
cox_model <- coxph(Surv(OS_Months, Vital_status) ~ NumericAge + Sex + Stages + Disease_specific + age_group, data = signature_df)

# Summarize the model
summary(cox_model)

cox_fit <- survfit(cox_model)

autoplot(cox_fit)


# Tidy the summary of the Cox model
tidy_cox <- tidy(cox_model)

# Print the tidy summary
print(tidy_cox)

# Add lower and upper confidence intervals for HR
tidy_cox$CI_lower <- exp(tidy_cox$estimate - 1.96 * tidy_cox$std.error)
tidy_cox$CI_upper <- exp(tidy_cox$estimate + 1.96 * tidy_cox$std.error)

# Create the forest plot
ggplot(tidy_cox, aes(y = term, x = estimate, xmin = CI_lower, xmax = CI_upper)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") + # Reference line at HR = 1
  geom_point() +  # Points for HRs
  geom_errorbarh(height = 0.2) +  # Horizontal error bars for CIs
  scale_x_continuous(name = "Hazard Ratio (95% CI)", trans = "log") + # Log-transformed x-axis
  theme_minimal() +
  theme(axis.title.y = element_blank()) # Hide y-axis title


## Median High Vs low risk
# Assuming 'extended_cox_model' is your final Cox model
# Calculate risk scores
risk_scores <- predict(extended_cox_model, type = "risk")

# Check the number of rows in the dataset and the length of risk_scores
nrow(signature_df)
length(risk_scores)

# Assuming the first 251 rows of signature_df were used in the Cox model
signature_df <- signature_df[1:97, ]

# Exclude observations from signature_df that were not used in the Cox model
signature_df <- signature_df[!is.na(risk_scores), ]

# Add risk scores to the dataset
signature_df$risk_score <- risk_scores

# Median split for high vs low risk
signature_df$risk_group <- ifelse(signature_df$risk_score > median(signature_df$risk_score), "High", "Low")

surv_obj <- Surv(signature_df$OS_Months, signature_df$Vital_status)
km_fit <- survfit(surv_obj ~ risk_group, data = signature_df)

signature_df$risk_group <- as.factor(signature_df$risk_group)

summary(km_fit, times = c(1,30,60,90*(1:10)))

# Plot survival curves
plot(km_fit, main = "Survival Curves for High vs. Low Risk Groups", xlab = "Time", ylab = "Survival Probability", col = 1:2)
legend("topright", legend = levels(signature_df$risk_group), col = 1:2, lty = 1)

# Log-rank test
survdiff(surv_obj ~ risk_group, data = signature_df)


# Making a copy of signature_df
signature_df_1 <- signature_df

# Assuming the first 13 columns are gene columns
gene_columns <- names(signature_df_1)[1:13]

# Create new columns for each gene to indicate risk category
for (gene in gene_columns) {
  new_col_name <- paste(gene, "Risk", sep = "_")
  signature_df_1 <- signature_df_1 %>%
    mutate(!!new_col_name := ifelse(get(gene) >= some_threshold, "High", "Low"))
}

# Define a function to categorize continuous variables
categorize_risk <- function(x) {
  median_val <- median(x, na.rm = TRUE)
  ifelse(x >= median_val, "High", "Low")
}

# Assuming signature_df_1 as dataframe
for (col_name in colnames(signature_df_1)) {
  # Skip for non-continuous variables or already existing risk categories
  if (col_name %in% c("Age", "Sex", "Stages", "Disease_specific", "age_group", "risk_group")) {
    next
  }

  # Continue with your operations for other columns
  # For example, creating a risk category for each gene
  risk_col_name <- paste0(col_name, "_Risk")
  # Assuming categorize_risk is a function you've defined to categorize risk
  signature_df_1[[risk_col_name]] <- categorize_risk(signature_df_1[[col_name]])
}


signature_df_1 <- signature_df_1 %>%
  mutate(
    Stages_Risk = ifelse(Stages %in% c("III/IV"), "High", "Low"),
    Disease_specific_Risk = ifelse(Disease_specific %in% c("Alive", "Dead -unk casue", "Dead-non OC", "Dead-oral ca"), "High", "Low"),
    age_group_Risk = ifelse(age_group %in% c("0-50", "51-65", "66-80", "80+"), "High", "Low")
  )

# Summary for Stages
table(signature_df_1$Stages, signature_df_1$Stages_Risk)

# Summary for Radiation Therapy
table(signature_df_1$Disease_specific, signature_df_1$Disease_specific_Risk)

# Summary for Age Group
table(signature_df_1$age_group, signature_df_1$age_group_Risk)

signature_df_1$Stages_Risk <- as.factor(signature_df_1$Stages_Risk)
signature_df_1$Disease_specific_Risk <- as.factor(signature_df_1$Disease_specific_Risk)
signature_df_1$age_group_Risk <- as.factor(signature_df_1$age_group_Risk)

levels(signature_df_1$Stages_Risk)
levels(signature_df_1$Disease_specific_Risk)
levels(signature_df_1$age_group_Risk)

## Linear model
linear_model_1 <- lm(Vital_status ~ risk_group + age_group + Sex + Stages +
                       ODC1_Risk + TUBB_Risk + ETV4_Risk + NUDT15_Risk + AK4_Risk + PPIB_Risk +
                       PLAU_Risk + LRRC59_Risk + CTSV_Risk +
                       HBEGF_Risk + PPIA_Risk + MAP4K4_Risk + AXL_Risk,
                     data = signature_df_1)
summary(linear_model_1)
forest_model(linear_model_1)


# Assuming columns 1 to 13 are your gene expression values
gene_expression_cols <- 1:13

# Checking for any negative or zero values in gene expression columns
if(any(signature_df_1[gene_expression_cols] <= 0)) {
  print("There are non-positive values in gene expression columns")
}


# Apply log2 transformation after adding a small constant
signature_df_1[gene_expression_cols] <- log2(signature_df_1[gene_expression_cols] + 1)

# Scale the log-transformed data
signature_df_1[gene_expression_cols] <- scale(signature_df_1[gene_expression_cols])

# Calculate the mean for each row
signature_df_1$Gene_Risk_Score <- rowMeans(signature_df_1[gene_expression_cols], na.rm = TRUE)

# Relabel the column
signature_df_1 <- dplyr::rename(signature_df_1, OS_Signature_risk_score = Gene_Risk_Score)


# make a copy of signatur_df_1
signature_df_3 <- signature_df_1

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
stages <- c("I/II", "III/IV")
sexes <- c("F", "M")
risk_groups <- c("High", "Low")


#  Age < 65 and High Risk
subset_data <- subset(signature_df_3, new_age_group == "<65" & risk_group == "High")


# Example Chi-square test using ODC1_Risk and new_age_group
if (all(c("ODC1_Risk", "new_age_group") %in% names(signature_df_3))) {
  # Create a contingency table
  contingency_table <- table(signature_df_3$new_age_group, signature_df_3$ODC1_Risk)

  # Print the contingency table for inspection
  print("Contingency Table for ODC1_Risk and new_age_group")
  print(contingency_table)

  # Check if each cell in the table has at least expected count 5
  if (all(contingency_table >= 5)) {
    # Perform Chi-square test
    test_result <- chisq.test(contingency_table)
    print(test_result)
  } else {
    print("Chi-square test not valid due to low expected counts in some cells.")
  }
} else {
  print("One or both variables (ODC1_Risk, new_age_group) do not exist in the dataset.")
}


# Check the data type of OS_Signature_risk_score
if ("OS_Signature_risk_score" %in% names(signature_df_3)) {
  data_type <- class(signature_df_3$OS_Signature_risk_score)
  print(paste("Data type of OS_Signature_risk_score:", data_type))
} else {
  print("OS_Signature_risk_score does not exist in the dataset.")
}


# Categorizing OS_Signature_risk_score
median_score <- median(signature_df_3$OS_Signature_risk_score, na.rm = TRUE)
signature_df_3$Categorized_Risk_Score <- ifelse(signature_df_3$OS_Signature_risk_score >= median_score, "High", "Low")

# Convert to factor
signature_df_3$Categorized_Risk_Score <- factor(signature_df_3$Categorized_Risk_Score)

# Create a contingency table
contingency_table <- table(signature_df_3$new_age_group, signature_df_3$Categorized_Risk_Score)

# Print the contingency table for inspection
print("Contingency Table for Categorized OS_Signature_risk_score and new_age_group")
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
    subset_data <- subset(signature_df_3, new_age_group == age & risk_group == risk)

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
  Test_Valid = logical(),
  Count = integer(),
  stringsAsFactors = FALSE
)

# Loop through each combination of age group and risk group
for (age in c("<65", ">65")) {
  for (risk in c("High", "Low")) {
    # Filter data for the specific combination
    subset_data <- subset(signature_df_3, new_age_group == age & risk_group == risk)

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
write.table(chi_square_results, file = "FHCRC_OS_signature_genes_risk_scores_age_group.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# Create a combined contingency table for new_age_group and risk_group
combined_contingency_table_age <- table(signature_df_3$new_age_group, signature_df_3$risk_group)

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
write.table(combined_chi_square_results_age, file = "FHCRC_OS_signature_genes_risk_scores_age_group_chi_square.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# make a copy of signatur_df_1
signature_df_4 <- signature_df_1

# Classify age_groups into <65 and >65
signature_df_4 <- signature_df_4 %>%
  mutate(new_age_group = ifelse(age_group %in% c("0-50", "51-65"), "<65",
                                ifelse(age_group %in% c("66-80", "80+"), ">65", NA)))

# Classify stages into I/II and III/IV
signature_df_4 <- signature_df_4 %>%
  mutate(new_stages = ifelse(Stages %in% c("I/II"), "Stage I/II",
                             ifelse(Stages %in% c("III/IV"), "Stage III/IV", NA)))


# Update the encoding for age_group to reflect the new classification
signature_df_4 <- signature_df_4 %>%
  mutate(age_group_encoded = ifelse(new_age_group == "<65", 0, 1)) # 0 for "<65", 1 for ">65"

# Encode Sex and Stages as numeric
signature_df_4$Sex_encoded <- ifelse(signature_df_4$Sex == "F", 0, 1)  # Example: 0 for Female, 1 for Male
signature_df_4$new_stages_encoded <- as.numeric(as.factor(signature_df_4$new_stages))


print(table(signature_df_4$new_stages))
print(table(signature_df_4$risk_group))

# Categorizing OS_Signature_risk_score
median_score <- median(signature_df_4$OS_Signature_risk_score, na.rm = TRUE)
signature_df_4$Categorized_Risk_Score <- ifelse(signature_df_4$OS_Signature_risk_score >= median_score, "High", "Low")

# Convert to factor
signature_df_4$Categorized_Risk_Score <- factor(signature_df_4$Categorized_Risk_Score)

# Create a contingency table
contingency_table <- table(signature_df_4$new_stages, signature_df_4$Categorized_Risk_Score)


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
    subset_data <- subset(signature_df_4, new_stages == stage & risk_group == risk)

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
write.table(chi_square_results_stages, file = "FHCRC_OS_signature_genes_risk_scores_stages.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# Create a combined contingency table for new_stages and risk_group
combined_contingency_table_stages <- table(signature_df_4$new_stages_encoded, signature_df_4$risk_group)

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
write.table(combined_chi_square_results_stages, file = "FHCRC_OS_signature_genes_risk_scores_stages_chi_square.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# make a copy of signatur_df_4
signature_df_5 <- signature_df_4


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
for (sex in unique(signature_df_5$Sex_encoded)) { # Use unique to get all sex groups
  for (risk in c("High", "Low")) {
    # Filter data for the specific combination
    subset_data <- subset(signature_df_5, Sex_encoded == sex & risk_group == risk)

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
write.table(chi_square_results_sex, file = "FHCRC_OS_signature_genes_risk_scores_Sex.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# Create a combined contingency table for Sex_encoded and risk_group
combined_contingency_table_sex <- table(signature_df_5$Sex_encoded, signature_df_5$risk_group)

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
write.table(combined_chi_square_results_sex, file = "FHCRC_OS_signature_genes_risk_scores_Sex_chi_square.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# make a copy of signatur_df_5
signature_df_8 <- signature_df_5

# Renaming columns
colnames(signature_df_8)[colnames(signature_df_8) == "new_age_group"] <- "age"
colnames(signature_df_8)[colnames(signature_df_8) == "new_stages"] <- "stages"


ggsurvplot(
  survfit(Surv(OS_Months, Vital_status) ~ stages + strata(risk_group), data = signature_df_8),
  data = signature_df_8,
  pval = NULL,
  risk.table = TRUE,
  xlim = c(0, 60),
  break.time.by = 12,
  xlab = "Overall Survival (months)"
)

ggsurvplot(
  survfit(Surv(OS_Months, Vital_status) ~ age + strata(risk_group), data = signature_df_8),
  data = signature_df_8,
  pval = NULL,
  risk.table = TRUE,
  xlim = c(0, 60),
  break.time.by = 12,
  xlab = "Overall Survival (months)"
)


ggsurvplot(
  survfit(Surv(OS_Months, Vital_status) ~ Sex + strata(risk_group), data = signature_df_8),
  data = signature_df_8,
  pval = NULL,
  risk.table = TRUE,
  xlim = c(0, 60),
  break.time.by = 12,
  xlab = "Overall Survival (months)"
)


## Linear model
linear_model_1 <- lm(Vital_status ~ risk_group + age + Sex + stages +
                       ODC1_Risk + TUBB_Risk + ETV4_Risk + NUDT15_Risk + AK4_Risk + PPIB_Risk +
                       PLAU_Risk + LRRC59_Risk + CTSV_Risk +
                       HBEGF_Risk + PPIA_Risk + MAP4K4_Risk + AXL_Risk,
                     data = signature_df_8)
summary(linear_model_1)
forest_model(linear_model_1)
