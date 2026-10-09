## =============================================================================
## Project : Prognostic gene-expression signature for radiotherapy outcome in HPV-negative HNSCC
## Module  : 04 | Clinical covariates: multivariable Cox, risk groups vs. clinical variables
## Script  : 12_FHCRC_PFS_t_test.R
## Purpose : FHCRC cohort, PFS signature: multivariable Cox models (signature genes + clinical
##           covariates), forest plots, median-split high/low risk groups, Kaplan-Meier curves
##           and t-tests of risk groups/risk scores against clinical variables.
## Input   : data/clinical_covariates/merged_df_FHCRC_PFS_for_analysis.txt
## Output  : results/04_clinical_covariate_analysis/FHCRC_PFS_*.txt and plots
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
file_name <- file.path(data_dir, "merged_df_FHCRC_PFS_for_analysis.txt")
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
genes_to_convert = c("ADAM8", "RAC2", "CSF2", "PPP1R18", "HBEGF", "ITGA5",
                     "LRRC59", "PLAU", "FOSL1", "PLAUR")

# Replace all periods with nothing and convert to numeric
for (gene in genes_to_convert) {
  signature_df[[gene]] <- as.numeric(gsub("\\.", "", signature_df[[gene]]))
}


# Convert Vital_status to a Surv object
surv_obj <- Surv(time = signature_df$OS_Months, event = signature_df$Vital_status)


# Example: Fit a Cox model with fewer variables
# Example: Fit a Cox model with fewer variables
cox_model_simple <- coxph(surv_obj ~ Age + Sex + Stages, data = signature_df)
summary(cox_model_simple)


# Example: Expanding the model to include gene variables
extended_cox_model <- coxph(surv_obj ~ Age + Sex + Stages + ADAM8 + RAC2 +
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


# Step 1: Clean the numerical data
numerical_cols <- c("ADAM8", "RAC2", "CSF2", "PPP1R18", "HBEGF", "ITGA5",
                    "LRRC59", "PLAU", "FOSL1", "PLAUR")

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

# Assuming the first 10 columns are gene columns
gene_columns <- names(signature_df_1)[1:10]

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
                       ADAM8_Risk + RAC2_Risk + CSF2_Risk + PPP1R18_Risk + HBEGF_Risk + ITGA5_Risk +
                       LRRC59_Risk + PLAU_Risk + FOSL1_Risk +
                       PLAUR_Risk,
                     data = signature_df_1)
summary(linear_model_1)
forest_model(linear_model_1)


# columns 1 to 10 are gene expression values
gene_expression_cols <- 1:10

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
signature_df_1 <- dplyr::rename(signature_df_1, PFS_Signature_risk_score = Gene_Risk_Score)


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
  high_risk_scores <- signature_df_3$PFS_Signature_risk_score[signature_df_3$new_age_group == age & signature_df_3$Categorized_Risk_Score == "High"]

  # Filter data for the specific age group and Low risk
  low_risk_scores <- signature_df_3$PFS_Signature_risk_score[signature_df_3$new_age_group == age & signature_df_3$Categorized_Risk_Score == "Low"]

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
write.table(t_test_results, file = "FHCRC_PFS_signature_genes_risk_scores_age_group_t_test.txt",
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
median_score <- median(signature_df_4$PFS_Signature_risk_score, na.rm = TRUE)
signature_df_4$Categorized_Risk_Score <- ifelse(signature_df_4$PFS_Signature_risk_score >= median_score, "High", "Low")


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
write.table(t_test_results_stages, file = "FHCRC_PFS_signature_genes_risk_scores_stages_t_test.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# make a copy of signatur_df_4
signature_df_5 <- signature_df_4

# Assuming the Categorized_Risk_Score column has been created
# Categorize the risk score based on the median score
median_score <- median(signature_df_5$PFS_Signature_risk_score, na.rm = TRUE)
signature_df_5$Categorized_Risk_Score <- ifelse(signature_df_5$PFS_Signature_risk_score >= median_score, "High", "Low")

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
write.table(t_test_results_sex, file = "FHCRC_PFS_signature_genes_risk_scores_Sex_t_test.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
