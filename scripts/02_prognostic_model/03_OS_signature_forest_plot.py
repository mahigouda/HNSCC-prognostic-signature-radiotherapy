# =============================================================================
# Project : Prognostic gene-expression signature for radiotherapy outcome in HPV-negative HNSCC
# Module  : 02 | Prognostic model: univariate Cox screening, Lasso-Cox feature selection, risk score
# Script  : 03_OS_signature_forest_plot.py
# Purpose : Re-derivation of the OS signature and per-gene forest plot (HR, events, log-rank p,
#           AIC, concordance).
# Input   : data/tcga/merged_data_final_threshold_0_5.txt
# Output  : results/02_prognostic_model/OS_signature_genes_forest_plot.{pdf,svg}
# =============================================================================
# Cells are delimited with '# %%' (Jupyter/VS Code percent format).

# %%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import os
from sklearn.model_selection import KFold
from lifelines import CoxPHFitter
from sklearn.model_selection import TimeSeriesSplit
from sklearn.feature_selection import SelectKBest, f_regression
from sksurv.util import Surv
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sklearn.model_selection import cross_val_score
from sklearn.linear_model import Lasso
from lifelines.utils import concordance_index
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sksurv.metrics import concordance_index_censored
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sklearn.metrics import get_scorer
from sklearn.metrics import make_scorer
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import pairwise_logrank_test
from sklearn.decomposition import PCA
from lifelines.utils import survival_table_from_events
import matplotlib.ticker as mticker
from tabulate import tabulate
from matplotlib.backends.backend_pdf import PdfPages
from lifelines.statistics import logrank_test
from sklearn.linear_model import LassoCV
from sklearn.linear_model import MultiTaskLassoCV
from sklearn.model_selection import train_test_split


# %%
from sksurv.compare import compare_survival
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.collections import LineCollection
import itertools
from matplotlib.gridspec import GridSpec
from itertools import product
import matplotlib.ticker as ticker


# %%
# ---- Paths -------------------------------------------------------------------
# Run from the repository root. Inputs are read from data/tcga/
# (see data/README.md); outputs are written to results/02_prognostic_model/.
from pathlib import Path
data_dir = Path.cwd() / "data" / "tcga"
out_dir = Path.cwd() / "results" / "02_prognostic_model"
out_dir.mkdir(parents=True, exist_ok=True)
os.chdir(out_dir)


# %%
df = pd.read_csv(data_dir / 'merged_data_final_threshold_0_5.txt', sep='\t')


# %%
# Identify the gene columns
gene_columns = list(df.columns[:202])


# %%
# Remove rows containing "#NAME?" in any column
df = df[~df.applymap(lambda x: x == "#NAME?").any(axis=1)]


# %%
# Remove periods and convert to float
for column in gene_columns:
    # Check if the column is of type string before applying string operations
    if df[column].dtype == 'object':
        df[column] = df[column].str.replace('.', '').astype(float)


# %%
# Remove periods and convert to float
df['OS_Months'] = df['OS_Months'].str.replace('.', '').astype(float)


# %%
# Divide the values in the 'OS_Months' column by 1 x 10^8
df['OS_Months'] = df['OS_Months'] / (1 * 10**8)


# %%
# Remove periods and convert to float
df['Progress_free_survival_Months'] = df['Progress_free_survival_Months'].str.replace('.', '').astype(float)
# Divide the values in the 'OS_Months' column by 1 x 10^8
df['Progress_free_survival_Months'] = df['Progress_free_survival_Months'] / (1 * 10**8)


# %%
# Remove periods and convert to float
df['Disease_specific_survival_Months'] = df['Disease_specific_survival_Months'].str.replace('.', '').astype(float)
# Divide the values in the 'OS_Months' column by 1 x 10^8
df['Disease_specific_survival_Months'] = df['Disease_specific_survival_Months'] / (1 * 10**8)


# %%
# Remove the decimal part and keep the integer part of the "OS_Months" column
df['OS_Months'] = df['OS_Months'].apply(lambda x: int(x) if not np.isnan(x) else x)

# Remove the decimal part and keep the integer part of the "Progress_free_survival_Months" column
df['Progress_free_survival_Months'] = df['Progress_free_survival_Months'].apply(lambda x: int(x)
                                                                                if not np.isnan(x) else x)

# Remove the decimal part and keep the integer part of the "Disease_specific_survival_Months" column
df['Disease_specific_survival_Months'] = df['Disease_specific_survival_Months'].apply(lambda x: int(x)
                                                                                if not np.isnan(x) else x)


# %%
# Remove rows with NaNs
df = df.dropna()


# %%
# Check the data types of each column
print(df.dtypes)


# %%
# Suppress SettingWithCopyWarning
pd.options.mode.chained_assignment = None

# Convert 'Vital_status' column to numeric representation
df['Vital_status'] = df['Vital_status'].map({'ALIVE': 0, 'DEAD': 1})

# Convert 'Radiation therapy' column to numeric representation
df['Radiation_therapy'] = df['Radiation_therapy'].map({'Yes': 1, 'No': 0})

# Convert 'Neoplasm_histologic_grade' column to numeric representation
df['Neoplasm_histologic_grade'] = df['Neoplasm_histologic_grade'].map({'GX': 0, 'G1': 1, 'G2':2, 'G3':3})


# %%
# Keep rows with 'HPV-' in the 'HPV_status' column
df = df[df['HPV_status'] == 'HPV-']


# %%
# Remove the 'HPV_status' column from df
df = df.drop('HPV_status', axis=1)


# %%
df['Disease_specific_survival_Months'] = df['Disease_specific_survival_Months'].astype(int)

df['Progress_free_survival_Months'] = df['Progress_free_survival_Months'].astype(int)

df['OS_Months'] = df['OS_Months'].astype(int)

df['Neoplasm_histologic_grade'] = df['Neoplasm_histologic_grade'].astype(int)


# %%
y = Surv.from_dataframe("Vital_status", "OS_Months", df)

y = pd.DataFrame({'event': df['Vital_status'], 'time': df['OS_Months']}).to_records(index=False)

# Forward feature selection using Lasso Cox regression
X = df[gene_columns]
y = pd.DataFrame({'event': df['Vital_status'], 'time': df['OS_Months']}).to_records(index=False)

tscv = TimeSeriesSplit(n_splits=5)
scores = []


# %%
# Create an empty list to store results
univariate_results = []


# %%
# Convert y_true to structured array
y_true_structured = np.array(list(zip(df['Vital_status'].astype(bool), df['OS_Months'])),
                             dtype=[('event', '?'), ('time', '<f8')])

for gene in df.columns[:-8]: # exclude last seven columns
    gene_X = df[[gene]]

    scaler = StandardScaler().fit(gene_X)
    gene_X_std = scaler.transform(gene_X)

    cph = CoxPHSurvivalAnalysis()

    try:
        cph.fit(gene_X_std, y_true_structured)

        hr = np.exp(cph.coef_[0])

        group_indicator = gene_X[gene].values > gene_X[gene].median()

        test_statistic, p_value = compare_survival(y_true_structured, group_indicator)
        ci_score = concordance_index_censored(y_true_structured['event'], y_true_structured['time'],
                                              cph.predict(gene_X_std))

        univariate_results.append([gene, ci_score[0], hr, p_value])
    except ValueError as e:
        print(f"Error for gene {gene}: {e}")
        continue


# %%
# Print the results
print("Gene\t\tCI Score\tHR\t\tp-value")
for result in univariate_results:
    print("{:<12}\t{:.3f}\t\t{:.3f}\t\t{:.3f}".format(result[0], result[1], result[2], result[3]))


# %%
# Open a new file for writing
with open("univariate_results_threshold_0.5.txt", "w") as f:

    # Write the column headers to the file
    f.write("Gene\tCI Score\tHR\tp-value\n")

    # Write the results to the file
    for result in univariate_results:
        f.write("{:<12}\t{:.3f}\t\t{:.3f}\t\t{:.3f}\n".format(result[0], result[1], result[2], result[3]))

# Close the file
f.close()


# %%
print(len(univariate_results))


# %%
# Sort the results by hazard ratio in descending order
univariate_results_sorted = sorted(univariate_results, key=lambda x: x[2], reverse=True)

# Select the genes with HR > 1 or HR < 1 and p-value < 0.05
HR_top_genes = []
p_value_threshold = 0.05
for result in univariate_results_sorted:
    hr = result[2]
    p_value = result[3]
    if (hr > 1 or hr < 1) and p_value < p_value_threshold:
        HR_top_genes.append(result[0])

# Print the selected genes
print("Genes with HR > 1 or HR < 1 and p-value < 0.05:")
print(HR_top_genes)


# %%
with open('final_univariate_cox_results_threshold_0.5.txt', 'w') as f:
    f.write("Gene\t\tCI Score\tHR\t\tp-value\n")
    for result in univariate_results_sorted:
        f.write("{:<12}\t{:.3f}\t\t{:.3f}\t\t{:.3f}\n".format(result[0], result[1], result[2], result[3]))


# %%
# Save the complete result table for the top genes to a text file
with open("result_table.txt", "w") as f:
    # Write the header to the text file
    f.write("Gene\t\tCI Score\tHR\t\tp-value\n")

    for result in univariate_results_sorted:
        if result[0] in HR_top_genes:
            f.write("{:<12}\t{:.3f}\t\t{:.3f}\t\t{:.3f}\n".format(result[0], result[1], result[2], result[3]))

# Print the result table to the console
with open("result_table.txt", "r") as f:
    print(f.read())


# %%
# Add the additional columns to the list of top genes
selected_columns = HR_top_genes + ['Age', 'Aneuploidy_score', 'Neoplasm_histologic_grade', 'Disease_specific_survival_Months',
                                   'OS_Months', 'Progress_free_survival_Months', 'Vital_status', 'Radiation_therapy']


# %%
# Filter the list to include only the column names that are in the DataFrame
filtered_selected_columns = [col for col in selected_columns if col in df.columns]

# Create a subset DataFrame with the selected columns
univariate_df = df.loc[:, filtered_selected_columns]


# %%
print(univariate_df.columns)


# %%
y = Surv.from_dataframe("Vital_status", "OS_Months", univariate_df)

y = pd.DataFrame({'event': df['Vital_status'], 'time': df['OS_Months']}).to_records(index=False)

# Identify the gene columns
gene_columns = list(univariate_df.columns[:14])
X = univariate_df[gene_columns]
y = pd.DataFrame({'event': df['Vital_status'], 'time': df['OS_Months']}).to_records(index=False)

tscv = TimeSeriesSplit(n_splits=5)
scores = []


# %%
# Function to calculate the concordance index for Lasso Cox regression
def cox_concordance_index(estimator, X, y):
    y_event, y_time = y['event'], y['time']
    risk_scores = -estimator.predict(X)
    return concordance_index(y_time, risk_scores, y_event)
# Create a pipeline with StandardScaler and CoxnetSurvivalAnalysis (Lasso Cox regression)
lasso_cox_pipeline = make_pipeline(StandardScaler(), CoxnetSurvivalAnalysis(alphas=[0.01], n_alphas=1,
                                                                            l1_ratio=1.0, fit_baseline_model=True))
# Convert 'Vital_status' to boolean
df['Vital_status'] = df['Vital_status'].astype(bool)
# Create the structured array
y = Surv.from_dataframe("Vital_status", "OS_Months", df)

# Create a pipeline with StandardScaler and CoxnetSurvivalAnalysis (Lasso Cox regression)
lasso_cox_pipeline = make_pipeline(StandardScaler(), CoxnetSurvivalAnalysis(alphas=[0.01],
                                                                            n_alphas=1, l1_ratio=1.0, fit_baseline_model=True))

lasso_cox = CoxnetSurvivalAnalysis(alphas=[0.01], n_alphas=1, l1_ratio=1.0, fit_baseline_model=True)


lasso_cox_scores = []


# %%
print(df['Vital_status'].unique())


# %%
for train_index, test_index in tscv.split(X):
    X_train, X_test = X.iloc[train_index], X.iloc[test_index]
    y_train, y_test = y[train_index], y[test_index]

    # Convert y_train and y_test to structured arrays
    y_train_structured = np.array(list(zip(y_train["Vital_status"].astype(bool), y_train["OS_Months"])),
                                  dtype=[('event', '?'), ('time', '<f8')])
    y_test_structured = np.array(list(zip(y_test["Vital_status"].astype(bool), y_test["OS_Months"])),
                                 dtype=[('event', '?'), ('time', '<f8')])

    # Standardize the data
    scaler = StandardScaler().fit(X_train)
    X_train_std = scaler.transform(X_train)
    X_test_std = scaler.transform(X_test)

    # Fit the Lasso Cox model
    lasso_cox.fit(X_train_std, y_train_structured)

    # Calculate the concordance index
    ci_score = concordance_index_censored(y_test_structured["event"], y_test_structured["time"], lasso_cox.predict(X_test_std))
    lasso_cox_scores.append(ci_score[0])


# %%
print("Lasso Cox regression concordance index scores:", lasso_cox_scores)


# %%
print("Mean Concordance Index:", np.mean(lasso_cox_scores))


# %%
print("Standard deviation of Concordance Index:", np.std(lasso_cox_scores))


# %%
# Standardize the data
scaler = StandardScaler().fit(X)
X_std = scaler.transform(X)


# %%
# Convert y to a structured array
y_structured = np.array(list(zip(y["Vital_status"].astype(bool), y["OS_Months"])),
                        dtype=[('event', '?'), ('time', '<f8')])


# %%
# Fit the Lasso Cox model on the entire dataset
lasso_cox.fit(X_std, y_structured)


# %%
# Extract the non-zero coefficients and the corresponding feature names
non_zero_coefficients = lasso_cox.coef_[lasso_cox.coef_ != 0]
significant_gene_indices = np.nonzero(lasso_cox.coef_)[0]
significant_genes = X.columns[significant_gene_indices]


# %%
# Print the significant genes and their corresponding coefficients
print("Significant genes and their coefficients:")
for gene, coef in zip(significant_genes, non_zero_coefficients):
    print(f"{gene}: {coef}")

# Store the significant genes in a list
signature_genes = list(significant_genes)
print("\nSignature genes:", signature_genes)

# Specify the name of the file you want to create
filename = "significant_genes.txt"


# %%
# Redirect the standard output of the script to the file
with open(filename, "w") as f:
    # Write the output of the code to the file
    f.write("Significant genes and their coefficients:\n")
    for gene, coef in zip(significant_genes, non_zero_coefficients):
        f.write(f"{gene}: {coef}\n")

# Close the file
f.close()


# %%
surv_sig_genes_univariable = univariate_df.copy()

surv_sig_genes_univariable = surv_sig_genes_univariable.drop(['CSF2'], axis=1)

signature_genes = ['ODC1', 'TUBB', 'ETV4', 'NUDT15', 'AK4', 'PPIB',
                   'PLAU', 'LRRC59', 'CTSV', 'HBEGF', 'PPIA', 'MAP4K4', 'AXL']

additional_columns = ['OS_Months', 'Vital_status']

signature_df = surv_sig_genes_univariable[signature_genes + additional_columns]


# %%
# Convert the 'Vital_status' and 'OS_Months' columns to a structured array
y_structured = np.array(list(zip(signature_df['Vital_status'].astype(bool), signature_df['OS_Months'])),
                        dtype=[('event', '?'), ('time', '<f8')])

# Identify the gene columns
gene_columns = signature_genes
X = signature_df[gene_columns]

# Standardize the data
scaler = StandardScaler().fit(X)
X_std = scaler.transform(X)

# Fit the Lasso Cox model
lasso_cox.fit(X_std, y_structured)


# %%
# Construct a DataFrame from your significant genes and coefficients
data = pd.DataFrame({'Gene': significant_genes, 'coef': non_zero_coefficients})


# %%
# Calculate Hazard Ratios
data['HR'] = np.exp(data['coef'])


# %%
# Initialize the Cox Proportional Hazard Model
cph = CoxPHFitter()

# Calculate additional metrics for each gene
for gene in significant_genes:
    # Fit a CoxPHFitter model for each gene
    cph.fit(signature_df[[gene, 'OS_Months', 'Vital_status']], duration_col='OS_Months', event_col='Vital_status')
    data.loc[data['Gene'] == gene, 'Events'] = signature_df['Vital_status'].sum()  # number of events
    data.loc[data['Gene'] == gene, 'AIC'] = cph.AIC_partial_  # Akaike Information Criterion

    # Calculate log-rank test
    high_risk = signature_df[signature_df[gene] >= signature_df[gene].median()]
    low_risk = signature_df[signature_df[gene] < signature_df[gene].median()]
    high_risk_p = high_risk.loc[high_risk['OS_Months'] <= 60, 'OS_Months']
    low_risk_p = low_risk.loc[low_risk['OS_Months'] <= 60, 'OS_Months']
    results = logrank_test(high_risk_p, low_risk_p, event_observed_A=high_risk.loc[high_risk['OS_Months'] <= 60, 'Vital_status'], event_observed_B=low_risk.loc[low_risk['OS_Months'] <= 60, 'Vital_status'])
    data.loc[data['Gene'] == gene, 'Log_rank_p'] = results.p_value


# %%
# Perform bootstrapping for coefficient estimation
n_bootstrap_samples = 1000
boot_coefs = []

for _ in range(n_bootstrap_samples):
    sample_X, sample_y = resample(X_std, y_structured)
    lasso_cox.fit(sample_X, sample_y)
    boot_coefs.append(lasso_cox.coef_)

ci_lower = np.percentile(boot_coefs, 2.5, axis=0)  # lower limit of the 95% CI
ci_upper = np.percentile(boot_coefs, 97.5, axis=0)  # upper limit of the 95% CI


# %%
# Calculate hazard ratios and their confidence intervals
log_hr = lasso_cox.coef_
hr = np.exp(log_hr)
ci_lower_hr = np.exp(ci_lower)
ci_upper_hr = np.exp(ci_upper)
hr = np.array(hr).reshape(-1)
ci_lower_hr = np.array(ci_lower_hr).reshape(-1)
ci_upper_hr = np.array(ci_upper_hr).reshape(-1)

# Calculate the absolute difference between the hazard ratios and the confidence intervals
xerr_lower = np.abs(hr - ci_lower_hr)
xerr_upper = np.abs(ci_upper_hr - hr)


# %%
xerr = np.maximum(xerr_lower, xerr_upper)


# %%
plt.figure(figsize=(9, 5))
# Plot the hazard ratios with error bars
plt.errorbar(x=hr, y=range(len(hr)), xerr=xerr, fmt='o', capsize=5)
plt.xscale('log')
plt.xlabel('Hazard Ratio')
plt.ylabel('Genes')
plt.title('Significant Genes')

# Specify the desired tick positions
tick_positions = [0.01, 0.1, 1, 2, -2]

# Set the x-axis tick positions and labels
plt.xticks(tick_positions, tick_positions)

# Set the y-axis tick labels to the gene names
plt.yticks(range(len(signature_genes)), signature_genes)

plt.xlim(0.1, 2)  # Set the x-axis limits to 0 and 2

# Add text annotation for log_rank_p at the top of the plot
plt.text(2.1, len(signature_genes), "Log-rank p-value", ha='left', va='center', fontsize=12)
plt.text(5.8, len(signature_genes), "CI", ha='left', va='center', fontsize=12)

# Add text annotations for events, AIC, and log_rank_p
for i, gene in enumerate(signature_genes):
    event_count = data.loc[data['Gene'] == gene, 'Events'].values[0]
    aic_value = data.loc[data['Gene'] == gene, 'AIC'].values[0]
    log_rank_p = data.loc[data['Gene'] == gene, 'Log_rank_p'].values[0]

    # Position the text annotations next to the corresponding gene
    plt.text(2.1, i, f"{log_rank_p:.4f}", ha='left', va='center', fontsize=12)

# Add text annotations for average events, AIC, and concordance
plt.text(0.03, 0.1, f'Events: {avg_events:.2f}', fontsize=12, transform=plt.gca().transAxes)
plt.text(0.03, 0.16, f'AIC: {avg_aic:.2f}', fontsize=12, transform=plt.gca().transAxes)
plt.text(0.03, 0.22, f'Concordance: {avg_concordance:.2f}', fontsize=12, transform=plt.gca().transAxes)

# Add text annotations for ci_lower_hr and ci_upper_hr
for i, gene in enumerate(signature_genes):
    # Position the text annotations next to the corresponding gene
    plt.text(4.9, i, f"({ci_lower_hr[i]:.2f} - {ci_upper_hr[i]:.2f})", ha='left', va='center', fontsize=12)

# Remove grid lines
plt.grid(False)

# Adjust layout
plt.tight_layout()

# Save the plot to a PDF file
with PdfPages("OS_signature_genes_forest_plot.pdf") as pdf:
    pdf.savefig()

# Save the plot to an SVG file
plt.savefig("OS_signature_genes_forest_plot.svg", format='svg')

plt.show()
