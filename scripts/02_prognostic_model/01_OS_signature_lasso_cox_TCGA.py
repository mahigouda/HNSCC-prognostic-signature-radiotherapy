# =============================================================================
# Project : Prognostic gene-expression signature for radiotherapy outcome in HPV-negative HNSCC
# Module  : 02 | Prognostic model: univariate Cox screening, Lasso-Cox feature selection, risk score
# Script  : 01_OS_signature_lasso_cox_TCGA.py
# Purpose : HPV-negative TCGA-HNSC: univariate Cox screening of candidate genes (C-index, HR,
#           log-rank p), Lasso-penalised Cox regression (scikit-survival Coxnet, time-series
#           cross-validation, concordance index) for signature selection, multivariable Cox risk
#           score, median split into high/low risk, Kaplan-Meier analysis stratified by
#           radiotherapy for OS, PFS and DSS.
# Input   : data/tcga/merged_data_final_threshold_0_5.txt,
#           data/tcga/signature_genes_merged_df_for_DSS.txt
# Output  : results/02_prognostic_model/ (univariate Cox tables, signature genes and
#           coefficients, KM plots, risk tables)
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


# %%
# Store the significant genes in a list
signature_genes = list(significant_genes)
print("\nSignature genes:", signature_genes)


# %%
# Specify the name of the file you want to create
filename = "significant_genes.txt"

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


# %%
surv_sig_genes_univariable = surv_sig_genes_univariable.drop(['CSF2'], axis=1)


# %%
signature_genes = ['ODC1', 'TUBB', 'ETV4', 'NUDT15', 'AK4', 'PPIB',
                   'PLAU', 'LRRC59', 'CTSV', 'HBEGF', 'PPIA', 'MAP4K4', 'AXL']

additional_columns = ['Age', 'Aneuploidy_score', 'Neoplasm_histologic_grade',
                                   'OS_Months', 'Progress_free_survival_Months', 'Vital_status', 'Radiation_therapy']

signature_df = surv_sig_genes_univariable[signature_genes + additional_columns]


# %%
signature_df = signature_df.drop('Aneuploidy_score', axis=1)
signature_df = signature_df.drop('Neoplasm_histologic_grade', axis=1)


# %%
# Create Cox proportional hazards model
cph = CoxPHFitter()
cph.fit(signature_df, duration_col='OS_Months', event_col='Vital_status')
risk_scores = cph.predict_partial_hazard(signature_df)

signature_df['Risk_scores'] = risk_scores
signature_df['Risk_group'] = np.where(signature_df['Risk_scores'] > np.median(signature_df['Risk_scores']), 1, 0)

# Define high_risk and low_risk masks
high_risk = signature_df['Risk_group'] == 1
low_risk = signature_df['Risk_group'] == 0

# Define radiation_therapy mask
rt = signature_df['Radiation_therapy'] == 1
no_rt = signature_df['Radiation_therapy'] == 0


# %%
# Calculate hazard ratios for radiation therapy and high-risk patients
rt_high_risk_hr = cph.predict_partial_hazard(signature_df.loc[high_risk & rt]).mean() / cph.predict_partial_hazard(signature_df.loc[low_risk & rt]).mean()

# Calculate hazard ratios for radiation therapy and low-risk patients
rt_low_risk_hr = cph.predict_partial_hazard(signature_df.loc[low_risk & rt]).mean() / cph.predict_partial_hazard(signature_df.loc[high_risk & rt]).mean()

# Calculate hazard ratios for no radiation therapy and high-risk patients
no_rt_high_risk_hr = cph.predict_partial_hazard(signature_df.loc[high_risk & no_rt]).mean() / cph.predict_partial_hazard(signature_df.loc[low_risk & no_rt]).mean()

# Calculate hazard ratios for no radiation therapy and low-risk patients
no_rt_low_risk_hr = cph.predict_partial_hazard(signature_df.loc[low_risk & no_rt]).mean() / cph.predict_partial_hazard(signature_df.loc[high_risk & no_rt]).mean()


# %%
# Print hazard ratios
print(f'Hazard ratio for radiation therapy: {rt_high_risk_hr:.3f}')
print(f'Hazard ratio for no radiation therapy: {rt_low_risk_hr:.3f}')
print(f'Hazard ratio for radiation therapy: {no_rt_high_risk_hr:.3f}')
print(f'Hazard ratio for no radiation therapy: {no_rt_low_risk_hr:.3f}')


# %%
kmf = KaplanMeierFitter()


# %%
# Calculate log-rank p-value for radiation therapy, high-risk group
T_rt_high_risk = signature_df.loc[(rt) & (high_risk) & (signature_df['OS_Months'] <= 60), 'OS_Months']
E_rt_high_risk = signature_df.loc[(rt) & (high_risk) & (signature_df['OS_Months'] <= 60), 'Vital_status']
T_rt_low_risk = signature_df.loc[(rt) & (low_risk) & (signature_df['OS_Months'] <= 60), 'OS_Months']
E_rt_low_risk = signature_df.loc[(rt) & (low_risk) & (signature_df['OS_Months'] <= 60), 'Vital_status']
results_rt_high_vs_low = logrank_test(T_rt_high_risk, T_rt_low_risk, event_observed_A=E_rt_high_risk, event_observed_B=E_rt_low_risk)
p_value_rt_high_vs_low = results_rt_high_vs_low.p_value


# %%
# Calculate log-rank p-value for no radiation therapy, high-risk group
T_no_rt_high_risk = signature_df.loc[(no_rt) & (high_risk) & (signature_df['OS_Months'] <= 60), 'OS_Months']
E_no_rt_high_risk = signature_df.loc[(no_rt) & (high_risk) & (signature_df['OS_Months'] <= 60), 'Vital_status']
T_no_rt_low_risk = signature_df.loc[(no_rt) & (low_risk) & (signature_df['OS_Months'] <= 60), 'OS_Months']
E_no_rt_low_risk = signature_df.loc[(no_rt) & (low_risk) & (signature_df['OS_Months'] <= 60), 'Vital_status']
results_no_rt_high_vs_low = logrank_test(T_no_rt_high_risk, T_no_rt_low_risk, event_observed_A=E_no_rt_high_risk, event_observed_B=E_no_rt_low_risk)
p_value_no_rt_high_vs_low = results_no_rt_high_vs_low.p_value


# %%
# Labels and line styles for the stratified groups
labels = {
    (1, 1): 'High Risk',
    (1, 0): 'High Risk',
    (0, 1): 'Low Risk',
    (0, 0): 'Low Risk'
}

colors = {
    (1, 1): 'red',
    (1, 0): 'darkred',
    (0, 1): 'blue',
    (0, 0): 'darkblue'
}

line_styles = {
    (1, 1): '-',
    (1, 0): '--',
    (0, 1): '-',
    (0, 0): '--'
}

# Set Arial font with size 7
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 7

# Increase figure size based on font size
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3))

# Create plot for radiation therapy
for risk_group in [1, 0]:
    mask = (signature_df['Risk_group'] == risk_group) & (signature_df['Radiation_therapy'] == 1)
    kmf.fit(signature_df.loc[mask, 'OS_Months'], event_observed=signature_df.loc[mask, 'Vital_status'], label=labels[(risk_group, 1)])
    kmf.plot(ax=ax1, ci_show=True, ci_alpha=0.1, color=colors[(risk_group, 1)], linewidth=3, linestyle=line_styles[(risk_group, 1)], marker='o', markersize=4)

ax1.set_title('Radiation Therapy', fontsize=7)
ax1.set_xlabel('Overall survival (months)', fontsize=7)
ax1.set_ylabel('Survival probability', fontsize=7)
ax1.legend(fontsize=7)
ax1.set_xlim(0, 60)
ax1.set_yticks([0, 0.25, 0.5, 0.75, 1])
ax1.set_xticks([0, 12, 24, 36, 48, 60])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_linewidth(2)
ax1.spines['left'].set_linewidth(2)
# Add hazard ratios and log-rank p-value to plot
ax1.text(0.5, 0.05, f'logrank p high vs low={p_value_rt_high_vs_low:.3f}', transform=ax1.transAxes, fontsize=7, ha='right', va='bottom')
ax1.text(0.5, 0.1, f'HR (high-risk): {rt_high_risk_hr:.2f}', transform=ax1.transAxes, fontsize=7, ha='right', va='bottom')
ax1.text(0.5, 0.15, f'HR (low-risk): {rt_low_risk_hr:.2f}', transform=ax1.transAxes, fontsize=7, ha='right', va='bottom')


# Create plot for no radiation therapy
for risk_group in [1, 0]:
    mask = (signature_df['Risk_group'] == risk_group) & (signature_df['Radiation_therapy'] == 0)
    kmf.fit(signature_df.loc[mask, 'OS_Months'], event_observed=signature_df.loc[mask, 'Vital_status'], label=labels[(risk_group, 0)])
    kmf.plot(ax=ax2, ci_show=True, ci_alpha=0.1, color=colors[(risk_group, 0)], linewidth=3, linestyle=line_styles[(risk_group, 0)], marker='o', markersize=4)

ax2.set_title('No Radiation Therapy', fontsize=7)
ax2.set_xlabel('Overall survival (months)', fontsize=7)
ax2.set_ylabel('Survival probability', fontsize=7)
ax2.legend(fontsize=7)
ax2.set_xlim(0, 60)
ax2.set_yticks([0, 0.25, 0.5, 0.75, 1])
ax2.set_xticks([0, 12, 24, 36, 48, 60])
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_linewidth(2)
ax2.spines['left'].set_linewidth(2)
ax2.text(0.5, 0.05, f'logrank p high vs low={p_value_no_rt_high_vs_low:.3f}', transform=ax2.transAxes, fontsize=7, ha='right', va='bottom')
ax2.text(0.5, 0.1, f'HR (high-risk): {no_rt_high_risk_hr:.3f}', transform=ax2.transAxes, fontsize=7, ha='right', va='bottom')
ax2.text(0.5, 0.15, f'HR (low-risk): {no_rt_low_risk_hr:.3f}', transform=ax2.transAxes, fontsize=7, ha='right', va='bottom')

# Save the plot to a PDF file
with PdfPages("KM_plots.pdf") as pdf:
    pdf.savefig(fig)

# Save the plot to an SVG file
plt.savefig("KM_plots.svg", format='svg')

plt.show()


# %%
# Define time points
time_points = [0, 10, 20, 30, 40, 50, 60]

# Create empty DataFrames to store counts
low_risk_counts = pd.DataFrame(index=time_points, columns=['Count', 'Median risk score', 'At risk'])
high_risk_counts = pd.DataFrame(index=time_points, columns=['Count', 'Median risk score', 'At risk'])

# Filter patients who received radiation therapy
signature_df_rt_yes = signature_df[signature_df['Radiation_therapy'] == 1]

for t in time_points:
    df = signature_df_rt_yes[signature_df_rt_yes['OS_Months'] >= t]
    low_risk_df = df[df['Risk_group'] == 0].sort_values(by='Risk_scores')
    high_risk_df = df[df['Risk_group'] == 1].sort_values(by='Risk_scores')
    low_risk_counts.loc[t, 'Count'] = len(low_risk_df)
    high_risk_counts.loc[t, 'Count'] = len(high_risk_df)
    low_risk_counts.loc[t, 'Median risk score'] = low_risk_df['Risk_scores'].median()
    high_risk_counts.loc[t, 'Median risk score'] = high_risk_df['Risk_scores'].median()
    low_risk_counts.loc[t, 'At risk'] = len(low_risk_df)
    high_risk_counts.loc[t, 'At risk'] = len(high_risk_df)

# Create the plot with A4 dimensions
fig, ax = plt.subplots(figsize=(3, 0.5))

# Plot the number of patients at risk for each time point
for index, row in low_risk_counts.iterrows():
    t = index
    ax.text(t, 0.2, f"{int(row['Count'])}", color='blue', fontsize=7, ha='center', va='bottom')

for index, row in high_risk_counts.iterrows():
    t = index
    ax.text(t, 1.2, f"{int(row['Count'])}", color='red', fontsize=7, ha='center', va='bottom')

## Customize plot
plt.xlabel('Time (months)', fontsize=7)
plt.ylabel('Risk group', fontsize=7)

## Set x-axis limits and ticks
plt.xlim(-5, 65)
plt.xticks(time_points, fontsize=7)

## Set y-axis limits, ticks, and tick labels
plt.ylim(-0.3, 1.7)
ax.set_yticks([0.2, 1.2])
ax.set_yticklabels(['Low risk', 'High risk'], color='black', fontsize=7)
ax.get_yticklabels()[0].set_color('blue')
ax.get_yticklabels()[1].set_color('red')

## Remove gridlines and frame
ax.grid(False)
ax.spines[['top', 'right']].set_visible(False)
ax.spines['left'].set_linewidth(0.5)  # Decrease line thickness of the y-axis
ax.spines['bottom'].set_linewidth(0.5)  # Decrease line thickness of the x-axis

## Scale up the plot dimensions to match A4 paper size
scale = 1.0
dpi = 300
fig.set_size_inches(fig.get_size_inches() * scale)

# Save the plot to a PDF file with A4 dimensions
with PdfPages("Colored_Risk_numbers_plot.pdf", keep_empty=False) as pdf:
    pdf.savefig(fig, dpi=dpi, bbox_inches='tight')

# Save the plot as an SVG file with A4 dimensions
fig.savefig('Colored_Risk_numbers_plot.svg', format='svg', dpi=dpi, bbox_inches='tight')

# Show the plot
plt.show()


# %%
print("Low Risk Counts:")
print(low_risk_counts)

print("\nHigh Risk Counts:")
print(high_risk_counts)


# %%
# Define time points
time_points = [0, 10, 20, 30, 40, 50, 60]

# Create empty DataFrames to store counts
low_risk_counts = pd.DataFrame(index=time_points, columns=['Count', 'Median risk score', 'At risk'])
high_risk_counts = pd.DataFrame(index=time_points, columns=['Count', 'Median risk score', 'At risk'])

# Filter patients who received radiation therapy
signature_df_rt_no = signature_df[signature_df['Radiation_therapy'] == 0]

for t in time_points:
    df = signature_df_rt_no[signature_df_rt_no['OS_Months'] >= t]
    low_risk_df = df[df['Risk_group'] == 0].sort_values(by='Risk_scores')
    high_risk_df = df[df['Risk_group'] == 1].sort_values(by='Risk_scores')
    low_risk_counts.loc[t, 'Count'] = len(low_risk_df)
    high_risk_counts.loc[t, 'Count'] = len(high_risk_df)
    low_risk_counts.loc[t, 'Median risk score'] = low_risk_df['Risk_scores'].median()
    high_risk_counts.loc[t, 'Median risk score'] = high_risk_df['Risk_scores'].median()
    low_risk_counts.loc[t, 'At risk'] = len(low_risk_df)
    high_risk_counts.loc[t, 'At risk'] = len(high_risk_df)

# Create the plot with A4 dimensions
fig, ax = plt.subplots(figsize=(3, 0.5))

# Plot the number of patients at risk for each time point
for index, row in low_risk_counts.iterrows():
    t = index
    ax.text(t, 0.2, f"{int(row['Count'])}", color='blue', fontsize=7, ha='center', va='bottom')

for index, row in high_risk_counts.iterrows():
    t = index
    ax.text(t, 1.2, f"{int(row['Count'])}", color='red', fontsize=7, ha='center', va='bottom')

## Customize plot
plt.xlabel('Time (months)', fontsize=7)
plt.ylabel('Risk group', fontsize=7)

## Set x-axis limits and ticks
plt.xlim(-5, 65)
plt.xticks(time_points, fontsize=7)

## Set y-axis limits, ticks, and tick labels
plt.ylim(-0.3, 1.7)
ax.set_yticks([0.2, 1.2])
ax.set_yticklabels(['Low risk', 'High risk'], color='black', fontsize=7)
ax.get_yticklabels()[0].set_color('blue')
ax.get_yticklabels()[1].set_color('red')

## Remove gridlines and frame
ax.grid(False)
ax.spines[['top', 'right']].set_visible(False)
ax.spines['left'].set_linewidth(0.5)  # Decrease line thickness of the y-axis
ax.spines['bottom'].set_linewidth(0.5)  # Decrease line thickness of the x-axis

## Scale up the plot dimensions to match A4 paper size
scale = 1.0
dpi = 300
fig.set_size_inches(fig.get_size_inches() * scale)

# Save the plot to a PDF file with A4 dimensions
with PdfPages("Colored_Risk_numbers_plot.pdf", keep_empty=False) as pdf:
    pdf.savefig(fig, dpi=dpi, bbox_inches='tight')

# Save the plot as an SVG file with A4 dimensions
fig.savefig('Colored_Risk_numbers_plot.svg', format='svg', dpi=dpi, bbox_inches='tight')

# Show the plot
plt.show()


# %%
fig, ax = plt.subplots(figsize=(5, 3))

# Calculate hazard ratios for high-risk Vs low risk
high_vs_low_hr = cph.predict_partial_hazard(signature_df[high_risk]).mean() / cph.predict_partial_hazard(signature_df[low_risk]).mean()

# Calculate log-rank p-value for high-risk and low-risk group
high_risk_p = signature_df.loc[(high_risk) & (signature_df['OS_Months'] <= 60), 'OS_Months']
low_risk_p = signature_df.loc[(low_risk) & (signature_df['OS_Months'] <= 60), 'OS_Months']
results_high_vs_low = logrank_test(high_risk_p, low_risk_p, event_observed_A=signature_df.loc[(high_risk) & (signature_df['OS_Months'] <= 60), 'Vital_status'], event_observed_B=signature_df.loc[(low_risk) & (signature_df['OS_Months'] <= 60), 'Vital_status'])
p_value_high_vs_low = results_high_vs_low.p_value

colors = {1: 'red', 0: 'blue'}
labels = {1: 'High risk', 0: 'Low risk'}

for risk_group in [1, 0]:
    mask = signature_df['Risk_group'] == risk_group
    kmf.fit(signature_df.loc[mask, 'OS_Months'], event_observed=signature_df.loc[mask, 'Vital_status'], label=labels[risk_group])
    kmf.plot(ax=ax, ci_show=True, ci_alpha=0.1, color=colors[risk_group], linewidth=3)

plt.title('Overall survival', fontsize=7)
plt.xlabel('Overall survival (months)', fontsize=7)
plt.ylabel('Survival probability', fontsize=7)
plt.legend(fontsize=7)

# Set x-axis limits to show only up to 60 months
plt.xlim(0, 60)

# Set custom y-axis tick marks
plt.yticks([0, 0.25, 0.5, 0.75, 1], fontsize=7)

# Add the hazard ratio and log-rank p-value
ax.text(0.1, 0.05, f'logrank p high vs low={p_value_high_vs_low:.3f}', transform=ax.transAxes, fontsize=7)
ax.text(0.1, 0.1, f'HR (high-risk): {high_vs_low_hr:.2f}', transform=ax.transAxes, fontsize=7)

# Save the plot to a PDF file
with PdfPages("KM_plot.pdf") as pdf:
    pdf.savefig(fig)

# Save the plot as an SVG file with A4 dimensions
fig.savefig('KM_plot.svg', format='svg', dpi=dpi, bbox_inches='tight')

plt.show()


# %%
# Fit a Cox Proportional Hazards Model
cph = CoxPHFitter(penalizer=0.1)
cph.fit(signature_df, duration_col='OS_Months', event_col='Vital_status')

# Compute risk scores and assign risk group
signature_df['Risk_scores'] = cph.predict_partial_hazard(signature_df)
signature_df['Risk_group'] = np.where(signature_df['Risk_scores'] >= np.median(signature_df['Risk_scores']), 1, 0)

# Use the entire signature_df DataFrame
signature_df_0 = signature_df

# Define time points
time_points = [0, 10, 20, 30, 40, 50, 60]

# Create empty DataFrames to store counts
low_risk_counts = pd.DataFrame(index=time_points, columns=['Count', 'Median risk score', 'At risk'])
high_risk_counts = pd.DataFrame(index=time_points, columns=['Count', 'Median risk score', 'At risk'])

# Loop over time points
for t in time_points:
    df = signature_df_0[signature_df_0['OS_Months'] >= t]
    low_risk_df = df[df['Risk_group'] == 0].sort_values(by='Risk_scores')
    high_risk_df = df[df['Risk_group'] == 1].sort_values(by='Risk_scores')
    low_risk_counts.loc[t, 'Count'] = len(low_risk_df)
    high_risk_counts.loc[t, 'Count'] = len(high_risk_df)
    low_risk_counts.loc[t, 'Median risk score'] = low_risk_df['Risk_scores'].median()
    high_risk_counts.loc[t, 'Median risk score'] = high_risk_df['Risk_scores'].median()
    low_risk_counts.loc[t, 'At risk'] = len(low_risk_df)
    high_risk_counts.loc[t, 'At risk'] = len(high_risk_df)


# Create the plot with A4 dimensions
fig, ax = plt.subplots(figsize=(3, 0.5))

# Plot the number of patients at risk for each time point
for index, row in low_risk_counts.iterrows():
    t = index
    ax.text(t, 0.2, f"{int(row['Count'])}", color='blue', fontsize=7, ha='center', va='bottom')

for index, row in high_risk_counts.iterrows():
    t = index
    ax.text(t, 1.2, f"{int(row['Count'])}", color='red', fontsize=7, ha='center', va='bottom')

## Customize plot
plt.xlabel('Time (months)', fontsize=7)
plt.ylabel('Risk group', fontsize=7)

## Set x-axis limits and ticks
plt.xlim(-5, 65)
plt.xticks(time_points, fontsize=7)

## Set y-axis limits, ticks, and tick labels
plt.ylim(-0.3, 1.7)
ax.set_yticks([0.2, 1.2])
ax.set_yticklabels(['Low risk', 'High risk'], color='black', fontsize=7)
ax.get_yticklabels()[0].set_color('blue')
ax.get_yticklabels()[1].set_color('red')

## Remove gridlines and frame
ax.grid(False)
ax.spines[['top', 'right']].set_visible(False)
ax.spines['left'].set_linewidth(0.5)  # Decrease line thickness of the y-axis
ax.spines['bottom'].set_linewidth(0.5)  # Decrease line thickness of the x-axis

## Scale up the plot dimensions to match A4 paper size
scale = 1.0
dpi = 300
fig.set_size_inches(fig.get_size_inches() * scale)

# Save the plot to a PDF file with A4 dimensions
with PdfPages("Colored_Risk_numbers_plot.pdf", keep_empty=False) as pdf:
    pdf.savefig(fig, dpi=dpi, bbox_inches='tight')

# Save the plot as an SVG file with A4 dimensions
fig.savefig('Colored_Risk_numbers_plot.svg', format='svg', dpi=dpi, bbox_inches='tight')

## Show the plot
plt.show()


# %%
signature_genes = ['ODC1', 'TUBB', 'ETV4', 'NUDT15', 'AK4', 'PPIB',
                   'PLAU', 'LRRC59', 'CTSV', 'HBEGF', 'PPIA', 'MAP4K4', 'AXL']

additional_columns = ['Age', 'Aneuploidy_score', 'Neoplasm_histologic_grade',
                                   'OS_Months', 'Progress_free_survival_Months', 'Vital_status', 'Radiation_therapy']

signature_df = surv_sig_genes_univariable[signature_genes + additional_columns]


# %%
# Create Cox proportional hazards model
cph = CoxPHFitter()
cph.fit(signature_df, duration_col='Progress_free_survival_Months', event_col='Vital_status')
risk_scores = cph.predict_partial_hazard(signature_df)

signature_df['Risk_scores'] = risk_scores
signature_df['Risk_group'] = np.where(signature_df['Risk_scores'] > np.median(signature_df['Risk_scores']), 1, 0)

# Define high_risk and low_risk masks
high_risk = signature_df['Risk_group'] == 1
low_risk = signature_df['Risk_group'] == 0

# Define radiation_therapy mask
rt = signature_df['Radiation_therapy'] == 1
no_rt = signature_df['Radiation_therapy'] == 0


# %%
# Calculate hazard ratios for radiation therapy and high-risk patients
rt_high_risk_hr = cph.predict_partial_hazard(signature_df.loc[high_risk & rt]).mean() / cph.predict_partial_hazard(signature_df.loc[low_risk & rt]).mean()

# Calculate hazard ratios for radiation therapy and low-risk patients
rt_low_risk_hr = cph.predict_partial_hazard(signature_df.loc[low_risk & rt]).mean() / cph.predict_partial_hazard(signature_df.loc[high_risk & rt]).mean()

# Calculate hazard ratios for no radiation therapy and high-risk patients
no_rt_high_risk_hr = cph.predict_partial_hazard(signature_df.loc[high_risk & no_rt]).mean() / cph.predict_partial_hazard(signature_df.loc[low_risk & no_rt]).mean()

# Calculate hazard ratios for no radiation therapy and low-risk patients
no_rt_low_risk_hr = cph.predict_partial_hazard(signature_df.loc[low_risk & no_rt]).mean() / cph.predict_partial_hazard(signature_df.loc[high_risk & no_rt]).mean()


# %%
# Print hazard ratios
print(f'Hazard ratio for radiation therapy: {rt_high_risk_hr:.3f}')
print(f'Hazard ratio for no radiation therapy: {rt_low_risk_hr:.3f}')
print(f'Hazard ratio for radiation therapy: {no_rt_high_risk_hr:.3f}')
print(f'Hazard ratio for no radiation therapy: {no_rt_low_risk_hr:.3f}')


# %%
kmf = KaplanMeierFitter()


# %%
# Calculate log-rank p-value for radiation therapy, high-risk group
T_rt_high_risk = signature_df.loc[(rt) & (high_risk) & (signature_df['Progress_free_survival_Months'] <= 60), 'Progress_free_survival_Months']
E_rt_high_risk = signature_df.loc[(rt) & (high_risk) & (signature_df['Progress_free_survival_Months'] <= 60), 'Vital_status']
T_rt_low_risk = signature_df.loc[(rt) & (low_risk) & (signature_df['Progress_free_survival_Months'] <= 60), 'Progress_free_survival_Months']
E_rt_low_risk = signature_df.loc[(rt) & (low_risk) & (signature_df['Progress_free_survival_Months'] <= 60), 'Vital_status']
results_rt_high_vs_low = logrank_test(T_rt_high_risk, T_rt_low_risk, event_observed_A=E_rt_high_risk, event_observed_B=E_rt_low_risk)
p_value_rt_high_vs_low = results_rt_high_vs_low.p_value


# %%
# Calculate log-rank p-value for no radiation therapy, high-risk group
T_no_rt_high_risk = signature_df.loc[(no_rt) & (high_risk) & (signature_df['Progress_free_survival_Months'] <= 60), 'Progress_free_survival_Months']
E_no_rt_high_risk = signature_df.loc[(no_rt) & (high_risk) & (signature_df['Progress_free_survival_Months'] <= 60), 'Vital_status']
T_no_rt_low_risk = signature_df.loc[(no_rt) & (low_risk) & (signature_df['Progress_free_survival_Months'] <= 60), 'Progress_free_survival_Months']
E_no_rt_low_risk = signature_df.loc[(no_rt) & (low_risk) & (signature_df['Progress_free_survival_Months'] <= 60), 'Vital_status']
results_no_rt_high_vs_low = logrank_test(T_no_rt_high_risk, T_no_rt_low_risk, event_observed_A=E_no_rt_high_risk, event_observed_B=E_no_rt_low_risk)
p_value_no_rt_high_vs_low = results_no_rt_high_vs_low.p_value


# %%
# Labels and line styles for the stratified groups
labels = {
    (1, 1): 'High Risk',
    (1, 0): 'High Risk',
    (0, 1): 'Low Risk',
    (0, 0): 'Low Risk'
}

colors = {
    (1, 1): 'red',
    (1, 0): 'darkred',
    (0, 1): 'blue',
    (0, 0): 'darkblue'
}

line_styles = {
    (1, 1): '-',
    (1, 0): '--',
    (0, 1): '-',
    (0, 0): '--'
}

# Set Arial font with size 7
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 7

# Increase figure size based on font size
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3))

# Create plot for radiation therapy
for risk_group in [1, 0]:
    mask = (signature_df['Risk_group'] == risk_group) & (signature_df['Radiation_therapy'] == 1)
    kmf.fit(signature_df.loc[mask, 'Progress_free_survival_Months'], event_observed=signature_df.loc[mask, 'Vital_status'], label=labels[(risk_group, 1)])
    kmf.plot(ax=ax1, ci_show=True, ci_alpha=0.1, color=colors[(risk_group, 1)], linewidth=3, linestyle=line_styles[(risk_group, 1)], marker='o', markersize=4)

ax1.set_title('Radiation Therapy', fontsize=7)
ax1.set_xlabel('Progress free survival (months)', fontsize=7)
ax1.set_ylabel('Survival probability', fontsize=7)
ax1.legend(fontsize=7)
ax1.set_xlim(0, 60)
ax1.set_yticks([0, 0.25, 0.5, 0.75, 1])
ax1.set_xticks([0, 12, 24, 36, 48, 60])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_linewidth(2)
ax1.spines['left'].set_linewidth(2)
# Add hazard ratios and log-rank p-value to plot
ax1.text(0.5, 0.05, f'logrank p high vs low={p_value_rt_high_vs_low:.3f}', transform=ax1.transAxes, fontsize=7, ha='right', va='bottom')
ax1.text(0.5, 0.1, f'HR (high-risk): {rt_high_risk_hr:.2f}', transform=ax1.transAxes, fontsize=7, ha='right', va='bottom')
ax1.text(0.5, 0.15, f'HR (low-risk): {rt_low_risk_hr:.2f}', transform=ax1.transAxes, fontsize=7, ha='right', va='bottom')


# Create plot for no radiation therapy
for risk_group in [1, 0]:
    mask = (signature_df['Risk_group'] == risk_group) & (signature_df['Radiation_therapy'] == 0)
    kmf.fit(signature_df.loc[mask, 'Progress_free_survival_Months'], event_observed=signature_df.loc[mask, 'Vital_status'], label=labels[(risk_group, 0)])
    kmf.plot(ax=ax2, ci_show=True, ci_alpha=0.1, color=colors[(risk_group, 0)], linewidth=3, linestyle=line_styles[(risk_group, 0)], marker='o', markersize=4)

ax2.set_title('No Radiation Therapy', fontsize=7)
ax2.set_xlabel('Progress free survival (months)', fontsize=7)
ax2.set_ylabel('Survival probability', fontsize=7)
ax2.legend(fontsize=7)
ax2.set_xlim(0, 60)
ax2.set_yticks([0, 0.25, 0.5, 0.75, 1])
ax2.set_xticks([0, 12, 24, 36, 48, 60])
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_linewidth(2)
ax2.spines['left'].set_linewidth(2)
ax2.text(0.5, 0.05, f'logrank p high vs low={p_value_no_rt_high_vs_low:.3f}', transform=ax2.transAxes, fontsize=7, ha='right', va='bottom')
ax2.text(0.5, 0.1, f'HR (high-risk): {no_rt_high_risk_hr:.3f}', transform=ax2.transAxes, fontsize=7, ha='right', va='bottom')
ax2.text(0.5, 0.15, f'HR (low-risk): {no_rt_low_risk_hr:.3f}', transform=ax2.transAxes, fontsize=7, ha='right', va='bottom')

# Save the plot to a PDF file
with PdfPages("KM_plots.pdf") as pdf:
    pdf.savefig(fig)

# Save the plot to an SVG file
plt.savefig("KM_plots.svg", format='svg')

plt.show()


# %%
# Define time points
time_points = [0, 10, 20, 30, 40, 50, 60]

# Create empty DataFrames to store counts
low_risk_counts = pd.DataFrame(index=time_points, columns=['Count', 'Median risk score', 'At risk'])
high_risk_counts = pd.DataFrame(index=time_points, columns=['Count', 'Median risk score', 'At risk'])

# Filter patients who received radiation therapy
signature_df_rt_yes = signature_df[signature_df['Radiation_therapy'] == 1]

for t in time_points:
    df = signature_df_rt_yes[signature_df_rt_yes['Progress_free_survival_Months'] >= t]
    low_risk_df = df[df['Risk_group'] == 0].sort_values(by='Risk_scores')
    high_risk_df = df[df['Risk_group'] == 1].sort_values(by='Risk_scores')
    low_risk_counts.loc[t, 'Count'] = len(low_risk_df)
    high_risk_counts.loc[t, 'Count'] = len(high_risk_df)
    low_risk_counts.loc[t, 'Median risk score'] = low_risk_df['Risk_scores'].median()
    high_risk_counts.loc[t, 'Median risk score'] = high_risk_df['Risk_scores'].median()
    low_risk_counts.loc[t, 'At risk'] = len(low_risk_df)
    high_risk_counts.loc[t, 'At risk'] = len(high_risk_df)

# Create the plot with A4 dimensions
fig, ax = plt.subplots(figsize=(3, 0.5))

# Plot the number of patients at risk for each time point
for index, row in low_risk_counts.iterrows():
    t = index
    ax.text(t, 0.2, f"{int(row['Count'])}", color='blue', fontsize=7, ha='center', va='bottom')

for index, row in high_risk_counts.iterrows():
    t = index
    ax.text(t, 1.2, f"{int(row['Count'])}", color='red', fontsize=7, ha='center', va='bottom')

## Customize plot
plt.xlabel('Time (months)', fontsize=7)
plt.ylabel('Risk group', fontsize=7)

## Set x-axis limits and ticks
plt.xlim(-5, 65)
plt.xticks(time_points, fontsize=7)

## Set y-axis limits, ticks, and tick labels
plt.ylim(-0.3, 1.7)
ax.set_yticks([0.2, 1.2])
ax.set_yticklabels(['Low risk', 'High risk'], color='black', fontsize=7)
ax.get_yticklabels()[0].set_color('blue')
ax.get_yticklabels()[1].set_color('red')

## Remove gridlines and frame
ax.grid(False)
ax.spines[['top', 'right']].set_visible(False)
ax.spines['left'].set_linewidth(0.5)  # Decrease line thickness of the y-axis
ax.spines['bottom'].set_linewidth(0.5)  # Decrease line thickness of the x-axis

## Scale up the plot dimensions to match A4 paper size
scale = 1.0
dpi = 300
fig.set_size_inches(fig.get_size_inches() * scale)

# Save the plot to a PDF file with A4 dimensions
with PdfPages("Colored_Risk_numbers_plot.pdf", keep_empty=False) as pdf:
    pdf.savefig(fig, dpi=dpi, bbox_inches='tight')

# Save the plot as an SVG file with A4 dimensions
fig.savefig('Colored_Risk_numbers_plot.svg', format='svg', dpi=dpi, bbox_inches='tight')

# Show the plot
plt.show()


# %%
# Define time points
time_points = [0, 10, 20, 30, 40, 50, 60]

# Create empty DataFrames to store counts
low_risk_counts = pd.DataFrame(index=time_points, columns=['Count', 'Median risk score', 'At risk'])
high_risk_counts = pd.DataFrame(index=time_points, columns=['Count', 'Median risk score', 'At risk'])

# Filter patients who received radiation therapy
signature_df_rt_no = signature_df[signature_df['Radiation_therapy'] == 0]

for t in time_points:
    df = signature_df_rt_no[signature_df_rt_no['Progress_free_survival_Months'] >= t]
    low_risk_df = df[df['Risk_group'] == 0].sort_values(by='Risk_scores')
    high_risk_df = df[df['Risk_group'] == 1].sort_values(by='Risk_scores')
    low_risk_counts.loc[t, 'Count'] = len(low_risk_df)
    high_risk_counts.loc[t, 'Count'] = len(high_risk_df)
    low_risk_counts.loc[t, 'Median risk score'] = low_risk_df['Risk_scores'].median()
    high_risk_counts.loc[t, 'Median risk score'] = high_risk_df['Risk_scores'].median()
    low_risk_counts.loc[t, 'At risk'] = len(low_risk_df)
    high_risk_counts.loc[t, 'At risk'] = len(high_risk_df)

# Create the plot with A4 dimensions
fig, ax = plt.subplots(figsize=(3, 0.5))

# Plot the number of patients at risk for each time point
for index, row in low_risk_counts.iterrows():
    t = index
    ax.text(t, 0.2, f"{int(row['Count'])}", color='blue', fontsize=7, ha='center', va='bottom')

for index, row in high_risk_counts.iterrows():
    t = index
    ax.text(t, 1.2, f"{int(row['Count'])}", color='red', fontsize=7, ha='center', va='bottom')

## Customize plot
plt.xlabel('Time (months)', fontsize=7)
plt.ylabel('Risk group', fontsize=7)

## Set x-axis limits and ticks
plt.xlim(-5, 65)
plt.xticks(time_points, fontsize=7)

## Set y-axis limits, ticks, and tick labels
plt.ylim(-0.3, 1.7)
ax.set_yticks([0.2, 1.2])
ax.set_yticklabels(['Low risk', 'High risk'], color='black', fontsize=7)
ax.get_yticklabels()[0].set_color('blue')
ax.get_yticklabels()[1].set_color('red')

## Remove gridlines and frame
ax.grid(False)
ax.spines[['top', 'right']].set_visible(False)
ax.spines['left'].set_linewidth(0.5)  # Decrease line thickness of the y-axis
ax.spines['bottom'].set_linewidth(0.5)  # Decrease line thickness of the x-axis

## Scale up the plot dimensions to match A4 paper size
scale = 1.0
dpi = 300
fig.set_size_inches(fig.get_size_inches() * scale)

# Save the plot to a PDF file with A4 dimensions
with PdfPages("Colored_Risk_numbers_plot.pdf", keep_empty=False) as pdf:
    pdf.savefig(fig, dpi=dpi, bbox_inches='tight')

# Save the plot as an SVG file with A4 dimensions
fig.savefig('Colored_Risk_numbers_plot.svg', format='svg', dpi=dpi, bbox_inches='tight')

# Show the plot
plt.show()


# %%
fig, ax = plt.subplots(figsize=(5, 3))

# Calculate hazard ratios for high-risk Vs low risk
high_vs_low_hr = cph.predict_partial_hazard(signature_df[high_risk]).mean() / cph.predict_partial_hazard(signature_df[low_risk]).mean()

# Calculate log-rank p-value for high-risk and low-risk group
high_risk_p = signature_df.loc[(high_risk) & (signature_df['Progress_free_survival_Months'] <= 60), 'Progress_free_survival_Months']
low_risk_p = signature_df.loc[(low_risk) & (signature_df['Progress_free_survival_Months'] <= 60), 'Progress_free_survival_Months']
results_high_vs_low = logrank_test(high_risk_p, low_risk_p, event_observed_A=signature_df.loc[(high_risk) & (signature_df['Progress_free_survival_Months'] <= 60), 'Vital_status'], event_observed_B=signature_df.loc[(low_risk) & (signature_df['Progress_free_survival_Months'] <= 60), 'Vital_status'])
p_value_high_vs_low = results_high_vs_low.p_value

colors = {1: 'red', 0: 'blue'}
labels = {1: 'High risk', 0: 'Low risk'}

for risk_group in [1, 0]:
    mask = signature_df['Risk_group'] == risk_group
    kmf.fit(signature_df.loc[mask, 'Progress_free_survival_Months'], event_observed=signature_df.loc[mask, 'Vital_status'], label=labels[risk_group])
    kmf.plot(ax=ax, ci_show=True, ci_alpha=0.1, color=colors[risk_group], linewidth=3)

plt.title('Progress free survival', fontsize=7)
plt.xlabel('Progress free survival (months)', fontsize=7)
plt.ylabel('Survival probability', fontsize=7)
plt.legend(fontsize=7)

# Set x-axis limits to show only up to 60 months
plt.xlim(0, 60)

# Set custom y-axis tick marks
plt.yticks([0, 0.25, 0.5, 0.75, 1], fontsize=7)

# Add the hazard ratio and log-rank p-value
ax.text(0.1, 0.05, f'logrank p high vs low={p_value_high_vs_low:.3f}', transform=ax.transAxes, fontsize=7)
ax.text(0.1, 0.1, f'HR (high-risk): {high_vs_low_hr:.2f}', transform=ax.transAxes, fontsize=7)

# Save the plot to a PDF file
with PdfPages("KM_plot.pdf") as pdf:
    pdf.savefig(fig)

# Save the plot as an SVG file with A4 dimensions
fig.savefig('KM_plot.svg', format='svg', dpi=dpi, bbox_inches='tight')

plt.show()


# %%
# Fit a Cox Proportional Hazards Model
cph = CoxPHFitter(penalizer=0.1)
cph.fit(signature_df, duration_col='Progress_free_survival_Months', event_col='Vital_status')

# Compute risk scores and assign risk group
signature_df['Risk_scores'] = cph.predict_partial_hazard(signature_df)
signature_df['Risk_group'] = np.where(signature_df['Risk_scores'] >= np.median(signature_df['Risk_scores']), 1, 0)

# Use the entire signature_df DataFrame
signature_df_0 = signature_df

# Define time points
time_points = [0, 10, 20, 30, 40, 50, 60]

# Create empty DataFrames to store counts
low_risk_counts = pd.DataFrame(index=time_points, columns=['Count', 'Median risk score', 'At risk'])
high_risk_counts = pd.DataFrame(index=time_points, columns=['Count', 'Median risk score', 'At risk'])

# Loop over time points
for t in time_points:
    df = signature_df_0[signature_df_0['Progress_free_survival_Months'] >= t]
    low_risk_df = df[df['Risk_group'] == 0].sort_values(by='Risk_scores')
    high_risk_df = df[df['Risk_group'] == 1].sort_values(by='Risk_scores')
    low_risk_counts.loc[t, 'Count'] = len(low_risk_df)
    high_risk_counts.loc[t, 'Count'] = len(high_risk_df)
    low_risk_counts.loc[t, 'Median risk score'] = low_risk_df['Risk_scores'].median()
    high_risk_counts.loc[t, 'Median risk score'] = high_risk_df['Risk_scores'].median()
    low_risk_counts.loc[t, 'At risk'] = len(low_risk_df)
    high_risk_counts.loc[t, 'At risk'] = len(high_risk_df)


# Create the plot with A4 dimensions
fig, ax = plt.subplots(figsize=(3, 0.5))

# Plot the number of patients at risk for each time point
for index, row in low_risk_counts.iterrows():
    t = index
    ax.text(t, 0.2, f"{int(row['Count'])}", color='blue', fontsize=7, ha='center', va='bottom')

for index, row in high_risk_counts.iterrows():
    t = index
    ax.text(t, 1.2, f"{int(row['Count'])}", color='red', fontsize=7, ha='center', va='bottom')

## Customize plot
plt.xlabel('Time (months)', fontsize=7)
plt.ylabel('Risk group', fontsize=7)

## Set x-axis limits and ticks
plt.xlim(-5, 65)
plt.xticks(time_points, fontsize=7)

## Set y-axis limits, ticks, and tick labels
plt.ylim(-0.3, 1.7)
ax.set_yticks([0.2, 1.2])
ax.set_yticklabels(['Low risk', 'High risk'], color='black', fontsize=7)
ax.get_yticklabels()[0].set_color('blue')
ax.get_yticklabels()[1].set_color('red')

## Remove gridlines and frame
ax.grid(False)
ax.spines[['top', 'right']].set_visible(False)
ax.spines['left'].set_linewidth(0.5)  # Decrease line thickness of the y-axis
ax.spines['bottom'].set_linewidth(0.5)  # Decrease line thickness of the x-axis

## Scale up the plot dimensions to match A4 paper size
scale = 1.0
dpi = 300
fig.set_size_inches(fig.get_size_inches() * scale)

# Save the plot to a PDF file with A4 dimensions
with PdfPages("Colored_Risk_numbers_plot.pdf", keep_empty=False) as pdf:
    pdf.savefig(fig, dpi=dpi, bbox_inches='tight')

# Save the plot as an SVG file with A4 dimensions
fig.savefig('Colored_Risk_numbers_plot.svg', format='svg', dpi=dpi, bbox_inches='tight')

## Show the plot
plt.show()


# %%
signature_df = pd.read_csv(data_dir / 'signature_genes_merged_df_for_DSS.txt', sep='\t')


# %%
# Remove rows with NaN values
signature_df.dropna(axis=0, inplace=True)


# %%
# Remove the "Vital_status" column
signature_df.pop("OS_Months")


# %%
# Create Cox proportional hazards model
cph = CoxPHFitter()
cph.fit(signature_df, duration_col='Disease_specific_survival_Months', event_col='Disease_specific_survival_status')
risk_scores = cph.predict_partial_hazard(signature_df)

signature_df['Risk_scores'] = risk_scores
signature_df['Risk_group'] = np.where(signature_df['Risk_scores'] > np.median(signature_df['Risk_scores']), 1, 0)

# Define high_risk and low_risk masks
high_risk = signature_df['Risk_group'] == 1
low_risk = signature_df['Risk_group'] == 0

# Define radiation_therapy mask
rt = signature_df['Radiation_therapy'] == 1
no_rt = signature_df['Radiation_therapy'] == 0


# %%
# Calculate hazard ratios for radiation therapy and high-risk patients
rt_high_risk_hr = cph.predict_partial_hazard(signature_df.loc[high_risk & rt]).mean() / cph.predict_partial_hazard(signature_df.loc[low_risk & rt]).mean()

# Calculate hazard ratios for radiation therapy and low-risk patients
rt_low_risk_hr = cph.predict_partial_hazard(signature_df.loc[low_risk & rt]).mean() / cph.predict_partial_hazard(signature_df.loc[high_risk & rt]).mean()

# Calculate hazard ratios for no radiation therapy and high-risk patients
no_rt_high_risk_hr = cph.predict_partial_hazard(signature_df.loc[high_risk & no_rt]).mean() / cph.predict_partial_hazard(signature_df.loc[low_risk & no_rt]).mean()

# Calculate hazard ratios for no radiation therapy and low-risk patients
no_rt_low_risk_hr = cph.predict_partial_hazard(signature_df.loc[low_risk & no_rt]).mean() / cph.predict_partial_hazard(signature_df.loc[high_risk & no_rt]).mean()


# %%
# Print hazard ratios
print(f'Hazard ratio for radiation therapy: {rt_high_risk_hr:.3f}')
print(f'Hazard ratio for no radiation therapy: {rt_low_risk_hr:.3f}')
print(f'Hazard ratio for radiation therapy: {no_rt_high_risk_hr:.3f}')
print(f'Hazard ratio for no radiation therapy: {no_rt_low_risk_hr:.3f}')


# %%
kmf = KaplanMeierFitter()


# %%
# Calculate log-rank p-value for radiation therapy, high-risk group
T_rt_high_risk = signature_df.loc[(rt) & (high_risk) & (signature_df['Disease_specific_survival_Months'] <= 60), 'Disease_specific_survival_Months']
E_rt_high_risk = signature_df.loc[(rt) & (high_risk) & (signature_df['Disease_specific_survival_Months'] <= 60), 'Disease_specific_survival_status']
T_rt_low_risk = signature_df.loc[(rt) & (low_risk) & (signature_df['Disease_specific_survival_Months'] <= 60), 'Disease_specific_survival_Months']
E_rt_low_risk = signature_df.loc[(rt) & (low_risk) & (signature_df['Disease_specific_survival_Months'] <= 60), 'Disease_specific_survival_status']
results_rt_high_vs_low = logrank_test(T_rt_high_risk, T_rt_low_risk, event_observed_A=E_rt_high_risk, event_observed_B=E_rt_low_risk)
p_value_rt_high_vs_low = results_rt_high_vs_low.p_value


# %%
# Calculate log-rank p-value for no radiation therapy, high-risk group
T_no_rt_high_risk = signature_df.loc[(no_rt) & (high_risk) & (signature_df['Disease_specific_survival_Months'] <= 60), 'Disease_specific_survival_Months']
E_no_rt_high_risk = signature_df.loc[(no_rt) & (high_risk) & (signature_df['Disease_specific_survival_Months'] <= 60), 'Disease_specific_survival_status']
T_no_rt_low_risk = signature_df.loc[(no_rt) & (low_risk) & (signature_df['Disease_specific_survival_Months'] <= 60), 'Disease_specific_survival_Months']
E_no_rt_low_risk = signature_df.loc[(no_rt) & (low_risk) & (signature_df['Disease_specific_survival_Months'] <= 60), 'Disease_specific_survival_status']
results_no_rt_high_vs_low = logrank_test(T_no_rt_high_risk, T_no_rt_low_risk, event_observed_A=E_no_rt_high_risk, event_observed_B=E_no_rt_low_risk)
p_value_no_rt_high_vs_low = results_no_rt_high_vs_low.p_value


# %%
# Labels and line styles for the stratified groups
labels = {
    (1, 1): 'High Risk',
    (1, 0): 'High Risk',
    (0, 1): 'Low Risk',
    (0, 0): 'Low Risk'
}

colors = {
    (1, 1): 'red',
    (1, 0): 'darkred',
    (0, 1): 'blue',
    (0, 0): 'darkblue'
}

line_styles = {
    (1, 1): '-',
    (1, 0): '--',
    (0, 1): '-',
    (0, 0): '--'
}

# Set Arial font with size 7
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 7

# Increase figure size based on font size
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3))

# Create plot for radiation therapy
for risk_group in [1, 0]:
    mask = (signature_df['Risk_group'] == risk_group) & (signature_df['Radiation_therapy'] == 1)
    kmf.fit(signature_df.loc[mask, 'Disease_specific_survival_Months'], event_observed=signature_df.loc[mask, 'Disease_specific_survival_status'], label=labels[(risk_group, 1)])
    kmf.plot(ax=ax1, ci_show=True, ci_alpha=0.1, color=colors[(risk_group, 1)], linewidth=3, linestyle=line_styles[(risk_group, 1)], marker='o', markersize=4)

ax1.set_title('Radiation Therapy', fontsize=7)
ax1.set_xlabel('Disease specific survival (months)', fontsize=7)
ax1.set_ylabel('Survival probability', fontsize=7)
ax1.legend(fontsize=7)
ax1.set_xlim(0, 60)
ax1.set_yticks([0, 0.25, 0.5, 0.75, 1])
ax1.set_xticks([0, 12, 24, 36, 48, 60])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_linewidth(2)
ax1.spines['left'].set_linewidth(2)
# Add hazard ratios and log-rank p-value to plot
ax1.text(0.5, 0.05, f'logrank p high vs low={p_value_rt_high_vs_low:.3f}', transform=ax1.transAxes, fontsize=7, ha='right', va='bottom')
ax1.text(0.5, 0.1, f'HR (high-risk): {rt_high_risk_hr:.2f}', transform=ax1.transAxes, fontsize=7, ha='right', va='bottom')
ax1.text(0.5, 0.15, f'HR (low-risk): {rt_low_risk_hr:.2f}', transform=ax1.transAxes, fontsize=7, ha='right', va='bottom')


# Create plot for no radiation therapy
for risk_group in [1, 0]:
    mask = (signature_df['Risk_group'] == risk_group) & (signature_df['Radiation_therapy'] == 0)
    kmf.fit(signature_df.loc[mask, 'Disease_specific_survival_Months'], event_observed=signature_df.loc[mask, 'Disease_specific_survival_status'], label=labels[(risk_group, 0)])
    kmf.plot(ax=ax2, ci_show=True, ci_alpha=0.1, color=colors[(risk_group, 0)], linewidth=3, linestyle=line_styles[(risk_group, 0)], marker='o', markersize=4)

ax2.set_title('No Radiation Therapy', fontsize=7)
ax2.set_xlabel('Disease specific survival (months)', fontsize=7)
ax2.set_ylabel('Survival probability', fontsize=7)
ax2.legend(fontsize=7)
ax2.set_xlim(0, 60)
ax2.set_yticks([0, 0.25, 0.5, 0.75, 1])
ax2.set_xticks([0, 12, 24, 36, 48, 60])
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_linewidth(2)
ax2.spines['left'].set_linewidth(2)
ax2.text(0.5, 0.05, f'logrank p high vs low={p_value_no_rt_high_vs_low:.3f}', transform=ax2.transAxes, fontsize=7, ha='right', va='bottom')
ax2.text(0.5, 0.1, f'HR (high-risk): {no_rt_high_risk_hr:.3f}', transform=ax2.transAxes, fontsize=7, ha='right', va='bottom')
ax2.text(0.5, 0.15, f'HR (low-risk): {no_rt_low_risk_hr:.3f}', transform=ax2.transAxes, fontsize=7, ha='right', va='bottom')

# Save the plot to a PDF file
with PdfPages("KM_plots.pdf") as pdf:
    pdf.savefig(fig)

# Save the plot to an SVG file
plt.savefig("KM_plots.svg", format='svg')

plt.show()


# %%
# Define time points
time_points = [0, 10, 20, 30, 40, 50, 60]

# Create empty DataFrames to store counts
low_risk_counts = pd.DataFrame(index=time_points, columns=['Count', 'Median risk score', 'At risk'])
high_risk_counts = pd.DataFrame(index=time_points, columns=['Count', 'Median risk score', 'At risk'])

# Filter patients who received radiation therapy
signature_df_rt_yes = signature_df[signature_df['Radiation_therapy'] == 1]

for t in time_points:
    df = signature_df_rt_yes[signature_df_rt_yes['Disease_specific_survival_Months'] >= t]
    low_risk_df = df[df['Risk_group'] == 0].sort_values(by='Risk_scores')
    high_risk_df = df[df['Risk_group'] == 1].sort_values(by='Risk_scores')
    low_risk_counts.loc[t, 'Count'] = len(low_risk_df)
    high_risk_counts.loc[t, 'Count'] = len(high_risk_df)
    low_risk_counts.loc[t, 'Median risk score'] = low_risk_df['Risk_scores'].median()
    high_risk_counts.loc[t, 'Median risk score'] = high_risk_df['Risk_scores'].median()
    low_risk_counts.loc[t, 'At risk'] = len(low_risk_df)
    high_risk_counts.loc[t, 'At risk'] = len(high_risk_df)

# Create the plot with A4 dimensions
fig, ax = plt.subplots(figsize=(3, 0.5))

# Plot the number of patients at risk for each time point
for index, row in low_risk_counts.iterrows():
    t = index
    ax.text(t, 0.2, f"{int(row['Count'])}", color='blue', fontsize=7, ha='center', va='bottom')

for index, row in high_risk_counts.iterrows():
    t = index
    ax.text(t, 1.2, f"{int(row['Count'])}", color='red', fontsize=7, ha='center', va='bottom')

## Customize plot
plt.xlabel('Time (months)', fontsize=7)
plt.ylabel('Risk group', fontsize=7)

## Set x-axis limits and ticks
plt.xlim(-5, 65)
plt.xticks(time_points, fontsize=7)

## Set y-axis limits, ticks, and tick labels
plt.ylim(-0.3, 1.7)
ax.set_yticks([0.2, 1.2])
ax.set_yticklabels(['Low risk', 'High risk'], color='black', fontsize=7)
ax.get_yticklabels()[0].set_color('blue')
ax.get_yticklabels()[1].set_color('red')

## Remove gridlines and frame
ax.grid(False)
ax.spines[['top', 'right']].set_visible(False)
ax.spines['left'].set_linewidth(0.5)  # Decrease line thickness of the y-axis
ax.spines['bottom'].set_linewidth(0.5)  # Decrease line thickness of the x-axis

## Scale up the plot dimensions to match A4 paper size
scale = 1.0
dpi = 300
fig.set_size_inches(fig.get_size_inches() * scale)

# Save the plot to a PDF file with A4 dimensions
with PdfPages("DSS_rt_Risk_numbers_plot.pdf", keep_empty=False) as pdf:
    pdf.savefig(fig, dpi=dpi, bbox_inches='tight')

# Save the plot as an SVG file with A4 dimensions
fig.savefig('DSS_rt_Risk_numbers_plot.svg', format='svg', dpi=dpi, bbox_inches='tight')

# Show the plot
plt.show()


# %%
# Define time points
time_points = [0, 10, 20, 30, 40, 50, 60]

# Create empty DataFrames to store counts
low_risk_counts = pd.DataFrame(index=time_points, columns=['Count', 'Median risk score', 'At risk'])
high_risk_counts = pd.DataFrame(index=time_points, columns=['Count', 'Median risk score', 'At risk'])

# Filter patients who received radiation therapy
signature_df_rt_no = signature_df[signature_df['Radiation_therapy'] == 0]

for t in time_points:
    df = signature_df_rt_no[signature_df_rt_no['Disease_specific_survival_Months'] >= t]
    low_risk_df = df[df['Risk_group'] == 0].sort_values(by='Risk_scores')
    high_risk_df = df[df['Risk_group'] == 1].sort_values(by='Risk_scores')
    low_risk_counts.loc[t, 'Count'] = len(low_risk_df)
    high_risk_counts.loc[t, 'Count'] = len(high_risk_df)
    low_risk_counts.loc[t, 'Median risk score'] = low_risk_df['Risk_scores'].median()
    high_risk_counts.loc[t, 'Median risk score'] = high_risk_df['Risk_scores'].median()
    low_risk_counts.loc[t, 'At risk'] = len(low_risk_df)
    high_risk_counts.loc[t, 'At risk'] = len(high_risk_df)

# Create the plot with A4 dimensions
fig, ax = plt.subplots(figsize=(3, 0.5))

# Plot the number of patients at risk for each time point
for index, row in low_risk_counts.iterrows():
    t = index
    ax.text(t, 0.2, f"{int(row['Count'])}", color='blue', fontsize=7, ha='center', va='bottom')

for index, row in high_risk_counts.iterrows():
    t = index
    ax.text(t, 1.2, f"{int(row['Count'])}", color='red', fontsize=7, ha='center', va='bottom')

## Customize plot
plt.xlabel('Time (months)', fontsize=7)
plt.ylabel('Risk group', fontsize=7)

## Set x-axis limits and ticks
plt.xlim(-5, 65)
plt.xticks(time_points, fontsize=7)

## Set y-axis limits, ticks, and tick labels
plt.ylim(-0.3, 1.7)
ax.set_yticks([0.2, 1.2])
ax.set_yticklabels(['Low risk', 'High risk'], color='black', fontsize=7)
ax.get_yticklabels()[0].set_color('blue')
ax.get_yticklabels()[1].set_color('red')

## Remove gridlines and frame
ax.grid(False)
ax.spines[['top', 'right']].set_visible(False)
ax.spines['left'].set_linewidth(0.5)  # Decrease line thickness of the y-axis
ax.spines['bottom'].set_linewidth(0.5)  # Decrease line thickness of the x-axis

## Scale up the plot dimensions to match A4 paper size
scale = 1.0
dpi = 300
fig.set_size_inches(fig.get_size_inches() * scale)

# Save the plot to a PDF file with A4 dimensions
with PdfPages("DSS_no_rt_Risk_numbers_plot.pdf", keep_empty=False) as pdf:
    pdf.savefig(fig, dpi=dpi, bbox_inches='tight')

# Save the plot as an SVG file with A4 dimensions
fig.savefig('DSS_no_rt_Risk_numbers_plot.svg', format='svg', dpi=dpi, bbox_inches='tight')

# Show the plot
plt.show()


# %%
fig, ax = plt.subplots(figsize=(5, 3))

# Calculate hazard ratios for high-risk Vs low risk
high_vs_low_hr = cph.predict_partial_hazard(signature_df[high_risk]).mean() / cph.predict_partial_hazard(signature_df[low_risk]).mean()

# Calculate log-rank p-value for high-risk and low-risk group
high_risk_p = signature_df.loc[(high_risk) & (signature_df['Disease_specific_survival_Months'] <= 60), 'Disease_specific_survival_Months']
low_risk_p = signature_df.loc[(low_risk) & (signature_df['Disease_specific_survival_Months'] <= 60), 'Disease_specific_survival_Months']
results_high_vs_low = logrank_test(high_risk_p, low_risk_p, event_observed_A=signature_df.loc[(high_risk) & (signature_df['Disease_specific_survival_Months'] <= 60), 'Disease_specific_survival_status'], event_observed_B=signature_df.loc[(low_risk) & (signature_df['Disease_specific_survival_Months'] <= 60), 'Disease_specific_survival_status'])
p_value_high_vs_low = results_high_vs_low.p_value

colors = {1: 'red', 0: 'blue'}
labels = {1: 'High risk', 0: 'Low risk'}

for risk_group in [1, 0]:
    mask = signature_df['Risk_group'] == risk_group
    kmf.fit(signature_df.loc[mask, 'Disease_specific_survival_Months'], event_observed=signature_df.loc[mask, 'Disease_specific_survival_status'], label=labels[risk_group])
    kmf.plot(ax=ax, ci_show=True, ci_alpha=0.1, color=colors[risk_group], linewidth=3)

plt.title('Disease specific survival', fontsize=7)
plt.xlabel('Disease specific survival (months)', fontsize=7)
plt.ylabel('Survival probability', fontsize=7)
plt.legend(fontsize=7)

# Set x-axis limits to show only up to 60 months
plt.xlim(0, 60)

# Set custom y-axis tick marks
plt.yticks([0, 0.25, 0.5, 0.75, 1], fontsize=7)

# Add the hazard ratio and log-rank p-value
ax.text(0.1, 0.05, f'logrank p high vs low={p_value_high_vs_low:.3f}', transform=ax.transAxes, fontsize=7)
ax.text(0.1, 0.1, f'HR (high-risk): {high_vs_low_hr:.2f}', transform=ax.transAxes, fontsize=7)

# Save the plot to a PDF file
with PdfPages("DSS_Overall_KM_plot.pdf") as pdf:
    pdf.savefig(fig)

# Save the plot as an SVG file with A4 dimensions
fig.savefig('DSS_Overall_KM_plot.svg', format='svg', dpi=dpi, bbox_inches='tight')

plt.show()


# %%
# Fit a Cox Proportional Hazards Model
cph = CoxPHFitter(penalizer=0.1)
cph.fit(signature_df, duration_col='Disease_specific_survival_Months', event_col='Disease_specific_survival_status')

# Compute risk scores and assign risk group
signature_df['Risk_scores'] = cph.predict_partial_hazard(signature_df)
signature_df['Risk_group'] = np.where(signature_df['Risk_scores'] >= np.median(signature_df['Risk_scores']), 1, 0)

# Use the entire signature_df DataFrame
signature_df_0 = signature_df

# Define time points
time_points = [0, 10, 20, 30, 40, 50, 60]

# Create empty DataFrames to store counts
low_risk_counts = pd.DataFrame(index=time_points, columns=['Count', 'Median risk score', 'At risk'])
high_risk_counts = pd.DataFrame(index=time_points, columns=['Count', 'Median risk score', 'At risk'])

# Loop over time points
for t in time_points:
    df = signature_df_0[signature_df_0['Disease_specific_survival_Months'] >= t]
    low_risk_df = df[df['Risk_group'] == 0].sort_values(by='Risk_scores')
    high_risk_df = df[df['Risk_group'] == 1].sort_values(by='Risk_scores')
    low_risk_counts.loc[t, 'Count'] = len(low_risk_df)
    high_risk_counts.loc[t, 'Count'] = len(high_risk_df)
    low_risk_counts.loc[t, 'Median risk score'] = low_risk_df['Risk_scores'].median()
    high_risk_counts.loc[t, 'Median risk score'] = high_risk_df['Risk_scores'].median()
    low_risk_counts.loc[t, 'At risk'] = len(low_risk_df)
    high_risk_counts.loc[t, 'At risk'] = len(high_risk_df)


# Create the plot with A4 dimensions
fig, ax = plt.subplots(figsize=(3, 0.5))

# Plot the number of patients at risk for each time point
for index, row in low_risk_counts.iterrows():
    t = index
    ax.text(t, 0.2, f"{int(row['Count'])}", color='blue', fontsize=7, ha='center', va='bottom')

for index, row in high_risk_counts.iterrows():
    t = index
    ax.text(t, 1.2, f"{int(row['Count'])}", color='red', fontsize=7, ha='center', va='bottom')

## Customize plot
plt.xlabel('Time (months)', fontsize=7)
plt.ylabel('Risk group', fontsize=7)

## Set x-axis limits and ticks
plt.xlim(-5, 65)
plt.xticks(time_points, fontsize=7)

## Set y-axis limits, ticks, and tick labels
plt.ylim(-0.3, 1.7)
ax.set_yticks([0.2, 1.2])
ax.set_yticklabels(['Low risk', 'High risk'], color='black', fontsize=7)
ax.get_yticklabels()[0].set_color('blue')
ax.get_yticklabels()[1].set_color('red')

## Remove gridlines and frame
ax.grid(False)
ax.spines[['top', 'right']].set_visible(False)
ax.spines['left'].set_linewidth(0.5)  # Decrease line thickness of the y-axis
ax.spines['bottom'].set_linewidth(0.5)  # Decrease line thickness of the x-axis

## Scale up the plot dimensions to match A4 paper size
scale = 1.0
dpi = 300
fig.set_size_inches(fig.get_size_inches() * scale)

# Save the plot to a PDF file with A4 dimensions
with PdfPages("DSS_Overall_Risk_numbers_plot.pdf", keep_empty=False) as pdf:
    pdf.savefig(fig, dpi=dpi, bbox_inches='tight')

# Save the plot as an SVG file with A4 dimensions
fig.savefig('DSS_Overall_Risk_numbers_plot.svg', format='svg', dpi=dpi, bbox_inches='tight')

## Show the plot
plt.show()
