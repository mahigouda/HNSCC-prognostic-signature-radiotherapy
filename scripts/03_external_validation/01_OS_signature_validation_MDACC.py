# =============================================================================
# Project : Prognostic gene-expression signature for radiotherapy outcome in HPV-negative HNSCC
# Module  : 03 | External validation of the signatures (MDACC, FHCRC cohorts)
# Script  : 01_OS_signature_validation_MDACC.py
# Purpose : Application of the 13-gene OS signature to the independent MDACC cohort:
#           collinearity check (correlation, VIF), risk scores from TCGA-derived coefficients
#           and cohort-refitted Cox models, median split, Kaplan-Meier/log-rank,
#           number-at-risk tables and ROC/AUC at 60 months.
# Input   : data/validation_cohorts/exp_mdacc.txt, clinical_mdacc_OSCC.txt
# Output  : results/03_external_validation/ (KM, risk-number and ROC plots)
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
from sklearn.metrics import roc_curve, auc


# %%
# ---- Paths -------------------------------------------------------------------
# Run from the repository root. Inputs are read from data/validation_cohorts/
# (see data/README.md); outputs are written to results/03_external_validation/.
from pathlib import Path
data_dir = Path.cwd() / "data" / "validation_cohorts"
out_dir = Path.cwd() / "results" / "03_external_validation"
out_dir.mkdir(parents=True, exist_ok=True)
os.chdir(out_dir)


# %%
exp_df_mdacc = pd.read_csv(data_dir / 'exp_mdacc.txt', sep='\t')


# %%
# Transpose exp_df_mdacc
exp_df_mdacc = exp_df_mdacc.T


# %%
clinical_df_mdacc = pd.read_csv(data_dir / 'clinical_mdacc_OSCC.txt', sep='\t')


# %%
clinical_df_mdacc.set_index(clinical_df_mdacc.columns[0], inplace=True)


# %%
merged_df = exp_df_mdacc.merge(clinical_df_mdacc, how='inner', left_index=True, right_index=True)


# %%
# convert weeks to months
merged_df['OS_Months'] = merged_df['OS_Months'] / 12


# %%
## Remove the decimal part and keep the integer part of the "OS_Months" column
merged_df['OS_Months'] = merged_df['OS_Months'].apply(lambda x: int(x) if not np.isnan(x) else x)


# %%
# Remove the "treatment" column
merged_df.pop("treatment")


# %%
# Remove rows with NaNs
merged_df = merged_df.dropna()


# %%
df = merged_df


# %%
# Convert 'Vital_status' column to numeric representation
df['Vital_status'] = df['Vital_status'].map({'Alive': 0, 'Dead_OSCC': 1})


# %%
signature_genes = ['ODC1', 'TUBB', 'ETV4', 'NUDT15', 'AK4', 'PPIB',
                   'PLAU', 'LRRC59', 'CTSV', 'HBEGF', 'PPIA', 'MAP4K4', 'AXL']

additional_columns = ['Age', 'OS_Months', 'Vital_status']

signature_df = df[signature_genes + additional_columns]


# %%
## Remove index column
signature_df = signature_df.reset_index(drop=True)


# %%
# Remove rows with NaNs
signature_df = signature_df.dropna()


# %%
coef_dict = {
    'ODC1': 0.1992315970490915,
    'TUBB': 0.0834969947009525,
    'ETV4': 0.09505665026324066,
    'NUDT15': 0.09398444726033238,
    'AK4': 0.037355656003152944,
    'PPIB': 0.07645324849618733,
    'PLAU': 0.08081896172888454,
    'LRRC59': 0.022636856915514706,
    'CTSV': 0.057842885465445026,
    'HBEGF': 0.05605996171272654,
    'PPIA': -0.09441807186273042,
    'MAP4K4': -0.08514377784224743,
    'AXL': -0.15789917202480502
}


# %%
# Calculate correlation matrix
corr_matrix = signature_df[signature_genes].corr().abs()

# Select upper triangle of correlation matrix
upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(np.bool))

# Find index of feature columns with correlation greater than 0.95 (this is just an example, you can choose your threshold)
to_drop = [column for column in upper.columns if any(upper[column] > 0.95)]

# Drop highly correlated features
reduced_df = signature_df.drop(signature_df[to_drop], axis=1)


# %%
from statsmodels.stats.outliers_influence import variance_inflation_factor

vif = pd.DataFrame()
vif["VIF Factor"] = [variance_inflation_factor(signature_df[signature_genes].values, i) for i in range(signature_df[signature_genes].shape[1])]
vif["features"] = signature_df[signature_genes].columns


# %%
# Drop highly correlated features from signature_genes
signature_genes = [gene for gene in signature_genes if gene not in to_drop]

def calculate_risk(row):
    return sum(row[gene] * coef_dict[gene] for gene in signature_genes)

signature_df['Risk_scores'] = signature_df[signature_genes].apply(calculate_risk, axis=1)
signature_df['Risk_group'] = np.where(signature_df['Risk_scores'] > np.median(signature_df['Risk_scores']), 1, 0)

# Define high_risk and low_risk masks
high_risk = signature_df['Risk_group'] == 1
low_risk = signature_df['Risk_group'] == 0

# Create Cox proportional hazards model
cph = CoxPHFitter()

# The 'fit' method fits the Cox proportional hazard regression model on your dataset
cph.fit(signature_df[signature_genes + ['OS_Months', 'Vital_status']], duration_col='OS_Months', event_col='Vital_status')

# The 'print_summary' method will give you a nice summary of the fit
cph.print_summary()


# %%
# Create a KaplanMeierFitter instance
kmf = KaplanMeierFitter()


# %%
fig, ax = plt.subplots(figsize=(5, 3))

# Calculate hazard ratios for high-risk Vs low risk
high_vs_low_hr = cph.predict_partial_hazard(signature_df[high_risk]).mean() / cph.predict_partial_hazard(signature_df[low_risk]).mean()

# Calculate log-rank p-value for high-risk and low-risk group
results_high_vs_low = logrank_test(signature_df[high_risk]['OS_Months'], signature_df[low_risk]['OS_Months'],
                                   event_observed_A=signature_df[high_risk]['Vital_status'],
                                   event_observed_B=signature_df[low_risk]['Vital_status'])
p_value_high_vs_low = results_high_vs_low.p_value

colors = {1: 'red', 0: 'blue'}
labels = {1: 'High risk', 0: 'Low risk'}

# Kaplan-Meier Survival Curve
for risk_group in [1, 0]:
    mask = signature_df['Risk_group'] == risk_group
    kmf.fit(signature_df.loc[mask, 'OS_Months'], event_observed=signature_df.loc[mask, 'Vital_status'], label=labels[risk_group])
    kmf.plot(ax=ax, ci_show=True, ci_alpha=0.1, color=colors[risk_group], linewidth=3)

plt.title('Overall survival', fontsize=7)
plt.xlabel('Overall survival (months)', fontsize=7)
plt.ylabel('Survival probability', fontsize=7)
plt.legend(fontsize=7)

# Set custom y-axis tick marks
plt.yticks([0, 0.25, 0.5, 0.75, 1], fontsize=7)

# Add the hazard ratio and log-rank p-value
ax.text(0.3, 0.05, f'logrank p high vs low={p_value_high_vs_low:.3f}', transform=ax.transAxes, fontsize=7)
ax.text(0.3, 0.1, f'HR (high-risk): {high_vs_low_hr:.2f}', transform=ax.transAxes, fontsize=7)

# Restrict x-axis to 60 months
plt.xlim(0, 60)

# Save the plot to a PDF file
with PdfPages("Signature_genes_MDACC_KM_plot.pdf") as pdf:
    pdf.savefig(fig)

dpi = 300
# Save the plot as an SVG file with A4 dimensions
fig.savefig('Signature_genes_MDACC_plot.svg', format='svg', dpi=dpi, bbox_inches='tight')

plt.show()


# %%
# Fit a Cox Proportional Hazards Model
cph = CoxPHFitter(penalizer=0.1)
cph.fit(signature_df, duration_col='OS_Months', event_col='Vital_status')

# Compute risk scores and assign risk group
signature_df['Risk_scores'] = cph.predict_partial_hazard(signature_df)
signature_df['Risk_group'] = np.where(signature_df['Risk_scores'] >= np.median(signature_df['Risk_scores']), 1, 0)

# Define time points
time_points = [0, 10, 20, 30, 40, 50, 60]

# Create empty DataFrames to store counts
low_risk_counts = pd.DataFrame(index=time_points, columns=['Count', 'Median risk score', 'At risk'])
high_risk_counts = pd.DataFrame(index=time_points, columns=['Count', 'Median risk score', 'At risk'])

# Loop over time points
for t in time_points:
    df = signature_df[signature_df['OS_Months'] >= t]
    low_risk_df = df[df['Risk_group'] == 0].sort_values(by='Risk_scores')
    high_risk_df = df[df['Risk_group'] == 1].sort_values(by='Risk_scores')
    low_risk_counts.loc[t, 'Count'] = len(low_risk_df)
    high_risk_counts.loc[t, 'Count'] = len(high_risk_df)
    low_risk_counts.loc[t, 'Median risk score'] = low_risk_df['Risk_scores'].median()
    high_risk_counts.loc[t, 'Median risk score'] = high_risk_df['Risk_scores'].median()
    low_risk_counts.loc[t, 'At risk'] = len(low_risk_df)
    high_risk_counts.loc[t, 'At risk'] = len(high_risk_df)

# Create the plot
fig, ax = plt.subplots(figsize=(3, 0.5))

# Plot the number of patients at risk for each time point
for index, row in low_risk_counts.iterrows():
    t = index
    ax.text(t, 0.2, f"{int(row['Count'])}", color='blue', fontsize=7, ha='center', va='bottom')

for index, row in high_risk_counts.iterrows():
    t = index
    ax.text(t, 1.2, f"{int(row['Count'])}", color='red', fontsize=7, ha='center', va='bottom')

# Customize plot
plt.xlabel('Time (months)', fontsize=7)
plt.ylabel('Risk group', fontsize=7)

# Set x-axis limits and ticks
plt.xlim(-10, 65)
plt.xticks(time_points, fontsize=7)

# Set y-axis limits, ticks, and tick labels
plt.ylim(-0.3, 1.7)
ax.set_yticks([0.2, 1.2])
ax.set_yticklabels(['Low risk', 'High risk'], color='black', fontsize=7)
ax.get_yticklabels()[0].set_color('blue')
ax.get_yticklabels()[1].set_color('red')

# Remove gridlines and frame
ax.grid(False)
ax.spines[['top', 'right']].set_visible(False)
ax.spines['left'].set_linewidth(0.5)  # Decrease line thickness of the y-axis
ax.spines['bottom'].set_linewidth(0.5)  # Decrease line thickness of the x-axis

# Scale up the plot dimensions to match A4 paper size
scale = 1.0
dpi = 300
fig.set_size_inches(fig.get_size_inches() * scale)

# Save the plot to a PDF file with A4 dimensions
with PdfPages("Signature_genes_MDACC_Risk_numbers_plot.pdf", keep_empty=False) as pdf:
    pdf.savefig(fig, dpi=dpi, bbox_inches='tight')

# Save the plot as an SVG file with A4 dimensions
fig.savefig('Signature_genes_MDACC_Risk_numbers_plot.svg', format='svg', dpi=dpi, bbox_inches='tight')

# Show the plot
plt.show()


# %%
# set a time point
time_point = 60

# create binary labels
y_true = (signature_df['OS_Months'] <= time_point) & (signature_df['Vital_status'] == 1)

# calculate the roc curve
fpr, tpr, thresholds = roc_curve(y_true, signature_df['Risk_scores'])

# calculate the auc
roc_auc = auc(fpr, tpr)

fig, ax = plt.subplots(figsize=(5, 3))

# plot the roc curve
ax.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
ax.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
ax.set_xlim([0.0, 1.0])
ax.set_ylim([0.0, 1.05])
ax.set_xlabel('False Positive Rate', fontsize=7)
ax.set_ylabel('True Positive Rate', fontsize=7)
ax.set_title('Receiver operating characteristic example', fontsize=7)
ax.legend(loc="lower right", fontsize=7)

# Save the plot to a PDF file
with PdfPages("ROC_plot.pdf") as pdf:
    pdf.savefig(fig)

dpi = 300
# Save the plot as an SVG file with A4 dimensions
fig.savefig('ROC_plot.svg', format='svg', dpi=dpi, bbox_inches='tight')

plt.show()


# %%
fig, ax = plt.subplots(figsize=(5, 3))

# Calculate hazard ratios for high-risk Vs low risk
high_vs_low_hr = cph.predict_partial_hazard(signature_df[high_risk]).mean() / cph.predict_partial_hazard(signature_df[low_risk]).mean()

# Calculate log-rank p-value for high-risk and low-risk group
results_high_vs_low = logrank_test(signature_df[high_risk]['OS_Months'], signature_df[low_risk]['OS_Months'],
                                   event_observed_A=signature_df[high_risk]['Vital_status'],
                                   event_observed_B=signature_df[low_risk]['Vital_status'])
p_value_high_vs_low = results_high_vs_low.p_value

colors = {1: 'red', 0: 'blue'}
labels = {1: 'High risk', 0: 'Low risk'}

# Kaplan-Meier Survival Curve
for risk_group in [1, 0]:
    mask = signature_df['Risk_group'] == risk_group
    kmf.fit(signature_df.loc[mask, 'OS_Months'], event_observed=signature_df.loc[mask, 'Vital_status'], label=labels[risk_group])
    kmf.plot(ax=ax, ci_show=True, ci_alpha=0.1, color=colors[risk_group], linewidth=3)

plt.title('Overall survival', fontsize=7)
plt.xlabel('Overall survival (months)', fontsize=7)
plt.ylabel('Survival probability', fontsize=7)
plt.legend(fontsize=7)

# Set custom y-axis tick marks
plt.yticks([0, 0.25, 0.5, 0.75, 1], fontsize=7)

# Add the hazard ratio and log-rank p-value
ax.text(0.1, 0.05, f'logrank p high vs low={p_value_high_vs_low:.3f}', transform=ax.transAxes, fontsize=7)
ax.text(0.1, 0.1, f'HR (high-risk): {high_vs_low_hr:.2f}', transform=ax.transAxes, fontsize=7)

# Save the plot to a PDF file
with PdfPages("Signature_genes_MDACC_KM_plot.pdf") as pdf:
    pdf.savefig(fig)

dpi = 300
# Save the plot as an SVG file with A4 dimensions
fig.savefig('Signature_genes_MDACC_plot.svg', format='svg', dpi=dpi, bbox_inches='tight')

plt.show()


# %%
# Fit a Cox Proportional Hazards Model
cph = CoxPHFitter(penalizer=0.1)
cph.fit(signature_df, duration_col='OS_Months', event_col='Vital_status')

# Compute risk scores and assign risk group
signature_df['Risk_scores'] = cph.predict_partial_hazard(signature_df)
signature_df['Risk_group'] = np.where(signature_df['Risk_scores'] >= np.median(signature_df['Risk_scores']), 1, 0)

# Define time points
time_points = [0, 40, 80, 120, 160, 200]

# Create empty DataFrames to store counts
low_risk_counts = pd.DataFrame(index=time_points, columns=['Count', 'Median risk score', 'At risk'])
high_risk_counts = pd.DataFrame(index=time_points, columns=['Count', 'Median risk score', 'At risk'])

# Loop over time points
for t in time_points:
    df = signature_df[signature_df['OS_Months'] >= t]
    low_risk_df = df[df['Risk_group'] == 0].sort_values(by='Risk_scores')
    high_risk_df = df[df['Risk_group'] == 1].sort_values(by='Risk_scores')
    low_risk_counts.loc[t, 'Count'] = len(low_risk_df)
    high_risk_counts.loc[t, 'Count'] = len(high_risk_df)
    low_risk_counts.loc[t, 'Median risk score'] = low_risk_df['Risk_scores'].median()
    high_risk_counts.loc[t, 'Median risk score'] = high_risk_df['Risk_scores'].median()
    low_risk_counts.loc[t, 'At risk'] = len(low_risk_df)
    high_risk_counts.loc[t, 'At risk'] = len(high_risk_df)

# Create the plot
fig, ax = plt.subplots(figsize=(3, 0.5))

# Plot the number of patients at risk for each time point
for index, row in low_risk_counts.iterrows():
    t = index
    ax.text(t, 0.2, f"{int(row['Count'])}", color='blue', fontsize=7, ha='center', va='bottom')

for index, row in high_risk_counts.iterrows():
    t = index
    ax.text(t, 1.2, f"{int(row['Count'])}", color='red', fontsize=7, ha='center', va='bottom')

# Customize plot
plt.xlabel('Time (months)', fontsize=7)
plt.ylabel('Risk group', fontsize=7)

# Set x-axis limits and ticks
plt.xlim(-10, 65)
plt.xticks(time_points, fontsize=7)

# Set y-axis limits, ticks, and tick labels
plt.ylim(-0.3, 1.7)
ax.set_yticks([0.2, 1.2])
ax.set_yticklabels(['Low risk', 'High risk'], color='black', fontsize=7)
ax.get_yticklabels()[0].set_color('blue')
ax.get_yticklabels()[1].set_color('red')

# Remove gridlines and frame
ax.grid(False)
ax.spines[['top', 'right']].set_visible(False)
ax.spines['left'].set_linewidth(0.5)  # Decrease line thickness of the y-axis
ax.spines['bottom'].set_linewidth(0.5)  # Decrease line thickness of the x-axis

# Scale up the plot dimensions to match A4 paper size
scale = 1.0
dpi = 300
fig.set_size_inches(fig.get_size_inches() * scale)

# Save the plot to a PDF file with A4 dimensions
with PdfPages("Signature_genes_MDACC_Risk_numbers_plot.pdf", keep_empty=False) as pdf:
    pdf.savefig(fig, dpi=dpi, bbox_inches='tight')

# Save the plot as an SVG file with A4 dimensions
fig.savefig('Signature_genes_MDACC_Risk_numbers_plot.svg', format='svg', dpi=dpi, bbox_inches='tight')

# Show the plot
plt.show()


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


# %%
kmf = KaplanMeierFitter()


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
with PdfPages("Signature_genes_MDACC_KM_plot.pdf") as pdf:
    pdf.savefig(fig)

dpi = 300
# Save the plot as an SVG file with A4 dimensions
fig.savefig('Signature_genes_MDACC_plot.svg', format='svg', dpi=dpi, bbox_inches='tight')

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
with PdfPages("Signature_genes_MDACC_Risk_numbers_plot.pdf", keep_empty=False) as pdf:
    pdf.savefig(fig, dpi=dpi, bbox_inches='tight')

dpi = 300
# Save the plot as an SVG file with A4 dimensions
fig.savefig('Signature_genes_MDACC_Risk_numbers_plot.svg', format='svg', dpi=dpi, bbox_inches='tight')

## Show the plot
plt.show()


# %%
# set a time point
time_point = 60

# create binary labels
y_true = (signature_df['OS_Months'] <= time_point) & (signature_df['Vital_status'] == 1)

# calculate the roc curve
fpr, tpr, thresholds = roc_curve(y_true, signature_df['Risk_scores'])

# calculate the auc
roc_auc = auc(fpr, tpr)

fig, ax = plt.subplots(figsize=(5, 3))

# plot the roc curve
ax.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
ax.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
ax.set_xlim([0.0, 1.0])
ax.set_ylim([0.0, 1.05])
ax.set_xlabel('False Positive Rate', fontsize=7)
ax.set_ylabel('True Positive Rate', fontsize=7)
ax.set_title('Receiver operating characteristic example', fontsize=7)
ax.legend(loc="lower right", fontsize=7)

# Save the plot to a PDF file
with PdfPages("ROC_plot.pdf") as pdf:
    pdf.savefig(fig)

dpi = 300
# Save the plot as an SVG file with A4 dimensions
fig.savefig('ROC_plot.svg', format='svg', dpi=dpi, bbox_inches='tight')

plt.show()
