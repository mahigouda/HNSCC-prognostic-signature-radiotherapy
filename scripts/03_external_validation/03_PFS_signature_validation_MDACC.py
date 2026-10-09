# =============================================================================
# Project : Prognostic gene-expression signature for radiotherapy outcome in HPV-negative HNSCC
# Module  : 03 | External validation of the signatures (MDACC, FHCRC cohorts)
# Script  : 03_PFS_signature_validation_MDACC.py
# Purpose : Evaluation of the 10 PFS signature genes in the independent MDACC cohort
#           (Cox model refitted in the cohort, no transfer of TCGA coefficients): risk
#           scores, median split, Kaplan-Meier/log-rank, number-at-risk table and ROC/AUC.
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

clinical_df_mdacc = pd.read_csv(data_dir / 'clinical_mdacc_OSCC.txt', sep='\t')

clinical_df_mdacc.set_index(clinical_df_mdacc.columns[0], inplace=True)

merged_df = exp_df_mdacc.merge(clinical_df_mdacc, how='inner', left_index=True, right_index=True)


# convert weeks to months
merged_df['OS_Months'] = merged_df['OS_Months'] / 12

## Remove the decimal part and keep the integer part of the "OS_Months" column
merged_df['OS_Months'] = merged_df['OS_Months'].apply(lambda x: int(x) if not np.isnan(x) else x)


# %%
# Remove rows with NaNs
merged_df = merged_df.dropna()


# %%
df = merged_df


# %%
# Convert 'Vital_status' column to numeric representation
df['Vital_status'] = df['Vital_status'].map({'Alive': 0, 'Dead_OSCC': 1})


# %%
signature_genes = ['ADAM8', 'RAC2', 'CSF2', 'PPP1R18', 'HBEGF', 'ITGA5', 'LRRC59', 'PLAU', 'FOSL1', 'PLAUR']

additional_columns = ['Age', 'OS_Months', 'Vital_status']

signature_df = df[signature_genes + additional_columns]


# %%
## Remove index column
signature_df = signature_df.reset_index(drop=True)


# %%
# Remove rows with NaNs
signature_df = signature_df.dropna()


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
ax.text(0.3, 0.05, f'logrank p high vs low={p_value_high_vs_low:.3f}', transform=ax.transAxes, fontsize=7)
ax.text(0.3, 0.1, f'HR (high-risk): {high_vs_low_hr:.2f}', transform=ax.transAxes, fontsize=7)

# Save the plot to a PDF file
with PdfPages("PFS_Signature_genes_MDACC_KM_plot.pdf") as pdf:
    pdf.savefig(fig)

dpi = 300
# Save the plot as an SVG file with A4 dimensions
fig.savefig('PFS_Signature_genes_MDACC_plot.svg', format='svg', dpi=dpi, bbox_inches='tight')

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
with PdfPages("PFS_Signature_genes_MDACC_Risk_numbers_plot.pdf", keep_empty=False) as pdf:
    pdf.savefig(fig, dpi=dpi, bbox_inches='tight')

dpi = 300
# Save the plot as an SVG file with A4 dimensions
fig.savefig('PFS_Signature_genes_MDACC_Risk_numbers_plot.svg', format='svg', dpi=dpi, bbox_inches='tight')

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
with PdfPages("PFS_signature_genes_ROC_plot.pdf") as pdf:
    pdf.savefig(fig)

dpi = 300
# Save the plot as an SVG file with A4 dimensions
fig.savefig('PFS_Signature_genes_ROC_plot.svg', format='svg', dpi=dpi, bbox_inches='tight')

plt.show()
