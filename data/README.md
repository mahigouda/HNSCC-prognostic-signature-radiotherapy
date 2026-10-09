# Data

No data are distributed with this repository. All datasets are public. The scripts expect the
layout below, relative to the repository root.

```
data/
├── tcga/
│   ├── expData.txt                                   # TCGA-HNSC gene expression (Ensembl_ID x samples)
│   ├── clinical.txt                                  # curated clinical table (OS, PFS, DSS, RT, HPV, ...)
│   ├── clinical_DSS.txt                              # as above, with DSS status
│   ├── clinical_xena_cbioportal.txt                  # extended clinical table (sex, stage, site, outcomes)
│   ├── merged_data_final_threshold_0_5.txt           # output of 01_data_preprocessing/01_*.R
│   ├── merged_data_final_threshold_0_5_noDSS.txt     # same gene set, earlier export without the DSS column
│   └── signature_genes_merged_df_for_DSS.txt         # signature genes + DSS endpoint
├── validation_cohorts/
│   ├── GSE42743&41613.Rdata, Mdacc_cohort.Rdata, Fhcrc_cohort.Rdata
│   ├── exp_mdacc.txt, exp_fhcrc.txt                  # output of 01_data_preprocessing/03_*.R
│   └── clinical_mdacc_OSCC.txt, clinical_fhcrc_OSCC.txt
├── clinical_covariates/
│   ├── TCGA_OS_signature_clinical_df.txt, TCGA_PFS_signature_clinical_df.txt
│   ├── merged_df_MDACC_for_analysis.txt, merged_df_MDACC_PFS_for_analysis.txt
│   └── merged_df_FHCRC_for_analysis.txt, merged_df_FHCRC_PFS_for_analysis.txt
└── single_cell/
    └── data_malignant.txt                            # GSE103322, malignant cells
```

## Sources

| Dataset | Description | Source |
|---|---|---|
| TCGA-HNSC | Primary tumour gene expression | GDC / UCSC Xena (TCGA-HNSC) |
| TCGA-HNSC clinical | Survival endpoints, HPV status, radiotherapy, stage and other covariates | cBioPortal, *Head and Neck Squamous Cell Carcinoma (TCGA, PanCancer Atlas)* (`hnsc_tcga_pan_can_atlas_2018`); UCSC Xena |
| MDACC cohort | Oral cavity SCC, microarray | GEO [GSE42743](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42743) |
| FHCRC cohort | Oral cavity SCC, microarray | GEO [GSE41613](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41613) |
| HNSCC single-cell | scRNA-seq of primary HNSCC (Puram *et al.* 2017) | GEO [GSE103322](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322) |

## Preparation steps not covered by the scripts

- **Clinical tables.** The clinical tables were exported from cBioPortal/UCSC Xena and reduced
  to the columns that the merging scripts rename. Column order matters because the scripts
  rename the last *n* columns.
- **Candidate genes.** The candidate gene lists in `01_data_preprocessing/01_*.R` are
  druggable genes that were differentially expressed in the study's cell-model RNA-seq data,
  at log2 fold-change thresholds of 0.5 and 1 (see the article).
- **Analysis tables.** The tables in `clinical_covariates/` and
  `signature_genes_merged_df_for_DSS.txt` combine the signature genes with the clinical
  variables. They were prepared from the merged tables (module 01) and the cohort tables.
- **Outputs used as inputs.** Outputs of module 01 that later modules use as inputs must be
  copied from `results/01_data_preprocessing/` to the matching `data/` folder.
