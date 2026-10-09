# Code curation notes

The scripts in `scripts/` were curated from the working analysis scripts and Jupyter notebooks
used for the article. The aim was to make them readable and portable **without changing the
analysis logic, parameters, gene lists or thresholds**.

## What was changed

1. **One version per analysis.** Where several versions of an analysis existed, only the
   final version was kept. These were the scripts shared with collaborators for the revision,
   and the latest notebook for each model. Byte-identical copies, earlier drafts and
   exploratory scratch scripts were not carried over.
2. **Standard header.** Each script starts with a header that lists its purpose, input and
   output.
3. **Portable paths.** Hard-coded local paths, `setwd()` and `os.chdir()` calls were
   replaced by a `Paths` block. Inputs are read from `data/<dataset>/` and outputs are written
   to `results/<module>/`.
4. **Neutral file names.** Dates and personal names were removed from file names. For
   example, a dated revision table became `TCGA_OS_signature_clinical_df.txt`.
5. **Duplicates.** The following were removed:
   - Repeated `library()` calls and duplicated Python imports (both no-ops).
   - A copied block at the end of the extended-clinical merging script that repeated the
     threshold-1 merge.
   - Intermediate versions of single-cell heatmaps that were immediately redrawn. Only the
     final version is kept.
6. **Non-functional drafts.** The following were removed:
   - A forward-selection draft (`regsubsets`) that referenced an undefined column.
   - A placeholder plot that used simulated p-values.
   - A notebook cell that was never executed.
7. **Notebooks.** The Jupyter notebooks were converted to `.py` scripts in jupytext "percent"
   format. Cell order and code are unchanged. Header cells, empty cells and display-only cells
   (a bare variable name) were removed.
8. **Comments.** Commented-out code was removed. Explanatory comments were kept.
