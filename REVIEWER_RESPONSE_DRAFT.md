# Draft response to reviewer

**Reviewer comment:** Ensure that all code and hyperparameters used to generate the models are made available in a public repository (e.g., GitHub) to ensure reproducibility, which is a standard requirement for computational biology papers today.

**Response:** We thank the reviewer for emphasizing computational reproducibility. We have organized and deposited the analysis code in a public GitHub repository. The repository includes the R code used for data preparation, differential-expression analysis, LASSO feature selection, and enrichment analysis, together with the Python notebooks used for feature-subset selection, XGBoost model training, class-weighted risk stratification, performance evaluation, bootstrap confidence intervals, individual-miRNA ROC analysis, and SHAP interpretation. The repository also documents the model hyperparameters, class definitions, random seeds, expected input files, output files, and execution order. The repository is available at: **[INSERT GITHUB URL]**.

Do not submit this response until the missing items listed in `CODE_AUDIT.md` have been added and the repository has been rerun successfully.
