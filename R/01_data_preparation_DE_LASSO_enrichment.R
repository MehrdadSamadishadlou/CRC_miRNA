# Colorectal cancer serum miRNA analysis
# Run this script from the repository root.

library(GEOquery)
library(limma)
library(dplyr)
library(glmnet)
library(clusterProfiler)
library(org.Hs.eg.db)

# Create output directories.
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/intermediate", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/lasso/all_crc", recursive = TRUE, showWarnings = FALSE)
dir.create("results/lasso/stage_0", recursive = TRUE, showWarnings = FALSE)
dir.create("results/lasso/stage_4", recursive = TRUE, showWarnings = FALSE)


#### Read sample information from GEO ####

geo_series <- getGEO(
  "GSE211692",
  GSEMatrix = TRUE,
  AnnotGPL = TRUE
)
geo_dataset <- geo_series[[1]]

geo_accession <- geo_dataset@phenoData@data[["geo_accession"]]
disease_status <- geo_dataset@phenoData@data[["disease state:ch1"]]
sample_id <- geo_dataset@phenoData@data[["title"]]

crc_gsm <- geo_accession[disease_status == "colorectal cancer"]
crc_sample_id <- sample_id[disease_status == "colorectal cancer"]
crc_age <- geo_dataset@phenoData@data[["age:ch1"]][disease_status == "colorectal cancer"]
crc_stage <- geo_dataset@phenoData@data[["Stage:ch1"]][disease_status == "colorectal cancer"]
crc_sex <- geo_dataset@phenoData@data[["Sex:ch1"]][disease_status == "colorectal cancer"]

crc_metadata <- data.frame(
  row.names = crc_sample_id,
  gsm = crc_gsm,
  age = crc_age,
  sex = crc_sex,
  stage = crc_stage
)

crc_metadata$age <- as.integer(crc_metadata$age)
crc_metadata$sex <- as.factor(crc_metadata$sex)
crc_metadata$stage <- as.factor(crc_metadata$stage)

write.csv(
  crc_metadata,
  "data/processed/CRC_all_data.csv",
  quote = FALSE
)

crc_metadata_complete <- crc_metadata[!is.na(crc_metadata$stage), ]

# Read the saved file in the same format used in the original analysis.
crc_metadata <- read.csv(
  "data/processed/CRC_all_data.csv",
  check.names = FALSE
)


#### Combine CRC and non-cancer expression matrices ####

healthy_expression <- read.table(
  "data/raw/no_cancer_expression.txt",
  check.names = FALSE
)

crc_expression <- read.table(
  "data/raw/CR_expression.txt",
  check.names = FALSE
)

crc_expression <- crc_expression[
  , colnames(crc_expression) %in% rownames(crc_metadata_complete)
]

all_expression <- cbind(healthy_expression, crc_expression)

write.csv(
  all_expression,
  "data/processed/Expression_all.csv",
  quote = FALSE
)

# train_set.csv and test_set.csv were created in Python in the original workflow.
# The code used to create these files was not included in the supplied files.


#### Quantile normalization of the training data ####

training_data <- read.csv(
  "data/processed/train_set.csv",
  check.names = FALSE
)

rownames(training_data) <- training_data[, 1]
training_data <- training_data[, -1]
training_status <- training_data[, 2566]
training_data <- training_data[, -2566]
training_data <- t(training_data)

training_data_normalized <- normalizeQuantiles(training_data)
training_data_normalized <- rbind(
  training_data_normalized,
  training_status
)

write.table(
  training_data_normalized,
  "data/processed/train_normal.txt",
  quote = FALSE
)


#### Differential expression: CRC versus non-cancer ####

training_data_normalized <- read.delim(
  "data/processed/train_normal.txt",
  sep = " ",
  check.names = FALSE
)

training_status <- as.character(training_data_normalized[2566, ])
training_expression <- training_data_normalized[-2566, ]
training_expression <- as.data.frame(
  lapply(training_expression, type.convert, as.is = TRUE)
)
rownames(training_expression) <- rownames(training_data_normalized)[1:2565]

status_factor <- factor(training_status)
design_matrix <- model.matrix(
  ~ status_factor + 0,
  as.data.frame(training_expression)
)
colnames(design_matrix) <- levels(status_factor)

# Check the number of samples in each group.
colSums(design_matrix)

limma_fit <- lmFit(training_expression, design_matrix)
contrast_matrix <- makeContrasts(CRC - Healthy, levels = design_matrix)
limma_fit <- contrasts.fit(limma_fit, contrast_matrix)
limma_fit <- eBayes(limma_fit, 0.01)
crc_vs_healthy_results <- topTable(
  limma_fit,
  adjust = "fdr",
  sort.by = "B",
  number = Inf
)

crc_vs_healthy_results <- subset(
  crc_vs_healthy_results,
  select = c("adj.P.Val", "logFC")
)

write.csv(
  crc_vs_healthy_results,
  "results/tables/train_CRC_Healthy_statistics.csv",
  row.names = TRUE,
  quote = FALSE
)

crc_upregulated <- subset(
  crc_vs_healthy_results,
  logFC > 1 & adj.P.Val < 0.05
)

healthy_upregulated <- subset(
  crc_vs_healthy_results,
  logFC < -1 & adj.P.Val < 0.05
)

differential_mirnas <- rbind(
  crc_upregulated,
  healthy_upregulated
)

write.csv(
  crc_upregulated,
  "results/tables/train_CRC_upregulated.csv",
  quote = FALSE,
  row.names = TRUE
)

write.csv(
  healthy_upregulated,
  "results/tables/train_Healthy_upregulated.csv",
  quote = FALSE,
  row.names = TRUE
)

write.csv(
  differential_mirnas,
  "results/tables/train_differential_miRNAs.csv",
  quote = FALSE,
  row.names = TRUE
)

training_differential_expression <- training_expression[
  rownames(differential_mirnas),
]

write.table(
  training_differential_expression,
  "results/intermediate/expression_differential_miRNAs_train.txt",
  quote = FALSE
)


#### Prepare stage-specific training matrices ####

crc_metadata <- read.csv(
  "data/processed/CRC_all_data.csv",
  check.names = FALSE
)

crc_metadata_complete <- crc_metadata[!is.na(crc_metadata$stage), ]
rownames(crc_metadata_complete) <- crc_metadata_complete[, 1]
crc_metadata_complete <- crc_metadata_complete[, -1]

training_crc_metadata <- crc_metadata_complete %>%
  filter(rownames(crc_metadata_complete) %in% colnames(training_expression))

stage_0_samples <- rownames(
  training_crc_metadata[
    training_crc_metadata$stage == 0 & !is.na(training_crc_metadata$stage),
  ]
)
stage_0_expression <- training_expression[, stage_0_samples]

stage_1_samples <- rownames(
  training_crc_metadata[
    training_crc_metadata$stage == 1 & !is.na(training_crc_metadata$stage),
  ]
)
stage_1_expression <- training_expression[, stage_1_samples]

stage_2_samples <- rownames(
  training_crc_metadata[
    training_crc_metadata$stage == 2 & !is.na(training_crc_metadata$stage),
  ]
)
stage_2_expression <- training_expression[, stage_2_samples]

stage_3_samples <- rownames(
  training_crc_metadata[
    training_crc_metadata$stage == 3 & !is.na(training_crc_metadata$stage),
  ]
)
stage_3_expression <- training_expression[, stage_3_samples]

stage_4_samples <- rownames(
  training_crc_metadata[
    training_crc_metadata$stage == 4 & !is.na(training_crc_metadata$stage),
  ]
)
stage_4_expression <- training_expression[, stage_4_samples]


#### Differential expression: stage 0 versus non-cancer ####

healthy_training_samples <- colnames(training_expression)[
  !(colnames(training_expression) %in% crc_metadata[, 1])
]
healthy_training_expression <- training_expression[, healthy_training_samples]

stage_comparison_expression <- cbind(
  stage_0_expression,
  healthy_training_expression
)

stage_comparison_status <- c(
  rep("Stage", dim(stage_0_expression)[2]),
  rep("Healthy", 4232)
)

status_factor <- factor(stage_comparison_status)
design_matrix <- model.matrix(
  ~ status_factor + 0,
  as.data.frame(stage_comparison_expression)
)
colnames(design_matrix) <- levels(status_factor)

# Check the number of samples in each group.
colSums(design_matrix)

limma_fit <- lmFit(stage_comparison_expression, design_matrix)
contrast_matrix <- makeContrasts(Stage - Healthy, levels = design_matrix)
limma_fit <- contrasts.fit(limma_fit, contrast_matrix)
limma_fit <- eBayes(limma_fit, 0.01)
stage_results <- topTable(
  limma_fit,
  adjust = "fdr",
  sort.by = "B",
  number = Inf
)

stage_results <- subset(
  stage_results,
  select = c("adj.P.Val", "logFC")
)

write.csv(
  stage_results,
  "results/lasso/stage_0/stage_statistics.csv",
  row.names = TRUE,
  quote = FALSE
)

stage_upregulated <- subset(
  stage_results,
  logFC > 1 & adj.P.Val < 0.05
)

stage_downregulated <- subset(
  stage_results,
  logFC < -1 & adj.P.Val < 0.05
)

stage_differential_mirnas <- rbind(
  stage_upregulated,
  stage_downregulated
)

write.csv(
  stage_upregulated,
  "results/lasso/stage_0/stage_upregulated.csv",
  quote = FALSE,
  row.names = TRUE
)

write.csv(
  stage_downregulated,
  "results/lasso/stage_0/stage_downregulated.csv",
  quote = FALSE,
  row.names = TRUE
)

write.csv(
  stage_differential_mirnas,
  "results/lasso/stage_0/stage_differential_miRNAs.csv",
  quote = FALSE,
  row.names = TRUE
)

stage_differential_expression <- stage_comparison_expression[
  rownames(stage_differential_mirnas),
]

write.table(
  stage_differential_expression,
  "results/lasso/stage_0/expression_differential_miRNAs.txt",
  quote = FALSE
)


#### LASSO: stage 4 versus non-cancer ####

stage_comparison_status <- c(
  rep("Stage", 90),
  rep("Healthy", 4232)
)

disease_group <- data.frame(disease = stage_comparison_status)
lasso_outcome <- model.matrix(~ disease - 1, data = disease_group)

stage_differential_expression <- read.table(
  "results/lasso/stage_4/expression_differential_miRNAs.txt",
  check.names = FALSE
)

lasso_predictors <- t(stage_differential_expression)

lasso_cv <- cv.glmnet(
  lasso_predictors,
  lasso_outcome,
  alpha = 1,
  family = "binomial",
  type.measure = "class",
  nfolds = 10
)

# Extract non-zero coefficients at lambda.1se.
coefficient_matrix <- as.matrix(
  coef(lasso_cv, s = "lambda.1se")
)

selected_mirnas <- rownames(coefficient_matrix)[
  coefficient_matrix[, 1] != 0
]
selected_mirnas <- selected_mirnas[selected_mirnas != "(Intercept)"]

print(selected_mirnas)

write.table(
  selected_mirnas,
  "results/lasso/stage_4/lasso_miRNAs.txt",
  quote = FALSE,
  row.names = FALSE
)

selected_mirna_expression <- stage_differential_expression[
  selected_mirnas,
]
selected_mirna_expression <- rbind(
  selected_mirna_expression,
  stage_comparison_status
)

write.table(
  selected_mirna_expression,
  "results/lasso/stage_4/expression_lasso_miRNAs.txt",
  quote = FALSE
)


#### Functional enrichment using TarBase targets ####

# This section produces enrichment result tables only. Plotting code was removed.
tarbase <- read.delim(
  "data/external/Homo_sapiens_TarBase-v9.tsv"
)

all_crc_lasso_mirnas <- read.table(
  "results/lasso/all_crc/lasso_miRNAs.txt",
  check.names = FALSE
)

tarbase_target_genes <- unique(
  tarbase %>%
    filter(mirna_name %in% all_crc_lasso_mirnas$V1) %>%
    pull(gene_name)
)

entrez_ids <- bitr(
  tarbase_target_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)$ENTREZID

kegg_results <- enrichKEGG(
  entrez_ids,
  pvalueCutoff = 0.05,
  organism = "hsa"
)
kegg_results <- kegg_results@result
kegg_results <- kegg_results %>% filter(p.adjust < 0.05)
write.table(
  kegg_results,
  "results/tables/KEGG_enrichment.tsv",
  quote = FALSE,
  sep = "\t"
)

cellular_component_results <- enrichGO(
  entrez_ids,
  org.Hs.eg.db,
  ont = "CC",
  pvalueCutoff = 0.05
)
cellular_component_results <- cellular_component_results@result
cellular_component_results <- cellular_component_results %>%
  filter(p.adjust < 0.05)
write.table(
  cellular_component_results,
  "results/tables/GO_cellular_component.tsv",
  quote = FALSE,
  sep = "\t"
)

biological_process_results <- enrichGO(
  entrez_ids,
  org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.05
)
biological_process_results <- biological_process_results@result
biological_process_results <- biological_process_results %>%
  filter(p.adjust < 0.05)
write.table(
  biological_process_results,
  "results/tables/GO_biological_process.tsv",
  quote = FALSE,
  sep = "\t"
)

molecular_function_results <- enrichGO(
  entrez_ids,
  org.Hs.eg.db,
  ont = "MF",
  pvalueCutoff = 0.05
)
molecular_function_results <- molecular_function_results@result
molecular_function_results <- molecular_function_results %>%
  filter(p.adjust < 0.05)
write.table(
  molecular_function_results,
  "results/tables/GO_molecular_function.tsv",
  quote = FALSE,
  sep = "\t"
)


#### Extract selected miRNAs for machine learning ####

training_set <- read.csv(
  "data/processed/train_set.csv",
  check.names = FALSE,
  row.names = 1
)

test_set <- read.csv(
  "data/processed/test_set.csv",
  check.names = FALSE,
  row.names = 1
)

candidate_mirnas <- c(
  "hsa-miR-663a", "hsa-miR-1233-5p", "hsa-miR-4783-3p",
  "hsa-miR-4730", "hsa-miR-125a-3p", "hsa-miR-1306-5p",
  "hsa-miR-1307-3p", "hsa-miR-134-3p", "hsa-miR-1587",
  "hsa-miR-3194-5p", "hsa-miR-4530", "hsa-miR-4701-5p",
  "hsa-miR-642b-3p", "hsa-miR-6792-3p", "hsa-miR-6800-5p",
  "hsa-miR-6860", "hsa-miR-6880-5p", "hsa-miR-6893-5p",
  "hsa-miR-8073", "hsa-miR-1273c", "hsa-miR-4723-5p",
  "hsa-miR-575", "hsa-miR-6772-5p", "hsa-miR-6784-5p"
)

selected_training_expression <- training_set[, candidate_mirnas]
selected_test_expression <- test_set[, candidate_mirnas]

stage_lookup <- setNames(crc_metadata$stage, crc_metadata[, 1])

selected_training_expression$stage <- stage_lookup[
  rownames(selected_training_expression)
]
selected_test_expression$stage <- stage_lookup[
  rownames(selected_test_expression)
]

selected_training_expression$stage[
  is.na(match(rownames(selected_training_expression), crc_metadata[, 1]))
] <- "Healthy"

selected_test_expression$stage[
  is.na(match(rownames(selected_test_expression), crc_metadata[, 1]))
] <- "Healthy"

write.table(
  selected_training_expression,
  "data/processed/train_exp_selected_mirs.txt",
  quote = FALSE
)

write.table(
  selected_test_expression,
  "data/processed/test_exp_selected_mirs.txt",
  quote = FALSE
)
