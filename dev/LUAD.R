library(TCGAbiolinks)
library(SummarizedExperiment)
library(org.Hs.eg.db)
library(AnnotationDbi)
devtools::load_all("C:/Users/berna/CancerSignComp")

# LUAD verisini yeniden hazırla
query_luad <- GDCquery(
  project       = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type   = "Primary Tumor"
)
luad_data <- GDCprepare(query_luad, directory = "C:/Users/berna/GDCdata")

# TPM matrix
expr_luad <- assay(luad_data, "tpm_unstrand")
ens_ids <- gsub("\\..*", "", rownames(expr_luad))
symbols <- mapIds(org.Hs.eg.db, keys = ens_ids, column = "SYMBOL",
                  keytype = "ENSEMBL", multiVals = "first")
valid <- !is.na(symbols) & !duplicated(symbols)
expr_sym_luad <- expr_luad[valid, ]
rownames(expr_sym_luad) <- symbols[valid]

# Klinik veri
clin_luad <- as.data.frame(colData(luad_data))
surv_time_luad  <- ifelse(!is.na(clin_luad$days_to_death),
                          clin_luad$days_to_death,
                          clin_luad$days_to_last_follow_up)
surv_event_luad <- ifelse(clin_luad$vital_status == "Dead", 1, 0)
names(surv_time_luad)  <- rownames(clin_luad)
names(surv_event_luad) <- rownames(clin_luad)
stage_named <- clin_luad$ajcc_pathologic_stage
names(stage_named) <- rownames(clin_luad)

# Filtrele
valid_luad  <- !is.na(surv_time_luad) & surv_time_luad > 0
common_luad <- intersect(colnames(expr_sym_luad), names(surv_time_luad)[valid_luad])
expr_sym_luad   <- expr_sym_luad[, common_luad]
surv_time_luad  <- surv_time_luad[common_luad]
surv_event_luad <- surv_event_luad[common_luad]
stage_named     <- stage_named[common_luad]

cat("LUAD hazır: n =", length(common_luad), "\n")
