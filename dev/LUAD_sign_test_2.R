# ============================================================
# Sarang_2025 İmzası Testi - CancerSignComp
# AGER, MGP, PECAM1, SLC2A1 (diagnostik → prognostik test)
# ============================================================

library(CancerSignComp)
library(SummarizedExperiment)
library(survival)

# ── 1. TCGA-LUAD verisini yükle ──────────────────────────────
load("C:/Users/berna/GDCdata/TCGA_LUAD_prepared.RData")
# luad_data değişken adın neyse onu kullan

# Ekspresyon matrisini al
expr_mat <- assay(luad_data, "fpkm_unstrand")  # veya "unstranded"

# ── 2. Klinik veriyi hazırla ─────────────────────────────────
clinical <- as.data.frame(colData(luad_data))

# Survival vektörü
os_time   <- as.numeric(clinical$days_to_death)
os_status <- as.numeric(clinical$vital_status == "Dead")

# Eksik ölüm tarihi → follow-up süresi ile doldur
na_idx <- is.na(os_time)
os_time[na_idx] <- as.numeric(clinical$days_to_last_follow_up[na_idx])

# Geçersiz satırları çıkar
valid <- !is.na(os_time) & !is.na(os_status) & os_time > 0
os_time   <- os_time[valid]
os_status <- os_status[valid]
expr_mat  <- expr_mat[, valid]
clinical  <- clinical[valid, ]

# İsimleri eşleştir
names(os_time)   <- colnames(expr_mat)
names(os_status) <- colnames(expr_mat)

# TNM stage
stage_raw <- clinical$ajcc_pathologic_stage
stage_num <- dplyr::case_when(
  grepl("Stage I$|Stage IA|Stage IB", stage_raw)   ~ 1,
  grepl("Stage II$|Stage IIA|Stage IIB", stage_raw) ~ 2,
  grepl("Stage III", stage_raw)                      ~ 3,
  grepl("Stage IV", stage_raw)                       ~ 4,
  TRUE ~ NA_real_
)
names(stage_num) <- colnames(expr_mat)

cat("Geçerli hasta sayısı:", sum(valid), "\n")
cat("Stage dağılımı:\n")
print(table(stage_num, useNA = "ifany"))

# ── 3. Sarang_2025 imzasını tanımla ─────────────────────────
# NOT: SLC2A1 upregulated → pozitif ağırlık
#      AGER, MGP, PECAM1 downregulated → negatif ağırlık (yüksek = iyi)
# Risk skoru için downregulated genleri ters işaretle

sarang_sig <- load_signature(
  genes   = c("AGER", "MGP", "PECAM1", "SLC2A1"),
  name    = "Sarang_2025",
  source  = "Biomolecules",
  cancer  = "LUAD",
  weights = c(AGER = -1, MGP = -1, PECAM1 = -1, SLC2A1 = 1)
)

print(sarang_sig)

# ── 4. Hangi genler matriste mevcut? ────────────────────────
target_genes <- c("AGER", "MGP", "PECAM1", "SLC2A1")
found <- target_genes[target_genes %in% rownames(expr_mat)]
missing <- target_genes[!target_genes %in% rownames(expr_mat)]

cat("\nBulunan genler:", paste(found, collapse = ", "), "\n")
cat("Eksik genler :", paste(missing, collapse = ", "), "\n")
