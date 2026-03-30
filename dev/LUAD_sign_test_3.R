# ── 7. Risk skoru hesapla ────────────────────────────────────
scored <- score_signature(
  signature   = sarang_sig,
  expr_matrix = expr_mat
)

cat("Risk skoru özeti:\n")
print(summary(scored$risk_scores))

# ── 8. TNM ile karşılaştır ───────────────────────────────────
survival_vec <- Surv(os_time, os_status)
names(survival_vec) <- colnames(expr_mat)

stage_matched <- stage_num[names(scored$risk_scores)]

comparison <- compare_to_staging(
  scored_signature = scored,
  staging          = stage_matched,
  survival         = survival_vec[names(scored$risk_scores)]
)

cat("\n=== SONUÇLAR ===\n")
print(comparison)

# ── 9. KM Plot ───────────────────────────────────────────────
plot_km(
  scored_signature = scored,
  survival         = survival_vec[names(scored$risk_scores)],
  title            = "Sarang_2025 (AGER, MGP, PECAM1, SLC2A1) - TCGA-LUAD"
)
