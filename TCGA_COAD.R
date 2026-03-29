# TCGA-COAD indir - daha küçük kohort (~500 sample)
query_coad <- GDCquery(
  project       = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type   = "Primary Tumor"
)

cat("Dosya sayısı:", nrow(getResults(query_coad)), "\n")
GDCdownload(query_coad, directory = "C:/Users/berna/GDCdata")
