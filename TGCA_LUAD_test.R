library(TCGAbiolinks)
library(SummarizedExperiment)

query_luad <- GDCquery(
  project       = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query_luad, method = "api", files.per.chunk = 10)
luad_data <- GDCprepare(query_luad)
