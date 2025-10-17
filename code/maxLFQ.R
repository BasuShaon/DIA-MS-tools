# Set CRAN mirror
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# install only if missing
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("iq", quietly = TRUE)) install.packages("iq")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) BiocManager::install("AnnotationDbi")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")

library(iq)
library(data.table)
library(AnnotationDbi)
library(org.Hs.eg.db)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 4) {
  setwd(args[1])
  file_path <- args[2]
  out_path  <- args[3]
  convert   <- as.logical(args[4])
} else {
  stop("Not enough arguments provided")
}

pasef <- fread(file_path, header = TRUE, sep = ",", dec = ".")

# Correct column names expected by iq::fast_MaxLFQ
colnames(pasef) <- c("peptide", "sample", "protein", "intensity")

# Run MaxLFQ
maxlfq <- iq::fast_MaxLFQ(pasef)
pasef_lfq <- as.data.frame(t(maxlfq$estimate))

# Optional UniProt → Gene symbol mapping
if (convert) {
  uniprot_ids <- colnames(pasef_lfq)
  gene_symbols <- mapIds(
    org.Hs.eg.db,
    keys     = uniprot_ids,
    column   = "SYMBOL",
    keytype  = "UNIPROT",
    multiVals = function(x) paste(x, collapse = ";")
  )
  colnames(pasef_lfq) <- ifelse(is.na(gene_symbols), colnames(pasef_lfq), gene_symbols)
}

# Write output (CSV)
write.csv2(pasef_lfq, out_path)
