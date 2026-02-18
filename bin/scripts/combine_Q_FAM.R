#!/usr/bin/env Rscript
## Usage: Rscript combine_Q_FAM.R file.fam file.Q [output.csv]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Need at least: <fam> <Q> [out]")
}
fam  <- args[1]
qmat <- args[2]
out  <- ifelse(length(args) >= 3, args[3], "combined_ancestry.csv")

# ----- read files -----------------------------------------------------------
fam_df <- read.table(fam, header = FALSE,
                     col.names = c("FID","IID","PID","MID","Sex","Pheno"),
                     stringsAsFactors = FALSE)

q_df   <- read.table(qmat, header = FALSE, stringsAsFactors = FALSE)
colnames(q_df) <- paste0("K", seq_len(ncol(q_df)))   # K1, K2, …

if (nrow(fam_df) != nrow(q_df)) {
  stop("Row mismatch: fam rows = ", nrow(fam_df),
       ", Q rows = ", nrow(q_df),
       ". Ensure they’re in the same individual order.")
}

# ----- combine & write ------------------------------------------------------
combined <- cbind(fam_df[c("FID","IID")], q_df)
write.csv(combined, file = out, row.names = FALSE, quote = FALSE)

cat("✅ Wrote", out, "with", nrow(combined), "rows and",
    ncol(combined), "columns\n")
