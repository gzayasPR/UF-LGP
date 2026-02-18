#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GWASTools)
  library(GENESIS)
  library(SNPRelate)
  library(data.table)
})

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript Selection_of_Unrelated_Script.R <NAME>\n",
       "  where <NAME> is the prefix for .bed/.bim/.fam files\n",
       "  Example: Rscript Selection_of_Unrelated_Script.R mydata")
}

pop.name <- args[1]

# Core file names
bed <- paste(pop.name, "bed", sep = ".")
bim <- paste(pop.name, "bim", sep = ".")
fam <- paste(pop.name, "fam", sep = ".")
gds <- paste(pop.name, "gds", sep = ".")

# Check that required files exist
if (!file.exists(bed)) stop("File not found: ", bed)
if (!file.exists(bim)) stop("File not found: ", bim)
if (!file.exists(fam)) stop("File not found: ", fam)

unrelated_list <- "Unrelated_list.txt"

# ------------------ Convert PLINK -> GDS ------------------
X_GDS <- snpgdsBED2GDS(bed.fn = bed, bim.fn = bim, fam.fn = fam, out.gdsfn = gds)

# ------------------ KING kinship ------------------
X_GDS4 <- snpgdsOpen(filename = gds)
X.king_kinship <- snpgdsIBDKING(
  X_GDS4,
  sample.id = NULL,
  snp.id = NULL,
  autosome.only = FALSE,
  type = "KING-robust",
  family.id = NULL,
  verbose = TRUE
)
saveRDS(X.king_kinship, file = "X.king_kinship.rds")

X.kinship_matrix <- X.king_kinship$kinship
colnames(X.kinship_matrix) <- X.king_kinship$sample.id
rownames(X.kinship_matrix) <- X.king_kinship$sample.id
snpgdsClose(X_GDS4)

# ------------------ GENESIS GenotypeData ------------------
X_geno <- GenotypeData(GdsGenotypeReader(filename = gds))
saveRDS(X_geno, file = "X_geno.rds")

# ------------------ PCAIR ------------------
X_pcair <- pcair(
  gdsobj = X_geno,
  kinobj = X.kinship_matrix,
  divobj = X.kinship_matrix
)
saveRDS(X_pcair, file = "X_pcair.rds")

# Extract PCs (up to PC20) and save to CSV with FID/IID
pcs_mat <- X_pcair$vectors
if (is.null(pcs_mat) || nrow(pcs_mat) == 0) stop("pcair returned no vectors.")
n_pc <- min(20L, ncol(pcs_mat))
pcs_dt <- as.data.table(pcs_mat[, 1:n_pc, drop = FALSE])
setnames(pcs_dt, paste0("PC", seq_len(n_pc)))

# Attach IIDs from rownames (or sample.id fallback)
rn <- rownames(pcs_mat)
if (!is.null(rn)) {
  pcs_dt[, IID := rn]
} else if (!is.null(X.king_kinship$sample.id)) {
  pcs_dt[, IID := as.character(X.king_kinship$sample.id)]
} else {
  stop("Cannot determine IID labels for PCs (no rownames and no sample.id).")
}

# Bring in FID from .fam
fam_dt <- fread(fam, header = FALSE)
# PLINK .fam: V1=FID, V2=IID
setnames(fam_dt, c("V1", "V2"), c("FID", "IID"))
pcs_dt <- merge(pcs_dt, fam_dt[, .(FID, IID)], by = "IID", all.x = TRUE)
setcolorder(pcs_dt, c("FID", "IID", paste0("PC", seq_len(n_pc))))

fwrite(pcs_dt, file = "PC_scores_PC1-20.csv")  # will contain up to PC20

# ------------------ PC-Relate using PC1 & PC2 ------------------
X.PC_matrix <- pcs_mat[, 1:2, drop = FALSE]
rownames(X.PC_matrix) <- X.king_kinship$sample.id

X_iterator <- GenotypeBlockIterator(X_geno, snpBlock = 2000)
X.kinship_adjusted <- pcrelate(
  X_iterator,
  pcs = X.PC_matrix,
  sample.block.size = 500
)
saveRDS(X.kinship_adjusted, file = "X.kinship_adjusted.rds")

# ------------------ Unrelated selection ------------------
X.unrelated_list <- pcairPartition(
  kinobj = X.kinship_matrix,
  divobj = X.kinship_matrix
)
saveRDS(X.unrelated_list, file = "X.unrelated_list.rds")

unrelsNO <- length(X.unrelated_list$unrels)
print(X.unrelated_list$unrels)
# Build Unrelated_list.txt with columns: IID FID (space-separated)
FAM.IDS <- fread(fam, header = FALSE)[, 1:2]
setnames(FAM.IDS, c("V1", "V2"), c("FID", "IID"))
x.unrels <- merge(
  FAM.IDS,
  data.table(IID = X.unrelated_list$unrels),
  by = "IID"
)[, .(FID, IID)]  # order: IID FID to match many downstream keep/merge conventions

# Print & write
print(x.unrels)
fwrite(x.unrels, file = "Unrelated_list.txt", sep = " ", quote = FALSE, col.names = FALSE)

message("Done. Wrote:")
message("  - X.king_kinship.rds")
message("  - X_geno.rds")
message("  - X_pcair.rds")
message("  - X.kinship_adjusted.rds")
message("  - X.unrelated_list.rds")
message("  - Unrelated_list.txt (IID FID)")
message("  - PC_scores_PC1-20.csv (FID, IID, PC1..PCk)")