#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(reshape2)
  library(ggpubr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat("Usage:\n",
      "  Rscript PCA.plots.R <PCA.csv> <admixture.csv> <Unrelated.txt> <out.png>\n",
      "Where:\n",
      "  <PCA.csv>       = CSV with columns: FID,IID,PC1..PCk\n",
      "  <admixture.csv> = CSV with columns: FID,IID,K1,K2,...\n",
      "  <Unrelated.txt> = Two columns (FID IID), header optional\n",
      "  <out.png>       = Output figure filename\n")
  quit(status = 1)
}

pca_file  <- args[1]
admix_csv <- args[2]
unrel_txt <- args[3]
outfile   <- args[4]

# ---------- Read PCA ----------
pcs <- fread(pca_file)
req_pcs <- c("FID","IID","PC1","PC2")
if (!all(req_pcs %in% names(pcs))) {
  stop("PCA CSV must contain at least: FID, IID, PC1, PC2")
}

# ---------- Read ADMIXTURE ----------
admix <- fread(admix_csv)
if (!all(c("FID","IID") %in% names(admix))) {
  stop("admixture.csv must have columns: FID, IID, and cluster columns (K1,K2,...)")
}
k_cols <- setdiff(names(admix), c("FID","IID"))
if (!length(k_cols)) stop("No cluster columns found in admixture.csv")
# Normalize
row_sums <- admix[, rowSums(.SD), .SDcols = k_cols]
if (any(abs(row_sums - 1) > 1e-6)) {
  admix[, (k_cols) := lapply(.SD, function(x) x / pmax(rowSums(.SD), .Machine$double.eps)), .SDcols = k_cols]
}

# ---------- Read unrelated list ----------
unrel <- fread(unrel_txt, header = FALSE)
if (ncol(unrel) < 2) stop("Unrelated list must have two columns: FID IID")
setnames(unrel, 1:2, c("FID","IID"))

# Flag unrelated
pcs[, Unrelated := FALSE]
pcs[unrel, on = .(FID,IID), Unrelated := TRUE]
# ---------- PCA Plot ----------
p_pca <- ggplot() +
  geom_point(
    data = pcs[Unrelated == FALSE & !is.na(FID)],
    aes(PC1, PC2, color = FID),
    size = 1.8, alpha = 0.85, shape = 16
  ) +
  geom_point(
    data = pcs[Unrelated == TRUE],
    aes(PC1, PC2),
    color = "black", size = 2.6, shape = 16
  ) +
  labs(title = "PCA (PC1 vs PC2)",
       x = "PC1", y = "PC2",
       color = "FID") +
  theme_bw(base_size = 13) +
  theme(legend.position = "right")

# ---------- ADMIXTURE Plot (Unrelated only, no faceting) ----------
admix_unrel <- merge(unrel, admix, by = c("FID","IID"), all.x = TRUE)
admix_unrel <- admix_unrel[complete.cases(admix_unrel[, ..k_cols])]

if (nrow(admix_unrel) == 0) {
  p_adm <- ggplot() + theme_void() + labs(title = "No unrelated in admixture CSV")
} else {
  # Order by K1,K2,... globally
  o <- do.call(order, as.list(admix_unrel[, ..k_cols]))
  admix_unrel <- admix_unrel[o]
  admix_unrel[, x := factor(IID, levels = unique(IID))]

  long <- melt(admix_unrel[, c("IID","x",k_cols), with = FALSE],
               id.vars = c("IID","x"),
               variable.name = "Cluster", value.name = "Prop")

  p_adm <- ggplot(long, aes(x = x, y = Prop, fill = Cluster)) +
    geom_bar(stat = "identity", width = 1) +
    scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
    labs(title = "ADMIXTURE — Unrelated only",
         x = "Individuals (unrelated)", y = "Ancestry proportion") +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}

# ---------- Combine ----------
combo <- ggarrange(p_pca, p_adm, ncol = 2, widths = c(1, 1.15))
ggsave(outfile, combo, width = 14, height = 6.5, dpi = 300)
message("Saved combined PCA + Admixture plot to: ", normalizePath(outfile))
