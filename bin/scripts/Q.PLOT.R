#!/usr/bin/env Rscript

# Usage:
#   Rscript plot_admixture_csv.R <input.csv> [out.png] [UNRELATED_FILE] [UNRELATED_OUT.png]
#
# Example:
#   Rscript plot_admixture_csv.R admixture.csv all.png unrelated.txt unrelated.png
#
# Input CSV must have columns: FID, IID, and one or more cluster columns (e.g., K1,K2,... or named clusters)

suppressPackageStartupMessages({
  library(data.table)
  library(reshape2)
  library(ggplot2)
  library(psych)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage:\n  Rscript plot_admixture_csv.R <input.csv> [out.png] [UNRELATED_FILE] [UNRELATED_OUT.png]\n")
  quit(status = 1)
}

infile     <- args[1]
outfile    <- ifelse(length(args) >= 2 && nzchar(args[2]), args[2], "admixture.png")
unrel_file <- ifelse(length(args) >= 3 && nzchar(args[3]), args[3], NA_character_)
unrel_out  <- ifelse(length(args) >= 4 && nzchar(args[4]), args[4], "admixture_unrelated.png")

# ---------- Read & validate main CSV ----------
dt <- fread(infile)
req <- c("FID", "IID")
if (!all(req %in% names(dt))) {
  stop("Input CSV must include columns: FID, IID, and one or more cluster columns.")
}
k_cols <- setdiff(names(dt), req)
if (length(k_cols) < 1) stop("No cluster columns found.")

# Normalize rows to sum to 1
row_sums <- dt[, rowSums(.SD), .SDcols = k_cols]
if (any(abs(row_sums - 1) > 1e-6)) {
  dt[, (k_cols) := lapply(.SD, function(x) x / pmax(rowSums(.SD), .Machine$double.eps)), .SDcols = k_cols]
}

# ---------- Plot helpers ----------
# Main plot: facet by FID and order within each FID by the first K column
make_plot_all <- function(dt_sub, out_path) {
  if (nrow(dt_sub) == 0) {
    warning("No rows to plot for ", out_path, "; skipping.")
    return(invisible(NULL))
  }

  # Order by first K within each FID (ascending)
  ord <- dt_sub[, .(IID, keyval = get(k_cols[1])), by = FID]
  ord[, rank := frank(keyval, ties.method = "first")]
  ord <- ord[order(FID, rank)]

  dt_sub[, x := paste0(FID, "::", IID)]
  ord[,   x := paste0(FID, "::", IID)]
  x_levels <- ord$x

  long <- melt(dt_sub[, c("FID", "IID", "x", k_cols), with = FALSE],
               id.vars = c("FID", "IID", "x"),
               variable.name = "Cluster",
               value.name = "Prop")
  long[, x := factor(x, levels = x_levels)]
  iid_labels <- function(v) sub("^.*::", "", v)

  p <- ggplot(long, aes(x = x, y = Prop, fill = Cluster)) +
    geom_bar(stat = "identity", width = 1) +
    facet_grid(~ FID, scales = "free_x", space = "free_x") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    scale_x_discrete(labels = iid_labels) +
    labs(x = "Individuals", y = "Ancestry proportion", title = "ADMIXTURE plot") +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text = element_text(face = "bold")
    )

  ggsave(out_path, p, width = 12+ 0.5*length(unique(long$FID)),
         height = 4 + 0.5*length(unique(long$FID)), dpi = 300)
  message("Plot saved to: ", normalizePath(out_path))
}

# Unrelated plot: NO faceting, global lexicographic order by all K columns (K1, then K2, …)
make_plot_unrelated <- function(dt_unrel, out_path) {
  if (nrow(dt_unrel) == 0) {
    warning("No unrelated rows to plot for ", out_path, "; skipping.")
    return(invisible(NULL))
  }

  # Global lexicographic order by all K columns (ascending)
  # Build an order data.table with only K columns to sort by
  ord_dt <- copy(dt_unrel)[, ..k_cols]
  # Order rows by K1, then K2, etc.
  o <- do.call(order, as.list(ord_dt))
  dt_unrel <- dt_unrel[o]

  # x = IID only (omit FID completely)
  dt_unrel[, x := IID]

  long <- melt(dt_unrel[, c("IID", "x", k_cols), with = FALSE],
               id.vars = c("IID", "x"),
               variable.name = "Cluster",
               value.name = "Prop")
  long[, x := factor(x, levels = unique(dt_unrel$x))]

  p <- ggplot(long, aes(x = x, y = Prop, fill = Cluster)) +
    geom_bar(stat = "identity", width = 1) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    labs(x = "Individuals (unrelated)", y = "Ancestry proportion",
         title = "ADMIXTURE plot — Unrelated only") +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_blank(),   # set to element_text(angle=90, vjust=0.5) if you want labels
      axis.ticks.x = element_blank()
    )

  ggsave(out_path, p, width = 12, height = 4, dpi = 300)
  message("Unrelated plot saved to: ", normalizePath(out_path))
}

# ---------- Main plot ----------
make_plot_all(copy(dt), outfile)

# ---------- Unrelated-only plot (no FID clustering/facets; sort by Ks globally) ----------
if (!is.na(unrel_file)) {
  unrel <- fread(unrel_file, header = FALSE)
  if (ncol(unrel) >= 2) setnames(unrel, 1:2, c("FID","IID"))
  # Join on FID+IID to keep only those present in the CSV
  dt_unrel <- merge(unrel, dt, by = c("FID","IID"), all.x = TRUE)
  dt_unrel <- dt_unrel[complete.cases(dt_unrel)]
  # Drop FID entirely for plotting/ordering
  dt_unrel[, FID := NULL]
  make_plot_unrelated(dt_unrel, unrel_out)
}

# ---------- Descriptive stats ----------
desc_all <- describe(dt[, ..k_cols])
fwrite(as.data.table(desc_all, keep.rownames = "Cluster"),
       "admixture_descriptive_stats_overall.csv")
