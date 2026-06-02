#!/usr/bin/env Rscript
# =============================================================================
# 03f_ohnolog_fc_pattern.R
# -----------------------------------------------------------------------------
# Purpose
#   Mechanistic follow-up to the C1 co-clustering analysis. For each ohnolog
#   pair, quantify three properties and ask whether cluster-discordant pairs
#   (the two copies in different SOM clusters) are predominantly:
#
#     (a) NOISE-driven       — one copy is silenced / near the TPM floor;
#                              its SOM assignment is effectively random.
#     (b) LEVEL-divergence   — both copies expressed with similar temporal or
#                              tissue pattern (high pattern_r) but different
#                              absolute levels (high mag_gap).
#     (c) PATTERN-divergence — both copies expressed at comparable levels but
#                              with genuinely different shapes (low pattern_r);
#                              classic neo-/sub-functionalisation signature.
#
#   Report per-map (DevMap + BodyMap), stratified by AORe vs LORe, and test
#   whether the mechanism mix differs between rediploidization classes.
#
# Metrics (per pair)
#   mag_gap       |log2((mean_tpm_A + 1) / (mean_tpm_B + 1))|      — absolute
#                 log2 fold-change on RAW TPM (the SOM input is row-scaled,
#                 so magnitude must come from raw TPM).
#   pattern_r     Pearson correlation of the two copies' row-scaled SOM-input
#                 profiles (sd_norm / sb_norm rows). This is literally the
#                 shape similarity the SOM sees.
#   min_mean_tpm  min(mean_tpm_A, mean_tpm_B). Noise-floor proxy; if the
#                 weaker copy is near zero TPM, its SOM assignment reflects
#                 sampling noise more than biology.
#   log10_min_tpm log10(min_mean_tpm + 1). Used in regression.
#
# Fixed mechanism thresholds (see header rationale; thresholds chosen so that
# bucket labels are biologically interpretable rather than distribution-shaped)
#   LOW_TPM       min_mean_tpm < 1 TPM  → near noise floor
#   HIGH_MAG      mag_gap      > 1      → >2-fold magnitude gap
#   HIGH_PATTERN  pattern_r    > 0.5    → clearly-similar shape
#
# Inputs
#   results/ohnolog/data/c1_coclustering_results.RData   (from 03)
#   data/original_norm_files/norm.RData                  (sd_norm, sb_norm)
#   data/devmap/salmon/rsem.merged.gene_tpm.tsv
#   data/bodymap/salmon/Atlantic_salmon_*_rsem.merged.gene_tpm.tsv
#
# Outputs
#   data/c1_fc_pattern_results.RData   (pair-level tables + test outputs)
#   panels/C1_{devmap,bodymap}_fc_pattern_hex.{pdf,png}
#   panels/C1_{devmap,bodymap}_fc_pattern_ridges.{pdf,png}
#   panels/C1_{devmap,bodymap}_mechanism_stacked.{pdf,png}
#   panels/C1_{devmap,bodymap}_logit_forest.{pdf,png}
#   panels/C1_{devmap,bodymap}_example_trajectories.{pdf,png}
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggridges)
  library(broom)
  library(scales)
  library(patchwork)
})

# --- Global font (Helvetica) ------------------------------------------------
font_family <- "Helvetica"
theme_set(theme_bw(base_size = 11, base_family = font_family))
update_geom_defaults("text",  list(family = font_family))
update_geom_defaults("label", list(family = font_family))

# --- Paths -------------------------------------------------------------------
# Generated intermediates from earlier ohnolog scripts (read) + this script's
# own output panels (written). All other inputs are loaded via the shared
# loaders / Salmobase, not from local data/ paths.
data_dir  <- "results/ohnolog/data"
outdir    <- "results/ohnolog/panels"

# --- Thresholds --------------------------------------------------------------
# LOW_STAGE_MAX: "near-silent" boundary on max condition-averaged TPM. The
# upstream SOM-inclusion filter is `max(replicate-averaged TPM across
# conditions) > 1` (see Data_Loading_and_Normalization_Functions.R::normalize
# line 100). Genes that just squeak past that gate (max < 3 TPM across all
# conditions) are operationally indistinguishable from noise; we treat them
# as "silenced" for mechanism classification.
LOW_STAGE_MAX <- 3     # copy's max condition-averaged TPM floor
HIGH_MAG      <- 1     # >2-fold magnitude gap (|log2 FC|)
HIGH_PATTERN  <- 0.5   # clearly-similar shape (Pearson r)

redip_cols <- c("AORe" = "#c0392b", "LORe" = "#2980b9")
same_cols  <- c("TRUE" = "#2ecc71", "FALSE" = "#7f8c8d")

# Two silencing categories now replace the old monolithic "Noise-driven"
# bucket:
#   Both silenced — both copies below the floor. The pair's SOM
#                   assignments are both unreliable; these pairs
#                   are low-confidence and uninterpretable.
#   One silenced  — exactly one copy below the floor. A real asymmetric
#                   silencing story — one paralog has lost expression.
mech_levels <- c("Both silenced", "One silenced", "Level-divergence",
                 "Pattern-divergence", "Near-boundary")
mech_cols <- c("Both silenced"      = "#95a5a6",   # medium grey, low-confidence
               "One silenced"       = "#34495e",   # dark slate, biological
               "Level-divergence"   = "#e67e22",   # orange
               "Pattern-divergence" = "#9b59b6",   # purple
               "Near-boundary"      = "#bdc3c7")   # light grey

# -----------------------------------------------------------------------------
# save_gg — paired PDF + PNG save, routed into
# panels/{pdf|png}/{devmap|bodymap}/{const|noconst}/ by filename pattern
# (see 03c for details). FC-pattern panels don't have a `_noconst` variant,
# so they all land in const/ by default.
save_gg <- function(fname_base, p, w, h, dpi = 300) {
  name    <- basename(fname_base)
  map     <- if (grepl("_bodymap_", name)) "bodymap" else "devmap"
  variant <- if (grepl("_noconst$", name)) "noconst" else "const"
  pdf_dir <- file.path(outdir, "pdf", map, variant)
  png_dir <- file.path(outdir, "png", map, variant)
  dir.create(pdf_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)
  # cairo_pdf so Unicode glyphs in any label render safely; also honours
  # theme's base_family.
  ggsave(file.path(pdf_dir, paste0(name, ".pdf")), p, width = w, height = h,
         device = cairo_pdf)
  ggsave(file.path(png_dir, paste0(name, ".png")), p, width = w, height = h, dpi = dpi)
  cat("  Saved", name, "→", map, "/", variant, "\n")
}

# -----------------------------------------------------------------------------
# Load per-pair table (from 03) and the row-scaled SOM-input matrices.
# -----------------------------------------------------------------------------
cat("Loading co-clustering results + norm matrices + gene symbols...\n")
load(file.path(data_dir, "c1_coclustering_results.RData"))   # written by script 03
load_deposited_norm()   # injects sd_norm, td_norm, sb_norm, tb_norm

# Optional: friendly gene symbols from biomaRt (fetch_gene_symbols.R). Fall
# back to the last 6 characters of the Ensembl ID when a symbol is missing.
if (file.exists("data/ohnolog/gene_symbols.RData")) {
  load("data/ohnolog/gene_symbols.RData")   # gene_symbols
} else {
  message("gene_symbols.RData not found — run fetch_gene_symbols.R first.")
  gene_symbols <- character(0)
}
gene_symbol_lookup <- function(gid) {
  sym <- gene_symbols[gid]
  ifelse(is.na(sym) | !nzchar(sym),
         paste0("…", substr(gid, nchar(gid) - 5, nchar(gid))),
         sym)
}

# -----------------------------------------------------------------------------
# DevMap raw TPM: single file, columns gene_id, transcript_id(s), then 42
# replicate columns. We keep only the numeric replicate columns and compute
# the per-gene mean TPM.
# -----------------------------------------------------------------------------
cat("Loading DevMap raw TPM (salmon) from Salmobase...\n")
# Absolute expression magnitude comes from the Salmobase RNA-seq TPM, not from
# the deposited normalised matrices (which carry expression SHAPE, not level).
# read_table_anywhere() downloads once and caches under cache/.
devmap_tpm_raw <- read_table_anywhere(paths$recompute$devmap_salmon)
devmap_tpm_mat <- devmap_tpm_raw %>%
  select(-`transcript_id(s)`) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()
devmap_mean_tpm <- rowMeans(devmap_tpm_mat)
devmap_max_tpm  <- apply(devmap_tpm_mat, 1, max)

# Per-condition (replicate-averaged) max TPM. Column names are
# "<stage>_R<N>"; strip the _R suffix to group and mean.
devmap_cond       <- sub("_R[0-9]+$", "", colnames(devmap_tpm_mat))
devmap_stage_avg  <- sapply(unique(devmap_cond), function(cond) {
  rowMeans(devmap_tpm_mat[, devmap_cond == cond, drop = FALSE])
})
devmap_stage_max_tpm <- apply(devmap_stage_avg, 1, max)

# -----------------------------------------------------------------------------
# BodyMap raw TPM: seven per-tissue files, each 69 390 × ~15 replicate columns.
# Load all, drop metadata columns, bind on `gene_id`, compute per-gene mean.
# -----------------------------------------------------------------------------
cat("Loading BodyMap raw TPM (salmon, 7 tissues) from Salmobase...\n")
# process_tissue_data() downloads the seven per-tissue Salmobase TPM files and
# concatenates them into one gene x 84-replicate matrix (rownames = gene_id) —
# the same loader 01_data_loading.Rmd uses. We then summarise per gene below.
bodymap_tpm_mat  <- as.matrix(process_tissue_data(
  "salmon", urls = paths$recompute$bodymap_salmon_urls))
bodymap_mean_tpm <- rowMeans(bodymap_tpm_mat)
bodymap_max_tpm  <- apply(bodymap_tpm_mat, 1, max)

# Per-condition (replicate-averaged) max TPM. BodyMap column names are of
# the form "AtlanticSalmon_RNA_Liver_Immature_Female_R1"; strip _R\d+.
bodymap_cond       <- sub("_R[0-9]+$", "", colnames(bodymap_tpm_mat))
bodymap_cond_avg   <- sapply(unique(bodymap_cond), function(cond) {
  rowMeans(bodymap_tpm_mat[, bodymap_cond == cond, drop = FALSE])
})
bodymap_stage_max_tpm <- apply(bodymap_cond_avg, 1, max)

# =============================================================================
# build_pair_metrics — compute (mag_gap, pattern_r, min_mean_tpm,
# log10_min_tpm, mech_bucket) for every pair in `detail_df`.
#
# Args:
#   detail_df     — $detail tibble from 03 (gene1, gene2, cluster1, cluster2,
#                   same_cluster, proxiPhylogeny).
#   norm_mat      — row-scaled SOM input (sd_norm or sb_norm).
#   mean_tpm_vec  — named numeric vector: gene_id → mean raw TPM.
#
# Returns: tibble with metrics + mechanism classification for discordant pairs.
# =============================================================================
build_pair_metrics <- function(detail_df, norm_mat, mean_tpm_vec,
                                max_tpm_vec, stage_max_tpm_vec) {
  df <- detail_df %>% as_tibble() %>%
    mutate(
      tpm_A         = mean_tpm_vec[gene1],
      tpm_B         = mean_tpm_vec[gene2],
      max_tpm_A     = max_tpm_vec[gene1],
      max_tpm_B     = max_tpm_vec[gene2],
      stage_max_A   = stage_max_tpm_vec[gene1],
      stage_max_B   = stage_max_tpm_vec[gene2]
    ) %>%
    filter(!is.na(tpm_A), !is.na(tpm_B),
           gene1 %in% rownames(norm_mat),
           gene2 %in% rownames(norm_mat))

  # Magnitude gap on raw TPM (absolute log2 FC; direction-agnostic).
  df <- df %>% mutate(
    mag_gap        = abs(log2((tpm_A + 1) / (tpm_B + 1))),
    min_mean_tpm   = pmin(tpm_A, tpm_B),
    log10_min_tpm  = log10(min_mean_tpm + 1)
  )

  # Pattern r: Pearson correlation of the two copies' row-scaled profiles.
  # Done in a tight loop (vectorising needs aligned subsetting; ~10k pairs
  # is fast enough without further optimisation).
  pr <- numeric(nrow(df))
  for (i in seq_len(nrow(df))) {
    pr[i] <- suppressWarnings(cor(norm_mat[df$gene1[i], ],
                                   norm_mat[df$gene2[i], ]))
  }
  df$pattern_r <- pr

  # Mechanism classification. Silencing is evaluated on max
  # condition-averaged TPM (matches the upstream SOM-inclusion filter),
  # split into Both / One silenced so the "one copy lost expression"
  # biological story is separated from the "both copies at noise, SOM
  # assignment random on both" low-confidence case.
  df <- df %>% mutate(
    mech_bucket = factor(case_when(
      stage_max_A < LOW_STAGE_MAX & stage_max_B < LOW_STAGE_MAX ~ "Both silenced",
      stage_max_A < LOW_STAGE_MAX | stage_max_B < LOW_STAGE_MAX ~ "One silenced",
      pattern_r   <= HIGH_PATTERN                               ~ "Pattern-divergence",
      mag_gap     >  HIGH_MAG  & pattern_r > HIGH_PATTERN       ~ "Level-divergence",
      TRUE                                                       ~ "Near-boundary"
    ), levels = mech_levels)
  )

  df
}

cat("Computing pair metrics...\n")
devmap_pair_metrics  <- build_pair_metrics(devmap_result$detail,  sd_norm,
                                           devmap_mean_tpm, devmap_max_tpm,
                                           devmap_stage_max_tpm)
bodymap_pair_metrics <- build_pair_metrics(bodymap_result$detail, sb_norm,
                                           bodymap_mean_tpm, bodymap_max_tpm,
                                           bodymap_stage_max_tpm)

cat("  DevMap pairs with complete metrics: ",  nrow(devmap_pair_metrics),  "\n")
cat("  BodyMap pairs with complete metrics:",  nrow(bodymap_pair_metrics), "\n")

# =============================================================================
# run_stats — Wilcoxon / KS / logistic regression / mechanism χ² for one map.
# Returns a list that is attached to the RData output; the panels below use
# the same input frames so no coupling to `run_stats` internals.
# =============================================================================
run_stats <- function(pair_df, map_label) {
  cat("\n=== ", map_label, " ===\n")

  # --- 1. Univariate comparisons: same_cluster=T vs F within each class ----
  metrics <- c("mag_gap", "pattern_r", "log10_min_tpm")
  uni <- expand_grid(class = c("AORe", "LORe"), metric = metrics) %>%
    rowwise() %>%
    mutate(test = list({
      sub <- pair_df %>% filter(proxiPhylogeny == class) %>%
        pull(all_of(metric))
      grp <- pair_df %>% filter(proxiPhylogeny == class) %>% pull(same_cluster)
      w <- suppressWarnings(wilcox.test(sub ~ grp))
      k <- suppressWarnings(ks.test(sub[grp], sub[!grp]))
      tibble(wilcox_p = w$p.value, ks_p = k$p.value,
             median_same = median(sub[grp],  na.rm = TRUE),
             median_diff = median(sub[!grp], na.rm = TRUE))
    })) %>%
    unnest(test) %>%
    ungroup() %>%
    mutate(wilcox_padj = p.adjust(wilcox_p, method = "BH"),
           ks_padj     = p.adjust(ks_p,     method = "BH"))
  cat("\n-- Univariate tests (Wilcoxon / KS, BH-adjusted) --\n")
  print(as.data.frame(uni))

  # --- 2. Multivariable logistic regression --------------------------------
  # Standardise predictors so OR is "per SD". Fit separately per class so
  # the mechanism question can be reported independently for AORe and LORe.
  fit_logit <- function(sub) {
    d <- sub %>% transmute(
      y       = as.integer(same_cluster),
      mag_gap = as.numeric(scale(mag_gap)),
      pat_r   = as.numeric(scale(pattern_r)),
      lmin    = as.numeric(scale(log10_min_tpm))
    )
    m <- glm(y ~ mag_gap + pat_r + lmin, data = d, family = binomial())
    broom::tidy(m, conf.int = TRUE, exponentiate = TRUE)
  }
  logit_aore <- pair_df %>% filter(proxiPhylogeny == "AORe") %>% fit_logit() %>%
    mutate(class = "AORe")
  logit_lore <- pair_df %>% filter(proxiPhylogeny == "LORe") %>% fit_logit() %>%
    mutate(class = "LORe")
  logit <- bind_rows(logit_aore, logit_lore) %>%
    filter(term != "(Intercept)")
  cat("\n-- Logistic regression (OR per SD) --\n")
  print(as.data.frame(logit))

  # --- 3. Mechanism bucket distribution among DISCORDANT pairs -------------
  # (co-clustered pairs by definition have mag_gap ~ 0 and high pattern_r,
  # so their mechanism label is uninformative; we report on discordants.)
  mech <- pair_df %>%
    filter(!same_cluster) %>%
    count(proxiPhylogeny, mech_bucket) %>%
    group_by(proxiPhylogeny) %>%
    mutate(pct = 100 * n / sum(n)) %>% ungroup()
  cat("\n-- Mechanism distribution (discordant pairs only) --\n")
  print(as.data.frame(mech))

  mech_mat <- mech %>% select(proxiPhylogeny, mech_bucket, n) %>%
    pivot_wider(names_from = mech_bucket, values_from = n, values_fill = 0) %>%
    column_to_rownames("proxiPhylogeny") %>% as.matrix()
  chi <- suppressWarnings(chisq.test(mech_mat))
  cat("Mechanism χ² (AORe vs LORe):  X2 =", round(chi$statistic, 2),
      "  df =", chi$parameter,
      "  p =", format(chi$p.value, digits = 3), "\n")

  list(univariate = uni, logit = logit, mechanism = mech, chisq = chi)
}

devmap_stats  <- run_stats(devmap_pair_metrics,  "DevMap")
bodymap_stats <- run_stats(bodymap_pair_metrics, "BodyMap")

# =============================================================================
# Panels
# =============================================================================
cat("\n=== Panels ===\n")

# --- 2D hex plot: mag_gap × pattern_r, split by proxiPhylogeny, shaded
# by same_cluster status (via facet). ----------------------------------------
# Using 2D hex with facet rather than ggExtra::ggMarginal so the facet grid
# shows all four (proxiPhylogeny × same_cluster) combinations cleanly.
make_hex_panel <- function(pair_df, fname) {
  p <- ggplot(pair_df,
              aes(x = mag_gap, y = pattern_r)) +
    geom_hex(bins = 55) +
    geom_hline(yintercept = HIGH_PATTERN, linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = HIGH_MAG,     linetype = "dashed", colour = "grey40") +
    scale_fill_viridis_c(option = "mako", trans = "log10",
                         name = "Pairs (log10)") +
    facet_grid(same_cluster ~ proxiPhylogeny,
               labeller = labeller(
                 same_cluster = c(`TRUE` = "Co-clustered",
                                  `FALSE` = "Cluster-discordant"),
                 proxiPhylogeny = c(AORe = "Early (AORe)", LORe = "Late (LORe)"))) +
    labs(x = expression("|log"[2]*"(mean TPM ratio between copies)|"),
         y = "Pattern correlation (Pearson, row-scaled profile)") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          legend.position  = "right")
  save_gg(fname, p, 8, 6)
}
make_hex_panel(devmap_pair_metrics,  file.path(outdir, "C1_devmap_fc_pattern_hex"))
make_hex_panel(bodymap_pair_metrics, file.path(outdir, "C1_bodymap_fc_pattern_hex"))

# --- Ridge plots: three metrics, split by same_cluster × proxiPhylogeny -----
make_ridge_panel <- function(pair_df, fname) {
  long <- pair_df %>%
    select(proxiPhylogeny, same_cluster, mag_gap, pattern_r, log10_min_tpm) %>%
    pivot_longer(c(mag_gap, pattern_r, log10_min_tpm),
                 names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric,
                            levels = c("mag_gap", "pattern_r", "log10_min_tpm"),
                            labels = c("|log2 mag. gap|",
                                       "Pattern r",
                                       "log10(min mean TPM + 1)")),
           facet = ifelse(same_cluster, "Co-clustered", "Cluster-discordant"))

  p <- ggplot(long, aes(x = value, y = facet, fill = proxiPhylogeny)) +
    geom_density_ridges(alpha = 0.55, scale = 0.95, rel_min_height = 0.01,
                        colour = NA) +
    scale_fill_manual(values = redip_cols,
                      labels = c("AORe" = "Early", "LORe" = "Late"),
                      name   = "Rediploidization") +
    facet_wrap(~ metric, scales = "free_x", nrow = 1) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position  = "bottom",
          panel.grid.minor = element_blank())
  save_gg(fname, p, 10, 4.5)
}
make_ridge_panel(devmap_pair_metrics,  file.path(outdir, "C1_devmap_fc_pattern_ridges"))
make_ridge_panel(bodymap_pair_metrics, file.path(outdir, "C1_bodymap_fc_pattern_ridges"))

# --- Stacked mechanism bars (discordant pairs only) -------------------------
make_mechanism_panel <- function(pair_df, fname) {
  d <- pair_df %>% filter(!same_cluster) %>%
    count(proxiPhylogeny, mech_bucket) %>%
    group_by(proxiPhylogeny) %>%
    mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
    mutate(proxiPhylogeny = factor(proxiPhylogeny,
             levels = c("AORe", "LORe"),
             labels = c("Early (AORe)", "Late (LORe)")))
  p <- ggplot(d, aes(x = proxiPhylogeny, y = pct, fill = mech_bucket)) +
    geom_col(width = 0.6, colour = "white") +
    geom_text(aes(label = sprintf("%.0f%%\nn=%d", pct, n)),
              position = position_stack(vjust = 0.5),
              size = 2.6, colour = "grey10") +
    scale_fill_manual(values = mech_cols, name = "Mechanism", drop = FALSE) +
    labs(x = NULL, y = "Cluster-discordant pairs (%)") +
    theme_bw(base_size = 10) +
    theme(panel.grid.major.x = element_blank(),
          legend.position    = "right")
  save_gg(fname, p, 6, 5)
}
make_mechanism_panel(devmap_pair_metrics,  file.path(outdir, "C1_devmap_mechanism_stacked"))
make_mechanism_panel(bodymap_pair_metrics, file.path(outdir, "C1_bodymap_mechanism_stacked"))

# --- Logistic regression forest ---------------------------------------------
make_forest_panel <- function(stats_obj, fname) {
  d <- stats_obj$logit %>%
    mutate(term = factor(term,
             levels = c("mag_gap", "pat_r", "lmin"),
             labels = c("|log2 mag. gap| (z)",
                        "Pattern r (z)",
                        "log10 min TPM (z)")),
           class = factor(class, levels = c("AORe", "LORe"),
                          labels = c("Early (AORe)", "Late (LORe)")))
  p <- ggplot(d, aes(x = estimate, y = term, colour = class)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high),
                    position = position_dodge(width = 0.5),
                    size = 0.35) +
    scale_colour_manual(values = c("Early (AORe)" = unname(redip_cols["AORe"]),
                                   "Late (LORe)"  = unname(redip_cols["LORe"])),
                        name = "Rediploidization") +
    scale_x_continuous(trans = "log2",
                       breaks = c(0.5, 0.75, 1, 1.5, 2, 3)) +
    labs(x = "Odds ratio (same cluster) per SD",
         y = NULL) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          legend.position  = "bottom")
  save_gg(fname, p, 6, 4)
}
make_forest_panel(devmap_stats,  file.path(outdir, "C1_devmap_logit_forest"))
make_forest_panel(bodymap_stats, file.path(outdir, "C1_bodymap_logit_forest"))

# --- Example trajectories: 3 pairs per mechanism bucket ---------------------
# For each mechanism bucket, pick the pair whose (mag_gap, pattern_r,
# log10_min_tpm) is closest to the bucket centroid — i.e. a "representative"
# pair. Three picks per bucket = 12 lines per figure.
make_example_trajectories <- function(pair_df, norm_mat, fname, n_per = 4) {
  # Pick n_per pairs per mechanism bucket whose (mag_gap, pattern_r,
  # log10_min_tpm) is closest to the bucket centroid (in standardised
  # space). These are "representative" pairs for each mechanism.
  rep_pairs <- pair_df %>% filter(!same_cluster) %>%
    group_by(mech_bucket) %>%
    filter(n() > 0) %>%
    mutate(
      z_mag  = as.numeric(scale(mag_gap)),
      z_pat  = as.numeric(scale(pattern_r)),
      z_tpm  = as.numeric(scale(log10_min_tpm)),
      dist   = sqrt(z_mag^2 + z_pat^2 + z_tpm^2)
    ) %>%
    slice_min(dist, n = n_per) %>% ungroup()

  # Friendly copy labels: "SYMBOL (max 17 TPM)". Falls back to the last 6
  # characters of the Ensembl ID when no external name is available.
  rep_pairs <- rep_pairs %>%
    mutate(
      sym_A   = gene_symbol_lookup(gene1),
      sym_B   = gene_symbol_lookup(gene2),
      label_A = sprintf("%s (max %.1f TPM)", sym_A, max_tpm_A),
      label_B = sprintf("%s (max %.1f TPM)", sym_B, max_tpm_B),
      pair_id = sprintf("%s / %s", sym_A, sym_B)
    )

  traj <- rep_pairs %>%
    rowwise() %>%
    mutate(prof = list(
             tibble(stage   = colnames(norm_mat),
                    value_A = norm_mat[gene1, ],
                    value_B = norm_mat[gene2, ])
           )) %>%
    unnest(prof) %>% ungroup() %>%
    pivot_longer(c(value_A, value_B), names_to = "copy_raw", values_to = "expr") %>%
    mutate(copy_label = ifelse(copy_raw == "value_A", label_A, label_B),
           stage      = factor(stage, levels = colnames(norm_mat)))

  # Legend shows the two copy labels per panel, which vary — easier to put
  # labels into the plot area. Use geom_text at the last stage for each copy.
  text_data <- traj %>%
    group_by(pair_id, mech_bucket, copy_raw, copy_label) %>%
    slice_max(as.integer(stage), n = 1, with_ties = FALSE) %>%
    ungroup()

  p <- ggplot(traj, aes(x = stage, y = expr, group = copy_raw, colour = copy_raw)) +
    geom_line(linewidth = 0.5) + geom_point(size = 0.8) +
    geom_text(data = text_data, aes(label = copy_label),
              hjust = 1, vjust = -0.4, size = 2.0, show.legend = FALSE) +
    scale_colour_manual(values = c("value_A" = "#e74c3c", "value_B" = "#3498db"),
                        labels = c("value_A" = "Copy A", "value_B" = "Copy B"),
                        name = NULL) +
    facet_wrap(mech_bucket ~ pair_id, ncol = n_per, scales = "free_y") +
    labs(x = NULL, y = "Row-scaled expression") +
    theme_bw(base_size = 8, base_family = font_family) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          strip.text  = element_text(size = 6),
          legend.position = "bottom",
          panel.grid.minor = element_blank())
  save_gg(fname, p, 14, 10)
}
make_example_trajectories(devmap_pair_metrics,  sd_norm,
                          file.path(outdir, "C1_devmap_example_trajectories"))
make_example_trajectories(bodymap_pair_metrics, sb_norm,
                          file.path(outdir, "C1_bodymap_example_trajectories"))

# -----------------------------------------------------------------------------
# Save pair-level frames + stat outputs so downstream scripts / manuscript text
# can pull numbers directly without re-deriving.
# -----------------------------------------------------------------------------
save(devmap_pair_metrics, bodymap_pair_metrics,
     devmap_stats,        bodymap_stats,
     LOW_STAGE_MAX, HIGH_MAG, HIGH_PATTERN,
     file = file.path(data_dir, "c1_fc_pattern_results.RData"))

cat("\nResults saved to data/c1_fc_pattern_results.RData\n")
cat("Done.\n")
