#!/usr/bin/env Rscript
# =============================================================================
# 03b_ohnolog_som_extended.R
# -----------------------------------------------------------------------------
# Purpose
#   Per-cluster, gene-level breakdown of AORe vs LORe ohnolog composition and
#   co-clustering behaviour. Complements 03 (which operates at the pair level
#   across the whole map) by asking:
#     1. Within each SOM cluster, what fraction of resident ohnolog genes are
#        AORe vs LORe? (tests non-uniform LORe enrichment across clusters).
#     2. Within each SOM cluster, of the ohnolog genes that live here, what
#        fraction have their partner ALSO in this cluster? (per-cluster
#        co-clustering rate — the metric shown in 03c's bar panels and in the
#        manuscript).
#
# Analysis chain (see 03 header for full list)
#   03  → c1_coclustering_results.RData   (pair-level)
#   03b → c1_extended_results.RData       (THIS SCRIPT — gene-level/per-cluster)
#   03c → reads c1_extended_results.RData to produce bar panels.
#
# Inputs (under data/ in this repository)
#   data/original_soms/soms.RData
#   data/original_norm_files/norm.RData          (needed here for cluster labels)
#   data/orthogroups/SalmonTroutOrthologs.tsv
#   data/ohnolog/som_info.RData                  (sb/sd_som_info + classifiers)
#
# Outputs → results/ohnolog/data/c1_extended_results.RData
#   devmap_gene_counts, bodymap_gene_counts
#       per-cluster AORe & LORe gene counts + pct_LORe + cluster label.
#   devmap_per_clust, bodymap_per_clust    (each: list(per_cluster, fisher))
#       $per_cluster  per-cluster co-clustering rate (pct_coclustered,
#                     n_genes, n_partner_here, proxiPhylogeny).
#       $fisher       per-cluster Fisher's Exact table for AORe vs LORe
#                     partner-here rates, BH-adjusted, clusters with n ≥ 20.
#   devmap_labels, bodymap_labels
#       cluster → ("Constitutive" | peak condition name) via label_clusters().
#
# Key analytical choices
#   * Cluster classification (Constitutive vs peak-condition) is INHERITED
#     from `sd_som_info` / `sb_som_info` in `data/ohnolog/som_info.RData`,
#     which are built by `som_file_info_d()` / `som_file_info_b()` (also
#     carried in that file). We do NOT re-derive it here — the 0.4 SD
#     threshold etc. live in the canonical function and re-implementing
#     them risks silent drift.
#   * Ohnolog filter uses nSsal==2 only (NOT nSsal==2 & nOmyk==2); see
#     rationale in the 03_ohnolog_som_coclustering.R header.
#   * Gene-level co-clustering is computed by exploding each pair into two
#     directed rows (gene, partner), so the per-cluster rate is symmetric in
#     the sense that both copies count towards whichever cluster they sit in.
#   * Per-cluster Fisher's tests only run where both AORe and LORe have
#     n_genes ≥ 20 (below that, the 2×2 is too sparse for a meaningful OR).
#     p-values are BH-adjusted across the surviving clusters.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# Sourced by 04_ohnolog_coclustering.Rmd with the working directory at
# transcriptome-analysis/, AFTER R/setup.R (so the shared loaders, `paths` and
# read_table_anywhere() are in scope).

# -----------------------------------------------------------------------------
# Load the deposited SOMs and their normalised matrices via the shared loaders,
# the orthogroup table from Salmobase, AND the precomputed cluster-info tables
# (sd_som_info, sb_som_info). The cluster-info tables are the canonical output
# of som_file_info_d()/_b() and are what the rest of the manuscript pipeline
# uses — we consume them directly so labels never drift.
# -----------------------------------------------------------------------------
cat("Loading data...\n")
load_deposited_soms()   # injects sd_som, td_som, sb_som, tb_som
load_deposited_norm()   # injects sd_norm, td_norm, sb_norm, tb_norm
load("data/ohnolog/som_info.RData")   # sb_som_info, sd_som_info, som_file_info_b, som_file_info_d

# Orthogroup table: always fetched from Salmobase (cached under cache/).
og <- read_table_anywhere(paths$active$orthogroups)

# Sanity: the SOM cluster indices in sd_som$som$unit.classif must be exactly
# 1..nrow(sd_som_info). Fails loudly if the upstream pipeline has drifted
# (e.g. refit SOM with a different grid size since sd_som_info was built).
stopifnot(max(sd_som$som$unit.classif) == nrow(sd_som_info),
          max(sb_som$som$unit.classif) == nrow(sb_som_info))

# Gene → cluster tibbles (same structure as 03; column renamed to `cluster`
# since this script doesn't need to carry devmap/bodymap disambiguation in
# the column name — the two maps are never joined together here).
devmap_clusters <- tibble(geneID = rownames(sd_norm),
                          cluster = sd_som$som$unit.classif)
bodymap_clusters <- tibble(geneID = rownames(sb_norm),
                           cluster = sb_som$som$unit.classif)

# =============================================================================
# Per-gene four-category classification for the composition panels (03c).
#
# Each Ssal gene falls into exactly one of:
#   Early      nSsal == 2 & proxiPhylogeny == "AORe"
#   Late       nSsal == 2 & proxiPhylogeny == "LORe"
#   Singleton  nSsal == 1 (any proxiPhylogeny)
#   Other      everything else — nSsal >= 3, nSsal == 2 with proxiPhylogeny
#              outside {AORe, LORe} (cryptic, NA, …), etc.
#
# The AORe/LORe bucket here uses the RAW orthogroup annotation (no subOG
# pair-consistency filter). That filter is specific to the pair-level
# co-clustering analysis; the composition plot is a per-gene histogram by
# cluster, so every annotated gene should contribute to exactly one bar.
# =============================================================================
ssal_category <- og %>%
  filter(spc == "Ssal") %>%
  transmute(
    geneID,
    category = case_when(
      nSsal == 2 & proxiPhylogeny == "AORe" ~ "Early",
      nSsal == 2 & proxiPhylogeny == "LORe" ~ "Late",
      nSsal == 1                            ~ "Singleton",
      TRUE                                   ~ "Other"
    )
  ) %>%
  distinct(geneID, .keep_all = TRUE)   # defensive: one row per geneID

cat("\nPer-gene category totals (Ssal, entire orthogroup table):\n")
print(table(ssal_category$category))

# Per-cluster composition: inner_join on geneID drops genes not in the SOM,
# then count (cluster, category) and compute within-cluster proportion.
build_composition <- function(cluster_df) {
  cluster_df %>%
    inner_join(ssal_category, by = "geneID") %>%
    count(cluster, category, name = "n") %>%
    group_by(cluster) %>%
    mutate(total = sum(n),
           pct   = 100 * n / total) %>%
    ungroup()
}
devmap_composition  <- build_composition(devmap_clusters)
bodymap_composition <- build_composition(bodymap_clusters)
cat("DevMap composition rows:",  nrow(devmap_composition),  "\n")
cat("BodyMap composition rows:", nrow(bodymap_composition), "\n")

# -----------------------------------------------------------------------------
# Thin adapters over the canonical cluster-info tables. Downstream code only
# reads `(cluster, label)`; we surface `label = "Constitutive"` for
# constitutive clusters and the upstream `class_max` (peak condition name)
# otherwise. The richer columns in *_som_info remain accessible via the
# loaded objects if needed.
# -----------------------------------------------------------------------------
devmap_labels <- sd_som_info %>%
  transmute(cluster = as.integer(class_number),
            label   = ifelse(as.logical(constitutive_class),
                             "Constitutive", as.character(class_max)))
bodymap_labels <- sb_som_info %>%
  transmute(cluster = as.integer(class_number),
            label   = ifelse(as.logical(constitutive_class),
                             "Constitutive", as.character(class_max)))

# -----------------------------------------------------------------------------
# Build ohnolog pair table. Exactly the same filtering as 03 — kept in this
# script rather than loaded from c1_coclustering_results.RData so that 03b
# can be run independently of 03 (both scripts read only from od_path).
# -----------------------------------------------------------------------------
ssal_genes <- og %>%
  filter(spc == "Ssal", nSsal == 2,
         proxiPhylogeny %in% c("AORe", "LORe")) %>%
  select(subOG, geneID, proxiPhylogeny)

# Keep only subOGs where both Salmon copies survived the filters (see 03).
pair_counts <- ssal_genes %>% count(subOG) %>% filter(n == 2)
ssal_genes <- ssal_genes %>% filter(subOG %in% pair_counts$subOG)

# Drop subOGs with inconsistent proxiPhylogeny across the two copies.
inconsistent <- ssal_genes %>%
  group_by(subOG) %>%
  summarise(n_phyl = n_distinct(proxiPhylogeny)) %>%
  filter(n_phyl > 1)
if (nrow(inconsistent) > 0) {
  cat("Removing", nrow(inconsistent), "pairs with inconsistent proxiPhylogeny\n")
  ssal_genes <- ssal_genes %>% filter(!subOG %in% inconsistent$subOG)
}

# Long → wide reshape (same as 03) — one row per pair.
ssal_pairs_wide <- ssal_genes %>%
  group_by(subOG) %>%
  mutate(copy = row_number()) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(subOG, proxiPhylogeny),
              names_from = copy, values_from = geneID,
              names_prefix = "gene")

# -----------------------------------------------------------------------------
# Explode back to one row per GENE (not per pair), with each gene carrying a
# pointer to its partner's geneID. This is the gene-level representation we
# need for the per-cluster "is my partner here?" test.
#
# The `bind_rows(…, …)` trick duplicates the pair table: first half takes
# gene1 as "me" and gene2 as "partner"; second half swaps. Result: each
# ohnolog gene appears exactly once as `geneID`, with `partnerID` filled.
# -----------------------------------------------------------------------------
gene_partner <- bind_rows(
  ssal_pairs_wide %>% select(subOG, proxiPhylogeny, geneID = gene1, partnerID = gene2),
  ssal_pairs_wide %>% select(subOG, proxiPhylogeny, geneID = gene2, partnerID = gene1)
)

cat("\nOhnolog genes:", nrow(gene_partner), "\n")
cat("By rediploidization:\n")
print(table(gene_partner$proxiPhylogeny))

# =============================================================================
# Singleton-inclusive composition gene pool.
#
# The composition analysis (Analysis 1) asks "of all salmon genes that can be
# classified Early or Late by proxiPhylogeny, how do they distribute across
# SOM clusters?". The relevant denominator is therefore every Ssal row with
# proxiPhylogeny in {AORe, LORe} that has a SOM cluster assignment, regardless
# of nSsal — this brings in salmon singletons (nSsal == 1) classified via
# rainbow trout retention as well as nSsal == 2 genes whose partner did not
# pass the SOM input filter.
#
# Co-clustering (Analysis 2) still requires both copies of an ohnolog pair in
# the SOM and so continues to use `gene_partner` and `result`.
# =============================================================================
ssal_classified <- og %>%
  filter(spc == "Ssal",
         proxiPhylogeny %in% c("AORe", "LORe")) %>%
  select(geneID, nSsal, proxiPhylogeny) %>%
  distinct(geneID, .keep_all = TRUE)

cat("\nAll Ssal genes classifiable as Early/Late (any nSsal):\n")
print(table(ssal_classified$proxiPhylogeny, ssal_classified$nSsal))

# =============================================================================
# run_gene_coclustering — two analyses per SOM map.
#
# Analysis 1: Gene counts per cluster (AORe vs LORe)
#   For each cluster, how many AORe and LORe ohnolog genes live there, and
#   what fraction are LORe? Chi-squared test on non-constitutive clusters
#   asks whether LORe proportion is uniform across clusters.
#
# Analysis 2: Gene-based co-clustering rate per cluster
#   For each cluster × class, what fraction of genes have their partner in
#   the same cluster? Per-cluster Fisher's Exact Test compares AORe vs LORe
#   rates; BH-adjusted across clusters with n ≥ 20 for both classes.
#
# Args:
#   gene_partner_df  gene-level table (geneID, partnerID, proxiPhylogeny).
#   cluster_df       geneID → cluster tibble.
#   cluster_labels   cluster → label (from label_clusters()).
#   label            log label for this map.
#
# Returns: list(gene_counts, coclust_rate, fisher).
# =============================================================================
run_gene_coclustering <- function(gene_partner_df, ssal_classified_df,
                                  cluster_df, cluster_labels, label) {
  cat("\n========================================\n")
  cat(label, "\n")
  cat("========================================\n")

  # comp_gene_pool = ALL Ssal genes with proxiPhylogeny in {AORe, LORe} that
  # have a SOM cluster assignment, regardless of nSsal. Used for ANALYSIS 1
  # (cluster composition). Includes singletons (nSsal==1) classified via
  # rainbow trout retention plus nSsal==2 ohnolog genes whose partner did
  # not pass the SOM input filter.
  comp_gene_pool <- ssal_classified_df %>%
    inner_join(cluster_df %>% rename(my_cluster = cluster), by = "geneID")

  # gene_pool = nSsal==2 ohnolog genes that landed in the SOM (any partner
  # status). result = subset where both copies are in the SOM. Used for
  # ANALYSIS 2 (co-clustering rate), which structurally requires the
  # partner to have a cluster assignment.
  gene_pool <- gene_partner_df %>%
    inner_join(cluster_df %>% rename(my_cluster = cluster), by = "geneID")
  result <- gene_pool %>%
    left_join(cluster_df %>% rename(partner_cluster = cluster), by = c("partnerID" = "geneID")) %>%
    filter(!is.na(partner_cluster)) %>%
    mutate(partner_here = my_cluster == partner_cluster)

  cat("Composition pool (singleton-inclusive, any nSsal):", nrow(comp_gene_pool), "\n")
  cat("  of which singletons (nSsal == 1):                ",
      sum(comp_gene_pool$nSsal == 1), "\n")
  cat("  of which two-copy ohnolog genes (nSsal == 2):    ",
      sum(comp_gene_pool$nSsal == 2), "\n")
  cat("Pair-only pool (nSsal == 2 ohnolog genes in SOM): ", nrow(gene_pool), "\n")
  cat("Both copies in SOM (used for co-clustering):      ", nrow(result), "\n")

  # --- ANALYSIS 1: Gene counts per cluster (AORe vs LORe), expanded pool ----
  cat("\n--- Gene counts per cluster (singleton-inclusive composition pool) ---\n")
  gene_counts <- comp_gene_pool %>%
    group_by(my_cluster, proxiPhylogeny) %>%
    summarise(n_genes = n(), .groups = "drop") %>%
    pivot_wider(names_from = proxiPhylogeny, values_from = n_genes, values_fill = 0) %>%
    rename(cluster = my_cluster) %>%
    left_join(cluster_labels %>% select(cluster, label), by = "cluster") %>%
    mutate(total      = AORe + LORe,
           pct_AORe   = round(100 * AORe / total, 1),
           pct_LORe   = round(100 * LORe / total, 1),
           el_ratio   = (AORe + 0.5) / (LORe + 0.5),
           # Wald CI on log-ratio with Haldane–Anscombe (+0.5) correction.
           se_log_r   = sqrt(1 / (AORe + 0.5) + 1 / (LORe + 0.5)),
           ratio_lo   = exp(log(el_ratio) - 1.96 * se_log_r),
           ratio_hi   = exp(log(el_ratio) + 1.96 * se_log_r)) %>%
    arrange(desc(total))

  print(as.data.frame(gene_counts))

  overall_aore <- sum(gene_counts$AORe)
  overall_lore <- sum(gene_counts$LORe)
  overall_ratio <- (overall_aore + 0.5) / (overall_lore + 0.5)
  cat(sprintf("Overall Early:Late counts: %d AORe / %d LORe (ratio %.2f)\n",
              overall_aore, overall_lore, overall_ratio))

  # Chi-squared: H0 = Early/Late ratio is the same in every cluster. Run on
  # both the full set (constitutive included — matches the figures the
  # manuscript cites) and the non-constitutive subset (kept for reference).
  if (nrow(gene_counts) > 1) {
    chi_full_mat <- gene_counts %>% select(AORe, LORe) %>% as.matrix()
    chi_full <- chisq.test(chi_full_mat)
    cat(sprintf("Chi-squared (all clusters):    chi2 = %.2f, df = %d, p = %s\n",
                chi_full$statistic, chi_full$parameter,
                format(chi_full$p.value, digits = 4)))
  }
  non_const <- gene_counts %>% filter(label != "Constitutive")
  if (nrow(non_const) > 1) {
    chi_mat <- non_const %>% select(AORe, LORe) %>% as.matrix()
    chi <- chisq.test(chi_mat)
    cat(sprintf("Chi-squared (non-constitutive): chi2 = %.2f, df = %d, p = %s\n",
                chi$statistic, chi$parameter,
                format(chi$p.value, digits = 4)))
  }

  # --- ANALYSIS 2: Gene-based co-clustering rate per cluster ----------------
  # For each cluster & class: of the ohnolog genes here, what fraction have
  # their partner in the same cluster? This is the metric plotted in the
  # per-cluster bar panels (03c) and cited in the manuscript E8 paragraph.
  cat("\n--- Gene-based co-clustering rate per cluster ---\n")
  coclust_rate <- result %>%
    group_by(my_cluster, proxiPhylogeny) %>%
    summarise(
      n_genes = n(),
      n_partner_here = sum(partner_here),
      .groups = "drop"
    ) %>%
    rename(cluster = my_cluster) %>%
    mutate(pct_coclustered = round(100 * n_partner_here / n_genes, 1)) %>%
    left_join(cluster_labels %>% select(cluster, label), by = "cluster") %>%
    arrange(cluster)

  print(as.data.frame(coclust_rate))

  # --- Per-cluster Fisher's (non-constitutive, both classes with n ≥ 20) ----
  # The n ≥ 20 floor is a pragmatic cutoff: below that, the 2×2 is sparse
  # enough that ORs blow up or collapse and the adjusted p isn't informative.
  cat("\n--- Per-cluster Fisher's Exact Tests (non-constitutive, n >= 20 per group) ---\n")
  clusters_to_test <- coclust_rate %>%
    filter(label != "Constitutive") %>%
    group_by(cluster) %>%
    filter(all(c("AORe", "LORe") %in% proxiPhylogeny)) %>%
    filter(all(n_genes >= 20)) %>%
    pull(cluster) %>% unique()

  fisher_results <- list()
  for (k in clusters_to_test) {
    a <- coclust_rate %>% filter(cluster == k, proxiPhylogeny == "AORe")
    l <- coclust_rate %>% filter(cluster == k, proxiPhylogeny == "LORe")
    # 2×2: rows AORe/LORe, cols PartnerHere/PartnerElsewhere.
    mat <- matrix(c(a$n_partner_here, a$n_genes - a$n_partner_here,
                     l$n_partner_here, l$n_genes - l$n_partner_here),
                  nrow = 2, byrow = TRUE,
                  dimnames = list(c("AORe", "LORe"), c("PartnerHere", "PartnerElsewhere")))
    ft <- fisher.test(mat)
    clust_label <- cluster_labels$label[cluster_labels$cluster == k]
    fisher_results[[as.character(k)]] <- tibble(
      cluster = k, label = clust_label,
      AORe_n = a$n_genes, AORe_pct = a$pct_coclustered,
      LORe_n = l$n_genes, LORe_pct = l$pct_coclustered,
      OR = round(ft$estimate, 3), p = ft$p.value
    )
  }

  if (length(fisher_results) > 0) {
    # BH across surviving clusters; star-annotate for at-a-glance reading.
    fisher_df <- bind_rows(fisher_results) %>%
      mutate(p_adj = p.adjust(p, method = "BH"),
             sig = case_when(p_adj < 0.001 ~ "***",
                             p_adj < 0.01 ~ "**",
                             p_adj < 0.05 ~ "*",
                             TRUE ~ "ns")) %>%
      arrange(p_adj)
    cat("\n")
    print(as.data.frame(fisher_df))
  } else {
    fisher_df <- tibble()
    cat("No clusters with sufficient data.\n")
  }

  list(gene_counts = gene_counts, coclust_rate = coclust_rate, fisher = fisher_df)
}

devmap_results <- run_gene_coclustering(
  gene_partner, ssal_classified, devmap_clusters, devmap_labels, "DevMap RNA-seq")

bodymap_results <- run_gene_coclustering(
  gene_partner, ssal_classified, bodymap_clusters, bodymap_labels, "BodyMap RNA-seq")

# =============================================================================
# Significance tables consumed by 03c overlays
#
# Three families of per-cluster tests, one stored per map:
#
# (1a) bar_tests — does a cluster's co-clustering rate deviate from the
#      map-wide mean for its class (AORe / LORe)?  One-sample proportion
#      test per (cluster × class), BH-adjusted within (class × map),
#      excluding Constitutive clusters from the mean calculation.
#
# (1b) ratio_tests — does the LORe fraction of a cluster's ohnolog pool
#      deviate from the genome-wide LORe fraction? Binomial test per
#      cluster (LORe vs AORe only), BH-adjusted across clusters within
#      the map.
#
# (2a) composition_tests — for the 4-category stacked panel, does each
#      (cluster × category) cell deviate from the genome-wide proportion
#      for that category? Binomial test; BH across cells within the map.
# =============================================================================

# (1a) — per-cluster co-clustering rate vs class mean.
# coclust_tbl already carries a `label` column from run_gene_coclustering
# (built from the joined cluster_labels), so we don't re-join — just
# filter on the existing column.
build_bar_tests <- function(coclust_tbl, cluster_labels) {
  non_const <- coclust_tbl %>%
    filter(label != "Constitutive")

  # Class-wide means over non-constitutive clusters only.
  class_means <- non_const %>%
    group_by(proxiPhylogeny) %>%
    summarise(
      total_n    = sum(n_genes),
      total_here = sum(n_partner_here),
      mean_rate  = total_here / total_n, .groups = "drop")

  tests <- coclust_tbl %>%
    left_join(class_means %>% select(proxiPhylogeny, mean_rate),
              by = "proxiPhylogeny") %>%
    rowwise() %>%
    mutate(pval = if (n_genes > 0) {
        suppressWarnings(prop.test(n_partner_here, n_genes, p = mean_rate)$p.value)
      } else NA_real_) %>%
    ungroup() %>%
    group_by(proxiPhylogeny) %>%
    mutate(padj = p.adjust(pval, method = "BH"),
           sig  = case_when(padj < 0.001 ~ "***",
                            padj < 0.01  ~ "**",
                            padj < 0.05  ~ "*",
                            TRUE         ~ ""),
           direction = ifelse(pct_coclustered > 100 * mean_rate, "up", "down")) %>%
    ungroup()

  list(class_means = class_means, tests = tests)
}

devmap_bar_tests  <- build_bar_tests(devmap_results$coclust_rate,  devmap_labels)
bodymap_bar_tests <- build_bar_tests(bodymap_results$coclust_rate, bodymap_labels)

# (1b) — per-cluster LORe:AORe ratio vs genome-wide fraction
build_ratio_tests <- function(gene_counts, cluster_labels) {
  all_AORe <- sum(gene_counts$AORe, na.rm = TRUE)
  all_LORe <- sum(gene_counts$LORe, na.rm = TRUE)
  p_lore <- all_LORe / (all_AORe + all_LORe)

  # gene_counts already has a `label` column (joined in run_gene_coclustering).
  tests <- gene_counts %>%
    mutate(total_aorelore = AORe + LORe) %>%
    rowwise() %>%
    mutate(pval = if (total_aorelore > 0) {
        suppressWarnings(binom.test(LORe, total_aorelore, p = p_lore)$p.value)
      } else NA_real_) %>%
    ungroup() %>%
    mutate(padj = p.adjust(pval, method = "BH"),
           sig  = case_when(padj < 0.001 ~ "***",
                            padj < 0.01  ~ "**",
                            padj < 0.05  ~ "*",
                            TRUE         ~ ""),
           enrichment = ifelse(
             (LORe / pmax(total_aorelore, 1)) > p_lore, "up", "down"))

  list(p_lore = p_lore, tests = tests)
}

devmap_ratio_tests  <- build_ratio_tests(devmap_results$gene_counts,  devmap_labels)
bodymap_ratio_tests <- build_ratio_tests(bodymap_results$gene_counts, bodymap_labels)

# (1c) — per-cluster Late/Early co-clustering RATE ratio with 95% CI.
# For each cluster, compute the ratio of LORe co-clustering rate to AORe
# co-clustering rate (with Haldane–Anscombe correction on counts), the
# Wald log-RR 95% CI, and a min-n gate. Significance is taken from the
# per-cluster Fisher's already in $fisher.
build_coclust_ratio_tests <- function(per_clust_list, cluster_labels, min_n = 10) {
  per_cluster <- per_clust_list$per_cluster

  wide <- per_cluster %>%
    pivot_wider(id_cols = c(cluster, label),
                names_from = proxiPhylogeny,
                values_from = c(n_genes, n_partner_here, pct_coclustered),
                values_fill = 0) %>%
    rename(n_AORe = n_genes_AORe,
           n_LORe = n_genes_LORe,
           here_AORe = n_partner_here_AORe,
           here_LORe = n_partner_here_LORe,
           rate_AORe = pct_coclustered_AORe,
           rate_LORe = pct_coclustered_LORe) %>%
    mutate(
      p_aore = (here_AORe + 0.5) / (n_AORe + 1),
      p_lore = (here_LORe + 0.5) / (n_LORe + 1),
      # Early-to-Late co-clustering rate ratio. Higher ratio = Early
      # ohnologs co-cluster at a rate closer to (or above) Late, i.e.
      # stronger regulatory constraint on Early ohnologs in this cluster.
      # Approaches 1 when Early matches Late. This direction recapitulates
      # the Fig. 3 constraint pattern: ratio rises with developmental
      # progression as Early ohnologs accumulate constraint.
      rr     = p_aore / p_lore,
      se_log_rr = sqrt((1 - p_aore) / (p_aore * (n_AORe + 1)) +
                       (1 - p_lore) / (p_lore * (n_LORe + 1))),
      rr_lo  = exp(log(rr) - 1.96 * se_log_rr),
      rr_hi  = exp(log(rr) + 1.96 * se_log_rr),
      pass_n = n_AORe >= min_n & n_LORe >= min_n
    ) %>%
    left_join(per_clust_list$fisher %>%
                select(cluster, p_adj_fisher = p_adj, sig_fisher = sig),
              by = "cluster")

  list(min_n = min_n, tests = wide)
}

devmap_per_clust <- list(per_cluster = devmap_results$coclust_rate,
                         fisher      = devmap_results$fisher)
bodymap_per_clust <- list(per_cluster = bodymap_results$coclust_rate,
                          fisher      = bodymap_results$fisher)

devmap_coclust_ratio_tests  <- build_coclust_ratio_tests(devmap_per_clust,  devmap_labels)
bodymap_coclust_ratio_tests <- build_coclust_ratio_tests(bodymap_per_clust, bodymap_labels)

# (2a) — per-cell composition vs map-wide category proportions.
#
# Baseline is the proportion of each category among the Ssal genes that
# ACTUALLY appear in this SOM — not the genome-wide proportion over the
# whole orthogroup table. This matters because the SOM input excludes
# low-expression genes, which over-represent the Singleton / Other
# categories. Testing against the full orthogroup distribution would
# therefore flag every cluster as "Other-depleted" — an artefact of the
# SOM input filter, not biology. Testing against the SOM-local pool
# instead asks "among the genes that made it into the SOM, does this
# cluster have a skewed composition?" — which is what we actually want.
build_composition_tests <- function(composition_df) {
  som_totals <- composition_df %>%
    group_by(category) %>% summarise(n_cat = sum(n), .groups = "drop")
  p_map <- setNames(som_totals$n_cat / sum(som_totals$n_cat),
                    as.character(som_totals$category))

  tests <- composition_df %>%
    rowwise() %>%
    mutate(p_cat = p_map[as.character(category)],
           pval = if (total > 0) {
             suppressWarnings(binom.test(n, total, p = p_cat)$p.value)
           } else NA_real_) %>%
    ungroup() %>%
    mutate(padj = p.adjust(pval, method = "BH"),
           sig  = case_when(padj < 0.001 ~ "***",
                            padj < 0.01  ~ "**",
                            padj < 0.05  ~ "*",
                            TRUE         ~ ""),
           direction = ifelse((n / pmax(total, 1)) > p_cat, "up", "down"))

  list(p_map = p_map, tests = tests)
}

devmap_composition_tests  <- build_composition_tests(devmap_composition)
bodymap_composition_tests <- build_composition_tests(bodymap_composition)

# -----------------------------------------------------------------------------
# Rename for compatibility with 03c's expected input variable names.
# `devmap_per_clust` / `bodymap_per_clust` wrap per_cluster + fisher into a
# single list; 03c indexes into $per_cluster and doesn't need $fisher.
# -----------------------------------------------------------------------------
devmap_gene_counts  <- devmap_results$gene_counts
bodymap_gene_counts <- bodymap_results$gene_counts

save(devmap_gene_counts, bodymap_gene_counts,
     devmap_per_clust, bodymap_per_clust,
     devmap_labels, bodymap_labels,
     ssal_category, ssal_classified,
     devmap_composition, bodymap_composition,
     devmap_bar_tests,            bodymap_bar_tests,
     devmap_ratio_tests,          bodymap_ratio_tests,
     devmap_composition_tests,    bodymap_composition_tests,
     devmap_coclust_ratio_tests,  bodymap_coclust_ratio_tests,
     file = "results/ohnolog/data/c1_extended_results.RData")

cat("\n\nResults saved to results/ohnolog/data/c1_extended_results.RData\n")
cat("Done.\n")
