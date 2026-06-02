#!/usr/bin/env Rscript
# =============================================================================
# 03_ohnolog_som_coclustering.R
# -----------------------------------------------------------------------------
# Purpose
#   Quantify whether AORe (Early rediploidization) vs LORe (Late rediploidization)
#   Salmo salar ohnolog pairs co-cluster in SOM expression clusters at different
#   rates, in both the DevMap (developmental, 16-unit) and BodyMap (adult tissue,
#   36-unit) SOMs. Addresses Reviewer 1 Comment 6 of the NEE revision.
#
# Analysis chain (run in this order)
#   03   → c1_coclustering_results.RData   (this script — pair-level tables)
#   03b  → c1_extended_results.RData       (per-cluster, gene-level metrics)
#   03c  → per-cluster bar-chart panels
#   03d  → Sankey/alluvial panels
#   03e  → UpSet panels
#   04   → stitches panels into supplementary composite figures
#
# Inputs (all under data/ in this repository)
#   data/original_soms/soms.RData          — sd_som, sb_som (+ td_som, tb_som)
#                                            kohonen SOM fits; we use
#                                            $som$unit.classif (gene → unit).
#   data/original_norm_files/norm.RData    — sd_norm, sb_norm matrices; row
#                                            names are gene IDs used here just
#                                            to recover the gene order that
#                                            matches `unit.classif`.
#   data/orthogroups/SalmonTroutOrthologs.tsv
#                                          — sub-orthogroup table with columns
#                                            spc, nSsal, nOmyk, proxiPhylogeny,
#                                            subOG, geneID.
#
# Output
#   results/ohnolog/data/c1_coclustering_results.RData
#     devmap_result, bodymap_result — each a list:
#       $summary      tibble: proxiPhylogeny × same/different_cluster totals
#       $fisher       htest:  AORe vs LORe 2×2 Fisher's Exact
#       $contingency  2×2 matrix: rows AORe/LORe, cols Same/Different
#       $detail       tibble: one row per pair with cluster1, cluster2,
#                             same_cluster, proxiPhylogeny — consumed by 03d
#                             (Sankey), 03e (UpSet), and 04 downstream.
#     devmap_per_cluster — per-SOM-cluster count of co-clustered pairs by class
#                          (sanity descriptive; figures are built from 03b).
#
# Ohnolog filter (why each gate matters)
#   spc == "Ssal"                   keep only Salmon rows of the orthogroup table
#   nSsal == 2                      require exactly two Salmon copies (the
#                                   ohnolog pair we are testing). Rainbow
#                                   Trout copy number (nOmyk) is intentionally
#                                   NOT constrained: rediploidization occurs
#                                   independently in each species post-Ss4R,
#                                   so the AORe/LORe timing of a Salmon pair
#                                   does not depend on how many Trout copies
#                                   survived. (Main-text Fig. 3 uses a 2:2
#                                   gate for a different question; replicating
#                                   it here would arbitrarily discard LORe
#                                   pairs.)
#   proxiPhylogeny %in% c(AORe,LORe) drop pairs with ambiguous timing.
#   subOG count == 2                removes any subOGs where filters above left
#                                   a single Salmon gene — we need a pair.
#   consistent proxiPhylogeny       both members of a pair must share the same
#                                   rediploidization label (a handful of rows
#                                   in the source table are inconsistent —
#                                   dropped here with a warning).
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# This script is sourced by 04_ohnolog_coclustering.Rmd with the working
# directory at transcriptome-analysis/, AFTER R/setup.R has run — so the shared
# loaders, `paths`, and read_table_anywhere() are all in scope. Generated
# outputs are written under results/ohnolog/.

# -----------------------------------------------------------------------------
# Load the deposited SOM unit assignments and normalised expression matrices
# via the shared loaders (99_helpers.R). These inject sd_som/td_som/sb_som/
# tb_som and sd_norm/td_norm/sb_norm/tb_norm into this environment — the exact
# same objects 02_soms.Rmd uses — so the ohnolog analysis can never drift from
# the rest of the pipeline. Here the matrices are used only for their row names:
# those row names are in the same order as `som$unit.classif`, the integer
# vector of per-gene SOM unit assignments.
# -----------------------------------------------------------------------------
cat("Loading SOM assignments...\n")
load_deposited_soms()   # injects sd_som, td_som, sb_som, tb_som
load_deposited_norm()   # injects sd_norm, td_norm, sb_norm, tb_norm

# Orthogroup table: always fetched from Salmobase (downloaded once, cached
# under cache/). Same source as every other notebook in this tree.
cat("Loading orthogroups...\n")
og <- read_table_anywhere(paths$active$orthogroups)

# -----------------------------------------------------------------------------
# Gene → SOM cluster lookup tables. `unit.classif[i]` is the SOM unit (1..K)
# that gene `rownames(<norm>)[i]` was mapped to. We zip them into tibbles for
# clean joins downstream.
# -----------------------------------------------------------------------------
devmap_clusters <- tibble(
  geneID = rownames(sd_norm),
  devmap_cluster = sd_som$som$unit.classif
)

bodymap_clusters <- tibble(
  geneID = rownames(sb_norm),
  bodymap_cluster = sb_som$som$unit.classif
)

# -----------------------------------------------------------------------------
# Extract 2:2 AORe/LORe Salmon pairs. See the "Ohnolog filter" note in the
# file header for the rationale behind each filter clause.
# -----------------------------------------------------------------------------
ssal_pairs <- og %>%
  filter(spc == "Ssal",
         nSsal == 2,
         proxiPhylogeny %in% c("AORe", "LORe")) %>%
  select(subOG, geneID, proxiPhylogeny)

# Sanity filter: a subOG that still has both copies after the species/n-copy
# filters above will have exactly n == 2 rows here. Anything else is dropped.
pair_counts <- ssal_pairs %>%
  count(subOG) %>%
  filter(n == 2)

ssal_pairs <- ssal_pairs %>%
  filter(subOG %in% pair_counts$subOG)

# A handful of subOGs in the source table carry different proxiPhylogeny
# labels across their two rows (upstream annotation inconsistency). Drop
# them rather than silently picking a winner.
inconsistent <- ssal_pairs %>%
  group_by(subOG) %>%
  summarise(n_phyl = n_distinct(proxiPhylogeny)) %>%
  filter(n_phyl > 1)
if (nrow(inconsistent) > 0) {
  cat("WARNING:", nrow(inconsistent), "pairs with inconsistent proxiPhylogeny - removing\n")
  ssal_pairs <- ssal_pairs %>% filter(!subOG %in% inconsistent$subOG)
}

# -----------------------------------------------------------------------------
# Reshape long → wide: two rows per subOG → one row per subOG with gene1/gene2.
# The `copy = row_number()` per-subOG counter is used as the pivot key. Which
# copy ends up in gene1 vs gene2 is arbitrary (row order within the table),
# but that's fine — co-clustering is a symmetric relation.
# -----------------------------------------------------------------------------
ssal_wide <- ssal_pairs %>%
  group_by(subOG) %>%
  mutate(copy = row_number()) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(subOG, proxiPhylogeny),
              names_from = copy,
              values_from = geneID,
              names_prefix = "gene")

cat("\nOhnolog pairs available:\n")
print(table(ssal_wide$proxiPhylogeny))

# =============================================================================
# run_coclustering — join ohnolog pairs to SOM cluster assignments and test
# whether AORe vs LORe co-clustering rates differ.
#
# Args:
#   pairs_df     tibble of pairs with columns subOG, proxiPhylogeny, gene1, gene2.
#   cluster_df   tibble mapping geneID → <cluster_col>.
#   cluster_col  name of the cluster column in cluster_df (string).
#   label        log label for this run ("DevMap RNA-seq" / "BodyMap RNA-seq").
#
# Returns: list(summary, fisher, contingency, detail). `detail` is the input
# `pairs_df` joined with cluster1/cluster2 plus a `same_cluster` indicator,
# and is the primary object consumed by the downstream figure scripts.
# =============================================================================
run_coclustering <- function(pairs_df, cluster_df, cluster_col, label) {
  cat("\n=== ", label, " ===\n")

  # Left-join twice: once on gene1 (→ cluster1), once on gene2 (→ cluster2).
  # The `rename(cluster1 = !!cluster_col)` step renames the joined-in column
  # so the two joins don't collide (both arrive named `<cluster_col>`).
  # `filter(!is.na(...))` drops pairs where either copy was absent from the
  # SOM's expression matrix (e.g. filtered out upstream for low expression).
  result <- pairs_df %>%
    left_join(cluster_df, by = c("gene1" = "geneID")) %>%
    rename(cluster1 = !!cluster_col) %>%
    left_join(cluster_df, by = c("gene2" = "geneID")) %>%
    rename(cluster2 = !!cluster_col) %>%
    filter(!is.na(cluster1), !is.na(cluster2)) %>%
    mutate(same_cluster = cluster1 == cluster2)

  cat("Pairs with both genes in SOM:", nrow(result), "\n")

  # Aggregate into 2×2 counts (AORe/LORe × Same/Different). The pivot_wider
  # creates columns named "FALSE" / "TRUE" which we then rename to descriptive
  # labels. `values_fill = 0` handles the edge case where one class has no
  # Same-cluster pairs at all.
  summary_tbl <- result %>%
    group_by(proxiPhylogeny, same_cluster) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = same_cluster,
                values_from = n,
                values_fill = 0) %>%
    rename(different_cluster = `FALSE`, same_cluster = `TRUE`) %>%
    mutate(total = same_cluster + different_cluster,
           pct_same = round(100 * same_cluster / total, 1))

  cat("\nCo-clustering summary:\n")
  print(as.data.frame(summary_tbl))

  # Fisher's Exact Test on the AORe-vs-LORe × Same-vs-Different 2×2. Exact
  # rather than chi-squared because the Same column is typically small enough
  # that asymptotic assumptions are borderline.
  aore <- summary_tbl %>% filter(proxiPhylogeny == "AORe")
  lore <- summary_tbl %>% filter(proxiPhylogeny == "LORe")

  fisher_mat <- matrix(c(
    aore$same_cluster, aore$different_cluster,
    lore$same_cluster, lore$different_cluster
  ), nrow = 2, byrow = TRUE,
  dimnames = list(c("AORe", "LORe"), c("Same", "Different")))

  cat("\nContingency table:\n")
  print(fisher_mat)

  ft <- fisher.test(fisher_mat)
  cat("\nFisher's Exact Test:\n")
  cat("  Odds ratio:", round(ft$estimate, 3), "\n")
  cat("  95% CI:", round(ft$conf.int[1], 3), "-", round(ft$conf.int[2], 3), "\n")
  cat("  p-value:", format(ft$p.value, digits = 4), "\n")

  return(list(summary = summary_tbl, fisher = ft, contingency = fisher_mat, detail = result))
}

# Run the two maps. The `detail` frames of both are the main downstream input
# to 03d/03e (Sankey/UpSet) and the source of the pair counts quoted in E8.
devmap_result <- run_coclustering(
  ssal_wide, devmap_clusters, "devmap_cluster", "DevMap RNA-seq")

bodymap_result <- run_coclustering(
  ssal_wide, bodymap_clusters, "bodymap_cluster", "BodyMap RNA-seq")

# -----------------------------------------------------------------------------
# Supplementary descriptive table: for DevMap, how many co-clustered pairs
# land in each SOM cluster, split by rediploidization class. Not used in any
# downstream figure (03b rebuilds this at the gene level, which is the metric
# the figures actually report), but kept for quick inspection.
# -----------------------------------------------------------------------------
cat("\n=== DevMap: Per-cluster breakdown ===\n")
devmap_per_cluster <- devmap_result$detail %>%
  filter(same_cluster) %>%
  group_by(proxiPhylogeny, cluster1) %>%
  summarise(n_coclustered = n(), .groups = "drop") %>%
  arrange(proxiPhylogeny, desc(n_coclustered))

cat("\nAORe co-clustered pairs per SOM cluster:\n")
print(as.data.frame(devmap_per_cluster %>% filter(proxiPhylogeny == "AORe")))
cat("\nLORe co-clustered pairs per SOM cluster:\n")
print(as.data.frame(devmap_per_cluster %>% filter(proxiPhylogeny == "LORe")))

# -----------------------------------------------------------------------------
# Save. Writes to results/ohnolog/data/ so the downstream figure scripts
# (03c/d/e) can load it from a known location.
# -----------------------------------------------------------------------------
output_path <- "results/ohnolog/data"
save(devmap_result, bodymap_result, devmap_per_cluster,
     file = file.path(output_path, "c1_coclustering_results.RData"))
cat("\nResults saved to", file.path(output_path, "c1_coclustering_results.RData"), "\n")
cat("\nDone.\n")
