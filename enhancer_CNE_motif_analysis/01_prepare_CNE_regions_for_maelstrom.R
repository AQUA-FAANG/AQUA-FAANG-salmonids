###############################################################
# CNE PREPROCESSING FOR GIMME MAELSTROM
#
# Purpose:
#   1. Merge nearby CNEs (<200 bp apart)
#   2. Prevent merging across biological classes
#   3. Expand short regions to a minimum width
#   4. Generate Maelstrom input files
#
# Analyses generated:
#   A. Evolutionary age
#   B. Conservation category
#   C. Age × Category interaction
###############################################################

# ============================================================
# LOAD LIBRARIES
# ============================================================

library(tidyverse)
library(purrr)

# ============================================================
# PARAMETERS
# ============================================================

MIN_WIDTH      <- 200
MERGE_DISTANCE <- 200

# ============================================================
# LOAD DATA
# ============================================================

all_CNE_salmon_coords <- read.table(
  "salmon_CNE_list_raw_input_gimme.txt",
  header = TRUE
)

# ============================================================
# GENERAL PROCESSING FUNCTION
# ============================================================
#
# Merges adjacent CNEs belonging to the same grouping
# variable and expands short regions to a minimum width.
#
# group_var can be:
#   age
#   category
#   age_category interaction
#
# ============================================================

process_cne_groups <- function(df,
                               group_var,
                               merge_distance = 200,
                               min_width = 200) {
  
  df <- df %>%
    mutate(
      chromosome = as.character(chromosome)
    ) %>%
    arrange(chromosome, pos_start)
  
  process_chr <- function(chr_df) {
    
    n <- nrow(chr_df)
    
    used <- rep(FALSE, n)
    
    output <- list()
    
    i <- 1
    
    while(i <= n) {
      
      if(used[i]) {
        i <- i + 1
        next
      }
      
      current_group <- chr_df[[group_var]][i]
      
      cluster <- i
      j <- i + 1
      
      # ------------------------------------------------------
      # Merge nearby regions belonging to same group
      # ------------------------------------------------------
      
      while(j <= n) {
        
        distance <- chr_df$pos_start[j] -
          chr_df$pos_end[cluster[length(cluster)]]
        
        same_group <- chr_df[[group_var]][j] ==
          current_group
        
        if(distance <= merge_distance & same_group) {
          
          cluster <- c(cluster, j)
          used[j] <- TRUE
          j <- j + 1
          
        } else {
          
          break
          
        }
      }
      
      start <- min(chr_df$pos_start[cluster])
      end   <- max(chr_df$pos_end[cluster])
      
      width <- end - start + 1
      
      # ------------------------------------------------------
      # Expand short regions to minimum width
      # ------------------------------------------------------
      
      if(width < min_width) {
        
        midpoint <- (start + end) / 2
        
        start <- round(midpoint - min_width / 2)
        end   <- round(midpoint + min_width / 2 - 1)
        
        start <- max(1, start)
      }
      
      output[[length(output) + 1]] <- data.frame(
        chromosome = chr_df$chromosome[i],
        start      = start,
        end        = end,
        group      = current_group,
        merged_n   = length(cluster)
      )
      
      i <- i + 1
    }
    
    bind_rows(output)
  }
  
  result <- df %>%
    group_split(chromosome) %>%
    map_dfr(process_chr)
  
  result %>%
    mutate(
      width = end - start + 1
    )
}

# ============================================================
# QC FUNCTION
# ============================================================

qc_summary <- function(result, label) {
  
  cat("\n")
  cat("========================================\n")
  cat(label, "\n")
  cat("========================================\n")
  
  cat("Regions:", nrow(result), "\n")
  
  cat(
    "Merged clusters:",
    sum(result$merged_n > 1),
    "\n"
  )
  
  cat(
    "Min width:",
    min(result$width),
    "\n"
  )
  
  cat(
    "Median width:",
    median(result$width),
    "\n"
  )
  
  cat(
    "Max width:",
    max(result$width),
    "\n"
  )
  
  cat(
    "Regions <200 bp:",
    sum(result$width < 200),
    "\n"
  )
  
  cat("========================================\n")
}

# ============================================================
# MAELSTROM EXPORT FUNCTION
# ============================================================

write_maelstrom_input <- function(result,
                                  outfile) {
  
  maelstrom_input <- result %>%
    transmute(
      region = paste0(
        chromosome,
        ":",
        start,
        "-",
        end
      ),
      cluster = group
    )
  
  write.table(
    maelstrom_input,
    outfile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

###############################################################
# ANALYSIS 1
# EVOLUTIONARY AGE
###############################################################

age_result <- process_cne_groups(
  df = all_CNE_salmon_coords,
  group_var = "age",
  merge_distance = MERGE_DISTANCE,
  min_width = MIN_WIDTH
)

qc_summary(
  age_result,
  "Evolutionary Age Analysis"
)

write_maelstrom_input(
  age_result,
  "maelstrom_input_age_200bp.txt"
)

###############################################################
# ANALYSIS 2
# CONSERVATION CATEGORY
###############################################################

category_result <- process_cne_groups(
  df = all_CNE_salmon_coords,
  group_var = "category",
  merge_distance = MERGE_DISTANCE,
  min_width = MIN_WIDTH
)

qc_summary(
  category_result,
  "Conservation Category Analysis"
)

write_maelstrom_input(
  category_result,
  "maelstrom_cons_category_200bp_input.txt"
)

###############################################################
# ANALYSIS 3
# AGE × CATEGORY INTERACTION
###############################################################

interaction_df <- all_CNE_salmon_coords %>%
  mutate(
    age_category = paste(
      category,
      age,
      sep = "_"
    )
  )

interaction_result <- process_cne_groups(
  df = interaction_df,
  group_var = "age_category",
  merge_distance = MERGE_DISTANCE,
  min_width = MIN_WIDTH
)

qc_summary(
  interaction_result,
  "Age × Category Interaction Analysis"
)

write_maelstrom_input(
  interaction_result,
  "maelstrom_age_category_interaction_200bp_input.txt"
)

