#
# After merging the bigwig files all the 0's were removed
# This script will fill all gaps with 0's.
#

library(rtracklayer)
library(GenomicRanges)
library(IRanges)
library(optparse)

# Define command line options
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input BigWig file path", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output BigWig file path", metavar="character")
)

# Parse command line options
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Check if input arguments are provided
if (is.null(opt$input) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Both input and output files must be supplied", call.=FALSE)
}

# Assign arguments to variables
input_bigwig_file <- opt$input
output_bigwig_file <- opt$output

# Print the provided arguments for confirmation
cat("Input BigWig file: ", input_bigwig_file, "\n")
cat("Output BigWig file: ", output_bigwig_file, "\n")

cat("Import bigwig...\n")
bw_data <- import(input_bigwig_file)

cat("Get gaps\n")
gaps_data <- gaps(bw_data) 
gaps_data <-  gaps_data[strand(gaps_data)=="*"] # ignore strand

cat("Assign score 0 to gaps\n")
mcols(gaps_data)$score <- 0

cat("Combine original ranges with gaps\n")
combined_data <- c(bw_data, gaps_data)

cat("Sort combined data by sequence names and ranges\n")
combined_data <- sort(combined_data)

cat("Export to BigWig file...\n")
export.bw(combined_data, output_bigwig_file)

cat("Done.\n")

