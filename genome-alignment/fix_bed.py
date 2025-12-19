# Convert coordinates in BED file
# BED file contains alignments performed on parts (blocks) of the genome
# The coordinates in the BED file are relative to the part of the genome that was aligned
# This script converts the coordinates to be relative to the entire genome
#
# The chromosome names in the BED contain block coordinates, i.e. chrom:block_start-block_end
# block coordinates are 1-based
#
# bed files have the following fields:
# chrom -- chrom:block_start-block_end
# start -- The start of the region in the block coordinates. This is a zero-based number.
# end -- The start of the region in the block coordinates. This is a zero-based number.
# rest of the fields -- anything
#
# We need to convert:
# chrom -- from chrom:block_start-block_end to chrom
# start -- block_start - 1 + start
# end -- block_start - 1 + end

import re

def fix_bed_coordinates(bed_file):
  # bed_file is a file objects
  # iterate through the bed file and convert the coordinates
  for line in bed_file:
    fields = line.strip().split("\t")
    chrom_block, start, end = fields[0:3]
    chrom,block = chrom_block.split(":")
    block_start = int(block.split("-")[0])
    if block_start==0: # dirty fix since the block start sometimes contains 0 
      block_start=1
    # convert the coordinates
    start = str(block_start - 1 + int(start))
    end = str(block_start - 1 + int(end))
    print("\t".join([chrom, start, end] + fields[3:]))
  return

def main():
  # parse arguments: bed_file
  import argparse
  import sys
  parser = argparse.ArgumentParser(description="Convert coordinates in bed file")
  parser.add_argument("bed_file", help="bed file")
  args = parser.parse_args()
  bed_file = args.bed_file
  # open the input file
  with open(bed_file) as input_stream:
    fix_bed_coordinates(input_stream)
  return

if __name__ == "__main__":
  main()