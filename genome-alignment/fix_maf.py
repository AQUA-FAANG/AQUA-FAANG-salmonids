# Convert coordinates in maf file
# Maf file contains alignments performed on parts (blocks) of the genome
# The coordinates in the maf file are relative to the part of the genome that was aligned
# This script converts the coordinates to be relative to the entire genome
#
# The chromosome names in the MAF contain block coordinates, i.e. chrom:block_start-block_end
#
# Also need chromosome sizes for each genome as input
# 
# MAF has the following tab separated columns (of the "s" lines):
# src -- The name of one of the source sequences for the alignment. formated as 'genome.chrom:block_start-block_end'
# start -- The start of the aligning region in the source sequence. This is a zero-based number. If 
#    the strand field is "-" then this is the start relative to the reverse-complemented source sequence 
# size -- The size of the aligning region in the source sequence. This number is equal to the number of 
#    non-dash characters in the alignment text field below.
# strand -- Either "+" or "-". If "-", then the alignment is to the reverse-complemented source.
# srcSize -- The size of the entire source sequence, not just the parts involved in the alignment.
# text -- The nucleotides in the alignment and any insertions (dashes) as well.
#
# We need to convert:
# src -- from genome.chrom:start-end to genome.chrom
# start -- if strand is "+" then block_start - 1 + start,
#       if strand is "-" then start + (chrom_size - block_end)
# srcSize -- set to chrom_size.

import re

def fix_maf_coordinates(maf_file, chrSizes, keep_subgenome=False):
  # maf_file is a file objects
  # chrSizes is a dictionary of dictionaries with the chromosome sizes for each genome
  # iterate through the maf file and convert the coordinates
  genomes = list(chrSizes.keys())
  for line in maf_file:
    if line.startswith("s"):
      # "s" line
      # split the line into fields
      fields = line.strip().split("\t")
      # extract the fields we need
      src, start, size, strand, srcSize, text = fields[1:]
      # extract the genome and chromosome from the src field
      # Since both the genome names chromosome can contain ".", we need to match against the list
      # of genomes to separate the genome name from the chromosome name
      # use regex to match any of the genomes and extract the genome name and chromosome name
      reMatch = re.match(f"({'|'.join(genomes)})(_[AB])?\\.([^:]+):([0-9]+)-([0-9]+)", src)
      if not reMatch:
        continue
      genome, subgenome, chrom, block_start, block_end  = reMatch.groups()
      # convert to integers
      block_start = int(block_start)
      block_end = int(block_end)
      if block_start==0:
        block_start = 1 # samtools faidx treats 0 as 1
      start = int(start)
      srcSize = int(srcSize)
      # convert the coordinates
      if strand == "+":
        start = block_start - 1 + start
      elif strand == "-":
        start = start + (chrSizes[genome][chrom] - block_end)
      # update the src field
      if keep_subgenome and (subgenome is not None):
        src = f"{genome}{subgenome}.{chrom}"
      else:
        src = f"{genome}.{chrom}"
      # update the start field
      start = str(start)
      # update the srcSize field
      srcSize = str(chrSizes[genome][chrom])
      # print the new line
      print("\t".join(fields[:1] + [src, start, size, strand, srcSize, text]))
    else:
      # not a "s" line, just print the line
      print(line, end="")
  return

def main():
  # parse arguments: maf_file chrSizes_file [--keep_subgenome]
  import argparse
  import sys
  parser = argparse.ArgumentParser(description="Convert coordinates in maf file")
  parser.add_argument("maf_file", help="maf file")
  parser.add_argument("chrSizes_file", help="chromosome sizes file")
  parser.add_argument("--keep_subgenome", action="store_true", help="keep the subgenome name (_A/_B) suffix in the chromosome name")
  args = parser.parse_args()
  maf_file = args.maf_file
  chrSizes_file = args.chrSizes_file
  keep_subgenome = args.keep_subgenome
  # read the chromosome sizes
  chrSizes = {}
  with open(chrSizes_file) as f:
    for line in f:
      genome, chrom, size = line.strip().split("\t")
      if genome not in chrSizes:
        chrSizes[genome] = {}
      chrSizes[genome][chrom] = int(size)
  # open the files
  with open(maf_file) as maf_file:
    fix_maf_coordinates(maf_file, chrSizes,keep_subgenome=keep_subgenome)
  return

if __name__ == "__main__":
  main()