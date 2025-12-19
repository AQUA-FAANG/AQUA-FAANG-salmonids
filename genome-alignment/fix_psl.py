# Convert coordinates in PSL file
# PSL file contains alignments performed on parts (blocks) of the genome
# The coordinates in the PSL file are relative to the part of the genome that was aligned
# This script converts the coordinates to be relative to the entire genome
#
# The chromosome names in the PSL contain block coordinates, i.e. chrom:block_start-block_end
# block coordinates are 1-based
#
# psl files have 21 fields, the ones we are interrested in are:
# 10: qName - Query sequence name.
# 12: qStart - Alignment start position in query.
# 13: qEnd - Alignment end position in query.
# 14: tName - Target sequence name.
# 16: tStart - Alignment start position in target.
# 17: tEnd - Alignment end position in target.
#
# We need to convert:
# q/tName -- from chrom:block_start-block_end to chrom
# q/tStart -- block_start - 1 + start
# q/tEnd -- block_start - 1 + end

import re

def fix_coords(chrom_block, start, end):
  chrom,block = chrom_block.split(":")
  block_start = int(block.split("-")[0])
  if block_start==0: # dirty fix since the block start sometimes contains 0 
    block_start=1
  # convert the coordinates
  start = str(block_start - 1 + int(start))
  end = str(block_start - 1 + int(end))
  return(chrom,start,end)

def fix_psl_coordinates(psl_file):
  # psl_file is a file objects
  # iterate through the psl file and convert the coordinates
  for line in psl_file:
    fields = line.strip().split("\t")
    # halLiftover can optionally add the bed name as an extra field at the beginning
    if len(fields)==22:
      name = [fields.pop(0)]
    else:
      name = []
    (qName,qStart,qEnd) = fix_coords(*[ fields[9], fields[11], fields[12] ])
    (tName,tStart,tEnd) = fix_coords(*[ fields[13], fields[15], fields[16] ])
    # fields 11 and 15 contain the chomosome size, which is not known... replacec with NA for now
    print("\t".join(name + fields[:9] + [qName, "NA", qStart, qEnd, tName, "NA", tStart, tEnd] + fields[17:]))
  return

def main():
  # parse arguments: psl_file
  import argparse
  import sys
  parser = argparse.ArgumentParser(description="Convert coordinates in psl file")
  parser.add_argument("psl_file", help="psl file")
  args = parser.parse_args()
  psl_file = args.psl_file
  # open the input file
  with open(psl_file) as input_stream:
    fix_psl_coordinates(input_stream)
  return

if __name__ == "__main__":
  main()
