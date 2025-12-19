# Convert coordinates in wiggle file
#
# Wiggle files have a header like this:
# fixedStep chrom=13:0-114417674 start=56512361 step=1
#


import re

def fix_coords(chrom_block, start):
  chrom,block = chrom_block.split(":")
  block_start = int(block.split("-")[0])
  if block_start==0: # dirty fix since the block start sometimes contains 0 
    block_start=1
  # convert the coordinates
  start = str(block_start - 1 + int(start))
  return(chrom,start)

def fix_wiggle_coordinates(wiggle_file):
  # wiggle_file is a file objects
  # iterate through the wiggle file and convert the coordinates
  for line in wiggle_file:
    # if header
    if "chrom=" in line:
      chrom_block, start = re.search(r"chrom=([^\s]+) start=(\d+)", line).groups()
      chrom, new_start = fix_coords(chrom_block, start)
      print(re.sub(r"chrom=[^\s]+", f"chrom={chrom}", re.sub(r"start=\d+", f"start={new_start}", line)), end="")
    else:
      print(line, end="")
  return

def main():
  # parse arguments: wiggle_file
  import argparse
  import sys
  parser = argparse.ArgumentParser(description="Convert coordinates in wiggle file")
  parser.add_argument("wiggle_file", help="Input wiggle file. Use 'stdin' to read from stdin")
  args = parser.parse_args()
  wiggle_file = args.wiggle_file
  # open the input file or read from stdin
  if wiggle_file == "stdin":
      input_stream = sys.stdin
  else:
      input_stream = open(wiggle_file)

  try:
      fix_wiggle_coordinates(input_stream)
  finally:
      if input_stream is not sys.stdin:
          input_stream.close()  
  return

if __name__ == "__main__":
  main()
