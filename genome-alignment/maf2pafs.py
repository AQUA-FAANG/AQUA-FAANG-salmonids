import re
import sys
import argparse


def makeCigar(seqQ, seqT, strandT):
  # iterate over each character in the strings seqQ and seqT (they are same length)
  # if both are gaps, skip
  # if seqQ is a gap, add a deletion
  # if seqT is a gap, add an insertion
  # if neither are gaps, add a match
  cigar = ""  # initialize the cigar string
  for i in range(len(seqQ)):
    if seqQ[i] == "-" and seqT[i] == "-":
      continue
    elif seqQ[i] == "-":
      cigar += "D"
    elif seqT[i] == "-":
      cigar += "I"
    else:
      cigar += "M"

  # if the target is on the reverse strand, reverse the cigar string
  if strandT == "-":
    cigar = cigar[::-1]
  # count recurring characters in the cigar string and write count + character
  # e.g. "MMMD" -> "3M1D"
  cigar = re.sub(r"(\w)\1*", lambda m: str(len(m.group(0))) + m.group(1), cigar)
  return cigar

def generate_paf(alnQ,alnT):
  """
  Generate a line of paf from a pair of converted maf lines
  
  The input alnQ and alnT are lists with:
    [genome, chrom, start, end, strand, srcSize, text]
  
  Output PAF format is tab separated with the following columns:
  1 	string 	Query sequence name
  2 	int 	Query sequence length
  3 	int 	Query start (0-based; BED-like; closed)
  4 	int 	Query end (0-based; BED-like; open)
  5 	char 	Relative strand: "+" or "-"
  6 	string 	Target sequence name
  7 	int 	Target sequence length
  8 	int 	Target start on original strand (0-based)
  9 	int 	Target end on original strand (0-based)
  10 	int 	Number of residue matches
  11 	int 	Alignment block length
  12 	int 	Mapping quality (0-255; 255 for missing)
  13  Tag 	Extra tags, in this case cigar in the format "cg:Z:<CIGAR>"
   
  For now we ignore columns 2,7,10,11,12 and set them to ""
  """
  # get the cigar string
  cigar = makeCigar(alnQ[6], alnT[6], alnT[4])
  # get the start and end positions of the alignment
  relStrand =  "+" if (alnQ[4]==alnT[4]) else "-"
  # generate the paf line
  return( "\t".join([alnQ[1], "", str(alnQ[2]), str(alnQ[3]), relStrand,
                     alnT[1], "", str(alnT[2]), str(alnT[3]), "", "", "", f"cg:Z:{cigar}"]))

def read_maf_block(maf_file):
  # read maf file and yield each block
  # maf_file is a file object
  block = []
  for line in maf_file:
    if line.startswith("a"): # "a" marks a new block
      if block:
        yield block
      block = []
    if line.startswith("s"): # "s" marks a row in the block
      block.append(line)
  yield block
  return

def get_outfiles(out_prefix,use_reference, remove_subgenome, genomes):
  """
  generate output file names for each genome combination
  """
  output_files = {}
  filenames2file = {}
  for g1 in genomes:
    if use_reference and g1 != genomes[0]:
      continue
    for g2 in genomes:
      if not g1+g2 in output_files:
        if remove_subgenome:
          # remove "_[AB]" from the end of the genome name
          g1_nosub = re.sub(r"_([AB])$", "", g1)
          g2_nosub = re.sub(r"_([AB])$", "", g2)
          filename = f'{out_prefix}{g1_nosub}_vs_{g2_nosub}.paf'
        else:
          filename = f'{out_prefix}{g1}_vs_{g2}.paf'
        # check if we have already opened a file with this name
        if filename not in filenames2file:
          filenames2file[filename] = open(filename, 'w')
        output_files[g1+g2] = filenames2file[filename]
  return output_files

def parse_maf_block(maf_block, genomes):
  """
  MAF file format (from https://genome.ucsc.edu/FAQ/FAQformat.html#format5)
  The maf_block only contains the "s" (sequence) lines of a single alignment block which is
  tab separated with the following columns:

  src -- The name of one of the source sequences for the alignment. formated as 'genome.chromosome'
  start -- The start of the aligning region in the source sequence. This is a zero-based number. If 
    the strand field is "-" then this is the start relative to the reverse-complemented source sequence 
    (see http://genomewiki.ucsc.edu/index.php/Coordinate_Transforms).
  size -- The size of the aligning region in the source sequence. This number is equal to the number of 
    non-dash characters in the alignment text field below.
  strand -- Either "+" or "-". If "-", then the alignment is to the reverse-complemented source.
  srcSize -- The size of the entire source sequence, not just the parts involved in the alignment.
  text -- The nucleotides (or amino acids) in the alignment and any insertions (dashes) as well. 
  
  The output should be a list of lists with the following fields:
  genome -- The name of one of the source sequences for the alignment.
  chrom -- The name of the chromosome of the source sequence.
  start -- The start of the aligning region in the source sequence (on forward strand). (0-based)
  end -- The end of the aligning region in the source sequence (on forward strand). (0-based)
  strand -- unchanged from the input
  srcSize -- unchanged from the input
  text -- unchanged from the input
  """
  # maf_block is a list of lists of strings
  outList = []
  for line in maf_block:
    (src,start,size,strand,srcSize,text) = line.strip().split("\t")[1:]
    #
    # "src" is the genome name and chromsome name separated by a dot
    # Since both the genome names chromosome can contain ".", we need to match against the list
    # of genomes to separate the genome name from the chromosome name
    # use regex to match any of the genomes and extract the genome name and chromosome name
    reMatch = re.match(f"({'|'.join(genomes)}).(.*)", src)
    if not reMatch:
      continue
    genome, chrom = reMatch.groups()
    #
    # convert the start and size to integers
    start = int(start)
    size = int(size)
    srcSize = int(srcSize)
    #
    if strand == "-":
      # need to calculate the start on the forward strand (count from the end of the chromosome)
      start = srcSize - (start + size) 
    end = start + size
    outList.append([genome, chrom, start, end, strand, srcSize, text])
  return(outList)


def getArguments( argv = None ):
  # parse arguments from command line
  epilog_text = """
There are two ways to use this script, either extract all pairwise combinations of given genomes or
use reference to extract all alignments from reference to all given genomes.

Example 1 (all pairwise combinations): 
  maf2pafs.py alignment.maf --out_prefix out/alignment_ Ssal Omyk

This will generate paf files for the following genome combinations:
  out/alignment_Ssal_vs_Ssal.paf
  out/alignment_Ssal_vs_Omyk.paf
  out/alignment_Omyk_vs_Ssal.paf
  out/alignment_Omyk_vs_Ssal.paf

Example 2 (alignment with reference): 
  maf2pafs.py alignment_with_subgenomes.maf --use_reference --out_prefix out/alignment_ Ssal_A Ssal_B Omyk_A Omyk_B

This will generate paf files for the following genome combinations:
  out/alignment_Ssal_A_vs_Ssal_A.paf (If duplications present in the maf file)
  out/alignment_Ssal_A_vs_Ssal_B.paf
  out/alignment_Ssal_A_vs_Omyk_A.paf
  out/alignment_Ssal_A_vs_Omyk_B.paf

Example 3 (alignment with reference, remove subgenomes): 
  maf2pafs.py alignment_with_subgenomes.maf --use_reference --remove_subgenome --out_prefix out/alignment_ Ssal_A Ssal_B Omyk_A Omyk_B

This will generate paf files for the following genome combinations:
  out/alignment_Ssal_vs_Ssal.paf (target genome is Ssal_A/B, and _A/B suffix is removed from the genome name)
  out/alignment_Ssal_vs_Omyk.paf (target genome is Omyk_A/B, and _A/B suffix is removed from the genome name)


"""
  parser = argparse.ArgumentParser(description="Convert maf file to multiple paf files",
                                   formatter_class=argparse.RawTextHelpFormatter,
                                   epilog=epilog_text)
  parser.add_argument('--out_prefix', type=str, default="", help='prefix for output files')
  parser.add_argument("--use_reference", action="store_true", help="Use the first genome in the list as the reference genome")
  parser.add_argument("--remove_subgenome", action="store_true", help="remove the subgenome name (_A/_B) suffix in the chromosome name")
  parser.add_argument('maf_file_path', type=str, help='path to maf file. Use - to read from stdin')
  parser.add_argument('genomes', type=str, nargs='+', help='list of genomes separated by spaces (minimum 1). Note that the genome names must match the prefixes (followed by .) in the maf file.')
  return(parser.parse_args( argv )) 

def main():
  # get the command line arguments
  #maf_file_path, out_prefix, genomes = getArguments(["Ssa03_prkag2b_fixed.maf", "--out_prefix", "test_", "Ssal", "Eluc"])

  # out_prefix, use_reference, remove_subgenome, maf_file_path, genomes 
  args = getArguments()
  # generate output file names for each genome combination
  output_files = get_outfiles(args.out_prefix, args.use_reference, args.remove_subgenome, args.genomes)

  block_count=0
  
  # read maf file
  if args.maf_file_path == '-':
    maf_file = sys.stdin
  else:
    maf_file = open(args.maf_file_path, "r")
  
  # iterate over each block in the maf file
  for maf_block in read_maf_block(maf_file):
    block = parse_maf_block(maf_block, args.genomes)
    block_count += 1
    print(f'block {block_count}. length: {len(block)}, genomes: {" ".join([b[0]for b in block])}')
    # for each pair of sequences in the block
    for iQ in range(len(block)):
      if args.use_reference and iQ != 0: # reference query is always first sequence in the block
        continue
      for iT in range(len(block)):
        if iQ == iT:
          continue
        # get the pairwise alignment
        alnQ = block[iQ]; alnT = block[iT]; 
        # generate the paf line
        paf = generate_paf(alnQ,alnT)
        # get the genome pair
        genomepair = block[iQ][0]+block[iT][0]
        # write the paf line to the output file
        output_files[genomepair].write(paf+"\n")
  # close the files
  maf_file.close()
  for f in output_files.values():
    f.close()

if __name__ == "__main__":
  main()
