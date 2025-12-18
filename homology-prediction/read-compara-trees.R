library(tidyverse)
library(ape)

# This script extracts a subset of species and converts gene trees from Ensembl 
# into tidytree::treedata format and generates OGtbl.
#
# Input:
# - bulk downloaded gene (protein) trees from Ensembl compara in .nhx.emf format
# - Species tree with internal node names (tips must match species in the gene trees)
# - List of species subset to extract
#
# Output:
# - spcTree.nw - pruned species tree
# - allTrees.RDS - list of tidytree::treedata
# - OGtbl.RDS - Ortholog groups for each internal node in the pruned species tree
#
#
# Each tree (or OG) is given a number corresponding to the order in the input file and is given
# a specified unique prefix.
# The tip labels of of the output trees contain the gene ID and species as a 4 letter suffix.
# The internal node labels indicates the MRCA of the clade (after pruning)
# The "D" NHX tag, which marks duplication nodes are inferred after very simple/strict criteria
# that assumes that the phylogeny is completely correct. This overwrites the tags in the original tree
# which were slightly different (probably more forgiving). The reason for the strictness is that it is 
# used in the OGtbl generation which needs to be consistent with the phylogeny. 
# (Note: not sure if the original duplication annotation would be correct after pruning anyway)



# load species tree
spcTree <- read.tree("data/download/ensembl106_PostSCORPiOs/SpeciesTree.nwk")

# load subset of species
spcs <- readLines("selected-species.txt")
spcs <- spcs[!grepl("^#",spcs)]
# prune and generate 4 letter species codes
spcTree <- keep.tip(spcTree,spcs)
spcTree$tip.label <- sub("(.)[^.]*\\.(...).*","\\1\\2",spcTree$tip.label)
# assert that all species codes are unique
if(any(duplicated(spcTree$tip.label))) stop("Could not produce unique 4 letter abbreviations of species")
# save to file
dir.create("data/Ens106s")
write.tree(spcTree,"data/Ens106s/spcTree.nw")

# plot(spcTree,use.edge.length = F,show.node.label = T)

# load some helper functions
source("R/makeOGtbl.R")
source("R/myReconcile.R")
source("R/myReadNHX.R")
source("R/read_emf_nhx_subtree.R")


TreesFile <- "data/download/Compara.106.protein_default.nhx.emf.gz"

# Prefix for naming the trees/OG
OGprefix <- "Ens106s" # Ensembl 106 salmonid species subset


# species names convertion table
species2spc <- tibble( spc=sub("(.)[^.]*\\.(...).*","\\1\\2",spcs), # 4-letter short form
                       species = sub("\\.","_",tolower(spcs)))      # Compara format






#### Load sub-trees from nhx emf file ####


rawTrees <- read_nhx_emf_subtree( nhx_emf_filename = TreesFile, species2spc)

# name the trees and remove missing trees
names(rawTrees) <- sprintf("%s%05i",OGprefix, seq_along(rawTrees))
rawTrees <- rawTrees[!sapply(rawTrees,is.null)]

#### Reconcile ####
cat( "reconcile gene-trees...\n")

# Adds species (leaf nodes) and MRCA (internal nodes) as S column of the metadata.
# Defines duplicates (D column) as any node where the MRCA is same as MRCA of one of the child nodes.
#
# Note that compara trees already have inferred duplications, but these are not as strict defined  and allows
# minor discrepancies in the phylogeny.

recon <- LCAreconcile_allTrees(lapply(rawTrees,function(x){x@phylo}),spcTree)
allTrees <- mapply( rawTree=rawTrees, rec=recon, FUN=function(rawTree,rec){
  rawTree@data[,colnames(rec)] <- rec
  rawTree
})

saveRDS(allTrees,"data/Ens106s/allTrees.RDS")



#### Generate OGtbl ####
cat( "Generate OGtbl...\n")

OGtbl <- makeOGtbl(trees = allTrees,spcTree)

saveRDS(OGtbl,"data/Ens106s/OGtbl.RDS")

cat( "Done.\n")

