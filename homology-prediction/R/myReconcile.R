library(tibble)

# This function reconciles a gene tree with a given species tree.
# It identifies duplication nodes in the gene tree.
#
# Lowest Common Ancestor (LCA) reconciliation works by mapping nodes of a gene 
# tree onto nodes of a known species tree. Each internal node in the gene tree 
# is assigned to the lowest common ancestor (LCA) of the species associated 
# with its descendant nodes. (note: LCA = MRCA)
#
# This is equivalent to the method described in:
#   Zmasek CM, Eddy SR. (2001). A simple algorithm to infer gene duplication 
#   and speciation events on a gene tree.
# 
#
# When a gene-tree node maps to the same species node as one of its descendants,
# this indicates a gene duplication event. Otherwise, it's considered a 
# speciation event.
#
# Inputs:
# - tree: Gene tree to reconcile.
# - tipSpc: Vector indicating the species for each tip in the gene tree.
# - mrcaMat: MRCA (Most Recent Common Ancestor) matrix computed from the species tree.
# - spcTreeLabels: Labels for all nodes (tips and internal) in the species tree.
#
# Outputs:
# - A tibble with two columns:
#     S: corresponding species-tree node label for each node in gene tree
#     D: duplication status ("Y" for duplication, "N" for speciation, NA for tips)
LCAreconcile <- function(tree, tipSpc, mrcaMat, spcTreeLabels){
  S <- c(match(tipSpc,spcTreeLabels),rep(NA,Nnode(tree)))
  D <- c(rep(NA_character_,Ntip(tree)),rep("N",Nnode(tree)))
  # postordered to iterate from tip to root
  e <- tree$edge[postorder(tree),]
  # Assert that e[i,1]==e[i+1,1] where i is 1,3,5,7...nrow(e). This should be the case if tree is binary
  if(!all(e[c(T,F),1]==e[c(F,T),1])){
    stop("Tree is not binary")
  }
  for( i in seq(1,nrow(e),by=2)){
    parent = e[i,1] # == e[i+1,1]
    child1 = e[i,2]
    child2 = e[i+1,2]
    S[parent] <- mrcaMat[S[child1],S[child2]]
    if(S[parent] %in% S[c(child1,child2)]){
      # either of the child nodes have the same species/ancestral node as their MRCA
      # this is a duplication node
      D[parent] <- "Y"
    }
  }
  tibble(S=spcTreeLabels[S],D)
}

LCAreconcile_allTrees <- function(trees, spcTree){
  mrcaMat <- mrca(spcTree,full = T)
  spcTreeLabels <- c(spcTree$tip.label,spcTree$node.label)
  lapply(trees, function(tree){
    if(Ntip(tree)==1){return(tibble(S=sub(".*(....)$","\\1",tree$tip.label),D=NA_character_))}
    tipSpc <- sub(".*(....)$","\\1",tree$tip.label)
    myReconcileInner(tree,tipSpc,mrcaMat,spcTreeLabels)
  })
}
