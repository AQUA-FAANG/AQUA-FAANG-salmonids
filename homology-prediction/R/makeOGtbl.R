#
# Functions to generate OGtbl from gene-trees
#

# Ortholog groups (OGs) are identified for each internal node in the species tree.
# An ortholog group is a clade in the gene tree with following properties:
# 1. Clade only contains genes from species descending from the internal node in the species tree
# 2. Root node in the clade must not be a "duplication" node, where "duplication" node is defined
#    as any node where the MRCA is same as MRCA of one of the child nodes.

#### Helper functions ####

getRootNode <- function(phylo){
  # get the node that has no edges pointing to it
  setdiff(phylo$edge[ ,1], phylo$edge[ ,2])
}


#' Get nodes between an ancestral node and a descendant node in a phylogeny
#'
#' @param anc_node Node number of the ancestral node
#' @param desc_node Node number of the descendant node
#' @param tree Phylogeny (ape::phylo object)
#'
#' @return vector of node numbers including desc_node, anc_node and all in between
#' @export
#'
#' @examples
getNodesBetween <- function(anc_node,desc_node,tree){
  retVal <- desc_node
  if(anc_node == desc_node) return(retVal)
  parent <- tree$edge[tree$edge[ ,2]==desc_node, 1]
  if(length(parent)==0) stop("desc_node is root!")
  while(parent != anc_node){
    retVal <- c(parent,retVal)
    parent <- tree$edge[tree$edge[ ,2]==parent, 1]
    if(length(parent)==0) stop("Root node reached! desc_node is not descendant of anc_node!")
  }
  return(c(anc_node,retVal))
}

# The OG matrix contains one row for each gene (i.e. tip node) in gene tree and one 
# column for each internal node in the species tree. 
# walk the gene-tree and species-tree and record the path 
#
# If not dup:
# * assign current and all earlier to new node
# If dup:
# * assign parent and all earlier to new (dup) (i.e. remove the desc_node)
#
getOGmatrix <- function(path, node, spcNode, tree, S, D, spcTree){
  isDupNode <- ifelse(is.na(D[node]), F, D[node]=="Y")
  # record the path
  
  # get current corresponding spcNode
  spcNode_new <- match(S[node],c(spcTree$tip.label,spcTree$node.label))
  # and all earlier spcNodes
  spcNodes <- getNodesBetween(anc_node = spcNode, desc_node = spcNode_new, spcTree )
  if( isDupNode ){
    spcNodes <- spcNodes[spcNodes!=spcNode_new] # remove the desc_node
  }
  spcNodes <- spcNodes[spcNodes>Ntip(spcTree)] # only include internal nodes
  spcNodesS <- c(spcTree$tip.label,spcTree$node.label)[spcNodes] # get name of node
  spcNodesS <- spcNodesS[is.na(path[spcNodesS])] # don't overwrite
  path[spcNodesS] <- node 
  spcNode <- spcNode_new # update spcNode
  if(node <= Ntip(tree)){ # node is tip
    return(matrix(path,nrow=1,dimnames = list(tree$tip.label[node],names(path))))
  }
  
  # recall for each child-node
  e <- tree$edge
  child_nodes <- e[e[,1]==node,2]
  
  retVal<- NULL
  for( child in child_nodes){
    # retVal <- c(retVal, Recall(path,node=child,spcNode,tree,S,D,spcTree))
    retVal <- rbind(retVal, Recall(path,node=child,spcNode,tree,S,D,spcTree))
  }
  return(retVal)
}


#### Main function ####


#' Create OGtbl from gene-trees
#'
#' @param trees named list of trees (treeio::treedata objects). 
#' @param spcTree The species tee
#'
#' @return
#' @export
#'
#' @examples
makeOGtbl <- function(trees,spcTree){
  path_blank <- setNames(rep(NA,Nnode(spcTree)),spcTree$node.label)
  spcNodeRoot <- getRootNode(spcTree)
  
  imap_dfr(trees,function(treedata,treeName){

    OGmat <- getOGmatrix(path=path_blank,
                         node = getRootNode(treedata@phylo),
                         spcNode = spcNodeRoot,
                         tree = treedata@phylo,
                         S = treedata@data$S,
                         D = treedata@data$D,
                         spcTree = spcTree)
    
    # add tree name as OG prefix and convert to dataframe
    OGmat[!is.na(OGmat)] <- paste0(treeName,".",OGmat[!is.na(OGmat)])
    as.data.frame(OGmat,stringsAsFactors = F) %>% 
      rownames_to_column(var="geneID") %>% 
      mutate( spc = sub(".*(....)$","\\1",geneID),
              geneID = sub("_....$","",geneID),
              OG=treeName ) %>% 
      select(OG, spc, geneID, everything()) # reorder columns
      
    
  })
}

