library(tidyverse)
library(treeio)
library(yaml)


#### Functions ####

#source("myReadNHX.R")

splitRecords <- function(x, sep){
  isSplitLines <- grepl(sep,x)
  if(isSplitLines[length(isSplitLines)]){
    cat("Split on last line of chunk!\n")
    # Keep the split line in the last record to prevent it from being merged with the next record
    isSplitLines[length(isSplitLines)] <- F
  }
  idx <- cumsum(isSplitLines)
  splitLines <- which(isSplitLines)
  split(x[-splitLines],idx[-splitLines])
}

getTreeFromRecord <- function(x){
  treeTxt <- x[grep("^DATA",x)+1]
  tmpTreeFile <- tempfile()
  writeLines(treeTxt,tmpTreeFile)
  # tree <- read.nhx(tmpTreeFile)
  tree <- myReadNHX(tmpTreeFile)
  unlink(tmpTreeFile)
  return(tree)  
}

getSeqMetaFromRecord <- function(x){
  suppressWarnings( # some entries have an extra column which generates warnings.. Just ignore those
    paste(x[grepl("^SEQ",x)],collapse = "\n") %>%
      read_delim(col_names = c("SEQ","species","protID","chrm","start","end","strand","geneID","symbol"),
                 delim = " ",show_col_types = FALSE) %>% 
      select(-SEQ)
  )
}

#### Main ####

read_nhx_emf_subtree <- function( nhx_emf_filename, species2spc){
  
  chunkSize = 5e4 # must be larger than the largest record
  
  con <- file(nhx_emf_filename,open = "r")
  
  leftOver <- character(0)
  endOfFile <- FALSE
  allTrees <- list()
  
  linesRead <- 0
  chunksRead <- 0
  recordsRead <- 0
  
  allTrees <- list()
  
  while(!endOfFile){
    txt <- readLines(con,n = chunkSize)
    chunksRead <- chunksRead + 1
    if(length(txt)==0){
      endOfFile <- TRUE
      if(length(leftOver)<2) break # last line of file
      treeEMFRecords <- leftOver
    } else{
      treeEMFRecords <- splitRecords(c(leftOver,txt),sep="^//")
  
      leftOver <- treeEMFRecords[[length(treeEMFRecords)]] # last record may be incomplete
      treeEMFRecords <- treeEMFRecords[-length(treeEMFRecords)]
    }
    
    cat("Reading chunk",chunksRead,"with",length(treeEMFRecords),"trees...\n")
    
    # read the trees
  
    trees <- lapply(treeEMFRecords, function(x){
      # For each record
      recordsRead <<- recordsRead + 1
      
      # read sequence metadata
      seqMeta <- 
        getSeqMetaFromRecord(x) %>% 
        inner_join(species2spc,by="species")
      
      # read tree
      if( nrow(seqMeta)>0 ){
        # cat(seqMeta %>% count(spc) %>% with(paste0(spc,":",n,collapse=", ")))
        # if(!all(species2spc$spc %in% seqMeta$spc))
        #   cat("    (Missing spc:",species2spc$spc[!(species2spc$spc %in% seqMeta$spc)],")")
        # cat("\nreading tree nr",recordsRead,"in chunk nr",chunksRead,"at line",linesRead,"with",nrow(seqMeta),"tips\n")
        fullTree <- getTreeFromRecord(x)
      } else {
        # cat("Skipping empty record nr",recordsRead,"at line",linesRead,"\n")
        linesRead <<- linesRead + length(x) + 1
        return(NULL)
      }
      linesRead <<- linesRead + length(x) + 1
      
      
      stopifnot(all(seqMeta$protID %in% fullTree@phylo$tip.label))
      
      tree <- drop.tip(fullTree,which(!(fullTree@phylo$tip.label %in% seqMeta$protID)))
      
      # set tip labels to be geneIDs with four letter species suffix
      tree@phylo$tip.label <- with(seqMeta, setNames(paste0(geneID,"_",spc),protID))[tree@phylo$tip.label]
      tree@phylo$node.label <- paste0("n",Ntip(tree) + 1:Nnode(tree))
      
      return(tree)
      
    })
    # saveRDS(trees,sprintf("treesPerChunk/trees_chunk%03i.RDS",chunksRead))
    
    allTrees <- c(allTrees, trees)
  }
  
  close(con)
  
  cat("Finished. Total number of trees read:",length(allTrees),"\n")
  
  
  names(allTrees) <- NULL
  
  return(allTrees)
}

read_nhx_emf <- function( nhx_emf_filename){
  
  chunkSize = 5e4 # must be larger than the largest record
  
  con <- file(nhx_emf_filename,open = "r")
  
  leftOver <- character(0)
  endOfFile <- FALSE
  allTrees <- list()
  
  linesRead <- 0
  chunksRead <- 0
  recordsRead <- 0
  
  allTrees <- list()
  
  while(!endOfFile){
    txt <- readLines(con,n = chunkSize)
    chunksRead <- chunksRead + 1
    if(length(txt)==0){
      endOfFile <- TRUE
      if(length(leftOver)<2) break # last line of file
      treeEMFRecords <- leftOver
    } else{
      treeEMFRecords <- splitRecords(c(leftOver,txt),sep="^//")
      
      leftOver <- treeEMFRecords[[length(treeEMFRecords)]] # last record may be incomplete
      treeEMFRecords <- treeEMFRecords[-length(treeEMFRecords)]
    }
    
    cat("Reading chunk",chunksRead,"with",length(treeEMFRecords),"trees...\n")
    
    # read the trees
    
    trees <- lapply(treeEMFRecords, function(x){
      # For each record
      recordsRead <<- recordsRead + 1
      # cat("Reading record",recordsRead,"\n")
      
      # read sequence metadata
      seqMeta <- getSeqMetaFromRecord(x) 
      
      tree <- getTreeFromRecord(x)
      linesRead <<- linesRead + length(x) + 1
      
      # assert that metadata matches tree...
      stopifnot(all(seqMeta$protID %in% tree@phylo$tip.label))
      
      # add metadata
      tree@data <- 
        tree@data %>%
        left_join(seqMeta, by=c("G" = "geneID"))
        
      return(tree)
      
    })
    # saveRDS(trees,sprintf("treesPerChunk/trees_chunk%03i.RDS",chunksRead))
    
    allTrees <- c(allTrees, trees)
  }
  
  close(con)
  
  cat("Finished. Total number of trees read:",length(allTrees),"\n")
  
  
  names(allTrees) <- NULL
  
  return(allTrees)
}


