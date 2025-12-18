library(treeio)

# Reads NHX trees a bit faster than the function provided in treeio (it is still slow though)
myReadNHX <- function (file) 
{
  treetext <- readLines(file, warn = FALSE)
  treetext <- treetext[treetext != ""]
  treetext <- treetext[treetext != " "]
  if (length(treetext) > 1) {
    treetext <- paste0(treetext, collapse = "")
  }
  treetext <- gsub(" ", "", treetext)
  # read phylogeny without NHX features
  phylo <- read.tree(text = treetext)
  nnode <- Nnode(phylo, internal.only = FALSE)
  nlab <- paste("X", 1:nnode, sep = "")
  # find the NHX features by matching the pattern:
  #   example "NODENAME:1.23e4[&&NHX blabla]", where both both the NODENAME and distance is optional.
  pattern <- "(\\w+)?(:?\\d*\\.?\\d*[Ee]?[\\+\\-]?\\d*)?\\[&&NHX.*?\\]"
  
  # tree2 <- treetext
  # for (i in 1:nnode) {
  #   tree2 <- sub(pattern, paste0( nlab[i], "\\2"), 
  #                tree2)
  # }
  
  # a faster way (Note it does not keep the edge lengths, but this should hopefully not affect the order of the nodes)
  tree2 <- character(length(nlab)*2+1)
  tree2[2*(1:length(nlab))] <- nlab
  tree2[1+2*(0:length(nlab))] <- regmatches(x = treetext,m = gregexpr(pattern,treetext),invert = T)[[1]]
  tree2 <- paste(tree2,collapse = "")
  
  phylo2 <- read.tree(text = tree2)
  node <- match(nlab, sub(".+(X\\d+)$", "\\1", c(phylo2$tip.label, 
                                                 phylo2$node.label)))
  nhx.matches <- gregexpr(pattern, treetext)
  matches <- nhx.matches[[1]]
  match.pos <- as.numeric(matches)
  if (length(match.pos) == 1 && (match.pos == -1)) {
    nhx_tags <- data.frame(node = 1:nnode)
  }
  else {
    match.len <- attr(matches, "match.length")
    nhx_str <- substring(treetext, match.pos, match.pos + 
                           match.len - 1)
    nhx_features <- gsub("^[^\\[]*", "", nhx_str) %>% gsub("\\[&&NHX:", 
                                                           "", .) %>% gsub("\\]", "", .)
    nhx_tags <- treeio:::get_nhx_feature(nhx_features)
    fields <- names(nhx_tags)
    for (i in ncol(nhx_tags)) {
      if (any(grepl("\\D+", nhx_tags[, i])) == FALSE) {
        nhx_tags[, i] <- as.numeric(nhx_tags[, i])
      }
    }
    nhx_tags$node <- as.integer(node)
  }
  nhx_tags <- nhx_tags[order(nhx_tags$node), ]
  new("treedata", file = treeio:::filename(file), phylo = phylo, data = as_tibble(nhx_tags))
}
