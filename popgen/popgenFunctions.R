# Copyright (c) 2017  Genome  Research  Ltd.
# Author: Alistair Dunham
# This  program  is free  software: you  can  redistribute  it and/or  modify  it  under
# the  terms  of the  GNU  General  Public  License  as  published  by the  Free  Software
# Foundation; either  version 3 of the  License , or (at your  option) any  later
# version.
# This  program  is  distributed  in the  hope  that it will be useful , but  WITHOUT
# ANY  WARRANTY; without  even  the  implied  warranty  of  MERCHANTABILITY  or  FITNESS
# FOR A PARTICULAR  PURPOSE. See  the  GNU  General  Public  License  for  more
# details.
# You  should  have  received a copy of the  GNU  General  Public  License  along  with
# this  program. If not , see <http :// www.gnu.org/licenses/>.

## Script containing functions supporting population genetic analysis
## mainly used for converting genotypes to genepop format for use with diveRsity package

## Function to convert a feature table to genepop format. Pops is a named vector of sample population pairs
## Deletion is always allele 2, wt allele 1
saveGenepop <- function(feats,pops,outFile,main="HGDP Populations GenomeStrip Genotypes"){
  sink(outFile)
  cat(main,"\n")
  cat(paste(rownames(feats),collapse=","),"\n")
  for (p in unique(pops)){
    cat("Pop","\n")
    for (i in names(pops)[pops == p]){
      cat(pops[i]," , ",paste(sapply(feats[,i],getGeno),collapse = " "),"\n")
    }
  }
  sink()
}

gens <- c("0101","0102","0202")
getGeno <- function(g){
  if (is.na(g)){
    return("0000")
  } else{
    return(gens[g + 1])
  }
}

###### Function to load pairwise table produced by diveRsity
readPairwiseTables <- function(x,nam = c("gst","Gst","GGst","D","Fst"),symm=TRUE){
  lines <- readLines(x)
  inds <- which(lines %in% nam)
  out <- list()
  out[["info"]] <- lines[1:(min(inds) - 2)]
  for (i in 1:length(nam)){
    if (i == length(nam)){
      tmp <- read.table(textConnection(lines[(inds[i] + 2):(length(lines) - 1)]),sep="\t",header=TRUE,row.names = 1)
    } else {
      tmp <- read.table(textConnection(lines[(inds[i] + 2):(inds[i+1] - 2)]),sep="\t",header=TRUE,row.names = 1)
    }
    rownames(tmp) <- gsub(",","",rownames(tmp))
    colnames(tmp) <- gsub("\\.","",colnames(tmp))
    if (symm){
      tmp[upper.tri(tmp)] <- t(tmp)[upper.tri(tmp)]
    }
    out[[nam[i]]] <- tmp
  }
  return(out)
}

#### Convert accession to sample name
accession2name <- function(x,pops){
  n <- pops[pops$sample_accession == x,"sample"]
  return(n[1])
}
