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

## Script to import data from filtered breakdancer and genomestrip calls for population genetic analysis
## Made separate to tidy up the other two scripts and avoid overlap
setwd('~/Project/')

############ Import Population Info ############
pops <- read.table("meta/hgdp.meta.april11.accessions.txt",header=TRUE,stringsAsFactors = FALSE)
p <- strsplit(pops$sample,".",fixed = TRUE)
pops$sample <- sapply(p,function(x){x[1]})
pops$library <- sapply(p,function(x){x[2]})
pops[grep("meyer2012$",pops$sample),"sample_accession"] <- grep("meyer2012$",pops$sample,value = TRUE)
pops[grep("sgdp",pops$source),"sample_accession"] <- make.names(pops[grep("sgdp",pops$source),"library"])
pops$sample <- gsub("_meyer2012","",pops$sample)
pops <- pops[order(pops$sample),]

regCols <- c("red","blue","orange","purple","green","yellow","pink")
names(regCols) <- c("EUROPE","EAST_ASIA","AFRICA","OCEANIA","CENTRAL_SOUTH_ASIA","MIDDLE_EAST","AMERICA")

populations <- unique(pops$population)
popRegs <- sapply(populations,function(x){pops[which(pops$population == x)[1],"region"]})
populations <- populations[order(popRegs)]
popCols <- c(colorRampPalette(c("darkorange2","orange"))(7),
             colorRampPalette(c("pink","hotpink"))(5),
             colorRampPalette(c("green","darkgreen"))(9),
             colorRampPalette(c("blue","lightblue"))(18),
             colorRampPalette(c("firebrick","red"))(8),
             colorRampPalette(c("yellow","goldenrod1"))(4),
             colorRampPalette(c("darkorchid4","darkorchid"))(2)
)
names(popCols) <- populations

popSummaries <- read.table("meta/hgdp_pops_BD.tsv",header=TRUE,stringsAsFactors = FALSE,fill = TRUE)

############ Import Feature Tables ############
## Binary
tenXfeats <- read.table("results/mergedTenXcalls_featureTable.tsv",header=TRUE,row.names = 1)
BDfeats <- read.table("results/mergedBDcalls_featureTable.tsv",header=TRUE,row.names = 1)

## With CN - also generate deletions number
GSfeats <- read.table("results/gs_cnv.featureTable.tsv",header=TRUE,row.names = 1)
GSfeats.delGenotype <- 2 - GSfeats
GSfeats.delGenotype[GSfeats.delGenotype < 0] <- 0

############ Import Deletions ############
BDdels <- read.table("results/mergedBDcalls.sorted.merged.bed")
colnames(BDdels) <- c("chr","start","end")
BDdels$len <- BDdels$end - BDdels$start

GSdels <- read.table("results/gs_cnv.bed")
colnames(GSdels) <- c("chr","start","end","ID")
GSdels$len <- GSdels$end - GSdels$start

############ Sample Colours ##############
samplePopColsBD <- sapply(colnames(BDfeats),function(x){p <- pops[pops$sample == x,"population"];return(unname(popCols[p[1]]))})
sampleRegColsBD <- sapply(colnames(BDfeats),function(x){p <- pops[pops$sample == x,"region"];return(unname(regCols[p[1]]))})

samplePopColsTenX <- sapply(colnames(tenXfeats),function(x){p <- pops[pops$sample == x,"population"];return(unname(popCols[p[1]]))})
sampleRegColsTenX <- sapply(colnames(tenXfeats),function(x){p <- pops[pops$sample == x,"region"];return(unname(regCols[p[1]]))})

samplePopColsGS <- sapply(colnames(GSfeats),function(x){p <- pops[pops$sample_accession == x,"population"];return(unname(popCols[p[1]]))})
sampleRegColsGS <- sapply(colnames(GSfeats),function(x){p <- pops[pops$sample_accession == x,"region"];return(unname(regCols[p[1]]))})

## Pop Subsets
popsBD <- pops[pops$sample %in% colnames(BDfeats),]
popsGS <- pops[pops$sample_accession %in% colnames(GSfeats),]
