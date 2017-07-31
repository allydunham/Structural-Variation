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

## Initiate
setwd('~/Project/')
library(data.table)

## Function to load breakdancer calls
bdLoad <- function(bd,tenx=NA,sou=NA){
  if (is.na(tenx)){
    bdCalls <- read.table(bd,stringsAsFactors = FALSE,blank.lines.skip = TRUE,comment.char = "#",fill = TRUE)
    colnames(bdCalls) <- c("Chr1","Pos1","Orientation1","Chr2","Pos2","Orientation2","Type","Size","Score","num_Reads","num_Reads_lib","bam")
    bdCalls$ID <- paste("BD",sou,1:dim(bdCalls)[1],sep = "_")
    bdCalls <- bdCalls[bdCalls[,7] == "DEL",]
    bdCalls$tenXmatches <- NA
    
  } else{
    bdCalls <- read.table(bd,header = TRUE,stringsAsFactors = FALSE)
    tenxCalls <- read.table(tenx,header = TRUE,stringsAsFactors = FALSE)
    calls <- sapply(strsplit(bdCalls$tenXmatches,':'),function(x){x[1]})
    bdCalls <- bdCalls[calls %in% tenxCalls$ID,]
    
  }
  ori1 <- strsplit(bdCalls$Orientation1,"[+-]")
  ori2 <- strsplit(bdCalls$Orientation2,"[+-]")
  bdCalls$plus1 <- as.numeric(sapply(ori1,function(x){x[1]}))
  bdCalls$minus1 <- as.numeric(sapply(ori1,function(x){x[2]}))
  bdCalls$plus2 <- as.numeric(sapply(ori2,function(x){x[1]}))
  bdCalls$minus2 <- as.numeric(sapply(ori2,function(x){x[2]}))
  return(bdCalls[,c("Chr1","Pos1","Orientation1","plus1","minus1","Chr2","Pos2","Orientation2","plus2","minus2","Type","Size","Score","num_Reads","num_Reads_lib","bam","ID","tenXmatches")])
}

## get hgdp ids with breakdancer output
ids <- grep("HGDP0*",dir("data/"),value = TRUE)

## import overlapping calls
bds_overlap <- lapply(ids[!ids == "HGDP01032"],function(x){
  lapply(grep("*BD_out.overlap",dir(paste0("data/",x)),value = TRUE),function(y){bdLoad(paste0("data/",x,"/",y),paste0("data/",x,"/tenx_rf_classified.tsv"))})
})

bds_overlap.all <- rbindlist(lapply(bds_overlap,rbindlist))


## Import all breakdancer calls
bds_all <- lapply(ids[!ids == "HGDP01032"],function(x){
  lapply(grep("*BD_out$",dir(paste0("data/",x)),value = TRUE),function(y){bdLoad(paste0("data/",x,"/",y),sou = gsub(".BD_out","",y))})
})

bds_all.all <- rbindlist(lapply(bds_all,rbindlist))

## Process breakdancer calls and extract those that don't overlap
bds_all.all$readRatio <- bds_all.all$num_Reads/((bds_all.all$plus1 + bds_all.all$minus1 + bds_all.all$plus2 + bds_all.all$minus2)/2)
bds_overlap.all$readRatio <- bds_overlap.all$num_Reads/((bds_overlap.all$plus1 + bds_overlap.all$minus1 + bds_overlap.all$plus2 + bds_overlap.all$minus2)/2)
bds_nonoverlap.all <- bds_all.all[!bds_all.all$ID %in% bds_overlap.all$ID,]

#### Compar statistics between overlap and non-overlap set
## Density of score for both cases - highly concentrated at high scores for calls in "gold-standard"
pdf("Figures/bdFiltering.pdf",10,14)
par(oma=c(4,1,4,1),xpd=FALSE)
layout(rbind(c(1,2),c(3,4),c(5,5)))
plot(density(bds_overlap.all$Score),col="blue",xlim=c(30,100),main = "A. Call Score",xlab = "Score")
lines(density(bds_nonoverlap.all$Score),col="red")
abline(v=c(93),lwd=1.5)

## Similarly for length. Also see LINE (and SINE?) peaks? - could be a result of call size difference between bd and longranger
hist(bds_nonoverlap.all$Size,breaks = 100000,freq = FALSE,col=rgb(1,0,0,1),xlim=c(0,100000),ylim=c(0,0.0004),
     main = "B. Length",xlab = "Length")
hist(bds_overlap.all$Size,freq = FALSE,col=rgb(0,0,1,0.6),add=TRUE)
abline(v=c(200,30000),lwd=1.5)

## Read Depths
hist(bds_overlap.all$num_Reads,col="blue",ylim=c(0,0.1),freq=FALSE,xlim = c(0,1000),
     main = "C. Supporting Reads",xlab = "Number of Reads")
hist(bds_nonoverlap.all$num_Reads,col="red",freq=FALSE,add=TRUE)
abline(v=c(0,50),lwd=1.5)

## estimated CN
hist(bds_nonoverlap.all$bam,col=rgb(1,0,0,1),freq=FALSE,xlim = c(-10,10),ylim = c(0,2),breaks = 10000,
     main = "D. Copy Number Estimate",xlab = "CN Estimate")
hist(bds_overlap.all$bam,col=rgb(0,0,1,0.6),freq=FALSE,add=TRUE)
abline(v=c(0,5),lwd=1.5)

## Read Ratio
hist(bds_nonoverlap.all$readRatio,col=rgb(1,0,0,0.75),main = "E. Supporting Read Ratio",xlab = "Read Ratio")
hist(bds_overlap.all$readRatio,col=rgb(0,0,1,0.75),add=TRUE)
abline(v=c(0.6,1),lwd=1.5)

legend(0.5,-8000,c("10X Overlap","No Overlap","Filter Boundary"),fill=c("blue","red",NA),
       horiz = TRUE,bty = "n",xpd=NA,xjust = 0.5,lty=c(NA,NA,1),border=c("blue","red",NA),
       x.intersp = c(-0.5,-0.5,0.5),lwd = c(NA,NA,2))
title("BreakDancer Filtering Criteria",outer=TRUE)
dev.off()

## oddities in breakpoint read mapping
plot(density(bds_overlap.all$plus1),col="blue",ylim = c(0,1))
lines(density(bds_overlap.all$plus2),col="red")
lines(density(bds_overlap.all$minus2),col="firebrick1")
lines(density(bds_overlap.all$minus1),col="cornflowerblue")

plot(density(bds_nonoverlap.all$plus1),col="blue",xlim=c(-0,200),ylim = c(0,1))
lines(density(bds_nonoverlap.all$plus2),col="red")
lines(density(bds_nonoverlap.all$minus2),col="firebrick1")
lines(density(bds_nonoverlap.all$minus1),col="cornflowerblue")

## Filter test - see how well filter criteria perform against overlap set
# Import genomic gap regions
gaps <- read.table("data/GRCh38_gaps.txt")
colnames(gaps) <- c("bin","chrom","start","end","ix","n","size","type","bridge")

# Initialise centromeres
centromeres <- as.data.frame(matrix(c('chr1',122026460,125184587,'chr2',92188146,94090557,'chr3',90772459,93655574,
                      'chr4',49708101,51743951,'chr5',46485901,50059807,'chr6',58553889,59829934,
                      'chr7',58169654,60828234,'chr8',44033745,45877265,'chr9',43236168,45518558,
                      'chr10',39686683,41593521,'chr11',51078349,54425074,'chr12',34769408,37185252,
                      'chr13',16000001,18051248,'chr14',16000001,18173523,'chr15',17000001,19725254,
                      'chr16',36311159,38280682,'chr17',22813680,26885980,'chr18',15460900,20861206,
                      'chr19',24498981,27190874,'chr20',26436233,30038348,'chr21',10864561,12915808,
                      'chr22',12954789,15054318,'chrX',58605580,62412542,'chrY',10316945,10544039),
                      ncol = 3,byrow = TRUE),stringsAsFactors = FALSE)

colnames(centromeres) <- c("chrom","start","end")
centromeres$start <- as.numeric(centromeres$start)
centromeres$end <- as.numeric(centromeres$end)

# function determining whether or not two regions overlap
overlap <- function(x,y){
  return(ifelse(x[1] <= y[2] & y[1] <= x[2],TRUE,FALSE))
}

# indeces of calls overlapping with centromeres and genomic gaps
centInd <- apply(bds_all.all,1,function(x){if (x[1] %in% centromeres$chrom){overlap(c(x[2],x[7]),centromeres[centromeres$chrom == x[1],c(2,3)])}else{FALSE}})
gapInd <- apply(bds_all.all,1,function(x){if (x[1] %in% gaps$chrom){any(apply(gaps[gaps$chrom == x[1],],1,function(y){overlap(c(x[2],x[7]),
                                                                                                                              c(y[3],y[4]))}))}else{FALSE}})

# Generate filtered set
bds_filtered <- bds_all.all[(bds_all.all$Size >= 200 &
                            bds_all.all$Size <= 30000) &
                            bds_all.all$Score >= 93 &
                            bds_all.all$num_Reads <= 50 &
                            ((bds_all.all$bam >= 0 &
                              bds_all.all$bam <= 5) |
                               is.na(bds_all.all$bam)) & 
                              bds_all.all$readRatio > 0.6,]
 

# Test overlap of filtered with overlap set
sum(bds_overlap.all$ID %in% bds_filtered$ID)/20003 ## 0.7839824
sum(bds_filtered$ID %in% bds_overlap.all$ID)/dim(bds_filtered)[1] ## 0.4838481


## Test against genome strip (filtered to same size)
gsOverlap <- read.table("results/gs_bd_intersect_dels",stringsAsFactors = FALSE)
colnames(gsOverlap) <- c("chrom1","start1","end1","chrom2","start2","end2","name2","tt","overlap")
gsOverlap <- gsOverlap[,c("chrom1","start1","end1","chrom2","start2","end2","name2","overlap")]
gsOverlap$size1 <- gsOverlap$end1 - gsOverlap$start1
gsOverlap$size2 <- gsOverlap$end2 - gsOverlap$start2
gsOverlap$prop1 <- gsOverlap$overlap/gsOverlap$size1
gsOverlap$prop2 <- gsOverlap$overlap/gsOverlap$size2

gsOverlap.filtered <- gsOverlap[ gsOverlap$size1 < 30000 & gsOverlap$size1 > 1000,]
gsOverlap.filtered <- gsOverlap.filtered[!duplicated(gsOverlap.filtered[,c("chrom1","start1","end1")]),]

sum(gsOverlap.filtered$prop1 >= 0.5)/dim(gsOverlap.filtered)[1]

hist(gsOverlap$size1,col=rgb(1,0,0,0.75))
hist(gsOverlap$size2,col=rgb(0,0,1,0.75),add=TRUE)

hist(gsOverlap.filtered[,"prop1"],col=rgb(0,0,1,0.75)) ## improves when sex chrom excluded !gsOverlap.filtered$chrom1 %in% c("chrX","chrY")

## good overlap of first match, when filtered to overlaping size range (1kb to 30kb)
