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

## Script containing code to perform all genomestrip population genetics and produce output plots

setwd('~/Project/')
source("scripts/importPopgen.R")

############ Import VEP ############
gsVEPimpacts <- read.table("results/gs_cnv.effects.vep.txt",header = TRUE,comment.char = '',stringsAsFactors = FALSE)
## Process consequences string and make columns
cons <- strsplit(gsVEPimpacts$Consequence,",")
t <- unique(unlist(cons))
for (i in t){
  gsVEPimpacts[,i] <- as.numeric(sapply(cons,function(x){any(i %in% x)}))
}
gsVEPimpacts <- gsVEPimpacts[,-7]

## Process additional info into columns 
sortExtras <- function(x){
  t <- strsplit(x,"=")
  n <- sapply(t,function(u){u[1]})
  v <- sapply(t,function(u){u[2]})
  names(v) <- n
  return(v)
}

ex <- strsplit(gsVEPimpacts$Extra,';')
ex <- lapply(ex,sortExtras)

t <- unique(unlist(sapply(ex,names)))
for (i in t){
  gsVEPimpacts[,i] <- sapply(ex,function(x){if (i %in% names(x)){x[i]}else{NA}})
}
gsVEPimpacts <- gsVEPimpacts[,-13]

## Assign call handles to each VEP line
gsVEPimpacts$call <- NA
for (i in 1:dim(GSdels)[1]){
  gsVEPimpacts[gsVEPimpacts$X.Uploaded_variation == paste(GSdels[i,1],GSdels[i,2],"deletion",sep = "_"),"call"] <- i
}

## Reduce vep table to those with an associated gene symbol and get genes
gsVEPimpactsRed <- gsVEPimpacts[!is.na(gsVEPimpacts$SYMBOL),]
genes <- unique(gsVEPimpactsRed$SYMBOL)

save(GSdels,GSfeats,GSfeats.delGenotype,pops,popsGS,popSummaries,gsVEPimpacts,gsVEPimpactsRed,genes,popCols,
     popRegs,populations,regCols,samplePopColsGS,sampleRegColsGS,file = "results/GSpopgenImport.Rdata")

############# POPGEN ##################
## (can start from here if data has been preloaded and saved)
setwd('~/Project/')
load("results/GSpopgenImport.Rdata")
source("scripts/popgenFunctions.R")

## Filter any unclear calls and those with no deletions
GSfeats.delGenotype.noNA <- GSfeats.delGenotype[rowSums(is.na(GSfeats.delGenotype)) == 0,]
GSfeats.delGenotype.noNA <- GSfeats.delGenotype.noNA[!rowSums(GSfeats.delGenotype.noNA) == 0,]
GSdels.noNA <- GSdels[which(!is.na(rowSums(GSfeats.delGenotype)) & rowSums(GSfeats.delGenotype) > 0),]

hom.freq <- rowSums(GSfeats == 0)/dim(GSfeats)[2]
het.freq <- rowSums(GSfeats == 1)/dim(GSfeats)[2]

###### dendogram ######
popClust <- hclust(dist(t(GSfeats.delGenotype.noNA)))
popDend <- as.dendrogram(popClust)

catNodePop <- function(x){
  if (is.leaf(x)){
    a <- attributes(x)
    p <- pops[pops$sample_accession == a$label,"population"]
    c <- unname(popCols[p[1]])
    attr(x, "nodePar") <- c(a$nodePar, col = c,pch=1)
    attr(x, "label") <- p[1]
  }
  x
}

popDend.popCols <- dendrapply(popDend,catNodePop)
plot(popDend.popCols,leaflab="none")

library(ape)
phylo <- as.phylo(popClust)
tipCols <- samplePopColsGS[phylo$tip.label]
phylo$tip.label <- sapply(phylo$tip.label,accession2name,pops=pops)

pdf("Figures/GSphylo.pdf",10,8)
par(oma=c(1,1,1,8),mar=c(1,1,4,4))
plot(phylo,type = "unrooted",use.edge.length = FALSE,tip.color = tipCols,
     lab4ut="axial",main="Heirarchical Clustering of GenomeSTRiP Derived Deletion Genotypes")
legend(50,20,xjust = 0.5,yjust=0.5, legend = c("Africa","America","Central and S. Asia",
                                                 "E. Asia","Europe","Middle East","Oceania"),
       col=c("darkorange","hotpink","green","blue","red","yellow","purple"),pch=20,xpd=NA,
       title = expression(bold("Sample Region")),bty="n")
dev.off()

############ PCA ############
pcGS <- prcomp(t(GSfeats.delGenotype.noNA),tol=10^-10)
par(mar=c(4,4,4,4),oma=c(1,1,1,1))
plot(pcGS$x[,1],pcGS$x[,2],pch=20,col=samplePopColsGS)

library(rgl)
plot3d(pcGS$x[,1],pcGS$x[,2],pcGS$x[,3],col=samplePopColsGS)

vars <- round(pcGS$sdev^2/sum(pcGS$sdev^2) * 100)
pdf("Figures/GSpca.pdf",10,14)
par(oma=c(4,4,4,4),mar=c(1,1,1,1),xpd=NA,cex.lab=1.5)
layout(rbind(c(1,2,3),c(0,4,5),c(7,0,6)))
plot(pcGS$x[,2],pcGS$x[,1],pch=20,col=samplePopColsGS,xlab = paste0("PC2 (",vars[2],"%)"),ylab = paste0("PC1 (",vars[1],"%)"))
plot(pcGS$x[,3],pcGS$x[,1],pch=20,col=samplePopColsGS,yaxt="n",xaxt="n",xlab = "",ylab = "")
plot(pcGS$x[,4],pcGS$x[,1],pch=20,col=samplePopColsGS,yaxt="n",xaxt="n",xlab = "",ylab = "")

plot(pcGS$x[,3],pcGS$x[,2],pch=20,col=samplePopColsGS,xlab = paste0("PC3 (",vars[3],"%)"),ylab = paste0("PC2 (",vars[2],"%)"))
plot(pcGS$x[,4],pcGS$x[,2],pch=20,col=samplePopColsGS,yaxt="n",xaxt="n",xlab = "",ylab = "")

plot(pcGS$x[,4],pcGS$x[,3],pch=20,col=samplePopColsGS,xlab = paste0("PC4 (",vars[4],"%)"),ylab = paste0("PC3 (",vars[3],"%)"))

plot(NA,axes=FALSE,ylim=c(0,1),xlim=c(0,1),xlab="",ylab = "")
legend(0.5,0.5,xjust = 0.5,yjust=0.5, legend = c("Africa","America","Central and S. Asia",
                                                 "E. Asia","Europe","Middle East","Oceania"),
       col=c("darkorange","hotpink","green","blue","red","yellow","purple"),pch=20,xpd=NA,
       title = expression(bold("Sample Region")),bty="n",cex=1.5)

title("Principal Component Analysis of GenomeSTRiP Deletion Calls",outer = TRUE,cex.main=2)
dev.off()
layout(1)


### By region
region <- "AFRICA"
ps <- popsGS[popsGS$region == region,"population"]
pcGS.reg <- prcomp(t(GSfeats.delGenotype.noNA[,popsGS[popsGS$region == region,"sample_accession"]]),tol=10^-10)
plot(pcGS.reg$x[,1],pcGS.reg$x[,2],pch=20,col=as.factor(ps))
legend("bottomright",legend = unique(ps),col=unique(as.factor(ps)),bty="n",pch=20)

pcGS.reg <- lapply(sort(unique(popsGS$region)),function(x){prcomp(t(GSfeats.delGenotype.noNA[,popsGS[popsGS$region == x,"sample_accession"]]),tol=10^-10)})
names(pcGS.reg) <- sort(unique(popsGS$region))
  
pdf("Figures/GSpca_regions1.pdf",8,12)
par(oma=c(1,1,3,1))
layout(rbind(c(1,1,1,2),c(3,3,3,4),c(5,5,5,6),c(7,7,7,8)))
nam <- c("A. Africa","B. America", "C. Oceania","tt")
reg <- c("CENTRAL_SOUTH_ASIA","EAST_ASIA","EUROPE","MIDDLE_EAST")
for (i in 1:length(reg)){
  r <- reg[i]
  ps <- popsGS[popsGS$region == r,"population"]
  vars = round(pcGS.reg[[r]]$sdev^2/sum(pcGS.reg[[r]]$sdev^2) * 100)
  par(mar=c(4,4,4,0))
  plot(pcGS.reg[[r]]$x[,1],pcGS.reg[[r]]$x[,2],pch=20,col=as.factor(ps),
       xlab = paste0("PC1 (",vars[1],"%)"),ylab = paste0("PC2 (",vars[2],"%)"),
       main = nam[i])
  par(mar=c(0,2,0,0))
  plot(NA,axes=FALSE,ylim=c(0,1),xlim=c(0,1),xlab="",ylab = "")
  if (r == "AFRICA"){
    legend(0,0.5,yjust = 0.5,legend = c("Mbuti","Biaka","Mandenka","Yoruba","San","Bantu (S. Africa)","Bantu (Kenya)"),col = as.factor(unique(ps)),bty="n",pch=20,xpd=NA,
           title = expression(bold("Population")))
  } else {
    legend(0,0.5,yjust = 0.5,legend = unique(ps),col = as.factor(unique(ps)),bty="n",pch=20,xpd=NA,
         title = expression(bold("Population")))
  }
}
title("Principal Component Analysis of GenomeSTRiP Deletion Genotypes by Region",outer = TRUE,cex.main=1.5)
mtext("Regions with Differentiation",outer = TRUE,cex=1.2,padj = 0.5)
dev.off()

pdf("Figures/GSpca_regions2.pdf",12,10)
par(mar=c(4,4,4,4),oma=c(1,1,5,1))
layout(rbind(c(1,2),c(3,4)))
nam <- c("A. Central & South Asia","B. East Asia","C. Europe","D. Middle East")
reg <- c("CENTRAL_SOUTH_ASIA","EAST_ASIA","EUROPE","MIDDLE_EAST")
for (i in 1:length(reg)){
  ps <- popsGS[popsGS$region == reg[i],"population"]
  r <- reg[i]
  vars = round(pcGS.reg[[r]]$sdev^2/sum(pcGS.reg[[r]]$sdev^2) * 100)
  plot(pcGS.reg[[r]]$x[,1],pcGS.reg[[r]]$x[,2],pch=20,col=popCols[ps],
       xlab = paste0("PC1 (",vars[1],"%)"),ylab = paste0("PC2 (",vars[2],"%)"),
       main = nam[i])
}
title("Principal Component Analysis of GenomeSTRiP Deletion Genotypes by Region",outer = TRUE,cex.main=1.5)
mtext("Regions with no Differentiation",outer = TRUE,cex=1.3)
dev.off()

### Other factors - no batch effects found
varGS <- pops[pops$sample_accession %in% colnames(GSfeats.delGenotype.noNA),c(1,9)]
varGS <- varGS[!duplicated(varGS),2]
plot(pcGS$x[,1],pcGS$x[,2],pch=20,col=as.factor(varGS))

############ Genetic Burden ############
### By bps
delInfo <- data.frame(total=colSums(GSfeats.delGenotype.noNA * GSdels.noNA$len),
                      region=sapply(colnames(GSfeats.delGenotype.noNA),function(x){r <- pops[pops$sample_accession == x,"region"];return(r[1])}),
                      population=sapply(colnames(GSfeats.delGenotype),function(x){r <- pops[pops$sample_accession == x,"population"];return(r[1])}),
                      regCol=sampleRegColsGS,
                      popCol=samplePopColsGS,stringsAsFactors = FALSE)

delInfo <- delInfo[order(delInfo$region,delInfo$population),]

pdf("Figures/GSburden_bp.pdf",14,10)
par(oma=c(0,0,0,18))
barplot(delInfo$total,border=delInfo$popCol,col=delInfo$popCol,main = "Total Number of Deleted Bases",xlab = "Population",ylab = "Base Pairs")
axis(1,c(0,cumsum(920/765*table(delInfo$region))),rep("",8))
text(cumsum(920/765*table(delInfo$region)) - 920/765*table(delInfo$region)/2,-250000,
     c("Africa","America","South Central Asia", "East Asia","Europe","Middle East","Oceania"),
     xpd=NA) 
legend(1200,3.5E6,names(popCols)[c(1:21,48:53)],col = popCols[c(1:21,48:53)],pch=20,bty = 'n',xpd=NA,yjust = 0.5)
legend(1420,3.5E6,names(popCols)[22:47],col = popCols[22:47],pch=20,bty = 'n',xpd=NA,yjust = 0.5)
text(1410,6.2E6,expression(bold("Population")),cex=1.5,xpd=NA)
dev.off()

### By count
delInfo.count <- data.frame(total=colSums(GSfeats.delGenotype.noNA),
                      region=sapply(colnames(GSfeats.delGenotype.noNA),function(x){r <- pops[pops$sample_accession == x,"region"];return(r[1])}),
                      population=sapply(colnames(GSfeats.delGenotype),function(x){r <- pops[pops$sample_accession == x,"population"];return(r[1])}),
                      regCol=sampleRegColsGS,
                      popCol=samplePopColsGS,stringsAsFactors = FALSE)

delInfo.count <- delInfo.count[order(delInfo.count$region,delInfo.count$population),]

delInfoHom <- data.frame(total=colSums((GSfeats.delGenotype.noNA == 2) * 2),
                      region=sapply(colnames(GSfeats.delGenotype.noNA),function(x){r <- pops[pops$sample_accession == x,"region"];return(r[1])}),
                      population=sapply(colnames(GSfeats.delGenotype),function(x){r <- pops[pops$sample_accession == x,"population"];return(r[1])}),
                      regCol=sampleRegColsGS,
                      popCol=samplePopColsGS,stringsAsFactors = FALSE)

delInfoHom <- delInfoHom[order(delInfoHom$region,delInfoHom$population),]

delInfoHet <- data.frame(total=colSums((GSfeats.delGenotype.noNA == 1)),
                         region=sapply(colnames(GSfeats.delGenotype.noNA),function(x){r <- pops[pops$sample_accession == x,"region"];return(r[1])}),
                         population=sapply(colnames(GSfeats.delGenotype),function(x){r <- pops[pops$sample_accession == x,"population"];return(r[1])}),
                         regCol=sampleRegColsGS,
                         popCol=samplePopColsGS,stringsAsFactors = FALSE)

delInfoHet <- delInfoHet[order(delInfoHet$region,delInfoHet$population),]

pdf("Figures/GSburden_count.pdf",18,12)
par(mar=c(5,5,3,1),oma=c(1,1,4,1))
layout(rbind(c(1,1,1,1,4),c(2,2,3,3,4)))
barplot(delInfo.count$total,border=delInfo.count$popCol,col=delInfo.count$popCol,
        main = "A. Total",xlab = "Population",ylab = "Count")
axis(1,c(0,cumsum(920/765*table(delInfo.count$region))),rep("",8))
text(cumsum(920/765*table(delInfo.count$region)) - 920/765*table(delInfo.count$region)/2,-50,
     c("Africa","America","South Central Asia", "East Asia","Europe","Middle East","Oceania"),
     xpd=NA) 

barplot(delInfoHom$total,border=delInfoHom$popCol,col=delInfoHom$popCol,ylim = c(0,350),
        main = "B. Homozygotes",xlab = "Population",ylab = "Count")
axis(1,c(0,cumsum(920/765*table(delInfoHom$region))),rep("",8))
text(cumsum(920/765*table(delInfoHom$region)) - 920/765*table(delInfoHom$region)/2,-15,
     c("Africa","America","South Central Asia", "East Asia","Europe","Middle East","Oceania"),
     xpd=NA) 

barplot(delInfoHet$total,border=delInfoHet$popCol,col=delInfoHet$popCol,
        main = "C. Heterozygotes",xlab = "Population",ylab = "Count")
axis(1,c(0,cumsum(920/765*table(delInfoHet$region))),rep("",8))
text(cumsum(920/765*table(delInfoHet$region)) - 920/765*table(delInfoHet$region)/2,-40,
     c("Africa","America","South Central Asia", "East Asia","Europe","Middle East","Oceania"),
     xpd=NA) 

par(mar=c(1,1,1,1))
plot(NA,axes=FALSE,ylim=c(0,1),xlim=c(0,1),xlab="",ylab = "")
legend(0,0.5,names(popCols)[c(1:21,48:53)],col = popCols[c(1:21,48:53)],pch=19,bty = 'n',xpd=NA,yjust = 0.5,cex=1.2)
legend(0.5,0.5,names(popCols)[22:47],col = popCols[22:47],pch=19,bty = 'n',xpd=NA,yjust = 0.5,cex=1.2)
text(0.4,0.73,expression(bold("Population")),cex=1.5)

title("Number of GenomeSTRiP Deletions in Different Individuals",outer = TRUE,cex.main=1.8)
dev.off()

## Other metrics
delInfoMet <- data.frame(total=colSums(GSfeats.delGenotype.noNA),
                          region=sapply(colnames(GSfeats.delGenotype.noNA),function(x){r <- pops[pops$sample_accession == x,"region"];return(r[1])}),
                          population=sapply(colnames(GSfeats.delGenotype),function(x){r <- pops[pops$sample_accession == x,"population"];return(r[1])}),
                          regCol=sampleRegColsGS,
                          popCol=samplePopColsGS,
                          metric=sapply(colnames(GSfeats.delGenotype.noNA),function(x){r <- pops[pops$sample_accession == x,"library_type"];return(r[1])}),
                          stringsAsFactors = FALSE)

delInfoMet <- delInfoMet[order(delInfoMet$region,delInfoMet$population),]

pdf("Figures/GSburden_libBatch.pdf",14,10)
par(oma=c(1,1,3,12),mar=c(5,5,1,1))
layout(rbind(c(1,1,1,2)))
co <- c("blue","red")
barplot(delInfoMet$total,border=co[as.factor(delInfoMet$metric)],
        col=co[as.factor(delInfoMet$metric)]
        ,xlab = "Population",ylab = "Number of Deleted Sites",ylim = c(0,1000))
axis(1,c(0,cumsum(920/765*table(delInfoMet$region))),rep("",8))
text(cumsum(920/765*table(delInfoMet$region)) - 920/765*table(delInfoMet$region)/2,-50,
     c("Africa","America","South Central Asia", "East Asia","Europe","Middle East","Oceania"),
     xpd=NA) 

boxplot(total~metric,col=co,data = delInfoMet,ylim=c(0,1000))
title("Effect of Library Prepartion Method on Number of Called Deletions Using GenomeSTRiP",outer = TRUE,cex.main=1.7)
legend(3,500,yjust = 0.5,legend = c("PCR","PCR Free"),fill = co,bty = "n",cex=1.5,xpd=NA)
dev.off()

######## Test difference in number of deletions ########
regions <- unique(delInfo.count$region)
tTests <- lapply(regions,function(x){
  lapply(regions,function(y){
    t.test(delInfo.count[delInfo.count$region == x,"total"],delInfo.count[delInfo.count$region == y,"total"])
  })
})
names(tTests) <- regions
for (i in tTests){
  names(i) <- regions
}

pVals <- sapply(tTests,function(x){sapply(x,function(y){y$p.value})})
rownames(pVals) <- regions


tTests.pops <- lapply(populations,function(x){
  lapply(populations,function(y){
    t.test(delInfo.count[delInfo.count$population == x,"total"],delInfo.count[delInfo.count$population == y,"total"])
  })
})
names(tTests.pops) <- populations
for (i in tTests.pops){
  names(i) <- populations
}

pVals.pops <- sapply(tTests.pops,function(x){sapply(x,function(y){y$p.value})})
rownames(pVals.pops) <- populations

diffPops <- pVals.pops
diffPops[,] <- 0
diffPops[which(pVals.pops < 0.01/53,arr.ind = TRUE)] <- 1

par(mar=c(8,8,5,5))
image(diffPops,axes=FALSE,xlab = "",ylab = "",x=0:53,y=0:53,col = c("red","black"))
axis(1,at=seq(0.5,52.5),labels = populations,las=2)
axis(2,at=seq(0.5,52.5),labels = populations,las=2)

############ Fst ############
getAlleleFreqs <- function(pop){
  s <- unique(popsGS[popsGS$population == pop,"sample_accession"])
  return(rowSums(GSfeats.delGenotype.noNA[,s[s %in% colnames(GSfeats.delGenotype.noNA)]])/(2*length(s)))
}

nPops <- length(populations)
popWeights <- sapply(populations,function(x){sum(popsGS$population == x)})
popWeights <- popWeights/sum(popWeights)


alleleFreqs <- sapply(populations,getAlleleFreqs)

p <- apply(alleleFreqs,1,function(x){sum(popWeights * x)})

sig <- sapply(1:dim(alleleFreqs)[1],function(i){sum(alleleFreqs[i,]^2*popWeights) - p[i]^2})

F_st <- sig/(p*(1-p))

pdf("Figures/GSfst_hist.pdf",12,10)
par(mar=c(5,5,4,4),oma=c(0,0,0,0))
hist(F_st,main = expression("Distribution of Per Loci F"[st]~"for GenomSTRiP Deletion Calls"),
     xlab = expression("F"[st]),col = "cornflowerblue",ylim = c(0,14000))
dev.off()

topFst <- F_st[order(F_st,decreasing = TRUE)[1:10]]
## 0.78 - GScall_46330 - private to san, modifier of IVD, involved in leucine synthesis
## 0.7 - GScall_46289 - oceania, lincRNA of unknown function, likely intron modifier
## 0.6 - GScall_10095 - papuan only, Argininosuccinate Synthetase 1 Pseudogene and linc00501 modifier
## 0.52 - GScall_53920 - mbutiPygmy,coding seq variants of SIGLEC14 and SIGLEC5 (cell lurface ig like lectins, expressed myeloid cells of the hemopoietic system, immune related)
## 0.5 - GScall_24264 - san, 3.5kb del on chr8, lose of exon in tcea1. "Necessary for efficient RNA polymerase II transcription elongation past template-encoded arresting sites"
## Association with HIV life cycle. Not CCR5 chemokine receptor which is found in some pops
## 0.45 - GScall_22544 - karitiana, 5' utr var of MGAM (maltase-glucoamylase), involved in starch metabolism
## 0.44 - GScall_19461 - bedouin (but low freq, 1%), STXBP5-AS1 intron variant, antisense RNA of STXBP5 (involved in synaptic vesicle fusion through memb)
## 0.44 - GScall_19462 - v. similar variant to above, also in bedouin
## 0.44 - GScall_35636 - bedouin, palestinian, mozabite,makrani low freq (<5%) and african high freq (about >40%) - intron var in DISC1FP1 lincRNA. Fusion protein with RNA is localised to mitochondria and inhibits oxidoreductase, mutations linked to psychiatric disease in other studies (https://www.ncbi.nlm.nih.gov/pubmed/25943690)
## 0.44 - GScall_46130 - bantuSA, stop lost and intron variants for OCA2 transcripts, melanocyt associated p protein, which may be involved in pigment and associated with albinoism - previously observed in SA pops!

getGenes <- function(x){
  return(unique(gsVEPimpactsRed[gsVEPimpactsRed$call == as.integer(strsplit(x,"_")[[1]][2]),"SYMBOL"]))
}

df_fst <- data.frame(call=names(F_st),
                     fst=F_st)

df_fst$genes <- lapply(names(F_st),getGenes)
df_fst$genes <- sapply(df_fst$genes,paste,collapse=",")
df_fst <- merge(df_fst,GSdels,by.x = "call",by.y = "ID")
df_fst <- df_fst[order(df_fst$fst,decreasing = TRUE),]
df_fst$highFreq <- df_fst$call %in% rownames(GSfeats.delGenotype.noNA)[rowSums(GSfeats.delGenotype.noNA > 0) > 700]

#### pairwise Fst ####
library(gplots)
pairwiseMetrics <- readPairwiseTables("results/GSpairwiseFst/pairwise.txt")
nam <- rownames(pairwiseMetrics$Fst)
nam[c(9,10,30,44,53)] <- c("Biaka","Mbuti","Bantu (SA)","Han (North)","Bantu (Kenya)")

pdf("Figures/fst_heatmap.pdf",20,16)
heatmap.2(as.matrix(pairwiseMetrics$Fst),symm = TRUE,col=colorRampPalette(c("blue","red"))(256),
          revC = TRUE,trace = "none",dendrogram = "row",symbreaks = FALSE,margins = c(12,14),
          symkey = FALSE,key.title = "Colour Key and Histogram",key.par = list(mar=c(4,4,4,14),cex.main=2,cex.lab=2),
          density.info = "none",
          key.xlab = expression("F"[st]),cexRow = 2,cexCol = 2,labRow = nam,labCol = nam,
          lmat = matrix(c(3,2,0,3,1,4),ncol = 2),lhei = c(0.5,6,1),lwid = c(1.5,6))
title(expression("Pairwise F"[st]~"Values Across All Populations"),cex.main=2.5)
dev.off()

diveRsityStats <- read.table("results/GSpairwiseFst/std_stats.txt",header=TRUE)
globalStats <- diveRsityStats[dim(diveRsityStats)[1],]
diveRsityStats <- diveRsityStats[1:(dim(diveRsityStats)[1]-1),]

hist(diveRsityStats$Fst)

###### Site Frequency ######
sitefreq <- rowSums(GSfeats.delGenotype.noNA)
impactCalls <- unique(gsVEPimpactsRed[gsVEPimpactsRed$IMPACT %in% c("MODERATE","HIGH"),"call"])

par(oma=c(0,0,0,0))
plot(density(sitefreq),col="blue",xlim=c(0,100))
lines(density(sitefreq[paste0("GScall_",impactCalls)],na.rm = TRUE),col="red")

pdf("Figures/GSsitefreq.pdf",12,10)
par(mar=c(5,5,3,1))
hist(log10(sitefreq[!names(sitefreq) %in% paste0("GScall_",impactCalls)]),xlim=c(0,4),ylim = c(0,3),
     col=rgb(0,0,1,0.5),border = rgb(0,0,1,0.5),freq=FALSE,xaxt="n",xlab = "Number of Occurances",
     main = "Distribution of the Number of Observations of Deletion Alleles")
axis(1,at=c(0,1,2,3,4),labels = c(1,10,100,1000,10000))
hist(log10(sitefreq[names(sitefreq) %in% paste0("GScall_",impactCalls)]),add=TRUE,
     col=rgb(1,0,0,0.5),border = rgb(1,0,0,0.5),freq=FALSE)
legend("right",legend = c("Low","High"),fill = c("blue","red"),bty="n",title = expression(bold("Impact")))
mtext("Using GenomeSTRiP Calls",padj = 0.5)
dev.off()

region <- "MIDDLE_EAST"
ps <- popsGS[popsGS$region == region,"sample_accession"]
sitefreq.reg <- rowSums(GSfeats.delGenotype.noNA[,ps])
hist(sitefreq.reg,breaks = 765)



