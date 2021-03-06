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

## Population names, as ordered by popCols
popNames <- c("Mbuti","Biaka","Mandenka","Yoruba","San","Bantu (S. Africa)","Bantu (Kenya)","Colombian","Surui","Maya","Karitiana","Pima","Brahui",
              "Balochi","Hazara","Makrani","Sindhi","Pathan","Kalash","Burusho","Uygur","Cambodian","Japanese","Han","Yakut","Tujia","Yi","Miao",
              "Oroqen","Daur","Mongola","Hezhen","Xibo","Han (North)","Dai","Lahu","She","Naxi","Tu","French","Sardinian","Orcadian","Russian",
              "Italian","Tuscan","Basque","Adygei","Druze","Bedouin","Palestinian","Mozabite","Melanesian","Papuan")


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
plot(pcGS$x[,7],pcGS$x[,8],pch=20,col=samplePopColsGS)

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
       fill=c("darkorange","hotpink","green","blue","red","yellow","purple"),xpd=NA,
       title = expression(bold("Sample Region")),bty="n",cex=1.8)

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
  
pdf("Figures/GSpca_regions1.pdf",10,12)
par(oma=c(1,2,3,1))
layout(rbind(c(1,1,1,2),c(3,3,3,4),c(5,5,5,6)))
nam <- c("A. Africa","B. America", "C. Oceania")
reg <- c("AFRICA","AMERICA","OCEANIA")
for (i in 1:length(reg)){
  r <- reg[i]
  ps <- popsGS[popsGS$region == r,"population"]
  vars = round(pcGS.reg[[r]]$sdev^2/sum(pcGS.reg[[r]]$sdev^2) * 100)
  par(mar=c(4,5,4,0))
  plot(pcGS.reg[[r]]$x[,1],pcGS.reg[[r]]$x[,2],pch=20,col=as.factor(ps),
       xlab = paste0("PC1 (",vars[1],"%)"),ylab = paste0("PC2 (",vars[2],"%)"),
       main = nam[i],cex.axis=1.5,cex.lab=1.5,cex.main=1.8)
  par(mar=c(0,2,0,0))
  plot(NA,axes=FALSE,ylim=c(0,1),xlim=c(0,1),xlab="",ylab = "")
  if (r == "AFRICA"){
    legend(0,0.5,yjust = 0.5,legend = c("Mbuti","Biaka","Mandenka","Yoruba","San","Bantu (S. Africa)","Bantu (Kenya)"),col = as.factor(unique(ps)),bty="n",pch=20,xpd=NA,
           title = expression(bold("Population")),cex=2)
  } else {
    legend(0,0.5,yjust = 0.5,legend = unique(ps),col = as.factor(unique(ps)),bty="n",pch=20,xpd=NA,
         title = expression(bold("Population")),cex=2)
  }
}
title("Principal Component Analysis of GenomeSTRiP Deletion Genotypes by Region",outer = TRUE,cex.main=2)
mtext("Regions with Differentiation",outer = TRUE,cex=1.5,padj = 0.6)
dev.off()

pdf("Figures/GSpca_regions2.pdf",12,10)
par(mar=c(4,5,4,4),oma=c(1,1,5,1))
layout(rbind(c(1,2),c(3,4)))
nam <- c("A. Central & South Asia","B. East Asia","C. Europe","D. Middle East")
reg <- c("CENTRAL_SOUTH_ASIA","EAST_ASIA","EUROPE","MIDDLE_EAST")
for (i in 1:length(reg)){
  ps <- popsGS[popsGS$region == reg[i],"population"]
  r <- reg[i]
  vars = round(pcGS.reg[[r]]$sdev^2/sum(pcGS.reg[[r]]$sdev^2) * 100)
  plot(pcGS.reg[[r]]$x[,1],pcGS.reg[[r]]$x[,2],pch=20,col=popCols[ps],
       xlab = paste0("PC1 (",vars[1],"%)"),ylab = paste0("PC2 (",vars[2],"%)"),
       main = nam[i],cex.axis=1.5,cex.lab=1.5,cex.main=1.8)
}
title("Principal Component Analysis of GenomeSTRiP Deletion Genotypes by Region",outer = TRUE,cex.main=2)
mtext("Regions with no Differentiation",outer = TRUE,cex=1.5)
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

pdf("Figures/GSburden_bp.pdf",15,10)
par(oma=c(0,0,0,18),mar=c(4,5,4,0))
barplot(delInfo$total,border=delInfo$popCol,col=delInfo$popCol,main = "Total Number of Deleted Bases",
        xlab = "Population",ylab = "Base Pairs",cex.lab=1.5,cex.main=2,cex.axis=1.5)
axis(1,c(0,cumsum(920/765*table(delInfo$region))),rep("",8))
text(cumsum(920/765*table(delInfo$region)) - 920/765*table(delInfo$region)/2,-250000,
     c("Africa","America","South Central Asia", "East Asia","Europe","Middle East","Oceania"),
     xpd=NA,cex=1.2) 
legend(1150,3.4E6,popNames[c(1:21,48:53)],col = popCols[c(1:21,48:53)],pch=20,bty = 'n',xpd=NA,yjust = 0.5,cex=1.2)
legend(1380,3.5E6,popNames[22:47],col = popCols[22:47],pch=20,bty = 'n',xpd=NA,yjust = 0.5,cex=1.2)
text(1340,6.6E6,expression(bold("Population")),cex=1.5,xpd=NA)
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
        main = "A. Total",xlab = "Population",ylab = "Count",cex.lab=1.5,cex.main=2,cex.axis=1.5)
axis(1,c(0,cumsum(920/765*table(delInfo.count$region))),rep("",8))
text(cumsum(920/765*table(delInfo.count$region)) - 920/765*table(delInfo.count$region)/2,-50,
     c("Africa","America","South Central Asia", "East Asia","Europe","Middle East","Oceania"),
     xpd=NA,cex=1.2) 

barplot(delInfoHom$total,border=delInfoHom$popCol,col=delInfoHom$popCol,ylim = c(0,350),
        main = "B. Homozygotes",xlab = "Population",ylab = "Count",cex.lab=1.5,cex.main=2,cex.axis=1.5)
axis(1,c(0,cumsum(920/765*table(delInfoHom$region))),rep("",8))
text(cumsum(920/765*table(delInfoHom$region)) - 920/765*table(delInfoHom$region)/2,-15,
     c("Africa","America","South Central Asia", "East Asia","Europe","Middle East","Oceania"),
     xpd=NA,cex=1.1) 

barplot(delInfoHet$total,border=delInfoHet$popCol,col=delInfoHet$popCol,
        main = "C. Heterozygotes",xlab = "Population",ylab = "Count",cex.lab=1.5,cex.main=2,cex.axis=1.5)
axis(1,c(0,cumsum(920/765*table(delInfoHet$region))),rep("",8))
text(cumsum(920/765*table(delInfoHet$region)) - 920/765*table(delInfoHet$region)/2,-40,
     c("Africa","America","South Central Asia", "East Asia","Europe","Middle East","Oceania"),
     xpd=NA,cex=1.1) 

par(mar=c(1,1,1,1))
plot(NA,axes=FALSE,ylim=c(0,1),xlim=c(0,1),xlab="",ylab = "")
legend(0,0.49,popNames[c(1:21,48:53)],col = popCols[c(1:21,48:53)],pch=19,bty = 'n',xpd=NA,yjust = 0.5,cex=1.5)
legend(0.5,0.5,popNames[22:47],col = popCols[22:47],pch=19,bty = 'n',xpd=NA,yjust = 0.5,cex=1.5)
text(0.4,0.78,expression(bold("Population")),cex=1.8)

title("Number of GenomeSTRiP Deletions in Different Individuals",outer = TRUE,cex.main=2.5)
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
par(oma=c(1,1,3,14),mar=c(5,5,1,1))
layout(rbind(c(1,1,1,2)))
co <- c("blue","red")
barplot(delInfoMet$total,border=co[as.factor(delInfoMet$metric)],
        col=co[as.factor(delInfoMet$metric)]
        ,xlab = "Population",ylab = "Number of Deleted Sites",ylim = c(0,1000),cex.lab=1.7,cex.main=2,cex.axis=1.7)
axis(1,c(0,cumsum(920/765*table(delInfoMet$region))),rep("",8))
text(cumsum(920/765*table(delInfoMet$region)) - 920/765*table(delInfoMet$region)/2,-20,
     c("Africa","America","South Central Asia", "East Asia","Europe","Middle East","Oceania"),
     xpd=NA,cex=1.5) 

boxplot(total~metric,col=co,data = delInfoMet,ylim=c(0,1000),cex.axis=1.5,cex.lab=1.5)
title("Effect of Library Prepartion Method on Number of Called Deletions Using GenomeSTRiP",outer = TRUE,cex.main=2)
legend(2.8,500,yjust = 0.5,legend = c("PCR","PCR-Free"),fill = co,bty = "n",cex=2,xpd=NA)
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
     xlab = expression("F"[st]),col = "cornflowerblue",ylim = c(0,14000),cex.main=2.2,cex.lab=1.5,cex.axis=1.5)
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

dis <- read.csv("results/60-targets-associated_diseases.csv",stringsAsFactors = FALSE)
colnames(dis) <- c("Disease","FullName","pVal","NumTargets","TheraputicArea","TopTargets","AllTargets")
dis$TopTargets <- lapply(dis$TopTargets,function(x){strsplit(x," ")[[1]]})
dis$AllTargets <- lapply(dis$AllTargets,function(x){strsplit(x," ")[[1]]})


#### pairwise Fst ####
library(gplots)
pairwiseMetrics <- readPairwiseTables("results/GSpairwiseFst/pairwise.txt")
nam <- rownames(pairwiseMetrics$Fst)
nam[c(9,10,30,44,53)] <- c("Biaka","Mbuti","Bantu (SA)","Han (North)","Bantu (Kenya)")

col <- gsub("yellow","gold",unname(sapply(rownames(pairwiseMetrics$Fst),function(x){regCols[popSummaries[popSummaries$Population == x,"Region"]]})))

f <- function(){legend(-0.2,1,xjust = 0.5,yjust=0.5, legend = c("Africa","America","Central and S. Asia",
                                                                 "E. Asia","Europe","Middle East","Oceania"),
                       fill=c("darkorange","hotpink","green","blue","red","gold","purple"),xpd=NA,
                       title = expression(bold("Sample Region")),bty="n",cex=2.5)
}

pdf("Figures/GSfst_heatmap.pdf",30,20)
heatmap.2(as.matrix(pairwiseMetrics$Fst),symm = TRUE,col=colorRampPalette(c("blue","red"))(256),
          revC = FALSE,trace = "none",dendrogram = "row",symbreaks = FALSE,margins = c(14,14),
          symkey = FALSE,key.title = "",key.par = list(mar=c(8,4,3,16),cex.lab=2.5,cex.axis=2.2),
          density.info = "none",key.xlab = expression("F"[st]),cexRow = 2.3,cexCol = 2.3,
          labRow = nam,labCol = nam,
          lmat = rbind(c(0,0,4),c(3,1,2),c(0,0,5)),lhei = c(0.5,6,1),lwid = c(1.7,0.1,6),
          RowSideColors = col,colCol = col,colRow = col,extrafun = f)
title(expression("Pairwise F"[st]~"Values Across All Populations"),cex.main=3)
dev.off()

diveRsityStats <- read.table("results/GSpairwiseFst/std_stats.txt",header=TRUE)
globalStats <- diveRsityStats[dim(diveRsityStats)[1],]
diveRsityStats <- diveRsityStats[1:(dim(diveRsityStats)[1]-1),]

hist(diveRsityStats$Fst)

## Check against distance
getCoords <- function(x){
  x <- gsub("[ \t\n\r\v\f]","",x)
  sp <- strsplit(x,",")[[1]]
  north <- as.numeric(gsub("[SN]","",sp[1]))
  if (grepl("S",sp[1])){
    north <- -1* north
  }
  east <- as.numeric(gsub("[EW]","",sp[2]))
  if (grepl("W",sp[2])){
    east <- -1* east
  }
  return(c(north,east))
}

## Degrees to radians
deg2rad <- function(x){pi * x / 180}

## Haversine function
hav <- function(x){(1 - cos(x))/2}

## Function for great circle distance
greatDist <- function(x1,x2,r = 6.3781366E6){
  h <- hav(x2[1] - x1[1]) + cos(x1[1])*cos(x2[1])*hav(x2[2] - x1[2])
  return(2 * r * asin(sqrt(h)))
}

loc <- read.table("meta/location",sep="\t",stringsAsFactors = FALSE)
colnames(loc) <- c("Country","Coords","Population")

t <- sapply(loc$Coords,getCoords)
loc$lat <- t[1,]
loc$long <- t[2,]

dists <- apply(loc,1,function(x1){apply(loc,1,function(x2){greatDist(deg2rad(as.numeric(x1[c(4,5)])),deg2rad(as.numeric(x2[c(4,5)])))})})
colnames(dists) <- gsub("-","",loc$Population)
rownames(dists) <- colnames(dists)
diag(dists) <- NA

dists <- dists[colnames(pairwiseMetrics$Fst),colnames(pairwiseMetrics$Fst)]

p <- xy.coords(dists/1000,as.matrix(pairwiseMetrics$Fst))
di <- p$x[!is.na(p$x)]
fst <- p$y[!is.na(p$y)]
fit <- lm(fst ~ di)
s <- summary(fit)

## Get colours for points of the same region (turned out to have no pattern)
tt <- gsub("HanNChina","Han-NChina",colnames(dists))
co <- sapply(tt,function(x){sapply(tt,function(y){if(popRegs[x] == popRegs[y]){return(regCols[popRegs[x]])}else{return("black")}})})

pdf("Figures/GSfst_dists.pdf",12,8)
par(mar=c(4,5,4,2),oma=c(0,0,0,0))
plot(di,fst,pch=20,col="cornflowerblue",
     xlab = "Distance (km)",ylab = expression("F"[st]),
     main = expression("Relationship between Distance and F"[st]),
     xlim = c(0,20000),ylim = c(0,0.4),bty="n", xaxs="i", yaxs="i",xpd=NA,
     cex.main=2,cex.lab=1.5,cex.axis=1.5)
abline(fit$coefficients[1],fit$coefficients[2],lwd=1.5)
abline(fit$coefficients[1] + qnorm(0.995) * s$coefficients[1,2],fit$coefficients[2] + qnorm(0.995) * s$coefficients[2,2],lwd=1.5,lty=2)
abline(fit$coefficients[1] - qnorm(0.995) * s$coefficients[1,2],fit$coefficients[2] - qnorm(0.995) * s$coefficients[2,2],lwd=1.5,lty=2)
dev.off()

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
     main = "Distribution of the Number of Observations of Deletion Alleles",cex.main=2,cex.axis=1.5,cex.lab=1.5)
axis(1,at=c(0,1,2,3,4),labels = c(1,10,100,1000,10000))
hist(log10(sitefreq[names(sitefreq) %in% paste0("GScall_",impactCalls)]),add=TRUE,
     col=rgb(1,0,0,0.5),border = rgb(1,0,0,0.5),freq=FALSE)
legend("right",legend = c("Low","High"),fill = c("blue","red"),bty="n",title = expression(bold("Impact")),cex=1.5)
mtext("Using GenomeSTRiP Calls",padj = 0.5,cex=1.5)
dev.off()

nonImp <- hist(log10(sitefreq[!names(sitefreq) %in% paste0("GScall_",impactCalls)]),plot=FALSE,breaks = seq(0,4,0.2))
imp <- hist(log10(sitefreq[impactCalls]),plot=FALSE,breaks = seq(0,4,0.2))

probs <- rbind(nonImp$counts/sum(nonImp$counts),imp$counts/sum(imp$counts))

pdf("Figures/GSsitefreq.pdf",12,10)
par(mar=c(5,5,6,6))
barplot(probs,beside = TRUE,space = c(0,0),col=c("blue","red"),ylim = c(0,0.5),xlim=c(0,4),width = 0.1,
        xlab = "Number of Occurances",ylab = "Frequency",
        main = "Distribution of the Frequency of Deletion Alleles",cex.main=2,cex.axis=1.5,cex.lab=1.5)
axis(1,at=seq(0,4,0.2),labels = rep("",21),tck=-0.01)
axis(1,at=c(0,1,2,3,4),labels = c(1,10,100,1000,10000),cex.axis=1.5)
legend(4.1,0.25,legend = c("Low","High"),fill = c("blue","red"),bty="n",title = expression(bold("Impact")),cex=1.5,xpd=NA)
mtext("Using GenomeSTRiP Calls",padj = 0,cex=1.5)
dev.off()

region <- "MIDDLE_EAST"
ps <- popsGS[popsGS$region == region,"sample_accession"]
sitefreq.reg <- rowSums(GSfeats.delGenotype.noNA[,ps])
hist(sitefreq.reg,breaks = 765)

### Number of variants against distance from Africa - not published
# Requires distance functions etc. loaded from above and delInfo.count
EAfrica.coords <- deg2rad(c(0,36))
loc$dist <- apply(loc,1,function(x){greatDist(deg2rad(as.numeric(x[c(4,5)])),EAfrica.coords)})
loc$count <- sapply(loc$Population,function(x){mean(delInfo.count[delInfo.count$population == x,"total"])})

delInfo.count$dist <- sapply(delInfo.count$population,function(x){loc[loc$Population == x,"dist"]})
plot(delInfo.count$dist,delInfo.count$total,col=delInfo.count$popCol,pch=20)

fit <- lm(total ~ dist,data = delInfo.count[!delInfo.count$region == "AFRICA",])

# Call:
#   lm(formula = total ~ dist, data = delInfo.count)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -78.58 -26.08  -5.12  19.09 441.42 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  5.715e+02  3.279e+00  174.27  < 2e-16 ***
#   dist        -3.335e-06  4.537e-07   -7.35 4.38e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 44.08 on 915 degrees of freedom
# Multiple R-squared:  0.05576,	Adjusted R-squared:  0.05472 <- not found to explain much variance in this data
# F-statistic: 54.03 on 1 and 915 DF,  p-value: 4.38e-13

