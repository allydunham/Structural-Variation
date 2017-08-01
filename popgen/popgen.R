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

## Script contining code to perform all population genetic analyses of breakdancer data and produce appropriate plots
setwd('~/Project/')
source("scripts/importPopgen.R")

############ Import VEP ############
VEPimpacts <- read.table("results/mergedBDcalls_effects.vep.txt",header = TRUE,comment.char = '',stringsAsFactors = FALSE)
cons <- strsplit(VEPimpacts$Consequence,",")
t <- unique(unlist(cons))
for (i in t){
  VEPimpacts[,i] <- as.numeric(sapply(cons,function(x){any(i %in% x)}))
}
VEPimpacts <- VEPimpacts[,-7]

sortExtras <- function(x){
  t <- strsplit(x,"=")
  n <- sapply(t,function(u){u[1]})
  v <- sapply(t,function(u){u[2]})
  names(v) <- n
  return(v)
}

ex <- strsplit(VEPimpacts$Extra,';')
ex <- lapply(ex,sortExtras)

t <- unique(unlist(sapply(ex,names)))
for (i in t){
  VEPimpacts[,i] <- sapply(ex,function(x){if (i %in% names(x)){x[i]}else{NA}})
}
VEPimpacts <- VEPimpacts[,-13]

VEPimpacts$call <- NA
for (i in 1:dim(BDdels)[1]){
  VEPimpacts[VEPimpacts$X.Uploaded_variation == paste(BDdels[i,1],BDdels[i,2],"deletion",sep = "_"),"call"] <- i
}

VEPimpactsRed <- VEPimpacts[!is.na(VEPimpacts$SYMBOL),]
genes <- unique(VEPimpactsRed$SYMBOL)

## Save data to binary here to save time on lengthy import process
save(BDdels,BDfeats,pops,popsBD,popSummaries,tenXfeats,VEPimpacts,VEPimpactsRed,genes,popCols,
     popRegs,populations,regCols,samplePopColsBD,samplePopColsTenX,sampleRegColsBD,
     sampleRegColsTenX,file = "results/popgenImport.Rdata")

############ POPGEN ANALYSIS #############
## (can start here with preloaded data)
setwd("~/Project")
load("results/popgenImport.Rdata")
source("scripts/popgenFunctions.R")

## Population names, as ordered by popCols
popNames <- c("Mbuti","Biaka","Mandenka","Yoruba","San","Bantu (S. Africa)","Bantu (Kenya)","Colombian","Surui","Maya","Karitiana","Pima","Brahui",
              "Balochi","Hazara","Makrani","Sindhi","Pathan","Kalash","Burusho","Uygur","Cambodian","Japanese","Han","Yakut","Tujia","Yi","Miao",
              "Oroqen","Daur","Mongola","Hezhen","Xibo","Han (North)","Dai","Lahu","She","Naxi","Tu","French","Sardinian","Orcadian","Russian",
              "Italian","Tuscan","Basque","Adygei","Druze","Bedouin","Palestinian","Mozabite","Melanesian","Papuan")

tenXfreqs <- rowSums(tenXfeats)/dim(tenXfeats)[2]
tenXprops <- colSums(tenXfeats)/dim(tenXfeats)[1]

BDfreqs <- rowSums(BDfeats)/dim(BDfeats)[2]
BDprops <- colSums(BDfeats)/dim(BDfeats)[1]

##ks test suggests proportions in each sig diff, san seem to have more

############ Hamming Distance ############
hamm <- sapply(1:765,function(x){sapply(1:765,function(y){
  sum(BDfeats[,x] == BDfeats[,y])/30884
})})
rownames(hamm) <- colnames(BDfeats)
colnames(hamm) <- colnames(BDfeats)

############ Clustering ############
## recapitulates a large amount of expected population structure - varies slightly on the distance measure
# binary dist seems most natural and creates a clear africa/rest of world split
popClust <- hclust(dist(t(BDfeats),"binary"))
popDend <- as.dendrogram(popClust)

catNodePop <- function(x){
  if (is.leaf(x)){
    a <- attributes(x)
    r <- pops[pops$sample == a$label,"population"]
    c <- unname(popCols[r[1]])
    attr(x, "nodePar") <- c(a$nodePar, lab.col = c)
    attr(x, "label") <- r[1]
  }
  x
}

popDend <- dendrapply(popDend,catNodePop)
plot(popDend)

library(ape)
phylo <- as.phylo(popClust)
tipCols <- samplePopColsBD[phylo$tip.label]
plot(phylo,type = "unrooted",use.edge.length = FALSE,tip.color = tipCols,lab4ut="axial",no.margin = TRUE)


## 10X results
popClust10X <- hclust(dist(t(tenXfeats),"binary"))
popDend10X <- as.dendrogram(popClust10X)
popDend10X <- dendrapply(popDend10X,catNodePop)
plot(popDend10X)


############ PCA ############
pcBD <- prcomp(t(BDfeats),tol=10^-10)
vars <- round(pcBD$sdev^2/sum(pcBD$sdev^2) * 100)
pdf("Figures/BDpca.pdf",10,14)
par(oma=c(4,4,4,4))
layout(rbind(c(1,2,3),c(0,4,5),c(7,0,6)))
par(mar=c(1,1,1,1),xpd=NA,cex.lab=1.5)
plot(pcBD$x[,2],pcBD$x[,1],pch=20,col=samplePopColsBD,xlab = paste0("PC2 (",vars[2],"%)"),ylab = paste0("PC1 (",vars[1],"%)"))
plot(pcBD$x[,3],pcBD$x[,1],pch=20,col=samplePopColsBD,yaxt="n",xaxt="n",xlab = "",ylab = "")
plot(pcBD$x[,4],pcBD$x[,1],pch=20,col=samplePopColsBD,yaxt="n",xaxt="n",xlab = "",ylab = "")

plot(pcBD$x[,3],pcBD$x[,2],pch=20,col=samplePopColsBD,xlab = paste0("PC3 (",vars[3],"%)"),ylab = paste0("PC2 (",vars[2],"%)"))
plot(pcBD$x[,4],pcBD$x[,2],pch=20,col=samplePopColsBD,yaxt="n",xaxt="n",xlab = "",ylab = "")

plot(pcBD$x[,4],pcBD$x[,3],pch=20,col=samplePopColsBD,xlab = paste0("PC4 (",vars[4],"%)"),ylab = paste0("PC3 (",vars[3],"%)"))

plot(NA,axes=FALSE,ylim=c(0,1),xlim=c(0,1),xlab="",ylab = "")
legend(0.5,0.5,xjust = 0.5,yjust=0.5, legend = c("Africa","America","Central and S. Asia",
                                    "E. Asia","Europe","Middle East","Oceania"),
       col=c("darkorange","hotpink","green","blue","red","yellow","purple"),pch=20,xpd=NA,
       title = expression(bold("Sample Region")),bty="n",cex=1.5)

title("Principal Component Analysis of BreakDancer Deletion Calls",outer = TRUE,cex.main=2)
dev.off()
layout(1)

library(rgl)
plot3d(pcBD$x[,1],pcBD$x[,2],pcBD$x[,3],col=samplePopColsBD)


# test other library stats - none found to be important
varBD <- pops[pops$sample %in% colnames(BDfeats),c(1,3)]
varBD <- varBD[!duplicated(varBD),2]
plot(pcBD$x[,1],pcBD$x[,2],pch=20,col=as.factor(varBD))


pcTenX <- prcomp(t(tenXfeats))
plot(pcTenX$x[,1],pcTenX$x[,2],pch=20,col=samplePopColsTenX)


############ Differential Analysis ############
Mode <- function(v) {
  u <- unique(v)
  u[which.max(tabulate(match(v, u)))]
}

nClusts <- 100
delClusts <- kmeans(as.matrix(BDfeats),nClusts)

BDfeats.reduced <- sapply(popNames,function(pop){sapply(1:nClusts,function(clu){Mode(as.matrix(BDfeats[which(delClusts$cluster == clu),
                                                                                             which(colnames(BDfeats) %in% pops[pops$population == pop,"sample"])
                                                                                                   ]))})})

heatmap(as.matrix(BDfeats.reduced),scale="none",col=colorRampPalette(c("blue","red"))(256))

###### Determine Zygosity ######
BDcurate <- read.table("data/randomSampleBD_manuallychecked.txt",fill=TRUE,stringsAsFactors = FALSE)
colnames(BDcurate) <- c("chrom1","pos1","map1","chr2","pos2","map2","check","zyg","type","length","score","reads","lib","CN","ID")
BDcurate <- BDcurate[BDcurate$check %in% c("T","T?") & grepl("hom|het",BDcurate$zyg) &BDcurate$type == "DEL",]

BDcurate$sample <- sapply(strsplit(BDcurate$ID,"[_.]"),function(x){x[4]})
BDcurate$reg <- sapply(BDcurate$sample,function(x){unique(pops[pops$sample == x,"region"])})
BDcurate$zyg <- gsub("[?()]","",BDcurate$zyg)
  
zyg <- sapply(unique(pops$region),function(x){t <- table(BDcurate[BDcurate$reg == x,"zyg"]); 1 + t["hom"]/sum(t)})

############ Genetic Burden ############
delInfo <- data.frame(total=colSums(BDfeats),
                      region=sapply(colnames(BDfeats),function(x){r <- pops[pops$sample == x,"region"];return(r[1])}),
                      population=sapply(colnames(BDfeats),function(x){r <- pops[pops$sample == x,"population"];return(r[1])}),
                      regCol=sampleRegColsBD,
                      popCol=samplePopColsBD,stringsAsFactors = FALSE)

delInfo <- delInfo[order(delInfo$region,delInfo$population),]

pdf("Figures/BDburden.pdf",16,10)
par(oma=c(1,1,1,17),mar=c(4,5,4,4))
barplot(delInfo$total,border=delInfo$popCol,col=delInfo$popCol,main = "Total Number of Deleted Sites",
        xlab = "Population",ylab = "Deleted Sites",cex.main=2,cex.axis=1.5,cex.lab=1.5)
axis(1,c(0,cumsum(920/765*table(delInfo$region))),rep("",8))
text(cumsum(920/765*table(delInfo$region)) - 920/765*table(delInfo$region)/2,-100,
     c("Africa","America","South Central Asia", "East Asia","Europe","Middle East","Oceania"),
     xpd=NA,cex=1.3) 
legend(1000,1150,popNames[c(1:21,48:53)],col = popCols[c(1:21,48:53)],pch=20,bty = 'n',xpd=NA,yjust = 0.5,cex=1.2)
legend(1200,1150,popNames[22:47],col = popCols[22:47],pch=20,bty = 'n',xpd=NA,yjust = 0.5,cex=1.2)
text(1170,2250,expression(bold("Population")),xpd=NA)
dev.off()

## Only common deletions
BDfeats.noSingles <- BDfeats[rowSums(BDfeats) > 100,]
n <- sapply(strsplit(rownames(BDfeats.noSingles),"_"),function(x){as.numeric(x[2])})
delInfo.noSingles <- data.frame(total=colSums(BDfeats.noSingles),
                      region=sapply(colnames(BDfeats.noSingles),function(x){r <- pops[pops$sample == x,"region"];return(r[1])}),
                      population=sapply(colnames(BDfeats.noSingles),function(x){r <- pops[pops$sample == x,"population"];return(r[1])}),
                      regCol=sampleRegColsBD,
                      popCol=samplePopColsBD,
                      stringsAsFactors = FALSE)

delInfo.noSingles <- delInfo.noSingles[order(delInfo.noSingles$region,delInfo.noSingles$population),]

par(oma=c(0,0,0,12))
barplot(delInfo.noSingles$total,border=delInfo.noSingles$popCol,col=delInfo.noSingles$popCol,
        main = "Number of Deleted Bases From Common Deletions",xlab = "Population",ylab = "Base Pairs")
axis(1,c(0,cumsum(920/765*table(delInfo.noSingles$region))),rep("",8))
text(cumsum(920/765*table(delInfo.noSingles$region)) - 920/765*table(delInfo.noSingles$region)/2,c(0.5,0.3,0.5,0.3,0.5,0.3,0.5)*-10^6,
     c("Africa","America","South Central Asia", "East Asia","Europe","Middle East","Oceania"),
     xpd=NA) 
legend(1000,2.5E6,names(popCols)[c(1:21,48:53)],col = popCols[c(1:21,48:53)],pch=20,bty = 'n',xpd=NA,yjust = 0.5,cex=0.8)
legend(1300,2.5E6,names(popCols)[22:47],col = popCols[22:47],pch=20,bty = 'n',xpd=NA,yjust = 0.5,cex=0.8)

## Personal deletions
BDfeats.singles <- BDfeats[rowSums(BDfeats) < 5,]
n <- sapply(strsplit(rownames(BDfeats.singles),"_"),function(x){as.numeric(x[2])})
delInfo.singles <- data.frame(total=colSums(BDfeats.singles),
                                region=sapply(colnames(BDfeats.singles),function(x){r <- pops[pops$sample == x,"region"];return(r[1])}),
                                population=sapply(colnames(BDfeats.singles),function(x){r <- pops[pops$sample == x,"population"];return(r[1])}),
                                regCol=sampleRegColsBD,
                                popCol=samplePopColsBD,
                                stringsAsFactors = FALSE)

delInfo.singles <- delInfo.singles[order(delInfo.singles$region,delInfo.singles$population),]

par(oma=c(0,0,0,12))
barplot(delInfo.singles$total,border=delInfo.singles$popCol,col=delInfo.singles$popCol,
        main = "Number of Deleted Bases From Personal Deletions",xlab = "Population",ylab = "Base Pairs")
axis(1,c(0,cumsum(920/765*table(delInfo.singles$region))),rep("",8))
text(cumsum(920/765*table(delInfo.singles$region)) - 920/765*table(delInfo.singles$region)/2,c(0.5,0.3,0.5,0.3,0.5,0.3,0.5)*-10^5,
     c("Africa","America","South Central Asia", "East Asia","Europe","Middle East","Oceania"),
     xpd=NA) 
legend(1000,2E5,names(popCols)[c(1:21,48:53)],col = popCols[c(1:21,48:53)],pch=20,bty = 'n',xpd=NA,yjust = 0.5,cex=0.8)
legend(1300,2E5,names(popCols)[22:47],col = popCols[22:47],pch=20,bty = 'n',xpd=NA,yjust = 0.5,cex=0.8)

## Other metrics
delInfo.met <- data.frame(total=colSums(BDfeats),
                          region=sapply(colnames(BDfeats),function(x){r <- pops[pops$sample == x,"region"];return(r[1])}),
                          population=sapply(colnames(BDfeats),function(x){r <- pops[pops$sample == x,"population"];return(r[1])}),
                          regCol=sampleRegColsBD,
                          popCol=samplePopColsBD,
                          metric=sapply(colnames(BDfeats),function(x){r <- pops[pops$sample == x,"library_type"];return(r[1])}),
                          stringsAsFactors = FALSE)

delInfo.met <- delInfo.met[order(delInfo.met$region,delInfo.met$population),]

pdf("Figures/BDburden_libBatch.pdf",14,10)
par(oma=c(1,1,4,14),mar=c(5,5,1,1))
layout(rbind(c(1,1,1,2)))
co <- c("blue","red")
barplot(delInfo.met$total,border=co[as.factor(delInfo.met$metric)],
        col=co[as.factor(delInfo.met$metric)],
        xlab = "Population",ylab = "Number of Deleted Sites",cex.axis=1.5,cex.lab=1.7)
axis(1,c(0,cumsum(920/765*table(delInfo.met$region))),rep("",8),cex.axis=1.5)
text(cumsum(920/765*table(delInfo.met$region)) - 920/765*table(delInfo.met$region)/2,-100,
     c("Africa","America","South Central Asia", "East Asia","Europe","Middle East","Oceania"),
     xpd=NA,cex=1.5) 

boxplot(total ~ metric,data = delInfo.met,col = co,cex.lab=1.7,cex.axis=1.5)
legend(2.8,1250,yjust = 0.5,legend = c("PCR","PCR Free"),fill=co,bty = "n",cex = 2,xpd=NA)
title("Effect of Library Prepartion Method on Number of Called Deletions Using BreakDancer",outer = TRUE,cex.main=2)
dev.off()
layout(1)

############ Fst ############
getAlleleFreqs <- function(pop,factor=2){
  s <- unique(popsBD[popsBD$population == pop,"sample"])
  return(factor*rowSums(BDfeats[,s[s %in% colnames(BDfeats)]])/(2*length(s)))
}

nPops <- length(populations)
popWeights <- sapply(populations,function(x){popSummaries[popSummaries$Population == x,"N"]})
popWeights <- popWeights/sum(popWeights)

alleleFreqs.het <- sapply(populations,getAlleleFreqs,factor=1)
alleleFreqs.hom <- sapply(populations,getAlleleFreqs,factor=2)

## use ratio from currated set
alleleFreqs.ratio <- sapply(populations,getAlleleFreqs,factor=1.53)

p.het <- apply(alleleFreqs.het,1,function(x){sum(popWeights * x)})
p.hom <- apply(alleleFreqs.hom,1,function(x){sum(popWeights * x)})
p.ratio <- apply(alleleFreqs.ratio,1,function(x){sum(popWeights * x)})

sig.het <- sapply(1:dim(alleleFreqs.het)[1],function(i){sum(alleleFreqs.het[i,]^2*popWeights) - p.het[i]^2})
sig.hom <- sapply(1:dim(alleleFreqs.hom)[1],function(i){sum(alleleFreqs.hom[i,]^2*popWeights) - p.hom[i]^2})
sig.ratio <- sapply(1:dim(alleleFreqs.ratio)[1],function(i){sum(alleleFreqs.ratio[i,]^2*popWeights) - p.ratio[i]^2})

F_st.het <- sig.het/(p.het*(1-p.het))
F_st.hom <- sig.hom/(p.hom*(1-p.hom))
F_st.ratio <- sig.ratio/(p.ratio*(1-p.ratio))

## F_st measured directly from population allele frequencies. Problems:
## Problems:
# Currently weighted by true population size (incompletely) rather than effective
# accounting for heterozygosity in frequency calculation is crude

hist(F_st.het)
hist(F_st.hom)

pdf("Figures/BDfstHist.pdf",12,8)
par(mar=c(4,4,4,4),oma=c(1,1,1,1),cex.lab=1.25)
hist(F_st.ratio,main = expression("Distribution of F"[st]~"for BreakDancer Deletion Calls"),
     xlab = expression("F"[st]),col = "cornflowerblue",cex.lab=1.25,cex.main=1.5)
dev.off()

getGenes <- function(x){
  return(unique(VEPimpactsRed[VEPimpactsRed$call == as.integer(strsplit(x,"_")[[1]][2]),"SYMBOL"]))
}

BDdels$ID <- paste0("del_",1:dim(BDdels)[1])

df_fst <- data.frame(call=names(F_st.het),
                     fst=F_st.het)
df_fst$genes <- lapply(names(F_st.het),getGenes)
df_fst$genes <- sapply(df_fst$genes,paste,collapse=",")
df_fst <- merge(df_fst,BDdels,by.x = "call",by.y = "ID")
df_fst <- df_fst[order(df_fst$fst,decreasing = TRUE),]


############ Site Frequency ############
sitefreq <- rowSums(BDfeats)
impactCalls <- unique(VEPimpactsRed[VEPimpactsRed$IMPACT %in% c("MODERATE","HIGH"),"call"])

pdf("Figures/BDsitefreq.pdf",12,10)
par(mar=c(5,5,3,1))
hist(log10(sitefreq[!names(sitefreq) %in% paste0("GScall_",impactCalls)]),xlim=c(0,3),ylim = c(0,3),
     col=rgb(0,0,1,0.5),border = rgb(0,0,1,0.5),freq=FALSE,xaxt="n",xlab = "Number of Occurances",
     main = "Distribution of the Number of Observations of Deletion Alleles",cex.main=2,cex.axis=1.5,cex.lab=1.5)
axis(1,at=c(0,1,2,3),labels = c(1,10,100,1000))
hist(log10(sitefreq[impactCalls]),add=TRUE,
     col=rgb(1,0,0,0.5),border = rgb(1,0,0,0.5),freq=FALSE)
legend("right",legend = c("Low","High"),fill = c("blue","red"),bty="n",title = expression(bold("Impact")),cex=1.5)
mtext("Using BreakDancer Calls",padj = 0.5,cex=1.5)
dev.off()

hist(sitefreq[impactCalls],breaks = 765,xlim = c(0,800))



############ Functional Analysis ############
t <- unique(VEPimpacts[grep("SLC",VEPimpacts$SYMBOL),"call"])
heatmap(alleleFreqs[t,],col=colorRampPalette(c("blue","white","red"))(256))

VEPimpacts.shared <- VEPimpactsRed[VEPimpactsRed$call %in% which(rowSums(BDfeats) > 30),]
sharedGenes <- unique(VEPimpacts.shared$SYMBOL)

VEPimpacts.sing <- VEPimpactsRed[duplicated(VEPimpactsRed[c("Location","BIOTYPE"),]),]

par(mar=c(4,16,1,1),oma=c(1,1,3,1))
barplot(table(VEPimpacts.sing$BIOTYPE),horiz = TRUE,las=1,
        names.arg = rev(c("Unprocessed Pseudogene","Unitary Pseudogene","Transcribed Unprocessed psudogene",
                      "Transcribed Unitary Pseudogene","Transcribed Processed Pseudogene","TR V Pseudogene", "TR V Gene",
                      "TR J Gene","TR D Gene","TR C Gene","TEC","sRNA","snRNA","snoRNA","Sense Overlapping",
                      "Sense Intronic", "scaRNA","rRNA","Ribozyme","retained Intron","Pseudogene","Protein Coding",
                      "Processed Transcript","Processed Pseudogene","Polymorphic Pseudogene","Nonsense Mediated Decay",
                      "Non Stop Decay","Misc RNA","miRNA","lincRNA","IG V Pseudogene","IG V Gene","IG J Pseudogene",
                      "IG J Gene","IG D Gene","IG C Pseudogene","IG C Gene","Bidirectional Promoter lncRNA",
                      "Antisense","3' Overlapping ncRNA")))

## Type of variants
par(oma=c(2,2,2,2),mar=c(13,4,4,4))
barplot(colSums(VEPimpacts[,13:32]),las=2)

############ Allele Freqs by gene ############
popSamples <- sapply(populations,function(x){popsBD[popsBD$population == x,"sample"]})
getGeneFreq <- function(gene,feat,vep){
  feats <- feat[paste0("del_",unique(vep[vep$SYMBOL == gene,"call"])),]
  return(sapply(populations,function(pop){sum(colSums(feats[,popSamples[[pop]]]) > 0)/length(popSamples[[pop]])}))
}

geneFreqs <- sapply(genes,getGeneFreq,feat=BDfeats,vep=VEPimpactsRed)

VEPimpacts.impact <- VEPimpacts.shared[VEPimpacts.shared$IMPACT %in% c("MODERATE","HIGH"),]
BDfeats.impact <- BDfeats[unique(VEPimpacts.impact$call),]
genes.impact <- unique(VEPimpacts.impact$SYMBOL)

geneFreqs.impact <- sapply(genes.impact,getGeneFreq,feat=BDfeats.impact,vep=VEPimpacts.impact)
heatmap(geneFreqs.impact[,colSums(geneFreqs.impact > 0.1) > 5 & colSums(geneFreqs.impact < 0.9) > 5],
        col=colorRampPalette(c("blue","red"))(256),scale="none")

##moderate/high impact deletions that are shared by at least 30 pops:
# GSTM1/2 in some non africans (detoxification inc of carcinogens etc)
# RNU7 snRNA in most non africans
# ORM1/2 in africans (protein acossiated with inflamation, potentially part of immuno supression)
# FOXP2 deletion found to be common in africa and uncommon outside (range of effects, only modifier, protein coding change and nonsense mediated decay most common)
# SLC24A2 deletion most common in Africa (also found elsewhere at lower rates) - light related Ca2+ exchanger in retinal/brain exchange (mostly intron variant)

heatmap(geneFreqs,col=colorRampPalette(c("blue","red"))(256),scale="none")

