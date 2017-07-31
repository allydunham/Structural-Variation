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

## Script containing code used to generate the various QC plots in the 10X quality control section
## ROC curves are found in the main combinedQC.R script
setwd('~/Project/')
source('scripts/QC_functions.R')

metrics <- c("Accuracy","Kappa","F1","FDR","FNR","FPR")
nam <- c("Random Forest","SVM","Boosted Tree","Logit","Adaboost","KNN","Ensemble")

### CV test ###
load("results/CVtests.Rdata")
cla <- list(rf=rf_CVtest,svm=svm_CVtest,boost=boost_CVtest,logit=logit_CVtest,adabost=adaboost_CVtest,knn=knn_CVtest,ens=ens_CVtest)
vals <- sapply(metrics,function(m){sapply(cla,function(y){y[m,]})})

pdf("Figures/CVtest_final.pdf",18,10)
par(mfrow=c(2,3),oma=c(2,2,3,2))
boxplot(vals[,"Accuracy"],main="Accuracy",names = nam,cex.main=1.5,ylim = c(0.8,1),col="cornflowerblue")
boxplot(vals[,"Kappa"],main="Kappa",names = nam,cex.main=1.5,ylim = c(0.2,0.8),col="cornflowerblue")
boxplot(vals[,"F1"],main="F1",names = nam,cex.main=1.5,ylim = c(0.8,1),col="cornflowerblue")
boxplot(vals[,"FDR"],main="FDR",names = nam,cex.main=1.5,ylim = c(0,0.2),col="cornflowerblue")
boxplot(vals[,"FNR"],main="FNR",names = nam,cex.main=1.5,ylim = c(0,0.3),col="cornflowerblue")
boxplot(vals[,"FPR"],main="FPR",names = nam,cex.main=1.5,ylim = c(0,1),col="cornflowerblue")
title("Classifier Assesment over 6 Metrics",outer = TRUE,cex.main=2.5)
dev.off()


### CV test - no Unknowns ###
load("results/CVtests_noUnknowns.Rdata")
cla <- list(rf=rf_CVtest,svm=svm_CVtest,boost=boost_CVtest,ada=adaboost_CVtest)
vals <- sapply(metrics,function(m){sapply(cla,function(y){y[m,]})})
vals <- vals[-4,]

pdf("Figures/CVtest_noUnknowns_final.pdf",18,10)
par(mfrow=c(2,3),oma=c(2,2,3,2))
boxplot(vals[,"Accuracy"],main="Accuracy",names = nam[1:3],cex.main=2.5,ylim = c(0.8,1),col="cornflowerblue",cex.labs=2.5,cex.axis=2)
boxplot(vals[,"Kappa"],main="Kappa",names = nam[1:3],cex.main=2.5,ylim = c(0.2,0.8),col="cornflowerblue",cex.labs=2.5,cex.axis=2)
boxplot(vals[,"F1"],main="F1",names = nam[1:3],cex.main=2.5,ylim = c(0.8,1),col="cornflowerblue",cex.labs=2.5,cex.axis=2)
boxplot(vals[,"FDR"],main="FDR",names = nam[1:3],cex.main=2.5,ylim = c(0,0.2),col="cornflowerblue",cex.labs=2.5,cex.axis=2)
boxplot(vals[,"FNR"],main="FNR",names = nam[1:3],cex.main=2.5,ylim = c(0,0.3),col="cornflowerblue",cex.labs=2.5,cex.axis=2)
boxplot(vals[,"FPR"],main="FPR",names = nam[1:3],cex.main=2.5,ylim = c(0,1),col="cornflowerblue",cex.labs=2.5,cex.axis=2)
title("Classifier Assesment over 6 Metrics with Ambiguous Training Data Excluded",outer = TRUE,cex.main=3)
dev.off()


### CV test - orig untuned###
## Must previously have loaded an adaboost result - it being a different size forces list into right format
load("results/CVtests_orig.Rdata")
cla <- list(rf=rf_CVtest,svm=svm_CVtest,ada=adaboost_CVtest)
vals <- sapply(metrics,function(m){sapply(cla,function(y){y[m,]})})
vals <- vals[-3,]

pdf("Figures/CVtest_orig_final.pdf",18,10)
par(mfrow=c(2,3),oma=c(2,2,3,2))
boxplot(vals[,"Accuracy"],main="Accuracy",names = nam[1:2],cex.main=2.5,ylim = c(0.8,1),col="cornflowerblue",cex.labs=2.5,cex.axis=2)
boxplot(vals[,"Kappa"],main="Kappa",names = nam[1:2],cex.main=2.5,ylim = c(0,0.5),col="cornflowerblue",cex.labs=2.5,cex.axis=2)
boxplot(vals[,"F1"],main="F1",names = nam[1:2],cex.main=2.5,ylim = c(0.8,1),col="cornflowerblue",cex.labs=2.5,cex.axis=2)
boxplot(vals[,"FDR"],main="FDR",names = nam[1:2],cex.main=2.5,ylim = c(0,0.2),col="cornflowerblue",cex.labs=2.5,cex.axis=2)
boxplot(vals[,"FNR"],main="FNR",names = nam[1:2],cex.main=2.5,ylim = c(0,0.3),col="cornflowerblue",cex.labs=2.5,cex.axis=2)
boxplot(vals[,"FPR"],main="FPR",names = nam[1:2],cex.main=2.5,ylim = c(0,1),col="cornflowerblue",cex.labs=2.5,cex.axis=2)
title("Initial Classifier Assesment over 6 Metrics Using the Raw Longranger Data Without Tuning",outer = TRUE,cex.main=3)
dev.off()

### CV test - orig tuned###
load("results/CVtests_orig_tuned.Rdata")
cla <- list(rf=rf_CVtest,svm=svm_CVtest,boost=boost_CVtest,logit=logit_CVtest,adabost=adaboost_CVtest,knn=knn_CVtest,ens=ens_CVtest)
vals <- sapply(metrics,function(m){sapply(cla,function(y){y[m,]})})

pdf("Figures/CVtest_orig_tuned_final.pdf",18,10)
par(mfrow=c(2,3),oma=c(2,2,3,2))
boxplot(vals[,"Accuracy"],main="Accuracy",names = nam,cex.main=1.5,ylim = c(0.6,1),col="cornflowerblue")
boxplot(vals[,"Kappa"],main="Kappa",names = nam,cex.main=1.5,ylim = c(0.2,0.6),col="cornflowerblue")
boxplot(vals[,"F1"],main="F1",names = nam,cex.main=1.5,ylim = c(0.8,1),col="cornflowerblue")
boxplot(vals[,"FDR"],main="FDR",names = nam,cex.main=1.5,ylim = c(0,0.2),col="cornflowerblue")
boxplot(vals[,"FNR"],main="FNR",names = nam,cex.main=1.5,ylim = c(0,0.5),col="cornflowerblue")
boxplot(vals[,"FPR"],main="FPR",names = nam,cex.main=1.5,ylim = c(0,1),col="cornflowerblue")
title("Initial Classifier Assesment over 6 Metrics Using the Raw Longranger Data After Tuning",outer = TRUE,cex.main=2.5)
dev.off()
