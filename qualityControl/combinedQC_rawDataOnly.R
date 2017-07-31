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

## Equivalent script to combinedQC.R but set up to only using the raw data from longranger rather than extra statistics
## It is expected sections will be run interactively during an analysis rather than performing a function as a whole


setwd('~/Project/')
source('scripts/QC_functions.R')

## Import deletions 
dat <- list(dat=c("data/HGDP00819/processedPadded.tsv","data/HGDP01029/processedPadded.tsv"),
            man=c("data/HGDP00819/BD_2113899_HGDP00819_manual.txt","data/HGDP01029/BD_2079013_HGDP01029_manual.txt"))

dels <- rbindlist(lapply(1:length(dat),function(x){loadDels(dat$dat[x],dat$man[x],FALSE)}))
dels$RepTotal <- rowSums(dels[,c("LINEs","SINEs","LTRs","LowComplexityRepeats","SimpleRepeats","OtherRepeats")])
  
testData <- dels[!is.na(dels$Manual_check),]
testData$cat <- as.integer(testData$Manual_check %in% c("T","T?"))

############################################ Test Classifiers ########################################

## Cross Validate Potential Classifiers
rf_CVtest <- classifier.crossValidate(function(x){randomForest(as.factor(cat) ~ Length + Ref + Quality + Chromosome + Genotype,
                                                               data = x,cutoff=c(0.25,0.75),replace=TRUE,ntree=2000)},testData)

svm_CVtest <- classifier.crossValidate(function(x){svm(as.factor(cat) ~ Length + Ref + Quality + Chromosome + Genotype,
                                                       data = x,type="C-classification",cost=300,class.weights=setNames(c(1,5),c(1,0)))},testData)

logit_CVtest <- classifier.crossValidate(function(x){train(as.factor(cat) ~ Length + Ref + Quality + Chromosome + Genotype,data = x,method='LogitBoost',metric='Kappa',trControl=trainControl('none'))},testData)

knn_CVtest <- classifier.crossValidate(function(x){train(as.factor(cat) ~ Length + Ref + Quality + Chromosome + Genotype,data = x,method='knn',metric='Kappa',trControl=trainControl('none'))},testData)

boost_CVtest <- classifier.crossValidate(function(x){train(as.factor(cat) ~ Length + Ref + Quality + Chromosome + Genotype,
                                                           data = x,method='gbm',metric='Kappa',verbose=FALSE,trControl=trainControl('none'))},testData)

adaboost_CVtest <- classifier.crossValidate(function(x){train(as.factor(cat) ~ Length + Ref + Quality + Chromosome + Genotype,
                                                              data = x,method='adaboost',metric='Kappa',verbose=FALSE,trControl=trainControl('cv',repeats=5))},testData,divs=5)

ens <- function(y){ensembleclassifier(list(rf=function(x){randomForest(as.factor(cat) ~ Length + Ref + Quality + Chromosome + Genotype,
                                                                       data = x,cutoff=c(0.25,0.75),replace=TRUE,ntree=2000)},
                        ada=function(x){train(as.factor(cat) ~ Length + Ref + Quality + Chromosome + Genotype,
                                              data = x,method='adaboost',metric='Kappa',verbose=FALSE,trControl=trainControl('cv',repeats=5))},
                        logit=function(x){train(as.factor(cat) ~ Length + Ref + Quality + Chromosome + Genotype,data = x,method='LogitBoost',metric='Kappa',trControl=trainControl('none'))},
                        boost=function(x){train(as.factor(cat) ~ Length + Ref + Quality + Chromosome + Genotype,
                                                data = x,method='gbm',metric='Kappa',verbose=FALSE,trControl=trainControl('none'))},
                        knn=function(x){train(as.factor(cat) ~ Length + Ref + Quality + Chromosome + Genotype,data = x,method='knn',metric='Kappa',trControl=trainControl('none'))}),
                        y)}
ens_CVtest <- classifier.crossValidate(ens,testData)

cla <- list(rf=rf_CVtest,boost=boost_CVtest,svm=svm_CVtest,logit=logit_CVtest,adabost=adaboost_CVtest,knn=knn_CVtest,ens=ens_CVtest)
plotCVtest(cla,nam=c("Random Forest","Boosted Tree","SVM","Logit","Adaboost","KNN","Ensemble"))

## Determine Potential variable lists
# Forwards
vars <- colnames(testData)[-c(2,3,4,9,10,19,22,26,28)]

rf_factors <- classifier.forwardSelection(function(x){randomForest(as.factor(cat) ~ .,data = x,cutoff=c(1/5,4/5),replace=TRUE,ntree=2000)},testData,c("Chromosome","Start","Stop","Length","Ref","Quality","Genotype"))

svm_factors <- classifier.forwardSelection(function(x){svm(as.factor(cat) ~ .,data = x,type="C-classification",cost=300,class.weights=setNames(c(1,5),c(1,0)))},testData,vars)

logit_factors <- classifier.forwardSelection(function(x){train(as.factor(cat) ~ .,data = x,method='LogitBoost',metric='Kappa',trControl=trainControl('none'))},testData,vars)

#Reverse
rf_factors_rev <- classifier.reverseSelection(function(x){randomForest(as.factor(cat) ~ .,data = x,cutoff=c(1/5,4/5),replace=TRUE,ntree=2000)},testData,vars)

svm_factors_rev <- classifier.reverseSelection(function(x){svm(as.factor(cat) ~ .,data = x,type="C-classification",cost=300,class.weights=setNames(c(1,5),c(1,0)))},testData,vars)

logit_factors_rev <- classifier.reverseSelection(function(x){train(as.factor(cat) ~ .,data = x,method='LogitBoost',metric='Kappa',trControl=trainControl('none'))},testData,vars)


## Boost checks its own factors
boost_factors <- train(as.factor(cat) ~ .,data = testData[,c(vars,"cat"),with=FALSE],method='gbm',metric='Kappa',verbose=FALSE,trControl=trainControl('cv',number = 10))

adaboost_factors <- train(as.factor(cat) ~ .,data = testData[,c(vars,"cat"),with=FALSE],method='adaboost',metric='Kappa',verbose=FALSE,trControl=trainControl('cv',number = 10))


## ROC curves
rf_roc <- sapply(seq(0.05,0.95,0.05),function(x){rowMeans(classifier.crossValidate(function(y){randomForest(as.factor(cat) ~ Length + Ref + Quality + GC + SINEs + int + coef1 + SDbins + CentBins + TeloBins + LINEs + LTRs,
                                                                                                   data = y,cutoff=c(x,1 - x),replace=TRUE,ntree=2000)},testData))})

svm_roc_weight <- sapply(seq(1,8,0.5),function(x){rowMeans(classifier.crossValidate(function(y){svm(as.factor(cat) ~ Length + GC + int + coef1 + SDbins + TeloBins + RepTotal,
                                                        data = y,type="C-classification",cost=300,class.weights=setNames(c(1,x),c(1,0)))},testData))})

svm_roc_cost <- sapply(seq(100,1000,50),function(x){rowMeans(classifier.crossValidate(function(y){svm(as.factor(cat) ~ Length + GC + int + coef1 + SDbins + TeloBins + RepTotal,
                                                                                                    data = y,type="C-classification",cost=x,class.weights=setNames(c(1,5),c(1,0)))},testData))})

## RF Weight ROC
pdf("figures/rf_weight_roc.pdf",width = 10,height = 6)
par(fig=c(0,0.8,0,1),mar=c(4,4,4,4),xpd=FALSE,oma=c(1,1,1,5))
rocCurve(rf_roc,pch=20,col=colorRampPalette(c("red","blue"))(19),main="Random Forrest ROC Curve for Changing Weights")
abline(h=0.85)
colourBar(c(0.9,1,0.1,0.9),labels = seq(0.05,0.95,0.1),main = "False Weight")
dev.off()

## SVM Weight ROC
pdf("figures/svm_weight_roc.pdf",width = 10,height = 6)
par(fig=c(0,0.8,0,1),mar=c(4,4,4,4),xpd=FALSE,oma=c(1,1,1,5))
rocCurve(svm_roc_weight,pch=20,col=colorRampPalette(c("red","blue"))(15),main="SVM ROC Curve for Changing Weights")
abline(h=0.85)
colourBar(c(0.9,1,0.1,0.9),labels = seq(1,8,1),main = "False Weight")
dev.off()

pdf("figures/svm_cost_roc.pdf",width = 10,height = 6)
par(fig=c(0,0.8,0,1),mar=c(4,4,4,4),xpd=FALSE,oma=c(1,1,1,5))
rocCurve(svm_roc_cost,pch=20,col=colorRampPalette(c("red","blue"))(19),main="SVM ROC Curve for Changing Cost")
abline(h=0.85)
colourBar(c(0.9,1,0.1,0.9),labels = seq(100,1000,100),main = "Cost")
dev.off()

############################################ GS and BD cross validation ########################################
############ BD ############
bd <- read.table("data/HGDP01029/bd_10x_overlap.tsv",stringsAsFactors = FALSE)
colnames(bd) <- c('chr','start','stop','id','tenx')
bd <- bd[!is.na(bd$tenx),]

bd$tenx <- lapply(strsplit(bd$tenx,','),function(x){sapply(strsplit(x,':'),function(y){y[1]})})
dels$BD <- FALSE
calls <- unlist(bd$tenx)
dels[dels$ID %in% calls & dels$Source == 'HGDP01029','BD'] <- TRUE

testDataBD <- dels[!is.na(dels$Manual_check) & dels$Source == 'HGDP01029',]

ggplot(data=testData, aes(Manual_check)) + 
  geom_bar(aes(fill=as.factor(BD)), position="fill")+ 
  ggtitle("Proportion of Overlapping Calls (BD/Longranger)")


############ GS ############
gs1029 <- read.table('data/HGDP01029/HGDP01029_dels.vcf.gz_GS_OLAP_percent',fill = TRUE,header = TRUE)
gs1029$Source <- 'HGDP01029'
gs1029$GS <- !gs1029$GS_START == 'NO'

gs819 <- read.table('data/HGDP00819/HGDP00819_dels.vcf.gz_GS_OLAP_percent',fill = TRUE,header = TRUE)
gs819$Source <- 'HGDP00819'
gs819$GS <- !gs819$GS_START == 'NO'

gsMer <- rbind(gs1029[,c('Source','X10X_CNV_ID','GS')],gs819[,c('Source','X10X_CNV_ID','GS')])

testDataGS <- merge(testData,gsMer,all.x = TRUE,by.x = c('Source','ID'),by.y=c('Source','X10X_CNV_ID'))
testDataGS <- testDataGS[!duplicated(testDataGS[,-22]),]

testDataGS <- testDataGS[!is.na(testDataGS$GS),]

ggplot(data=testDataGS, aes(Manual_check)) + 
  geom_bar(aes(fill=as.factor(GS)), position="fill") + 
  ggtitle("Proportion of Overlapping Calls (GS/Longranger)")


###### GS vs BD overlap ######
getOverlap <- function(x1,x2){
  if (x1[2] > x2[2]){
    if (x1[1] > x2[1]){
      return(x2[2] - x1[1])
    } else {
      return(x2[2] - x2[1])
    }
  } else {
    if (x2[1] > x1[1]){
      return(x1[2] - x2[1])
    } else {
      return(x1[2] - x1[1])
    }
  }
}

overlap <- read.table('data/HGDP01029/gs_bd_overlap.tsv',header=TRUE,colClasses = c("character","character","integer","integer","character","integer","integer"))

overlap$ov <- sapply(1:dim(overlap)[1],function(x){getOverlap(c(overlap[x,'BD_Start'],overlap[x,'BD_End']),
                                                              c(overlap[x,'GS_Start'],overlap[x,'GS_End']))})

overlap$GS_Prop <- overlap$ov / (overlap$GS_End - overlap$GS_Start)
overlap$BD_Prop <- overlap$ov / (overlap$BD_End - overlap$BD_Start)

hist(overlap$BD_Prop,col=rgb(1,0,0,0.75),ylim = c(0,1000))
hist(overlap$GS_Prop,col=rgb(0,0,1,0.75),add=T)
