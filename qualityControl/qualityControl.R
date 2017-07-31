#!/usr/bin/env Rscript
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

## Command line R tool for classifying 10X variant calls
## Supports classification of many files, model training and model saving/loading to increase speed over batches
## See qualityControl.R -h for more details

## Parse Arguments
library(argparser,quietly = TRUE)
p <- arg_parser("Script to train a classifier for quality control of 10X longranger deletion calls and perform classification.")
p <- add_argument(p,"--train",short='-t',nargs = Inf,type="character",
                       help = "Path(s) to training files, as processed by preProcess2.py. Or a single Rdata file containing a trained modelcalled 'model', if the first file is an Rdata file this option overrides training a new model.")

p <- add_argument(p,"--man",short='-m',nargs = Inf,type="character",
                  help = "Path(s) to classification files. List must be in equialent order to training files.")

p <- add_argument(p,"--class",short='-c',nargs = Inf,type="character",
                       help = "Path(s) to files containing unclassified calls, as processed by preProcess2.py")

p <- add_argument(p,"--save",short='-s',type = "character",
                       help = "Save training data in supplied Rdata file.")

p <- add_argument(p,"--all",short='-a',flag=TRUE,
                 help = "Return all calls with associated details, not only those classified as true.")

args <- parse_args(p)

## Import required libraries
library(data.table,quietly = TRUE)
suppressMessages(library(randomForest,quietly = TRUE))    

## Check arguments
if (any(is.na(args$train))){
  write("Error: No training files supplied",stderr())
  quit(save = 'no')
  
} else if (!((length(args$train) == length(args$man) & !any(is.na(args$man))) | grepl('.Rdata$',args$train[1]))){
  write("Error: incorrect training files supplied. Each file must have a matching classification file or a single Rdata file containing a trained model must be supplied.",stderr())
  quit(save = 'no')
}

if (!any(is.na(args$save)) & !grepl('.Rdata$',args$save)){
  args$save <- paste0(args$save,'.Rdata')
} else if (any(is.na(args$save)) & any(is.na(args$class))){
  write("No files to classify and no model output file suplied. Defaulting to save trained model as model.Rdata.",stderr())
  args$save <- "model.Rdata"
}

## Function Definitions
loadDels <- function(dat,cla=NA,drop=TRUE){
  dels <- read.table(dat,header=TRUE,
                     colClasses = c("character","factor","integer","integer","integer","factor","numeric","factor","factor",
                                    "integer","integer","numeric","integer","integer","integer","integer","integer","integer","integer"
                                    ,"character","character"))
  
  ## Process read depth
  dels$ReadDepth <- strsplit(as.character(dels$ReadDepth),",")
  dels$ReadDepth <- lapply(dels$ReadDepth,as.numeric)
  
  #poly(seq(0,1,length.out = x$Length + 200),2)
  fits <- apply(dels,1,function(x){w <- seq(0,1,length.out = x$Length + 200) - 0.5 ;lm(x$ReadDepth ~ I(w^2))})
  dels$int <- sapply(fits,function(x){x$coefficients[1]})
  dels$coef1 <- sapply(fits,function(x){x$coefficients[2]})
  
  dels <- dels[,c("ID","Chromosome","Start","Stop","Length","Ref","Quality","Genotype","CentromereDist","TelomereDist","GC",
                  "LINEs","SINEs","LTRs","LowComplexityRepeats","SimpleRepeats","OtherRepeats","SDs","Source","int","coef1","ReadDepth")]
  
  dels$Source <- as.factor(dels$Source)
  dels$CentBins <- as.factor(dels$CentromereDist < 0.5 * 10^6)
  dels$TeloBins <- as.factor(dels$TelomereDist < 0.5 * 10^6)
  dels$SDbins <- as.factor(dels$SDs > 0)
  
  ## Add annotations
  if (is.na(cla)){
    return(dels)
  } else {
    annot <- read.table(cla,header=TRUE)
    mer <- merge(dels,annot[,c("chr","vcf_start","Manual_check")],by.x = c("Chromosome","Start"),by.y=c("chr","vcf_start"),all.x = TRUE)
    
    if (drop){
      mer <- mer[!is.na(mer$Manual_check),]
    }
    return(mer)
  }
}

## Load files
# If the first training file is an Rdata file this is assumed to contain a trained model and is loaded
if (grepl('.Rdata$',args$train[1])){
  load(args$train[1])
  
} else {
  # Otherwise load and process training data and train a model
  train <- rbindlist(lapply(1:length(args$train),function(x){loadDels(args$train[x],args$man[x])}))
  
  # Categorise true and false data
  train$cat <- as.integer(train$Manual_check %in% c("T","T?"))
  
  # Train model
  model <- randomForest(as.factor(cat) ~ Length + Ref + Quality + GC + SINEs + int + coef1 + SDbins + CentBins + TeloBins + LINEs + LTRs,
                        data = train[!train$Manual_check == '??'],cutoff=c(1/5,4/5),replace=TRUE,ntree=2000)
  
  # Save model if required
  if (!any(is.na(args$save))){
    save(model,file = args$save)
  }
}

## Classify if required
if (!any(is.na(args$class))){
  # Load and process data to classify
  cla <- rbindlist(lapply(1:length(args$class),function(x){loadDels(args$class[x])}))
  
  # Classify calls
  cla$cat <- predict(model,newdata = cla)
  
  # Output results
  if (args$all){
    write.table(cla[,c("ID","Chromosome","Start","Stop","Length","Ref","Quality","Genotype",
                                   "CentromereDist","TelomereDist","GC","LINEs","SINEs","LTRs","LowComplexityRepeats",
                                   "SimpleRepeats","OtherRepeats","SDs","Source","int","coef1","cat"),with=FALSE],
                stdout(),quote=FALSE,row.names=FALSE,sep='\t')
  } else {
    write.table(cla[cla$cat == 1,c("Source","ID","Chromosome","Start","Stop"),with=FALSE],stdout(),quote=FALSE,row.names=FALSE,sep='\t')
  }
}
