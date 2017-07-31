# Structural Variation in Human Populations Code Repository
**Author**: Alistair S. Dunham  
**Affiliation**: University of Cambridge & Welcome Trust Sanger Institute  
**Email**: asd43@cam.ac.uk

This repository contains the code used in my MPhil project, titled "Structural Variation in Human Populations", for reference and future use. This project was carried out with the Tyler-Smith group at the Sanger Institute as part of the 'Computational Biology' MPhil course at DAMTP, University of Cambridge.

## Folders
* *qualityControl* - Scripts used for QC of the 10X variant calls
* *breakdancerFiltering* - Scripts used to determine and apply filtering to BreakDancer calls based on the 10X calls
* *popgen* - Scripts performing the population genetic analyses
* *misc* - Miscellaneous scripts used here but with more general function

A list of the different scripts and a more detailed description than in the main report is given here. For any that require command line arguments the best description is obtained using the appropriate *script.file -h* command.

## Quality Control Scripts
#### *combinedQC.R*
R script containing the code used to train and cross-validate teh various models. It was used interactively to train all the models (see example *combinedQC_rawDataOnly.R*), optimise them and determine which factors to use.

#### *combinedQC_rawDataOnly.R*
Equivalent script only using the raw longranger data and before tuning. It was used to train the original models.

#### *preProcess3.py*
Command line tool to process longranger output ready for classification (training and classifying). It can add any of the possible extra genomic statistics to a currently processed dataset or process one from scratch with the selected information. It can add:
* Telomere and Centromere distances (built-in, does not require extra information)
* Read depth (requires a read depth file processed by *regionDepths.py*)
* Segmental Duplications (requires UCSC formated SD table)
* Repetative Elements (requires UCSC formated repeat table)
* GC content (requires a fasta file of the regions, for instance made by *samtools faidx*)

#### *preProcess2.py*
Equivalent script to *preProcess3.py* but using the REST client (*misc/restClient.py*) rather than extracting GC from a fasta file of the regions.

#### *preProcess10xBatch.sh*
Example shell script to run all preprocessing analyses for longranger data, generating read depths, and a fasta file of the regions then running *preProcess3.py*. It also shows setup for running on a list of files in a batch job array, using an id table to direct each job to the right files. Directory structures have been made generic for clarity.

#### *QC_functions.R*
Supporting QC functions for classifiction, including functions for calculating all the metrics used (FDR, FPR, FNR, sensitivity, specificity, precision, accuracy, F1 and Kappa), generallised cross-validation and feature selection. It also includes the function for constructing the ensemble classifier and predicting with it.

#### *qc_plots.R*
Script to generate the QC plots seen in the report, apart from ROC curves which are made as part of *combinedQC.R*.

#### *qualityControl.R*
Command line R script implementing the final random forest model. It supports training, model saving/loading and classification, all over lists of file inputs.

#### *regionDepths.py*
Command line script to group region depth reading proced from bam read files using *samtools depth*, into lists of regions defined in an input bed file.

## Filtering Scripts
#### *bd_10X_overlap.py*
Command line tool to determine BreakDancer variant calls which overlap with 10X longranger calls, based on a customisable minimum overlap proportion. A binary search is used to narrow down the region of calls to check, which greatly increases speed but means inputs must be sorted. The BreakDancer calls are loaded into memory and then 10X calls are streamed, determining overlap sequentially to reduce memory requirements.

#### *bd_filtering.R*
Script containing the code used to determine filtering criteria for the BreakDancer calls, based on the set which overlapped with the high quality filtered 10X calls. The selected filtering criteria are then tested against the overlap set and the GENOMESTRiP results. The filtering plots seen in the report are also produced here.

#### *gs_bd_overlap.py*
Equivalent script to *bd_10X_overlap.py* but determining which GENOMESTRiP calls overlap with BreakDancer ones, using the same algorithm but adapted to the difeferent format. BreakDancer calls are loaded into memory since the GENOMESTRiP output tends to be a much larger file.

## Population genetics Scripts
#### *getFeatureTable.py*
#### *gsPopgen.R*
#### *gsProcess.py*
#### *importPopgen.R*
#### *merge10X2bed.py*
#### *mergeBD2bed.py*
#### *popgen.R*
#### *popgenFunctions.R*

## Miscellaneous Scripts
#### *colourbar.R*
R function to place a customisable colour bar on plots. This is placed at a specific position on the current plot using the *fig* argument to *par()*, meaning it is positioned directly on the plot canvas and does not account for *mfrow* or *layout* settings. The labels, title and colour pallate can all be customised.

#### *popSummary.py*
Simple script to condense a list of samples and their associated information into a table giving all populations that are found, how many samples they have and the sample sexes, source and sequencing method.

#### *restClient.py*
Simple REST client for querying the Ensembl database. It has sequence retrieval built in specifically alongside a general querying function, with automatic rate limiting to prevent server timeout.

#### *vcf2bed.py*
Simple script to convert a VCF file, for instance as produced by GENOMESTRiP, to bed format. Most columns are already available but the end position needs to be extracted from the VCF information field. Regions can optionally be padded, for
 instance to extract flanking sequence. Conversion from 1-based to 0-based is performed.
 
