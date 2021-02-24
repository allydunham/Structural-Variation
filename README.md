# Structural Variation in Human Populations Code Repository
**Author**: Alistair S. Dunham  
**Affiliation**: University of Cambridge & Welcome Trust Sanger Institute  
**Email**: asd43@cam.ac.uk

This repository contains the code used in my MPhil project, titled "Structural Variation in Human Populations", for reference and future use. This project was carried out with the Tyler-Smith group at the Sanger Institute as part of the 'Computational Biology' MPhil course at DAMTP, University of Cambridge. The work went on to be included in [Almarri et al. (2020)](https://doi.org/10.1016/j.cell.2020.05.024).

## Folders
* *qualityControl* - Scripts used for QC of the 10X variant calls
* *breakdancerFiltering* - Scripts used to determine and apply filtering to BreakDancer calls based on the 10X calls
* *popgen* - Scripts performing the population genetic analyses
* *misc* - Miscellaneous scripts used here but with more general function

A list of the different scripts and a more detailed description than in the main report is given here. For any that require command line arguments the best description is obtained using the appropriate *script.file -h* command. Scripts are split into sections and ordered by general usage point in the analysis, although primary scripts are listed before dependancies.

## Quality Control Scripts
#### *preProcess10xBatch.sh*
Example shell script to run all preprocessing analyses for longranger data, generating read depths, and a fasta file of the regions then running *preProcess3.py*. It also shows setup for running on a list of files in a batch job array, using an id table to direct each job to the right files. Directory structures have been made generic for clarity.

#### *preProcess3.py*
Command line tool to process longranger output ready for classification (training and classifying). It can add any of the possible extra genomic statistics to a currently processed dataset or process one from scratch with the selected information. It can add:
* Telomere and Centromere distances (built-in, does not require extra information)
* Read depth (requires a read depth file processed by *regionDepths.py*)
* Segmental Duplications (requires UCSC formated SD table)
* Repetative Elements (requires UCSC formated repeat table)
* GC content (requires a fasta file of the regions, for instance made by *samtools faidx*)

#### *preProcess2.py*
Equivalent script to *preProcess3.py* but using the REST client (*misc/restClient.py*) rather than extracting GC from a fasta file of the regions.

#### *regionDepths.py*
Command line script to group region depth reading proced from bam read files using *samtools depth*, into lists of regions defined in an input bed file.

#### *combinedQC.R*
R script containing the code used to train and cross-validate the various models. It was used interactively to train all the models (see example *combinedQC_rawDataOnly.R*), optimise them and determine which factors to use.

#### *combinedQC_rawDataOnly.R*
Equivalent script only using the raw longranger data and before tuning. It was used to train the original models.

#### *QC_functions.R*
Supporting QC functions for classifiction, including functions for calculating all the metrics used (FDR, FPR, FNR, sensitivity, specificity, precision, accuracy, F1 and Kappa), generallised cross-validation and feature selection. It also includes the function for constructing the ensemble classifier and predicting with it.

#### *qc_plots.R*
Script to generate the QC plots seen in the report, apart from ROC curves which are made as part of *combinedQC.R*.

#### *qualityControl.R*
Command line R script implementing the final random forest model. It supports training, model saving/loading and classification, all over lists of file inputs.



## Filtering Scripts
#### *bd_filtering.R*
Script containing the code used to determine filtering criteria for the BreakDancer calls, based on the set which overlapped with the high quality filtered 10X calls. The selected filtering criteria are then tested against the overlap set and the GenomeSTRiP results. The filtering plots seen in the report are also produced here.

#### *bd_10X_overlap.py*
Command line tool to determine BreakDancer variant calls which overlap with 10X longranger calls, based on a customisable minimum overlap proportion. A binary search is used to narrow down the region of calls to check, which greatly increases speed but means inputs must be sorted. The BreakDancer calls are loaded into memory and then 10X calls are streamed, determining overlap sequentially to reduce memory requirements.

#### *gs_bd_overlap.py*
Equivalent script to *bd_10X_overlap.py* but determining which GenomeSTRiP calls overlap with BreakDancer ones, using the same algorithm but adapted to the difeferent format. BreakDancer calls are loaded into memory since the GenomeSTRiP output tends to be a much larger file.

#### *bd_filter.py*
Command line tool to apply the selected filtering to a BreakDancer output file, with options to specify filters for:
* Length (upper and lower bounds)
* Quality Score (lower bound)
* Read Depth (upper and lower bounds)
* Supporting Read Ratio (lower bound)
* Copynumber (upper and lower bounds)
* Variant Type
* Chromosome (Restrict to the cannonical 1-22,X,Y)
* Centromere (Filter Ensembl GRCh38 centromere regions)
* Genomic Gaps (based on a supplied gap table)

## Population genetics Scripts
#### *merge10X2bed.py* and *mergeBD2bed.py*
Equivalent scripts which merge input 10X or BreakDancer output files into a single bed formatted file, with the ID field storing information about the source sample. This bed file can then be fed to *bedtools* *merge* and *cluster* to perform call merging and assigning each samples calls to their merged counterpart. This is not a very sophisticated method, not even checking for the degree of overlap, but it was found to perform very well accross all samples, with very few calls getting incorrectly merged. In general variants aligned reasonably well and did not overlap with others.

#### *getFeatureTable.py*
Script to convert merged and clustered BreakDancer and 10X call sets into a binary table of features. The table has rows for each calls and columns indicating whether the call was found in each sample.

#### *gsProcess.py*
Script to process GenomeSTRiP output, translating it into (optionally) a bed file, an information table and a copynumber table, with copynumber per call on each row in columns for each sample. The information table unpacks the VCF information field into a tsv file for easy analysis in R.

#### *popgen.R*
R script containing code to perform all population genetic analyses on the BreakDancer callset and produce the output plots. It performs:
* Hierarchical Clustering
* Principal Component Analysis
* Site Frequency
* Fst and Burden Estimation
* Functional Analysis

It also supports data saving and preloading to prevent the long import being done every time.

#### *gsPopgen.R*
Equivalent script to *popgen.R* but using the GenomeSTRiP calls. It performs the same set of analysis but extends them, using genotype information to calculate real Fst and burden values and partitions PCA by region.

#### *importPopgen.R*
Script to import the population genetic data (apart from VEP data which is imported in the appropriate script). It was separated out to tidy the two primary scripts.

#### *popgenFunctions.R*
Script containing supporting population genetic functions, primarily used to convert feature tables into GenePop format and to read in the pairwise statistics tables produced by *diveRsity*'s *diffCalc* function, used over *fastDivPart* because of the laters very high memory usage with this dataset (>90GB RAM without bootstrapping, *diffCalc* used <60GB with 10 bootstraps). Even so only 10 bootstraps were feasible in a reasonable time with the memory available.


## Miscellaneous Scripts
#### *colourbar.R*
R function to place a customisable colour bar on plots. This is placed at a specific position on the current plot using the *fig* argument to *par()*, meaning it is positioned directly on the plot canvas and does not account for *mfrow* or *layout* settings. The labels, title and colour pallate can all be customised.

#### *popSummary.py*
Simple script to condense a list of samples and their associated information into a table giving all populations that are found, how many samples they have and the sample sexes, source and sequencing method.

#### *restClient.py*
Simple REST client for querying the Ensembl database. It has sequence retrieval built in specifically alongside a general querying function, with automatic rate limiting to prevent server timeout.

#### *vcf2bed.py*
Simple script to convert a VCF file, for instance as produced by GenomeSTRiP, to bed format. Most columns are already available but the end position needs to be extracted from the VCF information field. Regions can optionally be padded, for instance to extract flanking sequence. Conversion from 1-based to 0-based is performed.
 
