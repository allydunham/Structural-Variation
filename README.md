# Structural Variation in Human Populations Code Repository
**Author**: Alistair S. Dunham  
**Affiliations**: University of Cambridge & Welcome Trust Sanger Institute  
**Email**: asd43@cam.ac.uk

This repository contains the code used in my MPhill project, titled "Structural Variation in Human Populations", for reference and future use.

## Folders
* *qualityControl* - Scripts used for QC of the 10X variant calls
* *breakdancerFiltering* - Scripts used to determine and apply filtering to BreakDancer calls based on the 10X calls
* *popgen* - Scripts performing the population genetic analyses
* *misc* - Miscellaneous scripts used here but with more general function

A list of the different scripts and a more detailed description than in the main report is given here. For any that require command line arguments the best description is obtained using the appropriate *script.file -h* command.

## Quality Control Scripts
#### combinedQC.R
#### combinedQC_rawDataOnly.R
#### preProcess3.py
#### preProcess2.py
#### preProcess10xBatch.sh
#### QC_functions.R
#### qc_plots.R
#### qualityControl.R
#### regionDepths.py

## Filtering Scripts
#### bd_10X_overlap.py
#### bd_filtering.R
#### gs_bd_overlap.py
#### merge10X2bed.py
#### mergeBD2bed.py

## Population genetics Scripts
#### getFeatureTable.py
#### gsPopgen.R
#### gsProcess.py
#### importPopgen.R
#### popgen.R
#### popgenFunctions.R

## Miscellaneous Scripts
#### colourbar.R
#### gs2bed.py
#### popSummary.py
Simple script to condense a list of samples and their associated information into a table giving all populations that are found, how many samples they have and the sample sexes, source and sequencing method.
#### restClient.py
#### vcf2bed.py

