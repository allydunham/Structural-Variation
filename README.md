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
R function to place a customisable colour bar on plots. This is placed at a specific position on the current plot using the *fig* argument to *par()*, meaning it is positioned directly on the plot canvas and does not account for *mfrow* or *layout* settings. The labels, title and colour pallate can all be customised.

#### popSummary.py
Simple script to condense a list of samples and their associated information into a table giving all populations that are found, how many samples they have and the sample sexes, source and sequencing method.
#### restClient.py

#### vcf2bed.py
Simple script to convert a VCF file, for instance as produced by GENOMESTRiP, to bed format. Most columns are already available but the end position needs to be extracted from the VCF information field.


