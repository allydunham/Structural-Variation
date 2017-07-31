#!/bin/bash
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

## Script to process a specific tenX callset ready for classification, submit as job array with index file to process many in a batch

## Input parameters - set to correct input/output directory structure
#main working directory
di=/main/working/directory

#location of index file (each row has sample name and location of its bam file)
sampleIndex=/dir/index/file/tenx_IDs.tsv

# indexed reference genome fasta file
genomeRef=/ref/genome/GRCh38.fa

# Directory containing all the processing scripts
scriptDir=/script/dir

# Reference line to use set by job array index - for batch processing
li=$LSB_JOBINDEX

## determine input statistics from job index

sample=$(awk -v var="$li" 'NR==var {print $1}' $sampleIndex)
tenXdir=$(awk -v var="$li" 'NR==var {print $2}' $sampleIndex)
echo "Sample $li"
echo $sample
echo $tenXdir
echo $di

## Calculate Read Depth
gunzip -c $tenXdir/outs/dels.vcf.gz > $di/$sample/dels.vcf
python3 $scriptDir/vcf2bed.py -p 100 $di/$sample/dels.vcf > $di/$sample/delsPadded.bed

samtools depth -b $di/$sample/delsPadded.bed $tenXdir/outs/phased_possorted_bam.bam | python3 $scriptDir/regionDepths.py $di/$sample/delsPadded.bed - > $di/$sample/delsReadDepthPadded

## Get GC content
python3 $scriptDir/vcf2bed.py $di/$sample/dels.vcf | awk -F "\t" '{ print $1":"$2+1"-"$3 }' - | xargs samtools faidx $genomeRef > $di/$sample/dels.fa

## Process data
python3 $scriptDir/preProcess3.py --source $sample --sd $di/genomicSuperDups.txt --reps $di/rmsk.txt --gc $di/$sample/dels.fa --depth $di/$sample/delsReadDepthPadded $di/$sample/dels.vcf > $di/$sample/tenx_dels_processedPadded.tsv

