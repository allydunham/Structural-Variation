#!/usr/bin/env python3
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

## Script to convert bedtools cluster output of merged breakdancer calls into a feature table
## which has columns for each sample indicating the presence of each deletion
import argparse
import fileinput
import re

## Read arguments
parser = argparse.ArgumentParser(description="Transform bedtools cluster output for deletion calls into a feature table of 'genotypes'.")
parser.add_argument('tenx',metavar='T',type=str,help="Bed file containing clustered deletion calls")
parser.add_argument('--bd','-b',action='store_true',help="Expect BreakDancer formatted IDs. Otherwise expect 10X formatted IDs.")
args = parser.parse_args()

## Determine function to use for setting sample ID depending on given source format
if args.bd:
	def getSample(x):
		"""Extract sample from BreakDancer formatted ID tags"""
		return(re.split("[_.]",x)[-2])
else:
	def getSample(x):
		"""Extract sample from 10X formatted ID tags"""
		return(x.split('.')[0])
	

## Extract each deletion call and its cluster number
dels = []
samples = set()
with fileinput.input(args.tenx) as bed:
	for li in bed:
		t = li.strip().split()
		s = getSample(t[3])
		n = int(t[4])
		samples.add(s)
		if len(dels) < n:
			dels.append(set([s]))
		
		else:
			dels[n - 1].add(s)


## Print feature table
samples = sorted(list(samples))
print("Deletion",*samples,sep='\t')
for n,delSamples in enumerate(dels):
	## generate feature string
	feats = [(1 if i in delSamples else 0) for i in samples]
	print('_'.join(["del",str(n + 1)]),*feats,sep='\t')