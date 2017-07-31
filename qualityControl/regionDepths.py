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

## Command line tool to group reads depths by region
## It uses the output of samtools depth command and a bed file containing regions of interest
import argparse
import sys
import numpy
import fileinput

## Parse arguments
parser = argparse.ArgumentParser(description="Assigns per base read depth values to their associated regions from an input bed file.")
parser.add_argument('bedFile',metavar='R',type=str,help="path to bed file containing regions of interest")
parser.add_argument('depthFile',metavar='D',type=str,help="path to file containing samtools depth output for the regions")
args = parser.parse_args()

## Initiate regions object
regions = {}
for i in range(1,23):
	regions[''.join(['chr',str(i)])] = []
regions['chrX'] = []
regions['chrY'] = []

## Import regions
with fileinput.input(args.bedFile) as bed:
	for li in bed:
		t=li.strip().split()
		regions[t[0]].append([int(t[1]),int(t[2]),t[3],[0]*(int(t[2])-int(t[1]))])

## Stream through depths and assign to regions
with fileinput.input(args.depthFile) as depth:
	for li in depth:
		t = li.strip().split()
		p = int(t[1])
		r = int(t[2])
		for i in regions[t[0]]:
			if p > i[0] and p <= i[1]:
				i[3][p - i[0] - 1] = r

## Print output as a 4 column bed file with a 5th column containing a comma separated list of per base depths
for k,v in regions.items():
	for i in v:
		print(k,i[0],i[1],i[2],','.join([str(x) for x in i[3]]))
