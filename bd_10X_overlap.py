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

## Command line tool to determine overlap between breakdancer and 10X longranger variant calls

# Import Modules
import argparse
import math
import fileinput

## Set up argument parsing
parser = argparse.ArgumentParser(description="Determines which deletion calls overlap between breakdancer and 10X.")
parser.add_argument('bdFile',metavar='B',type=str,help="path to breakdancer output")
parser.add_argument('tenXFile',metavar='T',type=str,help="path to processed 10X call file")
parser.add_argument('--prop','-p',type=float,default=0.8,help="Proportional overlap reqiuired for matches (default: 0.8)")
parser.add_argument('--source','-s',type=str,default='',help="Source genome for the calls (used as part of the BD IDs)")
parser.add_argument('--bdfilter','-b',action='store_true',help="Only output BD calls with overlaps")
args = parser.parse_args()

## Function to perform binary search on a sorted list of regions
## It returns a list of regions in the vicinity of the query region
def binaryRegSearch(regions,target,n=30,pad=1000):
	while len(regions) > 2*n:
		mid = math.ceil(len(regions)/2)
		if ((regions[mid][0] > target[0] - pad and regions[mid][0] < target[1] + pad) or
			(regions[mid][1] > target[0] - pad and regions[mid][1] < target[1] + pad) or
			(target[0] >= regions[mid][0] and target[1] <= regions[mid][1])):
			regions = regions[(mid-n):(mid+n)]
			break
		elif regions[mid][0] > target[0]:
			regions = regions[:mid]
		else:
			regions = regions[mid:]
	
	return(regions)

## Function determine the overlap between two regions, in absolute distance
def getOverlap(r1,r2):
	if r1[1] > r2[1]:
		if r1[0] > r2[0]:
			return(r2[1] - r1[0])
		else:
			return(r2[1] - r2[0])
	else:
		if r2[0] > r1[0]:
			return(r1[1] - r2[0])
		else:
			return(r1[1] - r1[0])

## Initiate Breakdancer calls dictionaries
bdCalls = {''.join(('chr',str(k))):[] for k in range(1,23)}
bdCalls['chrX'] = []
bdCalls['chrY'] = []

## Import breakdancer calls
n = 1
with fileinput.input(args.bdFile) as bd:
	for i in bd:
		if not i[0] == '#':
			t = i.strip().split()
			if t[0] == t[3] and t[0] in bdCalls.keys() and t[6] == 'DEL':
				bdCalls[t[0]].append([int(t[1]),int(t[4]),'.'.join(['BD',args.source,str(n)]),[],i.strip()])
				#Call Format: {chrom}[start,stop,ID,matches,call]
			
			n += 1

## Import 10X calls and check for overlap with breakdancer calls
with fileinput.input(args.tenXFile) as tenX:
	next(tenX)
	for i in tenX:
		t = i.strip().split()
		start = int(t[2])
		stop = int(t[3])
		ident = t[0] 
		r = binaryRegSearch(bdCalls[t[1]],(start,stop),n=50,pad=500)
		for j in r:
			ov = getOverlap((start,stop),(j[0],j[1]))
			if ov > args.prop * (stop - start) and ov > args.prop * (j[1] - j[0]):
				j[3].append([start,stop,ident])


## output results, either giving a tsv table of breakdancer calls with matches
## or a bed formated file giving each breakdancer calls and any overlaps
if args.bdfilter:
	print("Chr1", "Pos1", "Orientation1", "Chr2", "Pos2", "Orientation2", "Type",
		  "Size", "Score", "num_Reads", "num_Reads_lib", "bam","ID","tenXmatches",sep='\t')
	for k,v in bdCalls.items():
		for i in v:
			if len(i[3]) > 0:
				print(i[4],i[2],','.join([''.join([x[2],':',str(x[0]),'-',str(x[1])]) for x in i[3]]),sep='\t')
else:
	for k,v in bdCalls.items():
		for i in v:
			if len(i[3]) > 0:
				print(k,i[0],i[1],i[2],','.join([''.join([x[2],':',str(x[0]),'-',str(x[1])]) for x in i[3]]),sep='\t')
			else:
				print(k,i[0],i[1],i[2],'NA',sep='\t')

