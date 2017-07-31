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

## Command line tool to determine voerlap between breakdancer and genomestrip variant calls
import argparse
import math
import re
import fileinput

parser = argparse.ArgumentParser(description="Determines which deletion calls overlap between breakdancer and 10X.")
parser.add_argument('bdFile',metavar='B',type=str,help="path to breakdancer output")
parser.add_argument('gsFile',metavar='G',type=str,help="path to genome strip call file (VCF format)")
parser.add_argument('--prop','-p',type=float,default=0.8,help="Proportional overlap reqiuired for matches (default: 0.8)")
parser.add_argument('--source','-s',type=str,default='UNKNOWN',help="Source genome for the calls (used as part of the BD IDs)")
parser.add_argument('--upper','-u',type=int,default=10**10,help="Upper bound for considered deletion calls")
parser.add_argument('--lower','-l',type=int,default=0,help="Lower bound for considered deletion calls")
args = parser.parse_args()

## Function to perform a binary search on a sorted list of regions
## It returns all regions in the vicinity of the query
def binaryRegSearch(regions,target,n=30,pad=1000):
	while len(regions) > 2*n:
		mid = math.ceil(len(regions)/2)
		if (regions[mid][0] > target[0] - pad and regions[mid][0] < target[1] + pad) or (regions[mid][1] > target[0] - pad and regions[mid][1] < target[1] + pad) or (target[0] >= regions[mid][0] and target[1] <= regions[mid][1]):
			regions = regions[(mid-n):(mid+n)]
			break
		elif regions[mid][0] > target[0]:
			regions = regions[:mid]
		else:
			regions = regions[mid:]
	
	return(regions)

## Function to determine the overlap between two regions
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

## Read BreakDancer calls
bdCalls = {''.join(('chr',str(k))):[] for k in range(1,23)}
bdCalls['chrX'] = []
bdCalls['chrY'] = []

n = 1
with fileinput.input(args.bdFile) as bd:
	for i in bd:
		if not i[0] == '#':
			t = i.strip().split()
			l = int(t[4]) - int(t[1])
			if t[0] == t[3] and t[0] in bdCalls.keys() and t[6] == 'DEL' and l < args.upper and l > args.lower:
				bdCalls[t[0]].append([int(t[1]),int(t[4]),'_'.join(['BD',args.source,str(n)])]) #{chrom}[[start,stop,ID]]
			
			n += 1

## Read Genome Strip Calls
gsCalls = {''.join(('chr',str(k))):[] for k in range(1,23)}
gsCalls['chrX'] = []
gsCalls['chrY'] = []
with fileinput.input(args.gsFile) as gs:
	for i in gs:
		if not i[0] == '#':
			t = i.strip().split()
			end = re.split('[=;]',t[8])[1]
			l = int(end) - int(t[2])
			if t[1] in gsCalls.keys() and l < args.upper and l > args.lower:
				gsCalls[t[1]].append([int(t[2]),int(end),'_'.join([t[3],args.source])]) #{chrom}[[start,stop,ID]]

## Calculate overlapping reads and output
print('Source','BD_ID','BD_Start','BD_End','GS_ID','GS_Start','GS_End',sep='\t')
for chrom,v in bdCalls.items():
	for call in v:
		r = binaryRegSearch(gsCalls[chrom],call)
		for i in r:
			ov = getOverlap(call,i)
			if ov > args.prop * (call[1] - call[0]) or ov > args.prop * (i[1] - i[0]):
				print(args.source,call[2],call[0],call[1],i[2],i[0],i[1])

