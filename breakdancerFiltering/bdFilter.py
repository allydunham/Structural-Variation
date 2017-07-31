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

## Command line tool to filter breakdancer reads based on a set of input criteria

# Import modules
import argparse
import fileinput
import re

# Initialise argument parser
parser = argparse.ArgumentParser(description="Filter BreakDancer output based on the derived filtering criteria. Upper and lower bounds are input as -X lower upper. "
								 "Also adds a unique ID to each BD call, incorporating the sample source and a id number. "
								 "--centro and --gaps require --chrom, else non-standard chromosomes produce an errror.")
parser.add_argument('bdFile',metavar='B',type=str,help="Path to breakdancer output")
parser.add_argument('--length','-l',type=float,nargs=2,default=[200,30000],help="Lower and upper bounds for SV length")
parser.add_argument('--score','-s',type=float,default=93,help="Lower bound for SV score")
parser.add_argument('--depth','-d',type=float,nargs=2,default=[0,50],help="Lower and upper bounds for SV read depth")
parser.add_argument('--copynumber','-c',type=float,nargs=2,default=[0,5],help="Lower and upper bounds for SV Copy Number")
parser.add_argument('--type','-t',type=str,default='DEL',help="Deletion type. Must be one of INS, DEL, INV, ITX, CTX or Unknown.")
parser.add_argument('--source','-S',type=str,default='UNKNOWN',help="Sample the calls were made from.")
parser.add_argument('--chrom','-C',action='store_false',help="Filter to standard 24 chromosomes")
parser.add_argument('--centro','-e',action='store_false',help="Filter centromeres.")
parser.add_argument('--gaps','-g',type=str,default=False,help="Filter problematic regions such as telomeres contained in the supplied gap table (UCSC gap table format).")
parser.add_argument('--readRatio','-r',type=float,default=0.6,help="Minimum ratio between mapped reads (average of both ends) and supporting reads.")
args = parser.parse_args()

## Set up list of valid chromosomes
chroms = [''.join(['chr',str(x)]) for x in range(1,23)]
chroms.append('chrX')
chroms.append('chrY')

## Centromeres obtained from ensembl GRCh38 (June 2017)
centromeres = {'chr1':(122026460,125184587),'chr2':(92188146,94090557),'chr3':(90772459,93655574),'chr4':(49708101,51743951),'chr5':(46485901,50059807),'chr6':(58553889,59829934),
				   'chr7':(58169654,60828234),'chr8':(44033745,45877265),'chr9':(43236168,45518558),'chr10':(39686683,41593521),'chr11':(51078349,54425074),'chr12':(34769408,37185252),
				   'chr13':(16000001,18051248),'chr14':(16000001,18173523),'chr15':(17000001,19725254),'chr16':(36311159,38280682),'chr17':(22813680,26885980),'chr18':(15460900,20861206),
				   'chr19':(24498981,27190874),'chr20':(26436233,30038348),'chr21':(10864561,12915808),'chr22':(12954789,15054318),'chrX':(58605580,62412542),'chrY':(10316945,10544039)
				   }

## Function to calculate if two regions overlap
def overlap(r1,r2):
	"""Check if two regions overlap"""
	return (True if r1[0] <= r2[1] and r2[0] <= r1[1] else False)

## Function to check if a region overlaps with any in a list
def checkGaps(r,gaps):
	"""Check if a region overlaps with any regions in a list"""
	for i in gaps:
		if overlap(r,i):
			return True
	
	return False

## Determine supporting read ratio of breakdancer calls (reads that support call/all reads at locus)
def readRatio(t):
	"""Calculate supporting read ratio of a breakdancer variant call"""
	start = re.split("[+-]",t[2])
	end = re.split("[+-]",t[5])
	return 2*int(t[9])/sum([int(start[0]),int(start[1]),int(end[0]),int(end[1])])

## if requested import list of gap regions to filter
if args.gaps:
	ignore_gaps = False
	gaps = {c:[] for c in chroms}
	with fileinput.input(args.gaps) as g:
		for li in g:
			t = li.strip().split()
			if t[1] in chroms:
				gaps[t[1]].append((int(t[2]),int(t[3]),t[7]))
	
else:
	ignore_gaps = True

## Loop through reads and filter any that dont match the specified criteria
with fileinput.input(args.bdFile) as bd:
	n = 1
	for line in bd:
		if line[0] == '#':
			print(line.strip())
		
		else:
			t = line.strip().split()
			if (t[6] == args.type and
			int(t[8]) >= args.score and
			int(t[7]) <= args.length[1] and int(t[7]) >= args.length[0] and
			int(t[9]) <= args.depth[1] and int(t[9]) >= args.depth[0] and
			(t[11] == 'NA' or (float(t[11]) >= args.copynumber[0] and float(t[11]) <= args.copynumber[1])) and
			(args.chrom or t[0] in chroms) and
			(args.centro or not overlap((int(t[1]),int(t[4])),centromeres[t[0]])) and
			(ignore_gaps or not checkGaps((int(t[1]),int(t[4])),gaps[t[0]])) and
			readRatio(t) >= args.readRatio):
				
				print(line.strip(),'.'.join(['BD',args.source,str(n)]),sep='\t')
			
			n += 1

