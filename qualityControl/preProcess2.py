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

## Version 2 of the 10X preProcessing script to prepare variant calls for classification
## Uses the REST api to get genomes equences for GC content
## See preProcess2.py -h for full details
import sys
import restClient as rc
import os
import argparse
import math
import fileinput

## class representing deletion calls
class delCall:
	__slots__ = ('ID','chrom','start','stop','ref','qual','filt','len','genotype',
				 'gc','teloDist','centDist','reps','readDepth','sd')
	
	def __init__(self,chr,start,stop,ID,ref,qual,filt,geno):
		self.chrom = chr
		self.start = int(start)
		self.stop = int(stop)
		self.ID = ID
		self.ref = ref
		try:
			self.qual = float(qual)
		except:
			self.qual = 0
		
		self.filt = filt
		self.len = self.stop - self.start ## +1 or not?
		self.genotype = geno
		self.gc = -1
		
		self.teloDist = min(self.start, chrLen[self.chrom] - self.stop)
		self.centDist = -1
		
		self.reps = {'LINE':0,'SINE':0,'LTR':0,'Low_complexity':0,'Simple_repeat':0,'Other':0}
		self.readDepth = -1
		self.sd = 0

	def __str__(self):
		return('\t'.join([str(x) for x in (self.ID,self.chrom,self.start,self.stop,self.len,
						  self.ref,self.qual,self.filt,self.genotype,self.centDist,
						  self.teloDist,self.gc,self.reps['LINE'],self.reps['SINE'],
						  self.reps['LTR'],self.reps['Low_complexity'],
						  self.reps['Simple_repeat'],self.reps['Other'],self.sd,self.readDepth)]
						)
					)

	def getCentDist(self,cent):
		cent = cent[self.chrom]
		if cent[0] > self.stop:
			self.centDist = cent[0] - self.stop
		elif cent[0] < self.start and cent[1] > self.stop:
			self.centDist = 0
		else:
			self.centDist = self.start - cent[1]
	
	def getReps(self,candidates):
		for i in candidates:
			if (i[0] < self.stop and i[0] > self.start) or (i[1] < self.stop and i[1] > self.start) or (self.start >= i[0] and self.stop <= i[1]):
				if i[2] in ('LINE','SINE','LTR','Low_complexity','Simple_repeat'):
					self.reps[i[2]] += 1
				else:
					self.reps['Other'] += 1
		
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

## Parse arguments
parser = argparse.ArgumentParser(description="Script to add genomic information to 10X longranger deletion calls. Everything apart from the vcf file is optional and will be replaced by default output if not specified.")
parser.add_argument('inFile',metavar='V',type=str,help="path to input file containing deletion calls. Assumed to be in VCF format, can be a tsv table as produced by this script (to add additional values) using the --add flag")
parser.add_argument('--cen','-c',dest='cenFile',default=False,type=str,help="path to file containing centromere locations. Defaults to using those from ensembl GRCh38)")
parser.add_argument('--gc','-g',dest='fetchGC',action='store_true',help="Fetch GC%% for the deletion regions (+100bp flanking) from ensembl")
parser.add_argument('--reps','-r',dest='repsFile',type=str,default=False,help="File containing repetative regions (formated like that from UCSC)")
parser.add_argument('--depth','-d',dest='depthFile',type=str,default=False,help="File containing region per base coverage")
parser.add_argument('--add','-a',dest='add',action='store_true',help="Add (or recalculate) specified statistics from a processed file")
parser.add_argument('--rate','-v',dest='restRate',type=int,default=10,help="Set the per second rate of ensembl requests")
parser.add_argument('--sd','-s',dest='sdFile',type=str,default=False,help="File containing segmental duplications (formated like that from UCSC)")
parser.add_argument('--source','-S',dest='source',type=str,default=False,help="Name assigned to the call group (e.g. the genome they come from)")
parser.add_argument('--category','-C',dest='cat',type=str,default=False,help="File detailing manual category assignment")
args = parser.parse_args()

## Chromosome Lengths - retrieved manually from latest ensembl build fa.fai file
chrLen = {'chr1':248956422,"chr2":242193529,"chr3":198295559,"chr4":190214555,"chr5":181538259,"chr6":170805979,"chr7":159345973,"chr8":145138636,"chr9":138394717,"chr10":133797422,"chr11":135086622,
		  "chr12":133275309,"chr13":114364328,"chr14":107043718,"chr15":101991189,"chr16":90338345,"chr17":83257441,"chr18":80373285,"chr19":58617616,"chr20":64444167,"chr21":46709983,"chr22":50818468,
		  "chrX":156040895,"chrY":57227415
}

## import centromeres or use default (retrieved manually from ensembl)
if args.cenFile:
	cen = open(args.cenFile,"r").readlines()
	centromeres = {}
	for i in cen:
		t = i.strip().split()
		centromeres[t[0]] = (int(t[1]),int(t[2]))
else:
	centromeres = {'chr1':(122026460,125184587),'chr2':(92188146,94090557),'chr3':(90772459,93655574),'chr4':(49708101,51743951),'chr5':(46485901,50059807),'chr6':(58553889,59829934),
				   'chr7':(58169654,60828234),'chr8':(44033745,45877265),'chr9':(43236168,45518558),'chr10':(39686683,41593521),'chr11':(51078349,54425074),'chr12':(34769408,37185252),
				   'chr13':(16000001,18051248),'chr14':(16000001,18173523),'chr15':(17000001,19725254),'chr16':(36311159,38280682),'chr17':(22813680,26885980),'chr18':(15460900,20861206),
				   'chr19':(24498981,27190874),'chr20':(26436233,30038348),'chr21':(10864561,12915808),'chr22':(12954789,15054318),'chrX':(58605580,62412542),'chrY':(10316945,10544039)
				   }
## Check source
if args.source:
	source = args.source
	
elif args.add:
	with open(args.inFile,"r") as fi:
		next(fi)
		source = next(fi).strip().split()[20]
		
else:
	source = 'NA'

## Process deletions
dels = {}
with fileinput.input(args.inFile) as fi:
	if args.add:
		next(fi)
		for line in fi:
			t = line.strip().split()
			d = delCall(t[1],t[2],t[3],t[0],t[5],t[6],t[7],t[8])
			d.centDist = int(t[9])
			d.teloDist = int(t[10])
			d.gc = float(t[11])
			d.len = int(t[4])
			d.reps = {'LINE':int(t[12]),'SINE':int(t[13]),'LTR':int(t[14]),'Low_complexity':int(t[15]),'Simple_repeat':int(t[16]),'Other':int(t[17])}
			d.sd = int(t[18])
			d.readDepth = t[19]
			dels[d.ID] = d
			
	else:
		for line in fi:
			if line[0] == '#':
				continue
			
			t = line.strip().split("\t")
			info = {}
			for i in t[7].split(';'):
				pair = i.split('=')
				info[pair[0]] = pair[1]
				
			if t[4] == '<DEL>':
				end = info["END"]
				
				d = delCall(t[0],t[1],end,t[2],t[3],t[5],t[6],t[9])
				d.getCentDist(centromeres)
			
			else:
				end = next(fi).strip().split("\t")[1]
				d = delCall(t[0],t[1],end,t[2][:-2],t[3],t[5],t[6],t[9])
				d.getCentDist(centromeres)

			dels[d.ID] = d
			
# Get GC content
if args.fetchGC:
	rest = rc.RestClient(rate=args.restRate)
	for k,i in dels.items():
		s = rest.getSeq(chrom=i.chrom.strip("chr"),region=(i.start,i.stop),expand=(100,100))
		i.gc = (s.count('G') + s.count('C'))/len(s) * 100

## get repetative regions
if args.repsFile:
	## Import Reps
	repDict = {}
	for i in chrLen.keys():
		repDict[i] = []
	
	with fileinput.input(args.repsFile) as reps:
		for i in reps:
			t = i.strip().split()
			if t[5] in repDict.keys():
				repDict[t[5]].append((int(t[6]),int(t[7]),t[12])) #start stop strand name class family

	## Assign reps
	for k,i in dels.items():
		r = binaryRegSearch(repDict[i.chrom],(i.start,i.stop))
		i.getReps(r)

## Assign read depth based on an input file with read depth lists with tagged calls
if args.depthFile:
	with fileinput.input(args.depthFile) as depth:
		for li in depth:
			t = li.strip().split()
			dels[t[3]].readDepth = t[4]

## Calculate SD overlap based on an SD lookup file
if args.sdFile:
	with fileinput.input(args.sdFile) as sds:
		sdDict = {}
		for i in chrLen.keys():
			sdDict[i] = []
			
		for li in sds:
			t = li.strip().split()
			if t[1] in sdDict.keys():
				sdDict[t[1]].append((int(t[2]),int(t[3])))
		
	for k,i in dels.items():
		r = binaryRegSearch(sdDict[i.chrom],(i.start,i.stop),pad=3000)
		for sd in r:
			if (sd[0] >= i.start and sd[0] <= i.stop) or (sd[1] >= i.start and sd[1] <= i.stop) or (i.start >= sd[0] and i.stop <= sd[1]):
				i.sd += 1

## Output processed deletions
print('\t'.join(("ID","Chromosome","Start","Stop","Length",
		"Ref","Quality","Filter","Genotype","CentromereDist","TelomereDist","GC","LINEs","SINEs","LTRs","LowComplexityRepeats","SimpleRepeats","OtherRepeats","SDs","ReadDepth","Source")))

for k,i in dels.items():
	print(i,source,sep='\t')
