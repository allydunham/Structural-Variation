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

## Command line tool to process genomestrip output into a CN per sample/per call table, a bed file and a detailed information table 
import argparse
import fileinput

## Parse arguments
parser = argparse.ArgumentParser(description="Script to process GenomeStrip VCF output into a feature table, bed file and information table. At least one of -f, -b or -t must be given.")
parser.add_argument('vcfFile',metavar='V',type=str,help="path to input gs VCF file")
parser.add_argument('--feat','-f',dest='featOut',default=False,type=str,help="path to output feature table file")
parser.add_argument('--bed','-b',dest='bedOut',default=False,type=str,help="path to output bed file")
parser.add_argument('--tab','-t',dest='tabOut',default=False,type=str,help="path to output table file (tsv format)")
parser.add_argument('--chrom','-c',dest='filtChrom',action="store_false",help="Filter to the canonical 24 chromosomes")
args = parser.parse_args()


## Process and open output files
if not (args.featOut or args.bedOut or args.tabOut):
	print("Error: no output files supplied. Use at least one of -f, -b or -t",file=sys.stderr())
	sys.exit(1)

if args.featOut:
	featFile = open(args.featOut,"w+")
if args.bedOut:
	bedFile = open(args.bedOut,"w+")
if args.tabOut:
	tabFile = open(args.tabOut,"w+")
	print('id','type','chrom','start','startConf','end','endConf','len','ref','alt','qual','gc','altFreqEst','inbreedingCoef',
		  'callRate','alleles','CNdist','CNmin','CNmax','dupeOverlap','LODscore','novel',sep='\t',file=tabFile)

def getCN(x,test):
	if test(x):
		return(x.split(":")[1])
	else:
		return("NA")
	

##stream through VCF and process
n = 1
chroms = [''.join(['chr',str(x)]) for x in range(1,23)]
chroms.append('chrX')
chroms.append('chrY')

with fileinput.input(args.vcfFile) as fi:
	for li in fi:
		#process header lines
		if li[:2] == "##":
			pass
		
		elif li[0] == "#":
			t = li.strip().split()
			samples = t[9:]
			nSamps = len(samples)
			if args.featOut:
				print("call",*samples,sep="\t",file=featFile)
		
		#Process entry lines	
		else:
			t = li.strip().split()
			if args.filtChrom or t[0] in chroms:
				ident = '_'.join(['GScall',str(n)])
				n += 1
				
				## Extract information field if required
				if args.tabOut or args.bedOut:
					inf = {x[0]:x[1] for x in map(lambda x:x.split('='),t[7].split(';'))}
					
				
				if args.featOut:
					## Determine genotype of each sample
					if t[8][-2:] == "FT":
						gens = map(lambda x:getCN(x,lambda y: y[-4:] == "PASS"), t[9:])
					else:
						gens = map(lambda x:getCN(x,lambda y: float(y[-4:]) > 20), t[9:])
					
					## output table row
					print(ident,*gens,sep="\t",file=featFile)
				
				if args.bedOut:
					## Print bed entry
					print(t[0],int(t[1]) - 1,inf['END'],ident,sep='\t',file=bedFile)
					
				if args.tabOut:
					## Print information table
					for i in ('SVTYPE','CIPOS','END','CIEND','SVLEN','GCFRACTION','GLALTFREQ','GLINBREEDINGCOEFF',
							  'GSCALLRATE','GSCNALLELES','GSCNDIST','GSCNMIN','GSCNMAX','GSDUPLICATEOVERLAP',
							  'GSDUPLICATESCORE','NOVEL'):
						if not i in inf.keys():
							inf[i] = '.'
					
					print(ident,inf['SVTYPE'],t[0],int(t[1]) - 1,inf['CIPOS'],inf['END'],inf['CIEND'],inf['SVLEN'],t[3],t[4],t[5],inf['GCFRACTION'],
						  inf['GLALTFREQ'],inf['GLINBREEDINGCOEFF'],inf['GSCALLRATE'],inf['GSCNALLELES'],inf['GSCNDIST'],inf['GSCNMIN'],inf['GSCNMAX'],
						  inf['GSDUPLICATEOVERLAP'],inf['GSDUPLICATESCORE'],inf['NOVEL'],
						  sep='\t',file=tabFile)
