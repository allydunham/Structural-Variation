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

## Command line tool to convert VCF files to bed format (for instance for later use with bedtools or bedops)
## Currently set up for deletion calls only (including those spread over multiple liens)
import sys
import os
import argparse
import fileinput

parser = argparse.ArgumentParser(description="Converts VCF files to bed formated regions. Only chr, start, stop and name are preserved. currently only set up for <DEL> and BND calls")
parser.add_argument('vcfFile',metavar='V',type=str,help="path to VCF file")
parser.add_argument('-p','--pad',dest='padLen',type=int,default=0,help='Number of bases to extend deletion region by')
args = parser.parse_args()

## Strem through and extract relavent information for bed formet
with fileinput.input(args.vcfFile) as vcf:
	for line in vcf:
		if line[0] == '#':
			continue
		
		t = line.strip().split()
		info = {}
		for i in t[7].split(';'):
			pair = i.split('=')
			info[pair[0]] = pair[1]
			
		if t[4] == '<DEL>':
			end=info['END']
			name = t[2]
		else:
			end = next(vcf).strip().split("\t")[1]
			name = t[2][:-2]
		
		print(t[0],int(t[1]) - args.padLen,int(end) + args.padLen,name,sep='\t')
		
