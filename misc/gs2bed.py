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

## Command line tool to convert genomestrip output to bed format

import argparse
import fileinput
import sys

parser = argparse.ArgumentParser(description="Filter and Convert a GS output file to bed format.")
parser.add_argument('gsFile',metavar='G',type=str,help="path to genome strip call file (VCF format)")
parser.add_argument('--len','-l',type=float,nargs=2,default=[0,10**9],help="Range of variant sizes to accept")
parser.add_argument('--chrom','-c',action="store_false",help="Filter to canonical chromosome set")
args = parser.parse_args()

chroms = [''.join(['chr',str(x)]) for x in range(1,23)]
chroms.append('chrX')
chroms.append('chrY')

with fileinput.input(args.gsFile) as gs:
	for li in gs:
		if not li[0] == "#":
			t = li.strip().split()
			info = {i[0]:i[1] for i in map(lambda	x: x.split('='),t[7].split(';'))}
			length = int(info['END']) - int(t[1])
			if (length < args.len[1] and length > args.len[0] and (args.chrom or t[0] in chroms)):
				print(t[0],t[1],info['END'],t[2],sep='\t')
