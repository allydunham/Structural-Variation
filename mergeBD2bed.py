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

## Script to merge filtered BD call files to bed format for processing 
import argparse
import fileinput


parser = argparse.ArgumentParser(description="Merge filtered BreakDancer variant calls into a single bed file for genotyping. "
								 "If the ID flag is given an ID is expected in the column beyond the normal BD output, otherwise an id is generated.")
parser.add_argument('bd',metavar='B',type=str,help="File containing a list of BreakDancer variant call files (.BD_out format).")
parser.add_argument('--id','-i',action='store_true',help="Use IDs added to the BD_out file on filtering. Otherwise generate IDs enumerating only filtered calls.")
args = parser.parse_args()

bdFiles = []
with fileinput.input(args.bd) as fi:
	for i in fi:
		bdFiles.append(i.strip())


## Stream through each file and output bed formated versions
for fi in bdFiles:
	with fileinput.input(fi) as bd:
		if not args.id:
			f = fi.split('/')
			idBase = f[-1].strip('.BD_out')
			n = 1
		
		for li in bd:
			if not li[0] == '#':
				t = li.strip().split()
				if args.id:
					ID = t[12]
				else:
					ID = '.'.join(['BD',idBase,n])
					n += 1
				
				print(t[0],t[1],t[4],ID,sep='\t')