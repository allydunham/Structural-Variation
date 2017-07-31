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

# Script to merge 10X classification files to bed format for processing 

import argparse
import fileinput
import re

parser = argparse.ArgumentParser(description="Merge classified 10X deletion calls into a single bed file for genotyping. Each is given an id combining the source and the call number.")
parser.add_argument('tenx',metavar='T',type=str,help="File containing a list of classified 10X call paths (one per line)")
args = parser.parse_args()

tenXfiles = []
with fileinput.input(args.tenx) as fi:
	for i in fi:
		tenXfiles.append(i.strip())


## Stream through each file and output bed formated versions

for fi in tenXfiles:
	with fileinput.input(fi) as currentFile:
		next(currentFile)
		for li in currentFile:
			t = li.strip().split()
			print(t[2],t[3],t[4],'.'.join([t[0],t[1]]),sep='\t')
