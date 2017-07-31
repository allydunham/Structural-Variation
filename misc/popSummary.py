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

## Script to convert hgdp populations table into a summary for each population
## counts members of each population, their sex, source, sequence type etc.

import argparse
import fileinput
import re

parser = argparse.ArgumentParser(description="Summarises HGDP Population Table with statistics aboout each population.")
parser.add_argument('popFile',metavar='P',type=str,help="Path to population table")
parser.add_argument('--ids','-i',type=str,default=False,help="List of ids processed (only these will be included)")
args = parser.parse_args()

if args.ids:
	with fileinput.input(args.ids) as li:
		valid_IDs = []
		for i in li:
			valid_IDs.append(i.strip())
		
		def IDtest(x):
			return(True if x in valid_IDs else False)
else:
	def IDtest(x):
		return(True)

## Import Populations
pops = {}
with fileinput.input(args.popFile) as fi:
	next(fi)
	for i in fi:
		t = i.strip().split()
		ident = re.split('[\._]',t[0])[0]
		# Check exclude = 0 (is this required?)
		if IDtest(ident):
			if t[4] in pops.keys():
				pops[t[4]]["PCR"] += 1 if t[2] == 'PCR' else 0
				pops[t[4]]["PCRfree"] += 1 if t[2] == 'PCRfree' else 0
				pops[t[4]]["Sanger"] += 1 if t[1] == 'sanger' else 0
				pops[t[4]]["Meyer"] += 1 if t[1] == 'meyer2012' else 0
				if not ident in pops[t[4]]["libs"]:
					pops[t[4]]["Samples"] += 1
					pops[t[4]]["Male"] += 1 if t[8] == 'M' else 0
					pops[t[4]]["libs"].add(ident)
					
				pops[t[4]]["n"] += 1
				
			else:
				pops[t[4]] = {"Samples":1,"Region":t[5],
							  "Male":1 if t[8] == 'M' else 0,
							  "PCR":1 if t[2] == 'PCR' else 0,
							  "PCRfree":1 if t[2] == 'PCRfree' else 0,
							  "Sanger":1 if t[1] == 'sanger' else 0,
							  "Meyer":1 if t[1] == 'meyer2012' else 0,
							  "libs":set([ident]),"n":1}



## Print Results - currently only set of IDs, not librarys
print("Population","Region","N","Male","Female","PCR","PCRfree","Other","Sanger","Meyer","SGDP","IDs",sep="\t")
for k,v in pops.items():
	print(k,v["Region"],v["Samples"],v["Male"],v["Samples"] - v["Male"],
		  v["PCR"],v["PCRfree"],v["n"] - v["PCR"] - v["PCRfree"],
		  v["Sanger"],v["Meyer"],v["n"] - v["Sanger"] - v["Meyer"],
		  ','.join(v["libs"]),sep="\t")

