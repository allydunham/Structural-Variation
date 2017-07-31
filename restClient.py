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

## Based on the example client found at https://github.com/Ensembl/ensembl-rest/wiki/Example-Python-Client
"""Simple module containing a REST client to query the ensembl REST server.
Currently only implements sequence recall for a given region and fully user specified calls.
But it can easily be extended to support easier calling of other common requests in a similar manner to the getSeq function
"""


import sys
import requests
import time

class RestClient:
	"""REST client class.
	Methods are included to query a given REST server directly and for
	specific sequences.
	"""
	
	def __init__(self,server='http://rest.ensembl.org',rate = 15):
		self.server = server
		self.count = 0
		self.rate = rate
		self.lastTime = time.time()

	def restRequest(self,request,header={'Content-Type':'text/plain'}):
		"""Make a REST request, check to limit request rate and return the server resoponse"""
		
		self.count += 1
		## Check for self limiting
		if self.count >= self.rate:
			dt = time.time() - self.lastTime
			if dt < 1:
				time.sleep(1 - dt)
			self.lastTime = time.time()
			self.count = 0
		
		## Perform Rest Request
		r = requests.get(self.server+request,headers=header)
		
		## Check request suceeded
		if not r.ok:
			if 'Retry-After' in r.headers:
				time.sleep(float(r.headers['Retry-After'])+3)
				self.restRequest(request,header)
			else:
				r.raise_for_status()
				return(-1)
		else:
			return(r)
			
	def getSeq(self,chrom,region,species='human',strand=1,expand=(0,0),conType='text/plain'):
		"""Request the sequence of a specific genomic region"""
		
		req = '/sequence/region/{0}/{1}:{2}..{3}:{4}?expand_5prime={5};expand_3prime={6}'.format(species,chrom,region[0],region[1],strand,expand[0],expand[1])
		r = self.restRequest(req,header={'Content-Type':conType})
		if r == -1:
			return('<ERROR>')
		else:
			return(r.text)