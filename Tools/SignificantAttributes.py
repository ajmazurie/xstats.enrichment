#!/usr/bin/env python

__NAME = "SignificantAttributes"
__VERSION = "0.2"

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: Gathering options

from optparse import OptionParser
import sys

p = OptionParser(version = __NAME + " v" + __VERSION)

p.add_option("-m", "--mapping", dest = "attributes_file",
  help = "Mapping between attributes and population's entities (entity vs list of attributes)")

p.add_option("-q", "--query", dest = "query_file",
  help = "List of entities (query)")

p.add_option("--attributes-annotations", dest = "annotations_file",
  help = "Attributes' annotations ; optional") 

(p, a) = p.parse_args()

class Reader:
	def __init__ (self, file):
		try:
			self._source = open(file, "r")
		except:
			print " The file '%s' can't be opened."
			sys.exit(1)

	def __iter__ (self):
		return self

	# Return the next line available
	def next (self):

		while True:
			line = self._source.readline()

			if (line == ''):
				self._source.close()
				self._source = None
				raise StopIteration
			
			line = line.strip()
			if (line == '') or (line[0] == '#'):
				continue
			
			return line.split("	")

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: Loading population

G = 0  # number of entities in the population
C = {} # key: category, value: number of entities in this category

entity2categories = {}
g = 0

for data in Reader(p.attributes_file):

	if (data[0] == '@'):
		g = int(line[1:])
		continue
	
	entity, categories = data[0], data[1:]
	G += 1
	
	for category in categories:
		if C.has_key(category):
			c = C[category]
		else:
			c = 0

		C[category] = c + 1
	
	if (len(categories) > 0):
		entity2categories[entity] = categories

G = max(g, G)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: Loading query

N = 0	# size of the query

query2annotation = {} # key: entity in query, value: annotations (if any)
category2query = {}   # key: category, value: entities in query present in this category

for data in Reader(p.query_file):

	entity = data[0]
	N += 1

	if not entity2categories.has_key(entity):
		continue

	for category in entity2categories[entity]:
		if not category2query.has_key(category):
			category2query[category] = {}

		category2query[category][entity] = True

	if (len(data) > 1):
		query2annotation[entity] = data[1]

#:::::::::::::::::::::::::::::::::::::: Loading attributes annotations (if any)

category2annotation = {}

if p.annotations_file:
	for data in Reader(p.annotations_file):

		category, annotation = data[0], data[1]

		category2annotation[category] = annotation

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::: Computing enrichments

from math import log, exp

class FisherExactTest:

	# Reference: Lanczos, C. 'A precision approximation of the gamma function',
	# J. SIAM Numer. Anal., B, 1, 86-96, 1964." and http://www.matforsk.no/ola/fisher.htm 

	def _lngamma (self, z):
		x = 0
		x += 0.1659470187408462e-06 / (z + 7)
		x += 0.9934937113930748e-05 / (z + 6)
		x -= 0.1385710331296526 / (z + 5)
		x += 12.50734324009056 / (z + 4)
		x -= 176.6150291498386 / (z + 3)
		x += 771.3234287757674 / (z + 2)
		x -= 1259.139216722289 / (z + 1)
		x += 676.5203681218835 / (z)
		x += 0.9999999999995183

		return log(x) - 5.58106146679532777 - z + (z - 0.5) * log(z + 6.5)

	def _lnfactorial (self, n):
		if n <= 1: return 0
		return self._lngamma(n + 1)

	# Return the logarithm of combination of p elements between n
	#       / n \
	#  = ln \ p /
	def _lncombination (self, n, p):
		return self._lnfactorial(n) - self._lnfactorial(p) - self._lnfactorial(n - p)

	# Compute the P-value of the given selection of objects
	# Method: Fisher's exact test, with Gamma approximation of factorials
	def pvalue (self, k, n, C, G):
		v = 0
		for i in range(k):
			v = v + exp(
			  self._lncombination(C, i) +
			  self._lncombination(G - C, n - i) -
			  self._lncombination(G, n)
			 )
		
		v = 1.0 - v
		if v < 0: v = 0
		
		return v

	# Compute the enrichment of the given selection of objects
	def enrichment (self, k, n, C, G):
		return float((float(k) / C)) / (float(n) / G)

fisher = FisherExactTest()

categories = category2query.keys()
categories.sort()

print "Attributes		Enrichment	P-value	k	n	C	G	Entities	"

TAB = '	'

for category in categories:

	entities = category2query[category].keys()
	entities.sort()

	# compute p-value and enrichment
	K = len(entities)

	pvalue = fisher.pvalue(K, N, C[category], G)
	enrichment = fisher.enrichment(K, N, C[category], G)

	# display results
	cinfo = category2annotation.has_key(category) and category2annotation[category] or ''

	o = "%s	%s	%-.5f	%-6.3g	%s	%s	%s	%s" % (category, cinfo, enrichment, pvalue, K, N, C[category], G)

	for entity in entities:
		einfo = query2annotation.has_key(entity) and query2annotation[entity] or ''

		print o + "	%s	%s" % (entity, einfo)
