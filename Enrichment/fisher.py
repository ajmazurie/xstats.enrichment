#!/usr/bin/env python

# For people that are interested by the core of the SignificantAttributes
# script, here is the code that compute the enrichment and associated p-value

from math import log, exp

class EnrichmentEvaluator:

	# Reference: Lanczos, C. 'A precision approximation of the gamma function',
	# J. SIAM Numer. Anal., B, 1, 86-96, 1964."
	# and http://www.matforsk.no/ola/fisher.htm 
	#
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

	# return the logarithm of combination of p elements between n
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

e = EnrichmentEvaluator()

for c in range(100):
	print e.pvalue(c+1, 100, 100, 100),	" C =", c+1
