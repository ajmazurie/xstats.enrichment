
import enrichment_

def evaluate_subset (b, n, B, N):
	""" Evaluate the enrichment of a subset of objects (drawn without
	replacement from a larger population) in some attribute (i.e., value from
	a categorical variable).
	
	Given the following parameters,
		- **b**: number of objects in the subset with this attribute
		- **n**: total number of objects in the subset
		- **B**: number of objects in the whole population with this attribute
		- **N**: total number of objects in the whole population

	evaluate_subset() will compare, using the Fisher's exact test, the ratio of
	objects having this attribute in the subset with the ratio found in the whole
	population, and return a p-value (probability of getting **b** objects with
	the attribute in a sample of **n** distinct objects drawn without
	replacement from a population of **N** objects of which **B** have the
	attribute).

	For example, knowing that 1000 persons out of a population of 2000 are male,
	how significantly enriched is a subset of 10 persons in males if 8 of them
	are male?
	
	Return:
		- left, right and two-tailed p-values.
	
	.. note:: Adapted from WordHoard project (http://wordhoard.northwestern.edu),
		package edu.northwestern.at.utils.math.statistics.FishersExactTest

	.. note:: In 'Gene set enrichment analysis made simple' (Irizarry et al.,
		2009) the authors demonstrate how Khi-square is superior to previously
		published GSEA methods.
	"""
	if (b > n):
		raise ValueError("b > n")

	if (B > N):
		raise ValueError("B > N")

	if (b > B):
		raise ValueError("b > B")

	if (n > N):
		raise ValueError("n > N")

	um, lm = min(n, B), max(0, n + B - N)

	if (um == lm):
		return 1.0, 1.0, 1.0

	cutoff = enrichment_.hypergeometric_distribution(b, n, B, N)
	left_tail, right_tail, two_tailed = 0, 0, 0

	for i in range(lm, um + 1):
		p = enrichment_.hypergeometric_distribution(i, n, B, N)

		if (i <= b):
			left_tail += p

		if (i >= b):
			right_tail += p

		if (p <= cutoff):
			two_tailed += p

	left_tail = min(left_tail, 1)
	right_tail = min(right_tail, 1)
	two_tailed = min(two_tailed, 1)

	return left_tail, right_tail, two_tailed

def evaluate_list (occurrences, B = None, N = None, max_size = 1000, with_pivot = False):
	""" Evaluate the enrichment of the top of a ranked list of objects in some
	attribute (i.e., value from a categorical variable).
	
	Given the following parameters,
		- **occurrences**: vector of booleans (or binary integers) stating, for
		  all objects in a ranked list, if this object possess the attribute
		- **B**: total number of objects with this attribute in the population
		  (optional). Default: number of objects with this attribute as
		  reported by **occurrences**
		- **N**: total number of objects in the population (optional). Default:
		  number of objects described by **occurrences**

	evaluate_list() will calculate the probability of obtaining the observed
	density of objects with the attribute at the top of a ranked list of **N**
	objects under the null assumption that all repartitions of **B** successes
	in the list are equiprobable.

	Return:
		- p-value
		- pivot: position in the list chosen to distinguish between the top
		  and bottom (only returned if **with_pivot** is set to True)

	.. note:: Adapted from Eden et al., 2007 and GOrilla source code
	"""
	b = len(filter(lambda x: x, occurrences))

	if (B != None):
		B_ = min(B, b, max_size)
	else:
		B = b
		B_ = min(B, max_size)

	if (N != None):
		N_ = min(N, len(occurrences), max_size)
	else:
		N = len(occurrences)
		N_ = min(N, max_size)

	if (b > B):
		raise ValueError("b > B")

	if (B > N):
		raise ValueError("B > N")

	if (b == 0) or (B == 0):
		if (with_pivot):
			return None, None
		else:
			return None

	# calculate the minimum HGT for the vector 's' of
	# successes, considering all possible partitions
	min_hgt = 1
	pivot = 0
	b = 0
	for n in range(1, N_ - 1):
		if (occurrences[n - 1]):
			b += 1

		pvalue = fisher_exact_test(b, n, B, N)[1]
		if (pvalue <= min_hgt):
			min_hgt = pvalue
			pivot = n

	pvalue = float(enrichment_.mHG_pvalue(B, N, max_size, min_hgt))

	if (pvalue <= 0.0):
		pvalue = min_hgt * n

	if (with_pivot):
		return pvalue, pivot
	else:
		return pvalue
