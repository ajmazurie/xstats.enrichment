
import enrichment_

# Calculate the probability of getting b successes in a sample of n distinct
# objects drawn without replacement from a population of N objects with B
# successes. Return the left, right and two-tailed p-values.
#
# Adapted from WordHoard project - http://wordhoard.northwestern.edu
# (edu.northwestern.at.utils.math.statistics.FishersExactTest)
def fisher_exact_test (b, n, B, N):
	if (not ((b <= n) and (n <= N) and (B <= N) and (b <= B))):
		raise ValueError("Malformed contingency table: b is %s, n is %s, B is %s, N is %s" % (b, n, B, N))

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

# Calculate the probability of obtaining the observed density of successes
# at the top of a ranked list of N objects under the null assumption that all
# repartitions of B successes in the list are equiprobable. s is a ranked
# 'occurrence vector'; i.e., a vector that states, for all object in the ranked
# list, their success. If their is no success in the provided list, None
# is returned
#
# Adapted from the source code of GOrilla - Eden et al. (PMID 19192299)
def mHG (occurrences, B = None, N = None, max_size = 1000, with_pivot = False):
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

	if (not (b <= B <= N)):
		raise ValueError("Malformed contingency table: b is %s, B is %s, N is %s" % (b, B, N)) 

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
