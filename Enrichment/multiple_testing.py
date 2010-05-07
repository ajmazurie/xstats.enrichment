# Correct p-values for mutiple testing
# Equivalent in R: p.adjust ('stats' package)

# Correct a list of p-values using the Bonferroni adjustment
# Return a list of corrected p-values; null values are ignored.
# cf. http://en.wikipedia.org/wiki/Bonferroni_correction
def bonferroni_adjustment (p_values):
	n = len(filter(lambda x: x != None, p_values)) + 0.0

	adjusted_p_values = []
	for p_value in p_values:
		if (p_value == None):
			adjusted_p_values.append(None)
		else:
			adjusted_p_values.append(min(p_value * n, 1))

	return adjusted_p_values

# Correct a list of p-values using the Holm-Bonferroni adjustment
# Return a list of corrected p-values; null values are ignored.
# cf. http://en.wikipedia.org/wiki/Holm-Bonferroni_method
def holm_adjustment (p_values):
	# multiply p-values by a corrective factor, ignoring null entries
	n, c = len(filter(lambda x: x != None, p_values)), 0
	m = []

	adjusted_and_ranked_p_values = []
	for i, (i_, p_value) in enumerate(sorted(enumerate(p_values), lambda x, y: cmp(x[1], y[1]))):
		m.append(i_)
		if (p_value == None):
			adjusted_and_ranked_p_values.append(None)
		else:
			adjusted_and_ranked_p_values.append(min(p_value * (n - c), 1))
			c += 1

	# correct the p-values out of their proper order
	adjusted_p_values = [0 for c in range(i+1)]
	for i, p_value in enumerate(adjusted_and_ranked_p_values):
		if (p_value != None):
			adjusted_and_ranked_p_values[i] = max(filter(lambda x: x != None, adjusted_and_ranked_p_values[:i+1]))

		adjusted_p_values[m[i]] = adjusted_and_ranked_p_values[i]

	return adjusted_p_values

# Given a list of p-values, produce a ranked list of increasing q-values.
# A q-value in position k represents the false discovery rate, or expected
# proportion of false positives in the k first hypotheses.
# Assumes that the tests are independent or positively correlated.
#
# Uses the Benjamini-Hochberg algorithm
# cf. http://en.wikipedia.org/wiki/False_discovery_rate
def fdr (p_values, produce_ranking = False):
	# multiply p-values by a corrective factor, ignoring null entries
	n, c = len(filter(lambda x: x != None, p_values)), 0
	m = []

	q_values = []
	for i, (i_, p_value) in enumerate(sorted(enumerate(p_values), lambda x, y: cmp(y[1], x[1]))):
		m.append(i_)
		if (p_value == None):
			q_values.append(None)
		else:
			q_values.append(min(p_value * n / (n - c), 1))
			c += 1

	# correct the p-values out of their proper order
	fdr = [0 for c in range(i+1)]

	for i, q_value in enumerate(q_values):
		if (q_value != None):
			q_values[i] = min(filter(lambda x: x != None, q_values[:i+1]))

		if (produce_ranking):
			if (q_value == None):
				fp = None
			else:
				fp = int(round(q_values[i] * (n - i)))
			fdr[n-(i+1)] = (m[i], q_values[i], fp)
		else:
			fdr[m[i]] = q_values[i]

	return fdr
