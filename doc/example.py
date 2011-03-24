#!/usr/bin/env python

import Enrichment

# Analysis 1: how significant is it to have 10 objects out of 500
# that share a given annotation, knowing that 120 out of the 1800
# objects in the general population have this annotation?
l, r, t = Enrichment.fisher_exact_test(10, 500, 120, 1800)

# the left-tailed probability is the probability of having less
# than 10 objects out of 500 with this annotation:
print "left-tailed:", l # 5.25e-8

# the right-tailed probability is the probability of having more
# than 10 objects out of 500 with this annotation:
print "right-tailed:", r # 0.99

# the two-tailed probability is the probability of observing 10
# objects out of 500 with this annotation, plus the probabilities
# of observing even less likely proportions:
print "two-tailed:", t # 7.78e-8

# as a result, we demonstrate that finding 10 objects with this
# annotation is unexpected, as shown by the two-tailed p-value
# (significant at an error rate of 5%). To know exactly in which
# way the finding is unexpected, just look at the left- and right-
# tailed p-values. In this example the left-tailed p-value is very
# low, while the right-tailed p-value is almost 1. It means that
# finding 10 objects is unexpectedly low in regard of the general
# population.

# conversely, finding 100 objects in the selection of 500 with
# this property is unexpectedly high:
l, r, t = Enrichment.fisher_exact_test(100, 500, 120, 1800)

# the right-tailed p-value is very low, while the left-tailed is 1
print "left- and right-tailed:", l, r # 1, 1.32e-39


# Analysis 2: in a ranked list of 1000 objects, of which half of
# them share a given annotation, how significant is it to find 20
# objects with this annotation at the top of the list?

# we build an occurrence vector which, for each object in the list,
# contains either True or False to represent if this object have
# the annotation considered.

# we start with an homogeneous distribution:
occurrences = [True, False] * 500

# as expected, there is no significant enrichment at the top of the list:
print Enrichment.mHG(occurrences) # 0.99

# we now build a second occurrence vector, in which the first 20 objects
# all have the annotations:
occurrences = [True] * 20 + [False] * 20 + [True, False] * 480

# this time we found a significant enrichment at the top of the list,
# which in this case is determined as the 20 first entries:
p_value, pivot = Enrichment.mHG(occurrences, with_pivot = True)

print "p-value:", p_value # 3.67e-5
print "pivot:", pivot # 20
