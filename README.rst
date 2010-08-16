Enrichment
==========

``Enrichment`` is a Python_ library for conducting two types of enrichment analysis:

- given a set of annotated objects selected from a larger population, it can evaluate the statistical significance of finding a given annotation in the selection when compared to the population. E.g., it can answer the question of how significant is it to find 12 women in a group of 20 person knowing that half of the world population is female. This evaluation is performed by a `Fisher's exact test <http://en.wikipedia.org/wiki/Fisher's_exact_test>`_.

- given a ranked list of annotated objects, it can identify if a given annotation is significantly enriched at the top of the list. The method used, from [Eden2007a]_ and [Eden2007b]_, doesn't require a cut-off to decide what the 'top' of the list is, but rather evaluate all cut-offs and corrects for multiple testing.

Finally, ``Enrichment`` provides a set of methods to correct the resulting p-values for multiple testing, when multiple annotations are evaluated: `Bonferroni <http://en.wikipedia.org/wiki/Bonferroni_correction>`_ adjustment, `Holm-Bonferroni adjustment <http://en.wikipedia.org/wiki/Holm-Bonferroni_method>`_, and `Benjamini-Hochberg FDR <http://en.wikipedia.org/wiki/False_discovery_rate>`_ adjustment.

A typical application of enrichment analysis in bioinformatics is the search for any property that may be shared by genes or proteins selected by an experiment. E.g., which functional annotation do genes selected for their high expression level have in common? In a list of genes ranked by decreasing fold change, in which metabolic pathway the most differentially expressed genes tend to fall into?

.. [Eden2007a] Eden E, Lipson D, Yogev S and Yakhini Z. Motif discovery in ranked lists of DNA sequences. PLoS Computational Biology, 2007 Mar 23;3(3):e39
.. [Eden2007b] Eden E. Discovering motifs in ranked lists of DNA sequences. Research thesis, 2007 Jan

Contact
-------

Aurelien Mazurie, ajmazurie@oenone.net

Keywords
--------

Enrichment analysis, Bioinformatic, Python, Fisher's exact test

Getting started
---------------

- Download the latest version of the library from http://github/ajmazurie/Enrichment/downloads
- Unzip the downloaded file, and go in the resulting directory
- Run ``python setup.py install``. Alternatively, you can package the library by typing ``python setup.py bdist``, which will result in the creation of a file dist/Enrichment-xxx.tar.gz, with 'xxx' being the version number and the name of your platform. Installing the library is thus as simple as ``easy_install dist/Enrichment-xxx.tar.gz``

From then you only have to import ``Enrichment`` to use the library::

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

.. _Python: http://www.python.org/
