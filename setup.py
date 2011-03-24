#!/usr/bin/env python

try:
	from setuptools import setup, Extension, find_packages
except ImportError:
	from ez_setup import use_setuptools
	use_setuptools()
	from setuptools import setup, Extension, find_packages

#from distutils.core import setup, Extension

setup(
	name = "Enrichment",
	version = "1.0",
	description = "Statistics for enrichment analysis",
	long_description = open("README.rst").read(),
	url = "http://github.com/ajmazurie/Enrichment",
	license = open("LICENSE").read(),

	author = "Aurelien Mazurie",
	author_email = "ajmazurie@oenone.net",

	namespace_packages = ["xstats"],
	packages = find_packages("lib"),
	package_dir = {'': "lib"},
	ext_modules = [Extension("xstats.enrichment.enrichment_", ["lib/xstats/enrichment/enrichment_.c"])],
	scripts = ["tools/significant-attributes"],
)