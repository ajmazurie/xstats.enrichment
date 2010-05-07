#!/usr/bin/env python

from distutils.core import setup, Extension

setup(
	name = "Enrichment",
	version = "1.0b",
	author = "Aurelien Mazurie",
	author_email = "ajmazurie@oenone.net",
	url = "http://github.com/ajmazurie/Enrichment",
	license = "MIT/X11",

	packages = ["Enrichment"],
	package_dir = {"Enrichment": "Enrichment"},
	ext_modules = [Extension("Enrichment.enrichment_", ["Enrichment/enrichment_.c"])],
	scripts = ["Tools/significant-attributes"],
)