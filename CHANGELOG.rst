##########
Change Log
##########

All notable changes to this project are documented in this file.


==========
Unreleased
==========

Added
-----
- Suport ``hg38`` genome and add ``hg19`` and ``hg38`` blacklist bed
  files
- Add inputs and outputs description to the README file

Changed
-------

Fixed
-----
- Remove unused ``tfs`` parameter


================
1.0.0 2018-10-22
================

Added
-----
- Make ``crc`` a callable python function
- Set up testing infrastructure

Changed
-------
- Remove ``_CLIQUES_ALL.txt`` table from the outputs and report only
  the top 100 ranked cliques in the ``_CLIQUE_SCORES_DEGREE.txt`` table

Fix
---
- Remove the deprecated ``_CLIQUE_SCORES_VSA.txt`` table, the creation
  of which was very memory consuming
- Make the ``bam`` input work
- Fix assigning enhancers to missing genes
- Fix identifying gene column in activity file
- Make the ``mask_file`` input work
