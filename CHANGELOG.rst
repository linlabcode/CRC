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
- Add support for subpeak files in ``narrowPeak`` format
- Add new motifs for numerous TFs that already had at least one motif
  and add motifs for 23 TFs that didn't yet have a motif:

  -``ASCL1``
  -``BATF``
  -``BCL11A``
  -``BHLHB3``
  -``HES1``
  -``HES5``
  -``HES6``
  -``HSFY2``
  -``KLF2``
  -``MYF5``
  -``NEUROD1``
  -``NFYC``
  -``NOTCH1``
  -``PAX8``
  -``RAXL1``
  -``RBPJ``
  -``SIX5``
  -``TAL1``
  -``TP73``
  -``YAP1``
  -``YBX1``
  -``ZNF20``

Changed
-------

Fixed
-----
- Remove unused ``tfs`` parameter
- Fix ``chrom_path`` input formatting so it does not require a '/' at
  the end


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
