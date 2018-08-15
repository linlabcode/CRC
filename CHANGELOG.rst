##########
Change Log
##########

All notable changes to this project are documented in this file.


Added
-----
- Make ``crc`` a callable python function

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
