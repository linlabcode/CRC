===
CRC
===

Core transcriptional regulatory circuitry analysis

Dependencies
============

The CRC software uses the following dependencies:

- Bamliquidator_

- Samtools_

- FIMO_

- fasta-get-markov_

Both FIMO and fasta-get-markov can be found in `The MEME Suite`_.

.. _Bamliquidator: https://github.com/BradnerLab/pipeline/wiki/bamliquidator
.. _Samtools: http://www.htslib.org/
.. _FIMO: http://meme-suite.org/doc/fimo.html
.. _fasta-get-markov: http://meme-suite.org/doc/fasta-get-markov.html
.. _The MEME Suite: http://meme-suite.org/doc/install.html

Install
=======

.. code::

  pip install git+https://github.com/linlabcode/CRC.git


Usage
=====

As a command line tool::

  crc -e <enhacer_text_file> -g <genome_build> -s <peaks_bed_file> -c <chromosomes_dir> -o <out_dir> -n <name>

As a python library::

  import crc

  crc.crc(enhancers, genome_input, chrom_path, output, analysis_name, bam=None, subpeak_file=None,
          mask_file=None, activity_path=None, const_extension=100, number=1, motifs=None,
          config='')

Authors
=======

Charles Y. Lin
Donald R. Polaski
Jost V. Koren
Alexander J. Federation

