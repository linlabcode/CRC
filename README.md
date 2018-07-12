# CRC

Core transcriptional regulatory circuitry analysis

## Dependencies

The CRC software uses the following dependencies:

- ``Bamliquidator``

- ``Samtools``

- ``FIMO``

- ``fasta-get-markov``

Both ``FIMO`` and ``fasta-get-markov`` can be found in ``The MEME Suite``.

## Install

```
pip install git+https://github.com/linlabcode/CRC.git
```

## Usage

```
crc -e <enhacer_text_file> -g <genome_build> -s <peaks_bed_file> -c <chromosomes_dir> -o <out_dir> -n <name>
```

