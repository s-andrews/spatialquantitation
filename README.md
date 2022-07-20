# Spatial Transcriptomics Quantitation

This is a program which takes in a mapped BAM file from a 10X style 2D spatial transcriptomics dataset where the cell identity has been embedded into the read ID.  It produces a count matrix from a GTF file of gene annotations.

The mapping of reads to genes is done from the gene level, counting both exonic and intronic overlaps.  The association is directional, requiring that the read and gene are in the same orientiation.  The program uses the position of the match along with the sequence of the embedded UMI to deduplicate the counts.

## Expected Read ID structure
The program is designed to work with single-end mapped data where the read IDs have the following structure:

```
A00489.3.1166.25382.6464:23:20:23:02:CTGAGTAT
                               XXXXX YYYYYYYY
```

Where X is the cell ID (x,y coordinates for the spatial barcoding) and Y is the UMI sequence.

## Usage
```
$ ./quantitate_spatial_scrna.py  --help
usage: quantitate_spatial_scrna.py [-h] [--quiet] gtf_file bam_file outfile

Quantitate single cell spatial transcriptomic data

positional arguments:
  gtf_file    The GTF file of features
  bam_file    The BAM file to quantitate
  outfile     The output file to write to

optional arguments:
  -h, --help  show this help message and exit
  --quiet     Suppress progress messages
```
