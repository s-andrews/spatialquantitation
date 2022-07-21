# Multi-barcoding Single Cell Quantitation

This is a program which takes in a mapped BAM file from a 10X style transcriptomics dataset labelled with multiple rounds of barcodes, where the cell identity has been embedded into the read ID.  It produces a count matrix from a GTF file of gene annotations.

The mapping of reads to genes is done from the gene level, counting both exonic and intronic overlaps.  The association is directional, requiring that the read and gene are in the same orientiation.  The program uses the position of the match along with the sequence of the embedded UMI to deduplicate the counts.

## Expected Read ID structure
The program is designed to work with single-end mapped data where the read IDs have the following structure:

```
A00489.3.1166.25382.6464:23:20:23:02:CTGAGTAT
                         XXXXXXXXXXX YYYYYYYY
```

Where X is the cell ID and Y is the UMI sequence.

## Usage
The input must be a BAM file which is positionally sorted.  An index file is not required.

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

## Installation
Simply doing a ```git clone``` of this repository will get you the main python script for the quantitation.

The script requires python >= 3.7.

The only package requirement outside the standard library is [pysam](https://github.com/pysam-developers/pysam).  This can be installed with:

```python3 -m pip install pysam```

This package is only available on unix operating systems so this script will not work on windows.

## Output
The output of the script is a tab delimited file showing the genes (rows) vs cells (columns) and the raw overlap counts as data.

```
Gene    01:01   01:02   01:03   01:04   01:05   01:06   01:07   01:08
5S_rRNA 0       0       0       1       0       0       1       0
7SK     734     593     337     2205    663     1965    1209    933
A1BG    0       2       0       6       0       0       0       1
A1CF    0       0       0       0       0       0       0       0
A2M     0       0       0       1       2       0       0       1
A2M-AS1 0       0       0       0       0       0       0       0
A2ML1   12      4       16      17      11      10      0       3
A2MP1   6       3       3       15      7       4       6       9
A3GALT2 0       0       0       0       0       0       0       0
A4GALT  42      13      17      38      21      24      15      30
AAAS    9       3       5       33      2       38      4       14
AACS    39      5       8       16      7       24      5       8
AACSP1  0       0       0       0       0       0       0       0
AADAC   0       0       0       0       0       0       0       1
AADACL2 0       0       0       0       0       0       0       0
AADACL3 7       1       4       0       0       0       0       0
AADACL4 0       0       0       0       0       0       0       0
```

