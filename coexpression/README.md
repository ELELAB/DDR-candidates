# Introduction

This directory contains scripts to 1) calculate Spearman correlation between
expression profiles, from GTEX tissue data 2) process the calculated correlations
to keep only those with abs(r) >= 0.6, 3) further annotate the dataset with associated
GO terms

# Reproducing

## Requirements

in order to reproduce the results you will need:

  - a Linux system (we used Ubuntu Server 18.04)
  - R 3.5 with packages TCGAbiolinks (2.10.5), SummarizedExperiment (1.12.0), recount (1.8.2),
    TCGAutils (1.2.2), biomaRt (2.38.0)
  - Python 2.7.15 with pandas (0.24.2), Biopython (1.75), goatools (0.0.0)

## Running

from inside this folder, just run the following in this order:

1) Download the datasets from GTEX:

```
Rscript download.R
```

This downloads the GTEX transcription data in the form of .Rdata files.

2) calculate correlation values between pairs of genes for each tissue:

```
Rscript calc_corr.R
```

this calcualtes correlation between expression profiles, and save the result
as a matrix in a corr_TYPE.txt file, one for each TYPE of cancer. This is the
slowest step of the procedure and can take several hours.

3) perform filtering of the dataset, keeping pairs with correlation values >0.6
and keep those correlated genes found highly correlated in at least 10 tissues

```
python do_overlap_tissues_corrs.py
```

the output is a IID-like database (coexpression-interactions_0.6_min10.csv) as
a CSV file

4) annotate gene pairs by GO terms and filter out and keep only those gene pairs
in which at least one has a DNA-damage repair association terms, and other doesn't:

```
python do_annotate_coexpression.py
```

the output is the final data CSV file (coexpression-interactions_0.6_min10.csv.go_annotated.csv)


