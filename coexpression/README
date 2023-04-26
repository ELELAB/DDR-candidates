# Introduction

This directory contains scripts to 1) calculate Spearman correlation between
expression profiles, from GTEX tissue data 2) process the calculated correlations
to keep only those with abs(r) >= 0.6, 3) further annotate the dataset with associated
GO terms

# Reproducing

## Requirements

in order to reproduce the results you will need, in this folder:

  - R 3.5 with packages TCGAbiolinks 2.10.5, SummarizedExperiment, recount,
    TCGAutils, biomaRt
  - Python 2.7 with pandas, Biopython, goatools, 
  - a go.obo onthology definition file (the one we used was generated on 2018-10-24)
  - a goa_human.gaf human onthology file (the one we used was generated on 2018-05-22)

## Running

from inside this folder, just run the following in this order:

1) Download the datasets from GTEX:

```
Rscript download.R
```

2) calculate correlation values between pairs of genes for each tissue:

```
Rscript calc_corr.R
```

3) perform filtering of the dataset, keeping pairs with correlation values >0.6
and keep those correlated genes found highly correlated in at least 10 tissues

```
python do_overlap_tissues_corrs.py
```

the output is a IID-like database.

4) annotate gene pairs by GO terms and filter out and keep only those gene pairs
in which at least one has a DNA-damage repair association terms, and other doesn't:

```
python do_annotate_coexpression.py
```



