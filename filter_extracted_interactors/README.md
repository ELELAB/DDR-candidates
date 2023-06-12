# Introduction

This directory contain filtering for housekeeping genes and the Crapome database
of the 441 candidate gene pairs. 

# Reproducing

## Requirements

in order to reproduce the analysis you will need:

  - a Linux system (we used Ubuntu Server 18.04)
  - Python 2.7.15 with pandas (0.24.2), Biopython (1.75), goatools (0.0.0)

the run time should be in the order of minutes.

## Running

just run:

```
python do_filter.py
```

the output of this procedure is a filtered list of interactors in the form of
CSV files. For all of them, we filtered according to our housekeeping gene list
plus:

- nothing else (extracted_interactors_housekeeping.csv)
- Crapome frequently detected proteins in at least 10% of the experiments (extracted_interactors_housekeeping_crapome_10.csv)
- Crapome frequently detected proteins in at least 50% of the experiments (extracted_interactors_housekeeping_crapome_50.csv)
- Crapome all frequently detected proteins (extracted_interactors_housekeeping_crapome_all.csv)

