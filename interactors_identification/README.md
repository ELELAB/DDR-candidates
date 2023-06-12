# Introduction

This directory contains a script for the identification of potential DDR candidate
genes from IID and Gene Onthology.

# Reproducing

## Requirements

in order to reproduce the results you will need, in this folder:

  - a Linux system (we used Ubuntu Server 18.04)
  - Python 2.7.15 with pandas (0.24.2), Biopython (1.75), goatools (0.0.0)

## Running

from inside the folder, first, decompress the IID database file:

```
gunzip -d iid.human.2018-05.txt.gz
```

then just run the following:

```
python do_iid.py
```

the output is a CSV list of interactions (IID_interactions.csv). The runtime
is up to a few minutes