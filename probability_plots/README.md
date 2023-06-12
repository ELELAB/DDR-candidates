# Introduction

This directory contains script and data to generate plots for
paper figure 1 panels g-j,l-o.

# Reproducing

## Requirements

In order to redo the plotting, you will need:

  - a Linux system (We used a server with Ubuntu Server 18.04)
  - R 4.2
  - R packages igraph (1.4.2), extrafont (0.19), graphics, pracma (2.4.2)
  - the Arial font installed on your system (e.g. via the msttcorefonts 
    Debian/Ubuntu package)

## Running

from inside the folder, just run:

```
Rscript do.R
```

this will generate a set of pdf files, each containing a different panel.
The runtime is up to minutes.