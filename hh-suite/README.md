## Introduction

This directory contains an example script to filter multiple sequence alignments (MSAs)
generated in HHblits for probability plots shown in Figure 1f, 2f, 3b, and 4b. 

## Reproducing

### Requirements
In order to redo the MSA filtering, you will need:
  - a Linux system (We used a server with Ubuntu Server 18.04)
  - installation of the HH-suite3 package  (v3.3.0).
  - installation of Databases Uniclust30, pfam70_35, pdb70, scop70_1.45
  - an MSA of SPIDR_OB3 in a3m format (located in this directory), or produce
    it initially by running HHblits using the human SPIDR protein
    sequence (776-91) as a search query.

### Running
Once allocated to a defined folder, load the SPIDR_OB_3_2 alignment and follow
the instructions in the `SPIDR_OB3_filtering.txt` file.

