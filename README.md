# rabies_barcode_tech

This repository contains code used for "High-Complexity Barcoded Rabies Virus for Scalable Circuit Mapping Using Single-Cell and Single-Nucleus Sequencing" by Shin & Urbanek et al.

Each script is numbered based on its relative location in the analysis workflow. The main pipeline is divided into a set of scripts for analysis the diversity of our barcoded rabies virus (RVdG) libraries (scripts 1-6), and another set for analyzing experimental RVdG data (scripts 7-16). Additionally, there are some extraneous scripts (script 5, 13, 17-21) which aren't part of the standard RVdG workflow, but were used for analyses included within the manuscript.

Input and output files for each step of the pipeline are summarized below:

![alt text](https://github.com/MEUrbanek/rabies_barcode_tech/blob/main/images/pipeline_schematic.png)

Python modules and versions used for all analyses excluding 13_empty_droplet_modeling.ipynb are listed in the "requirements.txt" file housed in this repository. For 13_empty_droplet_modeling.ipynb, python modules and versions are listed in "empty_droplet_modeling_requirements.txt". All R packages and their versions are listed in the "renv.lock" file.
