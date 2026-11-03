# Enza-P circulating tumour DNA analyses

This repository contains the custom code used to analyse copy number profiles, somatic variants, and structural variants in the paper "Circulating tumour DNA fraction and androgen receptor genotype predict benefit for addition of 177Lu-PSMA-617 to enzalutamide in advanced prostate cancer in the ENZA-p trial". 

Access to patient sequencing data in the EnzaP trial is restricted but may be shared upon request (see Data availability section of paper).

## Pre-requisites
Our analysis software and scripts are written using the Julia and Rust programming languages. To use our code, you will first need to:

1. Install Julia version 1.8.5: [https://julialang.org/downloads/](https://julialang.org/downloads/oldreleases/)  
2. Install the Rust compiler and toolchain: https://www.rust-lang.org/tools/install  
3. Add the custom Julia modules in subfolder `julia_modules` included here to a folder on your system
4. Add your folder to the LOAD_PATH, and install required Julia packages from inside the Julia terminal:

```
push!(LOAD_PATH, "/path/to/your/folder")
using Pkg;
Pkg.add(split("Printf Random StatsBase DelimitedFiles Statistics Distributions CSV Plot2 PyPlot PDFIO FastaIO HypothesisTests BlackBoxOptim PyCall KernelDensity Interpolations Combinatorics Suppressor"))
```
5. Install Mutato: https://www.github.com/annalam/mutato
6. Install Breakfast: https://www.github.com/annalam/breakfast
7. Install ANNOVAR: https://annovar.openbioinformatics.org/en/latest/user-guide/download/
8. Install VarDictJava, GRIDSS, and GRIPSS in your conda environment

