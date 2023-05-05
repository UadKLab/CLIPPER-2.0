# CLIPPER 2.0
Peptide level annotation

Install **python 3.11.3**, and setup with *requirements.txt* through *pip*.

## Program init
- Init function

## Program preprocess class - tasks/functions:
- Infer software: use column headings
- Infer format
- Load data function: Load file to df, make a copy we can work on
- Sanitize/preprocess
- Initialize generic annotation
- Initialize degradomics pipeline annotation
- Pattern assessment and storing
- Read condition file if present, otherwise get channels/quant columns used
- Load databases here

## Entry class
- Threading or not
- Parse accession
- Parse sequence
- Get record
- Store general information (parse general)
- Parse cleavage (annotate cleavage environment)
- Parse protease cleavage uniprot
- Parse protease cleavage MEROPS
- Exopeptidase activity

## Statistics class - general functions
- Based on conditions, do general statistics: mean, CV, std
- Based on conditions, do percentile fold
- Ttest or ANOVA based on arguments and number of conditions
- Multiple testing correction

## Visualize class
- General statistics: Numbers about quantified, identified etc.
- CV plots
- Pie charts
- Heatmaps
- Clustermaps
- Volcano plots
- Fold plots
- Fold nterm
- Gallery
- Create logos (different class)

## Annotator utils functions
- initialize 
- parse arguments
- map dict ?
- write files

## To do list:
- Refactor code (made branch annotator_refactor with a folder for refactored code)
- Spectronaut support
- Dimethylation support
- Dimensionality reduction (PCA & UMAP)
- Multiple testing correction
- Protein atlas integration
- Visualize cleavages in sequences
- Visualize cleavages in structure
- Add secondary structure annotation
- Gprofiler - metascape GO enrichment
- Pathway annotation/significance
- Protease prediction (PSSM/GOtosubstrates)
- Fix condition bug in volcano plots
