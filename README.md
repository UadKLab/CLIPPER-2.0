# CLIPPER 2.0
Peptide level annotation

Install **python 3.11.3**, and setup with *requirements.txt* through *conda*.

Create a conda environment:

`conda create --n clipper python=3.11.3`

Install packages with pip:

`pip install -r requirements.txt`

Then install pymol with conda:

`conda install -c conda-forge pymol-open-source`

If successfully install, you should be able to run the following command:

`python clipper/run.py -h`

Example run:

`python clipper/run.py -i ..\tests\HUNTER_clean_100.xlsx -cf ..\tests\cond_HUNTER.txt -cs all -sig all -stat -spw -vis -clvis both`

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
- [ ] Refactor code (made branch annotator_refactor with a folder for refactored code)
    - [x] Rewritten for readability and modularity, but not refactored as described above
    - [x] Clean up destination folders and structured output better
    - [x] Clean up paths and use Path for relative paths
- [x] Spectronaut support
- [x] Dimethylation support
    - [x] Dimethylation support for Spectromine
- [x] Dimensionality reduction (PCA & UMAP)
    - [x] Added PCA
    - [x] Added UMAP
- [x] Multiple testing correction
- [x] Protein atlas integration
- [x] Visualize cleavages in sequences
    - [x] Added sequence plots
    - [x] Add colormap 
    - [x] Merge with structure visualization
- [x] Maybe only annotate MEROPS code in sequence visualization if within the peptides identified
- [x] Visualize cleavages in structure
    - [x] Downloaded Alphafold EBI database
    - [x] Plot structures
- [ ] Remove peptide duplicates before plotting on seq and structure
- [x] Add secondary structure annotation
- [x] Add solvent accessibility annotation
- [ ] Modify ss and sasa calculations to read each protein once
- [x] Gprofiler - metascape GO enrichment
- [ ] Pathway annotation/significance
- [ ] Protease prediction (PSSM/GOtosubstrates)
- [x] Fix condition bug in volcano plots
    - [x] Have not managed to recreate it yet, need to use input data with multiple (>4) conditions. Did not manage even with 6 conditions.
- [x] Static Uniprot entry lookup
    - [x] Will not be implemented