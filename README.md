# CLIPPER 2.0
Advanced peptide-level annotation software

---

## Getting Started

This document will guide you through the process of setting up and using CLIPPER 2.0, an advanced tool for peptide annotation and analysis of proteomics data utilizing various databases and visualization tools.

Please note that CLIPPER 2.0 is currently in beta, and we welcome any bug reports or feature requests at konka@dtu.dk or alemol@dtu.dk.

## Prerequisites

* Python 3.11.3
* [Graphvis](https://www.graphviz.org/download/), tested version is v7.1.0 (optional for enhanced layout options)

## Installation

1. Install Python 3.11.3 and create a conda environment:

```bash
conda create -n clipper python=3.11.3
```

2. Change to your newly created conda environament and install packages with *pip*:

```bash
pip install -r requirements.txt
```

3. To take advantage of the secondary structure calculation and structure plotting, install *pymol* with *conda*:

```bash
conda install -c conda-forge pymol-open-source
```

4. (Optional) To enable additional pathway layout features, install Graphviz (v7.1.0) separately. Then install *pygraphviz*:

```bash
conda install --channel conda-forge pygraphviz
```

Verify your installation by running the following commands:

```bash
cd clipper
python run.py -h
```

## Usage

You can use CLIPPER 2.0 by executing the 'run.py' script with specific arguments.

```bash
python run.py -i INFILE
```

This script will read an input file (INFILE) and provide peptide annotation and analysis. The command-line arguments allow you to control and customize the software's behavior.

For a complete list of available arguments and their description, run:

```bash
python run.py -h
```

### Output and available arguments:

```
usage: CLIPPER 2.0 [-h] -i INFILE [-it INFILE_TYPE] [-sw SOFTWARE] [-l LEVEL] [-dn] [-fn FILLNA] [-st SLEEPTIME] [-nx] [-nm] [-cs CALCSTRUCTURE] [-sc] [-cf CONDITIONFILE] [-stat] [-spw] [-sig SIGNIFICANCE] [-vis]
                   [-logo LOGO] [-psc PSEUDOCOUNTS] [-clvis CLEAVAGEVIS] [-enr] [-path] [-pf PROTEASEFILE] [-o OUTPUT_NAME] [-ot OUTFILE_TYPE] [-sep]

Peptide annotation and analysis of proteomics data utilizing databases and visulization tools

options:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        Input peptide group result file
  -it INFILE_TYPE, --infiletype INFILE_TYPE
                        File type of input file. Accepted values: excel/csv/infer
  -sw SOFTWARE, --software SOFTWARE
                        Software the data were analyzed with. Accepted values: pd/sm/infer
  -l LEVEL, --level LEVEL
                        Filtering on specified level of N-termini, or labelled N-termini. Accepted values: all/nterm/quant.
  -dn, --dropna         Flag to indicate whether to filter for empty quant rows
  -fn FILLNA, --fillna FILLNA
                        Value to fill empty quant rows or cells with
  -st SLEEPTIME, --sleeptime SLEEPTIME
                        Float that determine intervals between Biopython queries to Uniprot
  -nx, --noexo          Do not check for dipeptidase or aminopeptidase activity
  -nm, --nomerops       Do not check for cleavage site annotation in MEROPS database
  -cs CALCSTRUCTURE, --calcstructure CALCSTRUCTURE
                        Annotate cleavage site solvent accessibility and secondary structure for the significant peptides in the dataset using Pymol and Alphafold models Accepted values: all/sig.
  -sc, --singlecpu      Use a single process instead of threading for annotation
  -cf CONDITIONFILE, --conditionfile CONDITIONFILE
                        Map labels to conditions. Adds columns for fold change for each pairwise comparison, average and CV of each condition to the dataframe. Each line must start with the condition name followed by the     
                        channels used, separated by a single space
  -stat, --statistic    Performs statistical significance testing. Student T-test for two conditions, ANOVA from three or more
  -spw, --stat_pairwise
                        Performs statistical significance t-test for all conditions pairwise
  -sig SIGNIFICANCE, --significance SIGNIFICANCE
                        Performs fold change distribution significance check for all conditions pairwise Accepted values: all/nterm.
  -vis, --visualize     Draws various plots based on conditions passed from input and statistical tests
  -logo LOGO, --logo LOGO
                        Draws various logo plots based on all peptides or condition significant peptides. Values supported: [all|prob|pssm|shannon|kbl]
  -psc PSEUDOCOUNTS, --pseudocounts PSEUDOCOUNTS
                        Add pseudocounts to normalized matrix calculation
  -clvis CLEAVAGEVIS, --cleavagevisualization CLEAVAGEVIS
                        Whether to visualize significant cleavages in sequence or both structure and sequence Accepted values: seq/both.
  -enr, --enrichment    Draws heatmap plots based on conditions passed from input and statistical tests for enrichment of GO terms and KEGG pathways with the gProfiler API
  -path, --pathway      Draws pathway plots based on conditions passed from input and statistical tests and maps detected proteins and peptides
  -pf PROTEASEFILE, --proteasefile PROTEASEFILE
                        Protease MEROPS identifiers to predict activity of. Using weighted PSSM
  -o OUTPUT_NAME, --output_name OUTPUT_NAME
                        File name of output folder and annotated output file
  -ot OUTFILE_TYPE, --outfile_type OUTFILE_TYPE
                        File type of output file Accepted values: xlsx/csv/tsv/pkl/json
  -sep, --separate      Whether to merge or keep annotation as a separate file. False by default

Not extensively tested, this tool is still in beta version. Contact konka@dtu.dk for bug reports and requests.
```

## Input files

###Condition file (optional, but required for statistical tests and most visualizations)
The condition file is a text file where each line represents a condition. The first string on the line is the name of the condition, and the rest of the strings, space-separated, are columns corresponding to that condition.

Example of condition file format:
    
```	
Condition1 Column1a Column1b Column1c
Condition2 Column2a Column2b Column2c
...
```

###Protease file (optional, required for protease activity prediction)

The protease file is a text file containing one protease MEROPS code per line. These codes correspond to a list of proteins for which you want to predict cleavages.

Example of protease file format:

```	
MEROPSProteaseCode1
MEROPSProteaseCode2
...
```

## Examples

Here are some examples of how you can use CLIPPER 2.0:

1. Basic usage

```bash
python run.py -i ..\tests\HUNTER_clean_100.xlsx
```

2. Including pairwise statistical significance t-tests and fold change significance checks:

```bash
python run.py -i ..\tests\HUNTER_clean_100.xlsx -cf ..\tests\cond_HUNTER.txt -sig all -stat -spw
```

3. Adding visualizations like volcano plots, dimensionality reduction and heatmaps:

```bash
python run.py -i ..\tests\HUNTER_clean_100.xlsx -cf ..\tests\cond_HUNTER.txt -stat -spw -vis
```

4. Adding gene enrichment and pathway analysis and visualization:

```bash
python run.py -i ..\tests\HUNTER_clean_100.xlsx -cf ..\tests\cond_HUNTER.txt -sig all -stat -spw -vis -path -enr
```

5. Adding cleavage site solvent accessibility and secondary structure annotation:

```bash
python run.py -i ..\tests\HUNTER_clean_100.xlsx -cf ..\tests\cond_HUNTER.txt -cs all -sig all -stat -spw -vis
```

6. Adding both sequence and structural visualization of cleavage sites:

```bash
python run.py -i ..\tests\HUNTER_clean_100.xlsx -cf ..\tests\cond_HUNTER.txt -cs all -sig all -stat -spw -vis -clvis both
```

7. Predicting cleavages for specified proteins using the protease file:

```bash
python run.py -i ..\tests\HUNTER_clean_100.xlsx -cf ..\tests\cond_HUNTER.txt -stat -spw -pf ..\tests\proteases.txt
```

We hope you find CLIPPER 2.0 useful for your research. Feel free to contact us for any questions, bug reports, or feature requests.

# Development notes

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
- pathway enrichment
- secondary structure and solvent accessibility calculation

## To do list:
- [x] Refactor code (made branch annotator_refactor with a folder for refactored code)
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
- [x] Remove peptide duplicates before plotting on seq and structure
- [x] Add secondary structure annotation
- [x] Add solvent accessibility annotation
- [x] Modify ss and sasa calculations to read each protein once
- [x] Gprofiler - metascape GO enrichment
- [x] Pathway annotation/significance
- [x] Protease prediction (PSSM/GOtosubstrates)
- [x] Fix condition bug in volcano plots
    - [x] Have not managed to recreate it yet, need to use input data with multiple (>4) conditions. Did not manage even with 6 conditions.
- [x] Static Uniprot entry lookup
    - [x] Will not be implemented
- [ ] Update and prepare for publication
- [x] Write doc strings and prepare documentation
- [ ] Write tests
- [ ] Update application and website
- [ ] Figure out hosting
- [x] Clean up structure of folders
- [x] Rename folders and scripts to make more sense i.e. Clipper 2.0 instead of annotator for base folder, app_main.py instead of general.py
- [x] Merge with main
- [ ] Add support for other software like Fragpipe
- [ ] Use Path instead of os.path.join
- [x] Figure out way of pointing to Alphafold database that is device agnostic, and easy to change from user side
- [ ] Show column patterns in command line
- [ ] Make docs and examples for software export, column descriptions, etc.
- [x] Fix bug where model is not found in Alphafold database but still tries to read it
- [x] Delete tmp folder in the cmd result folder
- [ ] Change descriptions to reflect support for other software and new arguments
- [ ] Add C-terminomics support
- [ ] Add C-terminal exopeptidase activity check
- [ ] Add heatmap with cleavage site specificity to complement logos
- [ ] Add pycache to git ignore
- [ ] Create github issues instead of writing things to do here
- [ ] Add cleavage environemnt length as argument (right now is only 4)
- [ ] Why is log10 pvalue in the annotated dataframe negative instead of -log10?
- [ ] Add argument to control p-value and fold cutoff for logos
- [ ] Add argument to control p-value and fold cutoff for volcano plots

