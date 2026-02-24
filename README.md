# CLIPPER 2.0
[![CI (dev)](https://github.com/UadKLab/CLIPPER-2.0/actions/workflows/ci.yml/badge.svg?branch=dev)](https://github.com/UadKLab/CLIPPER-2.0/actions/workflows/ci.yml)

<p align="center">
  <img src="img/CLIPPER_logo.png">
</p>

Advanced peptide-level annotation software that can be used natively in a pipeline with Proteome Discoverer and Spectromine/Spectronaut. Many other MS software such as fragpipe can be used by simply modifying column names to fit those of either PD or spectronaut. We found MaxQuant to have certain limitations when it comes to degradomics data analysis, and therefore recommend Fragpipe as a free alternative.

---

## Getting Started

This document will guide you through the process of setting up and using CLIPPER 2.0, an advanced tool for peptide annotation and analysis of proteomics data utilizing various databases and visualization tools.

Please note that CLIPPER 2.0 is currently in beta, and we welcome any bug reports or feature requests at konka@dtu.dk or alemol@dtu.dk.

## Prerequisites

- Conda (Anaconda or Miniconda)
  - https://www.anaconda.com/
  - https://docs.conda.io/en/latest/miniconda.html
- [Graphviz](https://www.graphviz.org/download/) (optional, for better pathway layouts)
- [AlphaFold Swiss-Prot models](https://alphafold.ebi.ac.uk/download#swissprot-section) (optional, for structure annotation and 3D cleavage visualization)

## Installation

Run the following from a terminal:

1. (macOS only, if needed) install command line tools:

```bash
xcode-select --install
```

2. Create and activate environment:

```bash
conda create -n clipper python=3.11.3
conda activate clipper
```

3. Clone repository and install dependencies:

```bash
git clone https://github.com/UadKLab/CLIPPER-2.0.git
cd CLIPPER-2.0
pip install -r requirements.txt
```

4. Optional extras:

```bash
conda install -c conda-forge pymol-open-source
conda install -c conda-forge pycairo        # macOS users typically need this for PyMOL/graphics stack
conda install -c conda-forge pygraphviz
```

5. Verify installation:

```bash
python clipper/run.py -h
```

6. Optional setup for structure annotation:
- Download AlphaFold Swiss-Prot models.
- Set `alphafold_folder_name` in `clipper/bin/annutils.py` to your local model path.

7. Optional setup for email notifications in web GUI:
- Configure `clipper/data/credentials.json` with valid credentials used by `clipper/bin/mail.py`.

## Usage

### Command line

Run CLIPPER from repository root:

```bash
python clipper/run.py -i INFILE
```

For complete, current arguments:

```bash
python clipper/run.py -h
```

Note: The CLI help output is the source of truth for supported options. It may evolve over time.

### Web GUI (Flask)

Run locally:

```bash
python clipper/app.py
```

Then open the URL shown in terminal (typically `http://127.0.0.1:5000`).

GUI behavior:
- Upload input and optional condition/protease files.
- Job runs in the same process and redirects to a completion/error page.
- On success, the generated zip is downloaded automatically.
- On failure, the app redirects to `/JOB_ID/error` and shows a specific error message.

Troubleshooting:
- If job initialization fails (e.g., missing/invalid input file), the app flashes a message and returns to the index page.
- For runtime errors, check:
  - the error page message
  - the generated log file under `clipper/log/`

## Input files

### Condition file (optional, but required for statistical tests and most visualizations)
The condition file is a text file where each line represents a condition. The first string on the line is the name of the condition, and the rest of the strings are space-separated specific identifiers for the columns corresponding to that condition (An example is found in tests/...).

**File name, condition name, and column names may NOT contain the characters “:” or “/“ or “.” as this MAY result in error.**

Example of condition file format for a triplicate experiment with two conditions:
    
``` 
Condition1 Column1a Column1b Column1c
Condition2 Column2a Column2b Column2c
```

### Protease file (optional, required for protease activity prediction)

The protease file is a text file containing one protease MEROPS code per line. These codes correspond to a list of proteins for which you want to predict cleavages.

Example of protease file format:

``` 
MEROPSProteaseCode1
MEROPSProteaseCode2
...
```

## Examples

Here are some examples of how you can use CLIPPER 2.0 through the command line interface:

1. Basic usage

```bash
python clipper/run.py -i ../tests/HUNTER_clean_100.xlsx
```

2. Including pairwise statistical significance t-tests and fold change significance checks:

```bash
python clipper/run.py -i ../tests/HUNTER_clean_100.xlsx -cf ../tests/cond_HUNTER.txt -sig all -stat -spw
```

3. Adding visualizations like volcano plots, dimensionality reduction and heatmaps:

```bash
python clipper/run.py -i ../tests/HUNTER_clean_100.xlsx -cf ../tests/cond_HUNTER.txt -stat -spw -vis
```

4. Adding gene enrichment and pathway analysis and visualization:

```bash
python clipper/run.py -i ../tests/HUNTER_clean_100.xlsx -cf ../tests/cond_HUNTER.txt -sig all -stat -spw -vis -path -enr
```

5. Adding cleavage site solvent accessibility and secondary structure annotation:

```bash
python clipper/run.py -i ../tests/HUNTER_clean_100.xlsx -cf ../tests/cond_HUNTER.txt -cs all -sig all -stat -spw -vis
```

6. Adding both sequence and structural visualization of cleavage sites:

```bash
python clipper/run.py -i ../tests/HUNTER_clean_100.xlsx -cf ../tests/cond_HUNTER.txt -cs all -sig all -stat -spw -vis -clvis both
```

7. Predicting cleavages for specified proteins using the protease file:

```bash
python clipper/run.py -i ../tests/HUNTER_clean_100.xlsx -cf ../tests/cond_HUNTER.txt -stat -spw -pf ../tests/proteases.txt
```

## Description of output files and folders
Results are saved in a folder with the name of the input file and a timestamp (also saved as a zipped folder), unless an output folder name is specified. 

Depending on the arguments used, the output folder will contain the following files and folders:

1. **Annotated file**: A file containing the original input data with added columns for annotation and statistical tests. The file is saved in the format specified by the user (default is .xlsx with '_annot' in filename). The added columns to the original files and their descriptions are: 
  - **query_sequence**: Peptide sequence used for the annotation.
  - **query_accession**: Uniprot protein accession used for the annotation.
  - **name**: Gene name of the protein.
  - **full_sequence**: Full sequence of the protein.
  - **description**: Description of the protein.
  - **keywords**: keywords of the protein from Uniprot entry in the same column,  separated by '|'.
  - **go_codes**: Gene ontology codes of the protein from Uniprot entry in the same column, separated by '|'.
  - **go_names**: Gene ontology names of the protein from Uniprot entry in the same column, separated by '|'.
  - **proteoform_certainty%**: The certainty of the proteoform annotation, based on the number of proteins the peptide sequence is present.
  - **acc_length**: The number of residues in the protein.
  - **start_pep**: The location of the peptide in the protein, as the first residue.
  - **end_pep**: The location of the peptide in the protein, as the last residue.
  - **p1_position**: The cleavage site of the peptide in the protein, as the site of the residue before cleavage (p1 position).
  - **cleavage_site**: The cleavage site environment of the peptide in the protein, annotated as the 4 residues before the cleavage site, the location
  of the p1 residue in parenthesis, a period to indicate the cleavage site, the location of the p1' residue, and the full peptide identified.
  - **p4_p4prime**: The cleavage environment of the peptide as 4 residues before and after the cleavage site.
  - **nterm_annot**: Annotation of the cleavage event based on UniProt annotation, if available.
  - **protease_uniprot**: Proteases known to generate this cleavage site, based on UniProt annotation.
  - **protease_merops**: MEROPS code of proteases known to generate this cleavage site, based on MEROPS annotation.
  - **protease_merops_name**: Name of proteases known to generate this cleavage site, based on MEROPS annotation.
  - **ProteinAtlas_RNA tissue specific nTPM**: Normalized transcript per million (nTPM) expression of the protein in tissues as provided by ProteinAtlas,
  separated by a comma in the same column.
  - **ProteinAtlas_Chromosome**: Location of the protein in the chromosome as provided by ProteinAtlas.
  - **ProteinAtlas_Position**: Position of the protein in the chromosome as provided by ProteinAtlas.
  - **ProteinAtlas_Protein class**: Protein class of the protein as provided by ProteinAtlas.
  - **ProteinAtlas_Biological process**: Biological process of the protein as provided by ProteinAtlas.
  - **ProteinAtlas_Molecular function**: Molecular function of the protein as provided by ProteinAtlas.
  - **ProteinAtlas_Disease involvement**: Disease involvement of the protein as provided by ProteinAtlas.
  - **exopeptidase**: Annotation of the peptide as a potential exopeptidase substrate based on identification of upstream peptide in the same dataset.
  - **condition_mean**: The mean of the peptide abundance in the condition.
  - **condition_deviation**: The standard deviation of the peptide abundance in the condition.
  - **condition_cv**: The coefficient of variation of the peptide abundance in the condition.
  - **Fold change:**: The fold change of the peptide abundance between conditions. The first of the two conditions is the numerator, and the second is the denominator.
  - **Log2 fold change:**: The log2 fold change of the peptide abundance between conditions. The first of the two conditions is the numerator, and the second is the denominator.
  - **Independent T-test p-value:**: The p-value of the independent T-test between the two conditions.
  - **-Log10 Independent T-test p-value:**: The -log10 p-value of the independent T-test between the two conditions.
  - **Fold condition comparison significance**: The significance of the fold change between the two conditions, based on the distribution of peptide fold changes
  between conditions (5% tail). The first of the two conditions is the numerator, and the second is the denominator.
  - **predicted_protease_activity**: The predicted activity of the protease based on the peptide sequence and a MEROPS database PSSM for every protease provided by 
  the user. Proteases are separated by a pipe symbol in the same column. The score is the sum of log2 fold enrichment compared to the background frequency of
  amino acids across the p4-p4' positions of the peptide.

2. **Plots**: A number of folders are generated, containing plots specified.
  - **General plots**: 
    - General plot with numbers of identified peptides and proteins, and N-termini.
    - Clustermap of the peptide abundances across individual replicates/TMT channels.
    - Heatmap of the peptide abundances across individual replicates/TMT channels.
    - CV (coefficient of variation) of the peptide abundances across conditions.
    - PCA plot of the peptide abundances across conditions.
    - UMAP plot of the peptide abundances across conditions.
    - Gallery of significant peptides. If a specific peptide in a protein is significant, peptides of that protein and the protein quantification is plotted.
  - **Logo plots**: Sequence logo plots of the peptides. These logos are generated with the method the user specifies in the input arguments.
  - **Fold change plots**: Peptide fold change distribution plots across conditions.
  - **Volcano plots**: Volcano plots of the peptide fold change and p-values across conditions.
  - **Piechart plots**: Categories of peptides identified and their distribution in the dataset.
  - **Pathway plots**: Pathway plots of the proteins and peptides identified in the dataset.
  - **Enrichment plots**: Heatmap plots of the gene ontology terms, KEGG pathways, and other databases, enriched in the dataset as determined by gProfiler.

3. **Log file**: A log file containing information about the analysis, including the arguments used, and the time it took to complete the analysis. Please include this as 
an attachment if you contact us with bug reports, if possible.

## Extra tips

If you wish to plot using specific filter requirements, we recommend performing an initial annotation on the full dataset without visualizations, filter the peptides based on your preference (that might be multiple testing corrected pvalues), delete the rows containing peptides which does not satisfy your criteria, and run CLIPPER 2.0 again.

When repeatedly running the same file, for example if you wish to test different argument configurations, it is recommended to specify the output folder as "-o FOLDER_NAME". CLIPPER 2.0 will look in the output folder for previous Uniprot annotation data, and if present (if the annotation was run in the same folder previously), CLIPPER 2.0 will not fetch data from uniprot, but reuse the present annotation file, speeding up processing time a lot.

## Contact

<p align="center">
  <a href="https://www.dtu.dk/">
    <img src="img\DTU_logo.png" alt="DTU logo" width="100" height="140"/>
  </a>
  <a href="https://www.bioengineering.dtu.dk/research/research-sections/section-for-protein-science-and-biotherapeutics">
    <img src="img\PSB_logo.png" alt="PSB logo" width="140" height="140"/>
  </a>
</p>

We hope you find CLIPPER 2.0 useful for your research. Feel free to contact us for any questions, bug reports, or feature requests (mails konka@dtu.dk and alemol@dtu.dk).
