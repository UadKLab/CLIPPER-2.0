# CLIPPER 2.0

<p align="center">
  <img src="img/CLIPPER_logo.png">
</p>

Advanced peptide-level annotation software that can be used natively in a pipeline with Proteome Discoverer and Spectromine/Spectronaut. Many other MS software such as fragpipe can be used by simply modifying column names to fit those of either PD or spectronaut. We found MaxQuant to have certain limitations when it comes to degradomics data analysis, and therefore recommend Fragpipe as a free alternative.

---

## Getting Started

This document will guide you through the process of setting up and using CLIPPER 2.0, an advanced tool for peptide annotation and analysis of proteomics data utilizing various databases and visualization tools.

Please note that CLIPPER 2.0 is currently in beta, and we welcome any bug reports or feature requests at konka@dtu.dk or alemol@dtu.dk.

## Prerequisites

- Conda (Either Anaconda or Miniconda will work)
  - https://www.anaconda.com/
  - https://docs.conda.io/en/latest/miniconda.html
- [Graphvis](https://www.graphviz.org/download/), tested version is v7.1.0 (optional for enhanced layout options)
- [alphafold] (https://alphafold.ebi.ac.uk/download#swissprot-section) Download Swiss-Prot (CIF files) from alphafold if you wish to use structure annotation or protein plotting features.

## Installation

The following should be done in a command line interface, but the code below can be copy/pasted one at a time. If problems arise, feel free to contact konka@dtu.dk and alemol@dtu.dk.

1. **For Mac only:** Some users have experienced the need for xcode command line developer tools for the setup to succeed. If you experience any problems related to the installation run the following and retry what you were doing:

```bash
xcode-select --install
```

2. Open the Anaconda Prompt on windows or terminal on mac, and create a conda environment with Python version 3.11.3:

```bash
conda create -n clipper python=3.11.3
```

3. Activate the conda environment with:

```bash
conda activate clipper
```

4. (Optional) To take advantage of the secondary structure calculation and structure plotting, install *pymol* with *conda*:

```bash
conda install -c conda-forge pymol-open-source
```

5. **For Mac only:** Install pycairo from the conda repository

```bash
conda install -c conda-forge pycairo
```

6. (Optional) To enable additional (better) pathway layout features, install Graphviz (v7.1.0) separately (see prerequisites). Then install *pygraphviz* afterwards with:

```bash
conda install --channel conda-forge pygraphviz
```

7. Download the Clipper 2.0 application.

8. From the commandline navigate to the Clipper 2.0 application (the folder with the "*requirements.txt*" file) and install dependencies with *pip*:

```bash
pip install -r requirements.txt
```

Verify your installation by running the following commands:

```bash
cd clipper
python run.py -h
```

9. If you wish to do solvent accessibility calculations or plot peptides in sequence or 3D protein views, a local copy of the alphafold database is necessary (https://alphafold.ebi.ac.uk/download#proteomes-section, approx. 5.2 gb). Once downloaded, the path should be specified in the file "annutils.py" as alphafold_folder_name = r'INSERT_PATH_HERE' on line 24.

10. If you wish to set up a server with the browser GUI, and have emails sent with completed results, a valid gmail email and password combination should be specified in "clipper/data/credentials.json".

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
usage: CLIPPER 2.0 [-h] -i INFILE [-it INFILE_TYPE] [-pa] [-cf CONDITIONFILE] [-a ALPHA] [-sw SOFTWARE] [-l LEVEL] [-dn] [-fn FILLNA] [-st SLEEPTIME] [-nx] [-nm] [-cs CALCSTRUCTURE] [-c THREADINGCORES] [-stat]
                   [-spw] [-sig SIGNIFICANCE] [-mt] [-mtmethod MULTIPLETESTINGMETHOD] [-vis] [-logo LOGO] [-logo_fc LOGO_FC] [-cle CLEAVAGESITESIZE] [-vc_fc VOLCANO_FOLDCHANGE VOLCANO_FOLDCHANGE]
                   [-psc PSEUDOCOUNTS] [-clvis CLEAVAGEVIS] [-enr] [-path] [-pf PROTEASEFILE] [-o OUTPUT_NAME] [-ot OUTPUT_FILETYPE] [-sep] [-pv]

Peptide annotation and analysis of proteomics data utilizing databases and visulization tools

options:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        Input peptide group result file.
  -it INFILE_TYPE, --infiletype INFILE_TYPE
                        File type of input file. Accepted values: [excel|csv|infer].
  -pa, --preannotated   Flag, if given assume that file is preannotated, and will not annotate from uniprot, proteinatlas, or do exopeptidasecheck.
  -cf CONDITIONFILE, --conditionfile CONDITIONFILE
                        A .txt file which maps labels to conditions. Adds columns for fold change for each pairwise comparison, average and CV of each condition to the dataframe. Each line must start with the
                        condition name followed by the channels used, separated by a single space (see example files for examples).
  -a ALPHA, --alpha ALPHA
                        Float that sets the alpha value to be used for all relevant statistical thresholds and significance measurements. Default: 0.05.
  -sw SOFTWARE, --software SOFTWARE
                        Software the data were analyzed with. Accepted values: [pd|sm|infer].
  -l LEVEL, --level LEVEL
                        Filtering on specified level of N-termini, or labelled N-termini. Accepted values: [all|nterm|quant]. Default: all.
  -dn, --dropna         Flag to indicate whether to filter for empty quant rows.
  -fn FILLNA, --fillna FILLNA
                        Value to fill empty quant rows or cells with (must be a float). Must be an integer or float.
  -st SLEEPTIME, --sleeptime SLEEPTIME
                        Float that determine interval in seconds between Biopython queries to Uniprot. Default: 0.2.
  -nx, --noexo          Flag, if given do not check for dipeptidase or aminopeptidase activity.
  -nm, --nomerops       Flag, if given do not check for cleavage site annotation in MEROPS database.
  -cs CALCSTRUCTURE, --calcstructure CALCSTRUCTURE
                        Annotate cleavage site solvent accessibility and secondary structure for the significant peptides in the dataset using Pymol and Alphafold models Accepted values: [all|sig].
  -c THREADINGCORES, --cores THREADINGCORES
                        integer, Allows specifying the number of CPU cores used for gathering information from Uniprot. Accepted values are 'max' (using all cores) or an integer indicating the intended number of
                        cores to be used. Default is 'max'. The speed of Uniprot information fetching scales linearly with the number of cores used.
  -stat, --statistic    Flag, if given perform statistical significance testing. Student T-test for two conditions, ANOVA from three or more. In addition it provides multiple testing corrected p-values using the
                        Benjamini/Hochberg method.
  -spw, --stat_pairwise
                        Flag, if given perform statistical significance t-test for all conditions pairwise.
  -sig SIGNIFICANCE, --significance SIGNIFICANCE
                        Performs fold change distribution significance check for all conditions pairwise Accepted values: [all|nterm].
  -mt, --multipletesting
                        Flag, if given multiple testing corrected values are used to determine significance for all relevant procedures using the method defined in -mtmethod.
  -mtmethod MULTIPLETESTINGMETHOD, --multipletestingmethod MULTIPLETESTINGMETHOD
                        Sets the multiple testing method to be used. Accepted values: [None|bonferroni|sidak|holm-sidak|holm|simes-hochberg|hommel|fdr_bh|fdr_by|fdr_tsbh|fdr_tsbky]. For clarification see here:
                        https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html. Default: fdr_bh.
  -vis, --visualize     Flag, if given draws various plots based on conditions passed from input and statistical tests
  -logo LOGO, --logo LOGO
                        Draws various logo plots based on all peptides or condition significant peptides. Values supported: [all|prob|pssm|shannon|kbl].
  -logo_fc LOGO_FC, --logofoldchange LOGO_FC
                        Float, sets the peptide abundance fold change between conditions which is considered for the logo generation. Example: -logo_fc 3 --> Filters peptides with log2 fold change > 3 (high logo)
                        and < 3 (low logo). Requires that -logo and -stat is given. Default: 3.
  -cle CLEAVAGESITESIZE, --cleavageenvironment CLEAVAGESITESIZE
                        Integer, size of the cleavage environment on both sides of the cleavage site. Default: 4.
  -vc_fc VOLCANO_FOLDCHANGE VOLCANO_FOLDCHANGE, --volcanofoldchange VOLCANO_FOLDCHANGE VOLCANO_FOLDCHANGE
                        Float, sets the cutoff value for fold change in volcano plots. Example: -vc_fc 1.5, which colors peptides with log2 fold change > 1.5 and < -1.5. P-value threshold is set by -a. Default:
                        1.5.
  -psc PSEUDOCOUNTS, --pseudocounts PSEUDOCOUNTS
                        Flag, if given add pseudocounts to normalized matrix calculation for logo generation.
  -clvis CLEAVAGEVIS, --cleavagevisualization CLEAVAGEVIS
                        Whether to visualize significant cleavages in sequence or both structure and sequence. For a complete dataset this can take a VERY long time, especially if set to 'both'. When using this
                        argument please consider only adding proteins/peptides to the input file which you are really interested in (might be based on high fold change, pvalue, or prior knowledge). Accepted
                        values: [None|seq|both]. Default: None.
  -enr, --enrichment    Flag, if given draws heatmap plots based on conditions passed from input and statistical tests for enrichment of GO terms and KEGG pathways with the gProfiler API.
  -path, --pathway      Flag, if given draws pathway plots based on conditions passed from input and statistical tests and maps detected proteins and peptides.
  -pf PROTEASEFILE, --proteasefile PROTEASEFILE
                        A file with protease MEROPS identifiers to predict activity of. Using weighted PSSM. See examples for examples.
  -o OUTPUT_NAME, --output_name OUTPUT_NAME
                        File name of output folder and annotated output file. Defaults to a timestamp for the initialization of Clipper 2.0.
  -ot OUTPUT_FILETYPE, --output_filetype OUTPUT_FILETYPE
                        File type of output file Accepted values: [xlsx|csv|tsv|pkl|json].
  -sep, --separate      Flag, if given the output file will not contain the original columns of the inputfile (separate original and annotated columns).
  -pv, --pymol_verbose  Flag, if given the commandline will display info messages from pymol when used. Looks pretty ugly and unnecessary in this context, so muted by default

Not extensively tested, this tool is still in beta version. Contact konka@dtu.dk or alemol@dtu.dk for bug reports and requests.
```

### Run the GUI in local host instead of the command line interface
You can run the GUI in local host by running the following command using the terminal in the clipper folder:

```bash
python app.py
```

You should be able to access the GUI by opening a web browser and navigating to http:// ip address and the port number specified in the terminal. We recommend using Google Chrome for the best experience.

You will still be able to monitor the progress of the analysis in the terminal, and the results will be downloaded once the analysis is complete. The email feature will not work in the GUI unless you specified your own email as mentioned above. We are working on hosting the GUI on a server, and will update the documentation when this is available.

## Input files

### Condition file (optional, but required for statistical tests and most visualizations)
The condition file is a text file where each line represents a condition. The first string on the line is the name of the condition, and the rest of the strings are space-separated specific identifyers for the columns corresponding to that condition (An example is found in tests/...).

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

Here are some examples of how you can use CLIPPER 2.0 throught the command line interface:

1. Basic usage

```bash
python run.py -i ../tests/HUNTER_clean_100.xlsx
```

2. Including pairwise statistical significance t-tests and fold change significance checks:

```bash
python run.py -i ../tests/HUNTER_clean_100.xlsx -cf ..\tests\cond_HUNTER.txt -sig all -stat -spw
```

3. Adding visualizations like volcano plots, dimensionality reduction and heatmaps:

```bash
python run.py -i ../tests/HUNTER_clean_100.xlsx -cf ../tests/cond_HUNTER.txt -stat -spw -vis
```

4. Adding gene enrichment and pathway analysis and visualization:

```bash
python run.py -i ../tests/HUNTER_clean_100.xlsx -cf ../tests/cond_HUNTER.txt -sig all -stat -spw -vis -path -enr
```

5. Adding cleavage site solvent accessibility and secondary structure annotation:

```bash
python run.py -i ../tests/HUNTER_clean_100.xlsx -cf ../tests/cond_HUNTER.txt -cs all -sig all -stat -spw -vis
```

6. Adding both sequence and structural visualization of cleavage sites:

```bash
python run.py -i ../tests/HUNTER_clean_100.xlsx -cf ../tests/cond_HUNTER.txt -cs all -sig all -stat -spw -vis -clvis both
```

7. Predicting cleavages for specified proteins using the protease file:

```bash
python run.py -i ../tests/HUNTER_clean_100.xlsx -cf ../tests/cond_HUNTER.txt -stat -spw -pf ../tests/proteases.txt
```

## Description of the output
Results are saved in a folder with the name of the input file and a timestamp (also saved as a zipped folder), unless and output folder name is specified. 

Depending on the arguments used, the output folder will contain the following files and folders:

1. **Annotated file**: A file containing the original input data with added columns for annotation and statistical tests. The file is saved in the format specified by the user (default is .xlsx). The added columns to the original files and their descriptions are: 
  - **Uniprot annotation**: The Uniprot annotation of the protein the peptide is derived from.
  - **Protein Atlas annotation**: The Protein Atlas annotation of the protein the peptide is derived from.
2. **Plots**: A number of folders are generated, containing plots specified.
  - **General plots**: Volcano plots of the fold change and p-values of the peptides.
  - **Logo plots**: Sequence logo plots of the peptides.
  - **Heatmap plots**: Heatmap plots of the peptides.
  - **Pathway plots**: Pathway plots of the peptides.
  - **Enrichment plots**: Enrichment plots of the peptides.
  



## Extra tips

If you wish to plot using specific filter requirements, we recommend performing an initial annotation on the full dataset without visualizations, filter the peptides based on your preference (that might be multiple testing corrected pvalues), delete the rows containing peptides which does not satisfy your criteria, and run CLIPPER 2.0 again.

When repeatedly running the same file, for example if you wish to test different argument configurations, it is recommended to specify the output folder as "-o FOLDER_NAME". CLIPPER 2.0 will look in the output folder for previous Uniprot annotation data, and if present (if the annotation was run in the same folder previously), CLIPPER 2.0 will not fetch data from uniprot, but reuse the present annotation file, speeding up processing time a lot.

## Contact

<p align="center">
  <a href="https://www.dtu.dk/">
    <img src="img\DTU_logo.png" alt="DTU logo" width="400"/>
  </a>
  <a href="https://www.bioengineering.dtu.dk/research/research-sections/section-for-protein-science-and-biotherapeutics">
    <img src="img\PSB_logo.png" alt="PSB logo" width="400"/>
  </a>
</p>

We hope you find CLIPPER 2.0 useful for your research. Feel free to contact us for any questions, bug reports, or feature requests (mails konka@dtu.dk and alemol@dtu.dk).
