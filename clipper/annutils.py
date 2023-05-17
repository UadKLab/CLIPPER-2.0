import logging
import os
import re

import gzip
import tempfile

from argparse import ArgumentParser, HelpFormatter
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

import pymol

def initialize_logger(logfile):
    """Initializes the logger with a file handler and a console handler."""

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")

    file_handler = logging.FileHandler(logfile)
    file_handler.setFormatter(formatter)

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    logger.info(f"Annotator started at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.")

    return logger

def initialize_arguments():
    """Initializes a parser instance and adds arguments for annotator."""

    parser = ArgumentParser(
        prog="CLIPPER 2.0",
        description="Peptide annotation and analysis of proteomics data utilizing databases and visulization tools",
        epilog="Not extensively tested, this tool is still in beta version. Contact konka@dtu.dk for bug reports and requests.",
        formatter_class=HelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--infile",
        action="store",
        dest="infile",
        type=str,
        required=True,
        help="Input peptide group result file",
    )

    parser.add_argument(
    "-it",
    "--infiletype",
    default="infer",
    action="store",
    dest="infile_type",
    type=str,
    help="File type of input \
                        file. Accepted values: excel/csv/infer",
    )

    parser.add_argument(
        "-sw",
        "--software",
        default="infer",
        action="store",
        dest="software",
        type=str,
        help="Software the data were \
                            analyzed with. Accepted values: pd/sm/infer",
    )

    parser.add_argument(
        "-l",
        "--level",
        action="store",
        dest="level",
        default="all",
        type=str,
        help="Filtering on specified level of N-termini,\
                or labelled N-termini. Accepted values: all/nterm/quant.",
    )

    parser.add_argument(
        "-dn",
        "--dropna",
        action="store_true",
        dest="dropna",
        help="Flag to indicate whether to filter for empty quant rows",
    )

    parser.add_argument(
        "-fn",
        "--fillna",
        action="store",
        dest="fillna",
        help="Value to fill empty quant rows or cells with",
    )

    parser.add_argument(
        "-st",
        "--sleeptime",
        action="store",
        dest="sleeptime",
        default=0.2,
        type=float,
        help="Float that determine intervals between Biopython \
                            queries to Uniprot",
    )

    parser.add_argument(
        "-nx",
        "--noexo",
        action="store_true",
        dest="noexo",
        help="Do not check for dipeptidase or aminopeptidase activity",
    )

    parser.add_argument(
        "-nm",
        "--nomerops",
        action="store_true",
        dest="nomerops",
        help="Do not check for cleavage site annotation in MEROPS database",
    )

    parser.add_argument(
        "-cs",
        "--calcstructure",
        type=str,
        default=None,
        action="store",
        help="Annotate cleavage site solvent accessibility and secondary structure \
            for the significant peptides in the dataset using Pymol and Alphafold models \
                Accepted values: all/sig.",
    )

    parser.add_argument(
        "-sc",
        "--singlecpu",
        action="store_true",
        dest="singlecpu",
        help="Use a single process instead of threading for annotation",
    )

    parser.add_argument(
        "-cf",
        "--conditionfile",
        action="store",
        type=str,
        default=None,
        dest="conditionfile",
        help="Map labels to conditions. Adds columns for fold change for each pairwise \
        comparison, average and CV of each condition to the dataframe. Each line must \
        start with the condition name followed by the channels used, separated by a \
        single space",
    )

    parser.add_argument(
        "-stat",
        "--statistic",
        action="store_true",
        dest="stat",
        help="Performs statistical significance testing. Student T-test for two \
        conditions, ANOVA from three or more",
    )

    parser.add_argument(
        "-spw",
        "--stat_pairwise",
        action="store_true",
        dest="stat_pairwise",
        help="Performs statistical significance t-test for all conditions pairwise",
    )

    parser.add_argument(
        "-sig",
        "--significance",
        type=str,
        default=None,
        action="store",
        dest="significance",
        help="Performs fold change distribution significance check \
                    for all conditions pairwise \
                    Accepted values: all/nterm.",
    )

    parser.add_argument(
        "-vis",
        "--visualize",
        action="store_true",
        dest="visualize",
        help="Draws various plots based on conditions passed from input and statistical tests",
    )

    parser.add_argument(
        "-logo",
        "--logo",
        action="store",
        dest="logo",
        default=None,
        help="Draws various logo plots based on all peptides or condition significant peptides. Values supported: [all|prob|pssm|shannon|kbl]",
    )

    parser.add_argument(
        "-psc",
        "--pseudocounts",
        action="store",
        dest="pseudocounts",
        default=True,
        help="Add pseudocounts to normalized matrix calculation",
    )

    parser.add_argument(
        "-clvis",
        "--cleavagevisualization",
        type=str,
        default=None,
        action="store",
        dest="cleavagevis",
        help="Whether to visualize significant cleavages in \
                    sequence or both structure and sequence \
                    Accepted values: seq/both.",
    )

    parser.add_argument(
        "-enr",
        "--enrichment",
        action="store_true",
        dest="enrichment",
        help="Draws heatmap plots based on conditions passed from input and statistical tests \
                    for enrichment of GO terms and KEGG pathways with the gProfiler API",
    )

    parser.add_argument(
        "-o",
        "--output_name",
        action="store",
        dest="output_name",
        type=str,
        default=None,
        help="File name of output folder and annotated output file",
    )

    parser.add_argument(
        "-ot",
        "--outfile_type",
        action="store",
        dest="outfile_type",
        type=str,
        default="xlsx",
        help="File type of output \
                            file Accepted values: xlsx/csv/tsv/pkl/json",
    )

    parser.add_argument(
        "-sep",
        "--separate",
        action="store_true",
        dest="separate",
        help="Whether to merge or keep \
                            annotation as a separate file. False by default",
    )

    return vars(parser.parse_args())


def initialize(arguments=None):
    """Initializes the logger and the command line arguments."""

    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
    basefolder = Path.cwd().parent.absolute()
    logfile = basefolder / f"log/Annotator_{timestamp}.log"

    logger = initialize_logger(logfile)
    arguments = initialize_arguments() if arguments is None else arguments

    logger.info(f"Arguments: {arguments}")

    arguments["timestamp"] = timestamp
    arguments["logfile"] = logfile

    return arguments


def parse_arguments(parser):
    """Parses arguments, returns a dictionary with all key value argument
    pairs."""

    args = parser.parse_args()
    logging.info("Arguments parsed")
    arg_dict = vars(args)

    return arg_dict


def parse_sequence(seq: str):
    """Parses peptide sequence from annotated sequence input."""

    if not seq:
        return None

    match = re.search(r"\.([A-Z]+)\.", seq)
    if match:
        return match.group(1)
    else:
        return None


def parse_acc(acc_string: str):
    """Parses and selects the first accession for a list of uniprot ids."""

    if not acc_string:
        return None

    accs = acc_string.split(";")
    if len(accs) > 0:
        return accs[0]
    else:
        return None


def map_dict(annot_df_row: pd.core.series.Series, annot_dict: dict):
    """Maps values to annotation dataframe row."""

    for key in annot_dict:
        annot_df_row.loc[key] = annot_dict[key]

    return annot_df_row


def read_alphafold_accessions(accession_file: str):
    """Reads accessions from a file and returns a list of accessions."""

    with open(accession_file, "r") as f:
        accessions = f.read().splitlines()

    return accessions


def get_structure_properties(acc, position, env_length=4, available_models=None):
    """Computes the solvent accessible surface area of a peptide. Takes 
    the UniProt accession, the position of the cleavage in the protein, and
    a list of available models. Returns the solvent accessible surface area,
    around the cleavage site, or None if the model is not available."""

    if available_models is not None and acc not in available_models:
        logging.warning(f"Model for {acc} not available")
        return np.nan, np.nan
    
    alphafold_folder_name = r"\\ait-pdfs\services\BIO\Bio-Temp\Protease-Systems-Biology-temp\Kostas\CLIPPER\Datasets\Alphafold"
    alphafold_folder = Path(alphafold_folder_name)
    model_filename = f"AF-{acc}-F1-model_v4.cif.gz"
    model_path = alphafold_folder / model_filename
    
    # Load the model
    with tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as temp_cif:
        # Open the gzipped CIF file
        with gzip.open(model_path, 'rb') as f_in:
            # Write the decompressed contents to the temporary file
            temp_cif.write(f_in.read())
        
        temp_cif.flush()  # Ensure the file is written

        pymol.cmd.delete('all')
        pymol.cmd.load(temp_cif.name, acc)
        pymol.cmd.select('sel', f'resi {position - (env_length - 1)}-{position + env_length} and {acc}')

    # Compute the secondary structure
    pymol.cmd.dss(acc)
    ss_list = []
    pymol.cmd.iterate('sel and name ca', 'ss_list.append(ss)', space=locals())
    ss = ''.join(ss_list)

    # Compute the solvent accessible surface area
    pymol.cmd.set('dot_solvent', 1)
    sa= pymol.cmd.get_area('sel')

    return ss, sa


def save_figures(figures, folders):
    """Saves figures to output folder."""

    for k in figures:
        if figures[k] is not None:

            for i in figures[k]:
                if figures[k][i] is not None:

                    if "Fold" in k:
                        outfolder = folders["fold"]
                    elif "Piechart" in k:
                        outfolder = folders["piechart"]
                    elif "Logo" in k:
                        outfolder = folders["logo"]
                    elif "Volcano" in k:
                        outfolder = folders["volcano"]
                    elif "Enrichment" in k:
                        outfolder = folders["enrichment"]
                    else:
                        outfolder = folders["general"]
                    
                    try:
                        if k != "Clustermap" and not k.startswith("Logo") and not k.startswith("General") and not k.startswith("Piechart"):
                            figures[k][i].figure.savefig(outfolder / f"{k}_{i}.png", format="png", dpi=300)
                            figures[k][i].figure.savefig(outfolder / f"{k}_{i}.svg", format="svg")
                        else:
                            figures[k][i].savefig(outfolder / f"{k}_{i}.png", format="png", dpi=300)
                            figures[k][i].savefig(outfolder / f"{k}_{i}.svg", format="svg")
                    except Exception as e:
                        logging.info(f"Skipped {k, figures[k]} due to error: {str(e)}")
