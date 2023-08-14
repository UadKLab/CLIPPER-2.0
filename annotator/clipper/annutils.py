import logging
import os
import sys
import re
import string
import gzip
import tempfile
import platform

from argparse import ArgumentParser, HelpFormatter
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from tqdm import tqdm

import pymol
from reactome2py import analysis, content
import networkx as nx

if platform.system() == 'Darwin':
    alphafold_folder_name = r"/Volumes/Bio-1/PSBT-g-Protease-Systems-Biology/Kostas/CLIPPER/Datasets/Alphafold"
else:
    alphafold_folder_name = r"W:\Protease-Systems-Biology-temp\Kostas\CLIPPER\Datasets\Alphafold"

def format_seconds_to_time(s):
    hours, remainder = divmod(s, 3600)
    minutes, seconds = divmod(remainder, 60)
    return '{:02}:{:02}:{:02} hh:mm:ss'.format(int(hours), int(minutes), int(seconds))

def write_terminal_headers(text, length = 91):
    print("\n")
    print("".center(length, '*'))
    print(f"  {text}  ".center(length, '*'))
    print("".center(length, '*'))
    print(" ")

def initialize_logger(logfile):
    
    """
    Initializes the logger with a file handler and a console handler.
    
    Parameters:
    logfile (str): The name of the file to log output messages to.
    
    Returns:
    logger (logging.RootLogger): A logger with two handlers - one for console output and one for file output.
    """

    logger = logging.getLogger()
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    file_handler = logging.FileHandler(logfile)
    file_handler.setFormatter(formatter)
    file_handler.setLevel(logging.DEBUG)
    logger.addHandler(file_handler)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    logger.setLevel(logging.DEBUG)

    logger.info('')

    logger.info(f"Annotator started at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.")
    print("\n\n\n*******************************************************************************************")
    print("*******************************************************************************************")
    print("                  _____ _      _____ _____  _____  ______ _____          ___         ___  ")
    print("   _       ,/'   / ____| |    |_   _|  __ \|  __ \|  ____|  __ \        |__ \       / _ \ ")
    print("  (_).  ,/'     | |    | |      | | | |__) | |__) | |__  | |__) |          ) |     | | | |")
    print("   _  ::        | |    | |      | | |  ___/|  ___/|  __| |  _  /          / /      | | | |")
    print("  (_)'  `\.     | |____| |____ _| |_| |    | |    | |____| | \ \         / /_   _  | |_| |")
    print("           `\.   \_____|______|_____|_|    |_|    |______|_|  \_\       |____| (_)  \___/ ")
    print("\n*******************************************************************************************")
    print("******************  Welcome to CLIPPER 2.0, a degradomics data annotator  *****************")
    print("*******************************************************************************************\n")       

    return logger

def initialize_arguments():

    """
    Initializes a parser instance and adds arguments for the annotator.
    
    Returns:
    args (dict): A dictionary containing all of the command line arguments provided by the user.
    """

    parser = ArgumentParser(
        prog="CLIPPER 2.0",
        description="Peptide annotation and analysis of proteomics data utilizing databases and visulization tools",
        epilog="Not extensively tested, this tool is still in beta version. Contact konka@dtu.dk or alemol@dtu.dk for bug reports and requests.",
        formatter_class=HelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--infile",
        action="store",
        dest="infile",
        type=str,
        required=True,
        help="Input peptide group result file.",
    )

    parser.add_argument(
        "-it",
        "--infiletype",
        default="infer",
        action="store",
        dest="infile_type",
        type=str,
        help="File type of input file. Accepted values: [excel|csv|infer].",
    )

    parser.add_argument(
        "-pa",
        "--preannotated",
        action="store_true",
        dest="preannotated",
        help="Flag, if given assume that file is preannotated, and will not annotate from uniprot, proteinatlas, or do exopeptidasecheck.",
    )

    parser.add_argument(
        "-cf",
        "--conditionfile",
        action="store",
        type=str,
        default=None,
        dest="conditionfile",
        help="A .txt file which maps labels to conditions. Adds columns for fold change for each pairwise \
        comparison, average and CV of each condition to the dataframe. Each line must \
        start with the condition name followed by the channels used, separated by a \
        single space (see example files for examples).",
    )

    parser.add_argument(
        "-a",
        "--alpha",
        action="store",
        dest="alpha",
        default=0.05,
        type=float,
        help="Float that sets the alpha value to be used for all relevant statistical thresholds and significance measurements. Default: 0.05.",
    )


    parser.add_argument(
        "-sw",
        "--software",
        default="infer",
        action="store",
        dest="software",
        type=str,
        help="Software the data were \
                            analyzed with. Accepted values: [pd|sm|infer].",
    )

    parser.add_argument(
        "-l",
        "--level",
        action="store",
        dest="level",
        default="all",
        type=str,
        help="Filtering on specified level of N-termini,\
                or labelled N-termini. Accepted values: [all|nterm|quant]. Default: all.",
    )

    parser.add_argument(
        "-dn",
        "--dropna",
        action="store_true",
        dest="dropna",
        help="Flag to indicate whether to filter for empty quant rows.",
    )

    parser.add_argument(
        "-fn",
        "--fillna",
        action="store",
        dest="fillna",
        help="Value to fill empty quant rows or cells with (must be a float). Must be an integer or float.",
    )

    parser.add_argument(
        "-st",
        "--sleeptime",
        action="store",
        dest="sleeptime",
        default=0.2,
        type=float,
        help="Float that determine interval in seconds between Biopython \
        queries to Uniprot. Default: 0.2.",
    )

    parser.add_argument(
        "-nx",
        "--noexo",
        action="store_true",
        dest="noexo",
        help="Flag, if given do not check for dipeptidase or aminopeptidase activity.",
    )

    parser.add_argument(
        "-nm",
        "--nomerops",
        action="store_true",
        dest="nomerops",
        help="Flag, if given do not check for cleavage site annotation in MEROPS database.",
    )

    parser.add_argument(
        "-cs",
        "--calcstructure",
        type=str,
        dest="calcstructure",
        default=None,
        action="store",
        help="Annotate cleavage site solvent accessibility and secondary structure \
        for the significant peptides in the dataset using Pymol and Alphafold models \
        Accepted values: [all|sig].",
    )

    parser.add_argument(
        "-c",
        "--cores",
        action="store",
        dest="threadingcores",
        default=os.cpu_count(),
        help="integer, Allows specifying the number of CPU cores used for gathering information from Uniprot. Accepted values are 'max' (using all cores) or an integer indicating the intended number of cores to be used. Default is 'max'. The speed of Uniprot information fetching scales linearly with the number of cores used.",
    )

    parser.add_argument(
        "-stat",
        "--statistic",
        action="store_true",
        dest="stat",
        help="Flag, if given perform statistical significance testing. Student T-test for two \
        conditions, ANOVA from three or more. In addition it provides multiple testing corrected p-values using the Benjamini/Hochberg method."
    )

    parser.add_argument(
        "-spw",
        "--stat_pairwise",
        action="store_true",
        dest="stat_pairwise",
        help="Flag, if given perform statistical significance t-test for all conditions pairwise.",
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
                    Accepted values: [all|nterm].",
    )

    parser.add_argument(
        "-mt",
        "--multipletesting",
        action="store_true",
        dest="multipletesting",
        help="Flag, if given multiple testing corrected values are used to determine significance for all relevant procedures using the method defined in -mtmethod.",
    )

    parser.add_argument(
        "-mtmethod",
        "--multipletestingmethod",
        type=str,
        default='fdr_bh',
        action="store",
        dest="multipletestingmethod",
        help="Sets the multiple testing method to be used. Accepted values: [None|bonferroni|sidak|holm-sidak|holm|simes-hochberg|hommel|fdr_bh|fdr_by|fdr_tsbh|fdr_tsbky]. For clarification see here: https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html. Default: fdr_bh.",
    )

    parser.add_argument(
        "-vis",
        "--visualize",
        action="store_true",
        dest="visualize",
        help="Flag, if given draws various plots based on conditions passed from input and statistical tests",
    )

    parser.add_argument(
        "-logo",
        "--logo",
        action="store",
        dest="logo",
        default=None,
        help="Draws various logo plots based on all peptides or condition significant peptides. Values supported: [all|prob|pssm|shannon|kbl].",
    )

    parser.add_argument(
        "-logo_fc",
        "--logofoldchange",
        type=float,
        action="store",
        dest="logo_fc",
        default=3,
        help="Float, sets the peptide abundance fold change between conditions which is considered for the logo generation. Example: -logo_fc 3 --> Filters peptides with log2 fold change > 3 (high logo) and < 3 (low logo).  Requires that -logo and -stat is given. Default: 3.",
    )

    parser.add_argument(
        "-cle",
        "--cleavageenvironment",
        type=int,
        action="store",
        dest="cleavagesitesize",
        default=4,
        help="Integer, size of the cleavage environment on both sides of the cleavage site. Default: 4.",
    )

    parser.add_argument(
        "-vc_fc",
        "--volcanofoldchange",
        nargs=2,
        type=float,
        action="store",
        dest="volcano_foldchange",
        default=1.5,
        help="Float, sets the cutoff value for fold change in volcano plots. Example: -vc_fc 1.5, which colors peptides with log2 fold change > 1.5 and < -1.5. P-value threshold is set by -a. Default: 1.5.",
    )

    parser.add_argument(
        "-psc",
        "--pseudocounts",
        action="store",
        dest="pseudocounts",
        default=True,
        help="Flag, if given add pseudocounts to normalized matrix calculation for logo generation.",
    )

    parser.add_argument(
        "-clvis",
        "--cleavagevisualization",
        type=str,
        default=None,
        action="store",
        dest="cleavagevis",
        help="Whether to visualize significant cleavages in \
                    sequence or both structure and sequence. For a complete dataset this can take a VERY long time, especially if set to 'both'. When using this argument please consider only adding proteins/peptides to the input file which you are really interested in (might be based on high fold change, pvalue, or prior knowledge). \
                    Accepted values: [None|seq|both]. Default: None.",
    )

    parser.add_argument(
        "-enr",
        "--enrichment",
        action="store_true",
        dest="enrichment",
        help="Flag, if given draws heatmap plots based on conditions passed from input and statistical tests \
                    for enrichment of GO terms and KEGG pathways with the gProfiler API.",
    )

    parser.add_argument(
        "-path",
        "--pathway",
        action="store_true",
        dest="pathway",
        help="Flag, if given draws pathway plots based on conditions passed from input and statistical tests \
            and maps detected proteins and peptides.",
    )

    parser.add_argument(
        "-pf",
        "--proteasefile",
        action="store",
        type=str,
        default=None,
        dest="proteasefile",
        help="A file with protease MEROPS identifiers to predict activity of. Using weighted PSSM. See examples for examples.",
    )

    parser.add_argument(
        "-o",
        "--output_name",
        action="store",
        dest="output_name",
        type=str,
        default=None,
        help="File name of output folder and annotated output file. Defaults to a timestamp for the initialization of Clipper 2.0.",
    )

    parser.add_argument(
        "-ot",
        "--output_filetype",
        action="store",
        dest="output_filetype",
        type=str,
        default="xlsx",
        help="File type of output \
                            file Accepted values: [xlsx|csv|tsv|pkl|json].",
    )

    parser.add_argument(
        "-sep",
        "--separate",
        action="store_true",
        dest="separate",
        help="Flag, if given the output file will not contain the original columns of the inputfile (separate original and annotated columns).",
    )

    parser.add_argument(
        "-pv",
        "--pymol_verbose",
        action="store_true",
        dest="pymol_verbose",
        help="Flag, if given the commandline will display info messages from pymol when used. Looks pretty ugly and unnecessary in this context, so muted by default",
    )

    return vars(parser.parse_args())


def initialize(arguments):

    """
    Initializes the logger and the command line arguments.
    
    Parameters:
    arguments (dict, optional): Command line arguments parsed into a dictionary. Defaults to None.
    
    Returns:
    arguments (dict): Initialized arguments with logger and timestamp information.
    """

    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
    basefolder = Path(__file__).resolve().parent.parent
    print(basefolder)
    logfile = basefolder / f"log/Annotator_{timestamp}.log"

    logger = initialize_logger(logfile)
    logger.debug(f"Arguments: {arguments}")

    write_terminal_headers('INITIALIZING CLIPPER 2.0')

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S,%f")[0:-3]} [INFO] The arguments you have provided are:')
    for argument, value in arguments.items():
        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S,%f")[0:-3]} [INFO] - {argument}: {value}')
    print("")

    arguments["timestamp"] = timestamp
    arguments["logfile"] = logfile

    args_ok = True
    warnings = []

    if arguments['stat'] and not arguments["conditionfile"]:
        warning = 'Condition statistics was requested (-stat argument was supplied) but no condition file is given. Stats cannot be done, please add a condition file for propper function.'
        warnings.append(warning)
        args_ok = False

    if arguments["conditionfile"] and not arguments["stat"]:
        warning = 'You have supplied a condition file, but are not doing condition specific statistics. Consider adding the -stat flag to your query!'
        warnings.append(warning)
        args_ok = False
    
    if arguments["calcstructure"] == "sig" and not arguments["stat"]:
        warning = 'Calculating structure properties on significant peptides is not possible without doing statistics. Please add the -stat flag to the query!'
        warnings.append(warning)
        args_ok = False

    if arguments["pathway"] and not arguments["stat"]:
        warning = 'Pathway analysis cannot be performed as the -stat flag is not given. Please rerun your query with the -stat flag for pathway plots'
        warnings.append(warning)
        args_ok = False

    if arguments["stat"] and arguments["pathway"] and not arguments["stat_pairwise"]:
        warning = 'If you have provided only 2 conditions ignore this message. If you have provided more than 2 conditions, please add the -spw flag (for pairwise statistics) to your query in order to do pathway analysis!'
        warnings.append(warning)
        args_ok = False

    if arguments["enrichment"] and not arguments["stat"]:
        warning = 'Functional enrichment analysis cannot be performed as the -stat flag is not given. Please rerun your query with the -stat flag for enrichment plots'
        warnings.append(warning)
        args_ok = False

    if arguments["logo_fc"] and arguments["logo"] and not arguments["stat"]:
        warning = '-Logo_fc only has an effect if both -logo and -stat are given'
        warnings.append(warning)
        args_ok = False

    if arguments["cleavagevis"]:
        warning = 'Cleavage visualization takes a VERY long time on full datasets (in our experience around 0.5-5 seconds per peptide depending on your computer), especially when set to "both". Please consider using this argument only with a filtered input list containing only the peptides you wish to have plotted on a protein structure.'
        warnings.append(warning)
        args_ok = False

    if arguments["threadingcores"]:
        max_cores= os.cpu_count()
        cores = arguments["threadingcores"]

        if cores == 'max':
            cores = max_cores
        else:
            try:
                cores = int(cores)
                if cores > max_cores:
                    warning = f"The number of cores provided to --cores was {cores}, however only {max_cores} cores were detected. Falling back to using {max_cores} cores."
                    warnings.append(warning)
                    args_ok = False    
                    cores = max_cores
            except Exception as err:
                warning = f"An input error was registered related to the argument --cores (set as {cores}). {max_cores} cores were detected. Likely the argument is not an integer. Python error: {err}"
                warnings.append(warning)
                warning = f"Falling back to using {max_cores} cores for further processing"
                warnings.append(warning)
                args_ok = False    
                cores = max_cores
        arguments["threadingcores"] = cores

    if not args_ok:
        print("\n\n")
        print("".center(91, '*'))
        print("  ARGUMENT WARNINGS  ".center(91, '*'))
        print("".center(91, '*'))
        print("\n\n")
        for warning in warnings:
            logging.warning(warning)
        print("\n\n")
        print("".center(91, '*'))
        print("\n\n")

    return arguments

def parse_arguments(parser):

    """
    Parses arguments, returns a dictionary with all key value argument pairs.
    
    Parameters:
    parser (ArgumentParser): ArgumentParser object with defined command-line arguments.
    
    Returns:
    arg_dict (dict): Dictionary containing all the command-line arguments parsed.
    """

    args = parser.parse_args()
    logging.info("Arguments parsed")
    arg_dict = vars(args)

    return arg_dict

def parse_sequence(seq: str):

    """
    Parses peptide sequence from annotated sequence input.
    
    Parameters:
    seq (str): Annotated sequence input string.
    
    Returns:
    str: Peptide sequence.
    """

    if not seq:
        return None

    match = re.search(r"\.([A-Z]+)\.", seq)
    if match:
        return match.group(1)
    else:
        return None

def parse_acc(acc_string: str):

    """
    Parses and selects the first accession for a list of uniprot ids.
    
    Parameters:
    acc_string (str): String with list of uniprot ids.
    
    Returns:
    str: The first uniprot id in the list.
    """

    if not acc_string:
        return None

    accs = acc_string.split(";")
    if len(accs) > 0:
        return accs[0]
    else:
        return None

def map_dict(annot_df_row: pd.core.series.Series, annot_dict: dict):

    """
    Maps values to annotation dataframe row.
    
    Parameters:
    annot_df_row (pd.core.series.Series): DataFrame row where the values will be mapped.
    annot_dict (dict): Dictionary with values to be mapped.
    
    Returns:
    pd.core.series.Series: DataFrame row with mapped values.
    """

    for key in annot_dict:
        annot_df_row.loc[key] = annot_dict[key]

    return annot_df_row

def read_alphafold_accessions(accession_file: str):

    """
    Reads accessions from a file and returns a list of accessions.
    
    Parameters:
    accession_file (str): Name of the file containing accessions.
    
    Returns:
    list: List of accessions.
    """

    logging.info("Reading available AlphaFold models...")
    with open(accession_file, "r") as f:
        accessions = f.read().splitlines()
    logging.info(f"Read available AlphaFold models.")

    return accessions

def calculate_structure_properties(acc, cleavage_sites_indices, structure_properties, alphafold_folder_name, env_length = 4):

    # Load the model
    model_filename = f"AF-{acc}-F1-model_v4.cif.gz"
    alphafold_folder = Path(alphafold_folder_name)
    model_path = alphafold_folder / model_filename

    with tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as temp_cif:
        # Open the gzipped CIF file
        with gzip.open(model_path, 'rb') as f_in:
            # Write the decompressed contents to the temporary file
            temp_cif.write(f_in.read())

        temp_cif.flush()  # Ensure the file is written

        pymol.cmd.delete('all')
        pymol.cmd.load(temp_cif.name, acc)

        for index, cleavage_site in cleavage_sites_indices:
            pymol.cmd.select('sel', f'resi {cleavage_site - (env_length - 1)}-{cleavage_site + env_length} and {acc}')

            # Compute the secondary structure
            pymol.cmd.dss(acc)
            ss_list = []
            pymol.cmd.iterate('sel and name ca', 'ss_list.append(ss)', space=locals())
            ss = ''.join(ss_list)

            # Compute the solvent accessible surface area
            pymol.cmd.set('dot_solvent', 1)
            sa= pymol.cmd.get_area('sel')

            structure_properties[(acc, cleavage_site)] = (index, ss, sa)

    return structure_properties

@contextmanager
def stdout_redirected(to=os.devnull):
    '''
    Used to redirect stdout while running pymol commands. This removes the info messages from pymol in the terminal which obstructs tqdm bars.
    '''
    fd = sys.stdout.fileno()

    def _redirect_stdout(to):
        sys.stdout.close() # + implicit flush()
        os.dup2(to.fileno(), fd) # fd writes to 'to' file
        sys.stdout = os.fdopen(fd, 'w') # Python writes to fd

    with os.fdopen(os.dup(fd), 'w') as old_stdout:
        with open(to, 'w') as file:
            _redirect_stdout(to=file)
        try:
            yield # allow code to be run with the redirected stdout
        finally:
            _redirect_stdout(to=old_stdout) # restore stdout.
                                            # buffering and flags such as
                                            # CLOEXEC may be different

def get_structure_properties(acc_cleavage_sites, tmp_output_path, pymol_verbose, available_models=None):

    """
    Computes the solvent accessible surface area of a peptide.
    
    Parameters:
    acc_cleavage_sites (dict): A dictionary of UniProt accession and cleavage sites.
    env_length (int, optional): Length of the environment around the cleavage site.
    available_models (list, optional): A list of available models. 

    Returns:
    dict: A dictionary with the solvent accessible surface area around the cleavage site,
    or None if the model is not available.
    """

    logging.info("Calculating surface area and secondary structure...")

    structure_properties = {}

    # initialize the tmp file with an empty dict to not mess up the subprocess read/write system
    os.makedirs(os.path.dirname(tmp_output_path), exist_ok=True)
    with open(tmp_output_path, 'w') as f:
        f.write("{}")

    # if all pymol messages and output should be included in terminal set pymol_verbose = True
    warnings = []
    with tqdm(acc_cleavage_sites.items(), leave = 0) as t:
        for acc, cleavage_sites_indices in t:
            if available_models and acc in available_models:
                if pymol_verbose:
                    structure_properties = calculate_structure_properties(acc, cleavage_sites_indices, structure_properties, alphafold_folder_name)
                else:
                    with stdout_redirected():
                        structure_properties = calculate_structure_properties(acc, cleavage_sites_indices, structure_properties, alphafold_folder_name)
            else:
                warnings.append(f"Model for {acc} not available")
        with open(tmp_output_path, 'w') as f:
            f.write(str(structure_properties))
        elapsed = t.format_dict['elapsed']
        logging.info(f"Structure calculations took {format_seconds_to_time(elapsed)}")


    if warnings != []:
        logging.warning("Structure calculations were not possible for selected peptides, as alphafold models were not available")
        for warning in warnings:
            logging.warning(f'  - {warning}')


    

    return structure_properties


def map_accessions(accessions):

    """
    Map protein accessions to human gene names using the Reactome API.
    
    Parameters:
    accessions (list): A list of protein accessions.

    Returns:
    list: A list of human gene names corresponding to the input accessions.
    """

    accessions_str = ",".join(accessions)
    mapping = analysis.identifiers_mapping(ids=accessions_str, interactors=False, projection=True)

    # Extract gene names from the mapping result, ignoring entries without a mapping
    human_accessions = [entry['mapsTo'][0]['identifier'] for entry in mapping if entry['mapsTo']]

    return human_accessions
"""
def get_enriched_pathways(accs, cutoff=0.05):

#    ""
#    Get enriched pathways from Reactome using reactome2py.
#    
#    Parameters:
#    accs (list): A list of protein accessions.
#    cutoff (float, optional): The P-value cutoff for pathway enrichment.
#
#    Returns:
#    dict: A dictionary of enriched pathways and associated statistics.
#    ""

    # Use reactome2py to perform the analysis
    query = ",".join(accs)
    res = analysis.identifiers(ids=query)
    
    # Extract the token from the results
    token = res.get('summary', {}).get('token', None)
    if token is None:
        raise ValueError("Could not retrieve analysis token from Reactome")
    
    # Use the token to get the detailed results
    pathways = analysis.token(token, page_size='20', page='1', sort_by='ENTITIES_FDR', 
                              order='ASC', resource='TOTAL', p_value=cutoff, include_disease=True, 
                              min_entities=None, max_entities=None)
    
    return pathways
"""

def get_proteins_and_interactors(pathway_id):

    """
    Get the proteins and interactors for a pathway.
    
    Parameters:
    pathway_id (str): A Reactome pathway ID.

    Returns:
    tuple: A tuple with two lists, one of proteins and one of interaction maps.
    """
    
    path_res = content.participants_reference_entities(pathway_id)
    proteins = [p['identifier'] for p in path_res if p['identifier'][0] in list(string.ascii_uppercase) and len(p['identifier']) == 6]
    
    interactors_query = ','.join(proteins)
    inter_res = content.interactors_static_accs(accs=interactors_query)

    interaction_map = {}
    if inter_res is not None:
        for entry in inter_res['entities']:
            acc = entry['acc']
            if entry['count'] > 0:
                interactors = {inter['acc'] for inter in entry['interactors']}
                interactors_filtered = {inter for inter in interactors if inter in proteins}
            else:
                interactors_filtered = set()
            
            interaction_map[acc] = interactors_filtered

    return proteins, interaction_map

def construct_edgelists(subframe, interaction_map):

    """
    Constructs the edgelists for the protein and cleavage networks.
    
    Parameters:
    subframe (pd.DataFrame): A pandas DataFrame with sub-frame data.
    interaction_map (dict): An interaction map.

    Returns:
    tuple: A tuple with three elements - protein edgelist, cleavage edgelist, and cleavages.
    """

    protein_edgelist = []
    for source, targets in interaction_map.items():
        if len(targets) > 0:
            for target in targets:
                protein_edgelist.append((source, target))

    cleavage_edgelist = []
    cleavages = []
    for index in subframe.index:
        protein = subframe.loc[index, 'query_accession']

        # Only include cleavages from proteins in the pathway
        if protein in interaction_map:
            peptide = '-'.join([str(subframe.loc[index, 'start_pep']), str(subframe.loc[index, 'end_pep'])])
            position = ':'.join([protein, peptide])
            cleavage = ';'.join([position, subframe.loc[index, 'query_sequence']])
            cleavages.append(cleavage)
            cleavage_edgelist.append((protein, cleavage))
    
    return protein_edgelist, cleavage_edgelist, cleavages

def construct_network(protein_edgelist, cleavage_edgelist, proteins):

    """
    Constructs a network from the protein and cleavage edgelists.
    
    Parameters:
    protein_edgelist (list): A list of tuples representing edges between proteins.
    cleavage_edgelist (list): A list of tuples representing edges from proteins to cleavages.
    proteins (list): A list of protein identifiers.

    Returns:
    nx.DiGraph: A NetworkX directed graph object representing the network.
    """

    network = nx.DiGraph()
    # protein edges
    network.add_edges_from(protein_edgelist)
    # cleavage edges
    network.add_edges_from(cleavage_edgelist)
    # add proteins as nodes
    network.add_nodes_from(proteins)

    return network

def save_figures(figures, folders):
    
    """
    Saves figures to output folder.
    
    Parameters:
    figures (dict): A dictionary containing matplotlib figure objects.
    folders (dict): A dictionary containing output folder paths.

    Returns:
    None
    """

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