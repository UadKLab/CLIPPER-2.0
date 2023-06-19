import logging
import os
import re
import string
import gzip
import tempfile

from argparse import ArgumentParser, HelpFormatter
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from tqdm import tqdm

import pymol
from reactome2py import analysis, content
import networkx as nx


import subprocess



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
        "-path",
        "--pathway",
        action="store_true",
        dest="pathway",
        help="Draws pathway plots based on conditions passed from input and statistical tests \
            and maps detected proteins and peptides",
    )

    parser.add_argument(
        "-pf",
        "--proteasefile",
        action="store",
        type=str,
        default=None,
        dest="proteasefile",
        help="Protease MEROPS identifiers to predict activity of. Using weighted PSSM",
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

    parser.add_argument(
        "-pv",
        "--pymol_verbose",
        action="store",
        dest="pymol_verbose",
        type=bool,
        default=False,
        help="Whether to output all pymol warnings/information to terminal during run or keep quiet.",
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
    basefolder = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
    logfile = os.path.join(basefolder, f"log/Annotator_{timestamp}.log")

    logger = initialize_logger(logfile)
    logger.debug(f"Arguments: {arguments}")

    write_terminal_headers('INITIALIZING CLIPPER 2.0')

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S,%f")[0:-3]} [INFO] The arguments you have provided are:')
    for argument, value in arguments.items():
        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S,%f")[0:-3]} [INFO] - {argument}: {value}')
    print("")

    arguments["timestamp"] = timestamp
    arguments["logfile"] = logfile

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

    logging.info(" - Reading available AlphaFold models...")
    with open(accession_file, "r") as f:
        accessions = f.read().splitlines()
    logging.info(f" - Read available AlphaFold models.")

    return accessions

def calculate_structure_properties(acc, cleavage_sites_indices, structure_properties, alphafold_folder, env_length = 4):

    # Load the model
    model_filename = f"AF-{acc}-F1-model_v4.cif.gz"
    model_path = os.path.join(alphafold_folder, model_filename)

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

    logging.info(" - Calculating surface area and secondary structure...")

    structure_properties = {}
    alphafold_folder_name = r"/Volumes/Bio-Temp/Protease-Systems-Biology-temp/Kostas/CLIPPER/Datasets/Alphafold"
    alphafold_folder = Path(alphafold_folder_name)

    # initialize the tmp file with an empty dict to not mess up the subprocess read/write system
    os.makedirs(os.path.dirname(tmp_output_path), exist_ok=True)
    with open(tmp_output_path, 'w') as f:
        f.write("{}")

    # if all pymol messages and output should be included in terminal set pymol_verbose = True
    if pymol_verbose:
        with tqdm(acc_cleavage_sites.items(), leave = 0) as t:
            for acc, cleavage_sites_indices in t:
                if available_models is not None and acc not in available_models:
                    logging.warning(f"Model for {acc} not available")
                structure_properties = calculate_structure_properties(acc, cleavage_sites_indices, structure_properties, alphafold_folder)
            with open(tmp_output_path, 'w') as f:
                f.write(str(structure_properties))
            elapsed = t.format_dict['elapsed']
            logging.info(f"Structure calculations took {format_seconds_to_time(elapsed)}")

    # otherwise run the function as a subprocess to choke the otherwise very persistent output
    else:
        with tqdm(acc_cleavage_sites.items(), leave = 0) as t:
            for acc, cleavage_sites_indices in t:
                if available_models is not None and acc not in available_models:
                    logging.warning(f"Model for {acc} not available")
                subprocess.run(['python', 'pymol_subprocess_ss.py', '-tf', str(tmp_output_path), '-acc', str(acc), '-af', str(alphafold_folder), '-csi', str(cleavage_sites_indices)], stdout=subprocess.DEVNULL)
            elapsed = t.format_dict['elapsed']
            logging.info(f"Structure calculations took {format_seconds_to_time(elapsed)}")
    

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

def get_enriched_pathways(accs, cutoff=0.05):

    """
    Get enriched pathways from Reactome using reactome2py.
    
    Parameters:
    accs (list): A list of protein accessions.
    cutoff (float, optional): The P-value cutoff for pathway enrichment.

    Returns:
    dict: A dictionary of enriched pathways and associated statistics.
    """

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