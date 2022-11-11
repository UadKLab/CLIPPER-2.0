import logging
import os
import re
from argparse import ArgumentParser, HelpFormatter
from datetime import datetime


def initialize(arguments=None):
    """Initializes a parser instance and adds arguments for annotator."""

    timestamp = datetime.now()
    ftimestamp = timestamp.strftime(format="%d%m%y%H%M%S")
    logfile = os.path.join(
        os.path.dirname(os.getcwd()), f"log/Annotator_{ftimestamp}.log"
    )
    logging.basicConfig(filename=logfile, filemode="w", level=logging.INFO)

    formatted_timestamp = timestamp.strftime(format="%A %B %d %Y, %H:%M:%S")
    logging.info(f"Annotator started, {formatted_timestamp}")

    # add arguments if running from command line
    if arguments is None:

        parser = ArgumentParser(
            prog="TAILS annnotator",
            description="Peptide annotation and analysis of proteomics data \
                            utilizing Uniprot and MEROPS databases",
            epilog="Not extensively tested, this tool is still in \
                                        beta version. Contact konka@dtu.dk \
                                        for bug reports and requests.",
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

        logging.info("Arguments added")

        arguments = parse_arguments(parser)
        arguments["timestamp"] = ftimestamp
        arguments["logfile"] = logfile
        logging.info(f"Arguments returned with values:{arguments}")

        return arguments

    # else, just return passed arguments
    arguments["timestamp"] = ftimestamp
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

    match = re.search(r"\.([A-Z]+)\.", seq)
    if match:
        match_seq = match[1]

    return match_seq


def parse_acc(acc_string: str):
    """Parses and selects the first accession for a list of uniprot ids."""

    accs = acc_string.split(";")
    acc = accs[0]

    return acc


def map_dict(annot_df_row, annot_dict: dict):
    """Maps values to annotation dataframe row."""

    for key in annot_dict:
        annot_df_row.loc[key] = annot_dict[key]

    return annot_df_row
