import concurrent.futures
import logging
import os
import shutil
import sys
import warnings
from itertools import combinations, permutations

import annutils
import numpy as np
import pandas as pd
from entry import Entry
from logo import generate_logos
from scipy.stats import f_oneway, rv_histogram
from statsmodels.stats.weightstats import ttest_ind
from tqdm import tqdm
from visualize import Visualizer

# global variables
result_folder_name = "results"
data_folder_name = "data"
annotation_prefix = "_annot."


class Annotator:
    """Annotator class for processing and analyzing peptide proteomics data. All
    arguments are passed as a dictionary to the class constructor.

    Attributes:
        conditions (dict): A dictionary of conditions and their respective files.
        annot (pd.DataFrame): The input data as a Pandas DataFrame.
        outfolder (str): The path to the output folder.
        outfile_type (str): The type of output file to generate.
        logfolder (str): The path to the log folder.
        logfile (str): The path to the log file.
        figures (dict): A dictionary of generated figures.
        separate (bool): A flag to indicate if annotation data should be separate.
        pseudocounts (float): The pseudocount value for sequence logo generation.
        logo (str): The type of logo to generate.
        stat (bool): A flag to indicate if statistical calculations should be performed.
        df (pd.DataFrame): The dataframe to store processed data.
    """

    def __init__(self, args):
        """Initialize the Annotator class, and set the attributes from the arguments."""

        # input attributes
        self.infile_type = args["infile_type"]
        self.infile = args["infile"]
        self.software = args["software"]

        # filtering and sanitizing
        self.level = args["level"]
        self.dropna = args["dropna"]

        # annotation attributes
        self.separate = args["separate"]
        self.sleeptime = args["sleeptime"]
        self.noexo = args["noexo"]
        self.nomerops = args["nomerops"]

        self.conditionfile = args["conditionfile"]
        self.stat = args["stat"]
        self.significance = args["significance"]
        
        self.logo = args["logo"]
        self.pseudocounts = args["pseudocounts"]

        # logging and file handling
        self.timestamp = args["timestamp"]
        self.logfile = args["logfile"]
        self.outfile_type = args["outfile_type"]
        self.outname = args["output_name"]

        self.conditions = None
        self.figures = {}

        self.basefolder = os.path.dirname(os.getcwd())
        self.resultfolder = os.path.join(self.basefolder, result_folder_name)
        self.datafolder = os.path.join(self.basefolder, data_folder_name)

        if self.outname:
            self.outfolder = os.path.join(self.resultfolder, self.outname)
            self.outname = self.outname + annotation_prefix + self.outfile_type
        else:
            self.outfolder = os.path.join(self.resultfolder, self.timestamp)
            self.outname = self.infile.rsplit("/", 1)[-1].rsplit("\\", 1)[-1].rsplit(".", 1)[0] + annotation_prefix + self.outfile_type

        print("\nAnnotator initialized\n")
        logging.info("Initialization successful")

    def prepare(self):
        """Control of arguments from user input."""

        if self.infile_type == "infer":
            if self.infile.endswith(".csv"):
                self.infile_type = "csv"
            if self.infile.endswith(".xlsx") or self.infile.endswith(".xls"):
                self.infile_type = "excel"

        if self.outfile_type not in ["xlsx", "csv", "tsv", "pkl", "json"]:
            logging.critical(f"Invalid output extension: {self.outfile_type}. Exiting with code 1.")
            raise TypeError("Please select a valid output format and try again.")

        if self.infile_type not in ["excel", "csv"]:
            logging.critical(f"Invalid input extension: {self.infile_type}. Exiting with code 1")
            raise TypeError(
                "Please select a valid input format and try again. Accepted \
                values: excel/csv, excel by default"
            )

        # check if condition file exists before doing any work
        if self.conditionfile is not None:
            if not os.path.exists(self.conditionfile):
                logging.critical("Exit with code 6")
                raise TypeError(
                    "Invalid condition file. Check if the path exists, and try again."
                )
                logging.critical("Exit with code 6")

        print("Reading file...\n")
        self.read_file()
        logging.info("Read dataframe")
        logging.info(f"Read input with {len(self.df)} peptides")
        print(f"Read input with {len(self.df)} peptides")

        if self.software == "infer":
            try:
                self.df.loc[0, "Master Protein Accessions"]
                self.software = "pd"
                logging.info("Input software is pd")
            except KeyError:
                try:
                    self.df.loc[0, "PG.ProteinAccessions"]
                    self.software = "sm"
                    logging.info("Input software is sm")
                except KeyError:
                    logging.critical("Invalid input. Exiting with code 4.")
                    raise TypeError(
                        f"Invalid input. Please make sure input format is correct \
                        and contains accession and sequence columns with default \
                        names, and try again."
                    )
                except:
                    logging.critical("Invalid input. Exiting with code 5.")
                    raise TypeError("Invalid input")

        #get column patterns to be used for indexing
        self.patterns = self.get_patterns()
        logging.info("Successfully generated indexing patterns")
        logging.info("Format check complete")

        if self.level != "all":
            self.filter_df()
            logging.info(f"Filtered dataframe, {len(self.df)} peptides remaining")
            print(f"Filtered dataframe, {len(self.df)} peptides remaining\n")

        logging.info("Sanitizing")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.sanitize()
        logging.info(f"Sanitized dataframe, {len(self.df)} peptides remaining")
        print(f"Sanitized dataframe, {len(self.df)} peptides remaining\n")

        # create folder where results will be stored after all sanity checks
        os.mkdir(self.outfolder)

    def sanitize(self):
        """Checks for empty columns in sequence and accession rows of input
        dataframe.

        Also, checks for any non-standard amino acids in the sequences
        and removes any such rows.
        """

        col_acc = self.patterns['acc']
        col_seq = self.patterns['seq']
        pattern = self.patterns['amino']

        invalid_seq = self.df[self.df[col_seq].isna()]
        if len(invalid_seq) > 0:
            for ind in invalid_seq.index:
                logging.info(f"Empty sequence row {ind+1}")
            self.df = self.df.dropna(subset=[col_seq]).reset_index(drop=True)

        invalid_acc = self.df[self.df[col_acc].isna()]
        if len(invalid_acc) > 0:
            for ind in invalid_acc.index:
                logging.info(f"Empty accession row {ind+1}")
            self.df = self.df.dropna(subset=[col_acc]).reset_index(drop=True)

        invalid_alphabet = self.df[self.df[col_seq].str.contains(pattern)]
        if len(invalid_alphabet) > 0:
            for ind in invalid_alphabet.index:
                logging.info(
                    f"Invalid sequence character row {ind+1}, {invalid_alphabet.loc[ind, col_seq]}"
                )
            self.df = self.df.drop(index=invalid_alphabet.index).reset_index(drop=True)

    def read_MEROPS(self):
        """Reads the two files containing MEROPS data for protease
        annotation."""

        logging.info("Reading MEROPS data..")

        self.merops = pd.read_csv(os.path.join(self.datafolder, "cleavage.csv"))
        self.merops_name = pd.read_csv(os.path.join(self.datafolder, "protein_name.csv"))
        self.merops_name = self.merops_name[self.merops_name.type == "real"]

    def initialize_annotation(self, length: int):
        """Initialize a dataframe with empty annotation columns, same size as
        input df."""

        self.annot = pd.DataFrame(
            columns=[
                "query_sequence",
                "query_accession",
                "name",
                "full_sequence",
                "description",
                "keywords",
                "go_codes",
                "go_names",
                "proteoform_certainty%",
                "start_pep",
                "end_pep",
                "p1_position",
                "cleavage_site",
                "p4_p4prime",
                "nterm_annot",
                "protease_uniprot",
                "protease_merops_code",
                "protease_merops_name",
            ],
            index=range(length),
        )

        if self.nomerops:
            self.annot.drop(
                ["protease_merops_code", "protease_merops_name"], axis=1, inplace=True
            )

    def read_file(self):
        """Reads input file, returns dataframe object."""

        if self.infile_type == "csv":
            logging.info("Input is csv")

            try:
                df = pd.read_csv(self.infile, sep=";")
            except OSError as err:
                logging.critical("Could not read file, OSError")
                logging.critical(err)
                logging.critical("Exit with code 2")
                print(f"ERROR: {err}")
                raise TypeError(
                    f"Could not read file. make sure path {self.infile} is correct, \
                    and try again."
                )
            except Exception as err:
                logging.critical("Could not read file, Error")
                logging.critical(err)
                logging.critical("Exit with code 2")
                print(f"ERROR: {err}")
                raise TypeError(
                    f"Could not read file. make sure file type {self.infile_type} \
                    is correct, and try again."
                )    

        elif self.infile_type == "excel":
            logging.info("Input is excel")

            try:
                df = pd.read_excel(self.infile, engine="openpyxl")
            except OSError as err:
                logging.critical("Could not read file, OSError")
                logging.critical(err)
                logging.critical("Exit with code 2")
                print(f"ERROR: {err}")
                raise TypeError(
                    f"Could not read file. Make sure path {self.infile} is correct, \
                    and try again."
                )
            except Exception as err:
                logging.critical("Could not read file, Error")
                logging.critical(err)
                logging.critical("Exit with code 2")
                print(f"ERROR: {err}")
                raise TypeError(
                    f"Could not read file. Make sure file type {self.infile_type} is \
                    correct, and try again."
                )
                
        self.df = df

    def get_patterns(self):
        """Returns column patterns to be used for indexing."""

        patterns = {}

        if self.software == 'sm':
            try:
                annutils.parse_acc(self.df.loc[0, "PG.ProteinAccessions"])
                self.df.loc[0, "PEP.StrippedSequence"]
                
                patterns['acc'] = "PG.ProteinAccessions"
                patterns['seq'] = "PEP.StrippedSequence"
                patterns['mod'] = "P.MoleculeID"
                patterns['amino'] = "B|J|O|U|X|Z"
                patterns['label'] = r"\[TMT"
                patterns['nterm'] = r"N-?ter"
                patterns['nterm_label'] = r"TMT.*_Nter"
                patterns['lysine_label'] = r"K\[TMT.{0,3}_Lys\]"
            except:
                logging.critical("Invalid input. Exiting with code 4.")
                raise TypeError(
                    f"Invalid input. Please make sure input format {self.software} \
                    is correct and contains sequence columns with default names, \
                    and try again."
                )
            
            pat = r"PEP\.TMT"

            quant = self.df.columns[self.df.columns.str.contains(pat=pat)]
            if len(quant) > 0:
                patterns['quant'] = pat
                try:
                    columns = self.df.columns[self.df.columns.str.contains(patterns['quant'])]
                    for col in columns:
                        self.df[col] = (
                            self.df[col].astype(str).str.replace(",", ".").astype(float)
                        )
                except:
                    logging.critical("Invalid input. Exiting with code 4.")
                    raise TypeError(
                        f"Invalid input. Could not convert values to float. Make sure \
                        input format {self.software} is correct and there are no string \
                        literals in quant columns, and try again."
                    )
            else:
                patterns['quant'] = None

        elif self.software == 'pd':
            patterns['acc'] = "Master Protein Accessions"
            patterns['mod'] = "Modifications"
            patterns['nterm'] = r"\[N-Term\]"

            try:
                annutils.parse_acc(self.df.loc[0, "Master Protein Accessions"])
                self.df.loc[0, "Sequence"]
                patterns['seq'] = "Sequence"
                patterns['amino'] = "B|J|O|U|X|Z"
            except KeyError:
                try:
                    annutils.parse_acc(self.df.loc[0, "Master Protein Accessions"])
                    annutils.parse_sequence(self.df.loc[0, "Annotated Sequence"])
                    patterns['seq'] = "Annotated Sequence"
                    patterns['amino'] = "\.[A-Z]*(B|J|O|U|X|Z)[A-Z]*\."
                except KeyError:
                    logging.critical("Invalid input. Exiting with code 4.")
                    raise TypeError(
                        f"Invalid input. Please make sure input \
                        format {self.software} is correct and contains \
                        sequence columns with default names, and try again."
                    )

            pat_scale = r'Abundances \(?Scaled\)?.*[0-9]{3}'
            pat_norm = r'Abundances \(?Normalized\)?.*[0-9]{3}'
            pat_grouped = r'Abundances \(Grouped\)?.*[0-9]{3}'
            pat_raw = r'Abundance .*[0-9]{3}'
            pat_other = r'Abundance.*F[0-9]+:'
            
            quant_scale = self.df.columns[self.df.columns.str.contains(pat=pat_scale)]
            quant_norm = self.df.columns[self.df.columns.str.contains(pat=pat_norm)]
            quant_grouped = self.df.columns[self.df.columns.str.contains(pat=pat_grouped)]
            quant_raw = self.df.columns[self.df.columns.str.contains(pat=pat_raw)]
            quant_other = self.df.columns[self.df.columns.str.contains(pat=pat_other)]

            if len(quant_scale) > 0:
                patterns['quant'] = pat_scale
            elif len(quant_norm) > 0:
                patterns['quant'] = pat_norm
            elif len(quant_grouped) > 0:
                patterns['quant'] = pat_grouped
            elif len(quant_raw) > 0:
                patterns['quant'] = pat_raw
            elif len(quant_other) > 0:
                patterns['quant'] = pat_other
            else:
                patterns['quant'] = None

            pat_label = r"TMT"
            pat_alt_label = r"Dimethyl"

            quant_tmt = len(self.df[self.df['Modifications'].str.contains(pat=pat_label, na=False)])
            quant_dimethyl = len(self.df[self.df['Modifications'].str.contains(pat=pat_alt_label, na=False)])

            if quant_tmt == 0:
                if quant_dimethyl == 0:
                    patterns['label'] = None
                    patterns['nterm_label'] = None
                else:
                    patterns['label'] = pat_alt_label
                    patterns['nterm_label'] = r"Dimethyl \[N-Term\]"
                    patterns['lysine_label'] = r"Dimethyl \[K"

            else:
                patterns['label'] = pat_label
                patterns['nterm_label'] = r"TMT.* \[N-Term\]"
                patterns['lysine_label'] = r"TMT.{0,5} \[K"

        return patterns

    def proteoform_check(self):
        """Computes probability of proteoforms based on master accession
        column"""
        
        for row in self.df.index:
            accs = self.df.loc[row, self.patterns['acc']]
            proteoforms = [i.strip() for i in accs.split(';')]

            self.annot.loc[row, "proteoform_certainty%"] = 1/len(proteoforms) * 100

    def read_condition_file(self):
        """If a condition file was supplied, this function parses and stores
        information about conditions and corresponding channels in class
        variables."""

        with open(self.conditionfile, "r") as fh:
            lines = fh.readlines()

        conditions = {}
        for line in lines:
            line_list = line.strip().split()
            cond, channels = line_list[0], line_list[1:]
            conditions[cond] = channels

        self.conditions = conditions

    def general_conditions(self):
        """General statistics for the conditions supplied.

        Fold change, mean, std and CV. Fold change is reported for each
        pairwise comparison of input conditions
        """

        if self.conditions is None:
            self.conditions = {
                "all": ["126", "127", "128", "129", "130", "131", "132", "133", "134"]
            }

        for condition in self.conditions:

            column_mean = f"{condition}_mean"
            column_std = f"{condition}_deviation"
            column_cv = f"{condition}_CV"

            cols = self.df.columns[
                self.df.columns.str.contains("|".join(self.conditions[condition]))
            ]

            # remove any columns that have been accidentally picked up by the pattern
            cols = cols[cols.str.contains(self.patterns['quant'])]

            # initialize columns and fill only if quant columns were found
            if len(cols) > 0:
                for col in (column_mean, column_std, column_cv):
                    self.annot[col] = np.nan

                for loc in self.df.index:
                    mean = self.df.loc[loc][cols].astype(float).mean()
                    std = self.df.loc[loc][cols].astype(float).std()
                    self.annot.loc[loc, column_mean] = mean
                    self.annot.loc[loc, column_std] = std
                    self.annot.loc[loc, column_cv] = std / mean

        if len(self.conditions) > 1:
            conditions_iter = permutations(self.conditions.keys(), 2)
            for pair in tqdm(conditions_iter):

                column_name = f"Fold_change: {pair[0]}/{pair[1]}"
                column_log = f"Log2_fold_change: {pair[0]}/{pair[1]}"
                self.annot[column_name] = np.nan
                self.annot[column_log] = np.nan
                cols0 = self.df.columns[
                    self.df.columns.str.contains("|".join(self.conditions[pair[0]]))
                ]
                cols1 = self.df.columns[
                    self.df.columns.str.contains("|".join(self.conditions[pair[1]]))
                ]

                # remove any columns that have been accidentally picked up by the pattern
                cols0 = cols0[cols0.str.contains(self.patterns['quant'])]
                cols1 = cols1[cols1.str.contains(self.patterns['quant'])]

                for loc in self.df.index:
                    mean0 = self.df.loc[loc][cols0].mean()
                    mean1 = self.df.loc[loc][cols1].mean()

                    try:
                        self.annot.loc[loc, column_name] = mean0 / mean1
                        self.annot.loc[loc, column_log] = np.log2(mean0 / mean1)
                    except ZeroDivisionError:
                        continue
                    except:
                        logging.info(
                            f"Condition invalid \
                            assignment, {loc, pair[0], pair[1], mean0, mean1}"
                        )

    def percentile_fold(self, percentile):
        """Checks fold change distribution and marks rows above a certain
        percentile."""

        left_cutoff = percentile
        right_cutoff = 1 - percentile

        if len(self.conditions) > 1:
            conditions_iter = permutations(self.conditions.keys(), 2)

            if self.significance == 'all':

                for pair in conditions_iter:
                    self.annot[f"Fold {pair[0]}/{pair[1]} significance"] = np.nan
                    column = f"Fold_change: {pair[0]}/{pair[1]}"

                    hist = np.histogram(self.annot[column].dropna())
                    hist_dist = rv_histogram(hist)

                    for loc in self.annot.index:
                        cd = hist_dist.cdf(self.annot.loc[loc][column])
                        if cd > right_cutoff:
                            self.annot.loc[
                                loc, f"Fold {pair[0]}/{pair[1]} significance"
                            ] = "significant high"
                        elif cd < left_cutoff:
                            self.annot.loc[
                                loc, f"Fold {pair[0]}/{pair[1]} significance"
                            ] = "significant low"

            elif self.significance == 'nterm':

                subframe_internal = self.annot[self.annot["nterm_annot"] == "Internal"]
                subframe_natural = self.annot[self.annot["nterm_annot"] != "Internal"]
                
                for pair in conditions_iter:
                    self.annot[f"Fold {pair[0]}/{pair[1]} significance"] = np.nan
                    column = f"Fold_change: {pair[0]}/{pair[1]}"

                    hist = np.histogram(subframe_internal[column].dropna())
                    hist_dist = rv_histogram(hist)

                    for loc in subframe_internal.index:
                        cd = hist_dist.cdf(self.annot.loc[loc][column])
                        if cd > right_cutoff:
                            self.annot.loc[
                                loc, f"Fold {pair[0]}/{pair[1]} significance"
                            ] = "significant high"

                    hist = np.histogram(subframe_natural[column].dropna())
                    hist_dist = rv_histogram(hist)

                    for loc in subframe_natural.index:
                        cd = hist_dist.cdf(self.annot.loc[loc][column])
                        if cd < left_cutoff:
                            self.annot.loc[
                                loc, f"Fold {pair[0]}/{pair[1]} significance"
                            ] = "significant low"
                
            else:
                print(f"Significance invalid argument {self.significance}, skipping")
                logging.info(
                            f"Significance invalid argument {self.significance}, skipping"
                        )

    def condition_statistics(self, pairwise):
        """Perform a ttest or ANOVA statistical significance tests."""

        if pairwise:
            conditions_iter = combinations(self.conditions.keys(), 2)

            for pair in tqdm(conditions_iter):
                column_name = f"Ttest: {pair[0]}_{pair[1]}"
                column_log = f"Log10_ttest: {pair[0]}_{pair[1]}"
                self.annot[column_name] = np.nan
                self.annot[column_log] = np.nan

                cols0 = self.df.columns[
                    self.df.columns.str.contains("|".join(self.conditions[pair[0]]))
                ]
                cols1 = self.df.columns[
                    self.df.columns.str.contains("|".join(self.conditions[pair[1]]))
                ]

                # remove any columns that have been accidentally picked up by the pattern
                cols0 = cols0[cols0.str.contains(self.patterns['quant'])]
                cols1 = cols1[cols1.str.contains(self.patterns['quant'])]

                for loc in self.df.index:
                    values0 = list(self.df.loc[loc][cols0].values)
                    values1 = list(self.df.loc[loc][cols1].values)

                    result = ttest_ind(values0, values1)
                    self.annot.loc[loc, column_name] = result[1]
                    self.annot.loc[loc, column_log] = np.log10(result[1])

        else:
            if len(self.conditions) == 2:
                conditions = list(self.conditions.keys())
                column_name = f"Ttest: {conditions[0]}_{conditions[1]}"
                column_log = f"Log10_ttest: {conditions[0]}_{conditions[1]}"
                self.annot[column_name] = np.nan
                self.annot[column_log] = np.nan

                for loc in tqdm(self.df.index):
                    values = []

                    for condition in self.conditions:
                        cols = self.df.columns[
                            self.df.columns.str.contains(
                                "|".join(self.conditions[condition])
                            )
                        ]

                        # remove any columns that have been accidentally picked up by the pattern
                        cols = cols[cols.str.contains(self.patterns['quant'])]

                        col_values = list(self.df.loc[loc][cols].values)
                        values.append(col_values)

                    result = ttest_ind(*values)
                    self.annot.loc[loc, column_name] = result[1]
                    self.annot.loc[loc, column_log] = np.log10(result[1])

            elif len(self.conditions) > 2:
                
                self.annot["ANOVA"] = np.nan
                self.annot["Log10_ANOVA"] = np.nan

                for loc in tqdm(self.df.index):
                    values = []

                    for condition in self.conditions:
                        cols = self.df.columns[
                            self.df.columns.str.contains(
                                "|".join(self.conditions[condition])
                            )
                        ]

                        # remove any columns that have been accidentally picked up by the pattern
                        cols = cols[cols.str.contains(self.patterns['quant'])]

                        col_values = list(self.df.loc[loc][cols].values)
                        values.append(col_values)

                    result = f_oneway(*values)
                    self.annot.loc[loc, "ANOVA_pval"] = result.pvalue
                    self.annot.loc[loc, "Log10_ANOVA_pval"] = np.log10(result.pvalue)

    def filter_df(self):
        """Uses --level to remove peptides not in desired level."""

        if self.level == "nterm":
            logging.info("Level is N-term")

            try:
                self.df = self.df[
                    self.df[self.patterns['mod']].str.contains(pat=self.patterns['nterm'], na=False)
                ].reset_index(drop=True)
            except:
                logging.critical("Exit with code 3")
                raise TypeError(
                    f"Could not filter dataframe. Make sure software \
                    {self.software} is correct, and try again."
                )

        elif self.level == "quant":
            logging.info("Level is quant N-term")

            try:
                self.df = self.df[
                    self.df[self.patterns['mod']].str.contains(pat=self.patterns['nterm_label'], na=False)
                ].reset_index(drop=True)
            except:
                logging.critical("Exit with code 3")
                raise TypeError(
                    f"Could not filter dataframe. Make sure software \
                    {self.software} is correct, and try again."
                )

        else:
            logging.warning("Unrecognized level argument. Falling back to all")
            warnings.warn("Unrecognized level argument. Falling back to all...")

        if self.dropna:
            columns = self.df.columns[
                self.df.columns.str.contains(pat=self.patterns['quant'])
            ]
            self.df = self.df.dropna(subset=columns)

    def threaded_annotate(self):
        """Annotation with multiple threading.

        Fold change, mean, std and CV. Fold change is reported for each
        pairwise comparison of input conditions
        """

        length = len(self.df)
        batch_length = os.cpu_count()

        if self.nomerops is False:
            self.read_MEROPS()
            logging.info("Read MEROPS data")

        self.initialize_annotation(length)
        logging.info("Initialized annotation dataframe..")

        for i in tqdm(range(0, length, batch_length)):
            batch = list(range(i, i + batch_length))

            with concurrent.futures.ThreadPoolExecutor(
                max_workers=batch_length
            ) as executor:
                executor.map(self.entry_annotate, batch)

    def entry_annotate(self, loc):
        """Single entry annotation function used with multiple threading."""

        acc = annutils.parse_acc(self.df.loc[loc, self.patterns['acc']])
        if self.patterns['seq'] == 'Annotated Sequence':
            seq = annutils.parse_sequence(self.df.loc[loc, self.patterns['seq']])
        else:
            seq = self.df.loc[loc, self.patterns['seq']]

        ent = Entry(acc, seq)
        ent.get_record(self.sleeptime)

        if ent.record is not None:
            ent.parse_general()
            ent.parse_cleavage()

            if ent.cleavage_site is not None:
                ent.parse_protease()
                if self.nomerops is False:
                    ent.merops_protease(self.merops, self.merops_name)

            self.annot.loc[loc] = annutils.map_dict(self.annot.loc[loc], ent.annot)

        else:
            self.annot.loc[loc]["name"] = "HTTPError, not found"

    def annotate(self):
        """Main function that calls all other functions apart from exopeptidase
        and write_file."""

        length = len(self.df)
        if self.nomerops is False:
            self.read_MEROPS()
            logging.info("Read MEROPS data")

        self.initialize_annotation(length)
        logging.info("Initialized annotation dataframe..")

        for loc in tqdm(range(length)):

            acc = annutils.parse_acc(self.df.loc[loc, self.patterns['acc']])
            if self.patterns['seq'] == 'Annotated Sequence':
                seq = annutils.parse_sequence(self.df.loc[loc, self.patterns['seq']])
            else:
                seq = self.df.loc[loc, self.patterns['seq']]

            ent = Entry(acc, seq)
            ent.get_record(self.sleeptime)

            if ent.record is not None:
                ent.parse_general()
                ent.parse_cleavage()

                if ent.cleavage_site is not None:
                    ent.parse_protease()
                    if self.nomerops is False:
                        ent.merops_protease(self.merops, self.merops_name)

                self.annot.loc[loc] = annutils.map_dict(self.annot.loc[loc], ent.annot)

            else:
                self.annot.loc[loc]["name"] = "HTTPError, not found"

    def exopeptidase(self):
        """Annotate dipeptidase and aminopeptidase activity by checking
        sequences for rugging patterns."""

        cleared = set()
        self.annot["exopeptidase"] = np.nan
        sequences = self.annot["query_sequence"].dropna()
        sequences.sort_values(
            key=lambda x: x.str.len(), kind="mergesort", ascending=False, inplace=True
        )

        for seq in tqdm(sequences):

            if seq not in cleared:

                rag_flag = False
                cleared.add(seq)

                char_match = seq[-5:]
                matching_peptides = sequences[sequences.str.endswith(char_match)]

                if len(matching_peptides) > 1:
                    compare = seq
                    lpep = len(seq)

                    for ind in matching_peptides.index:

                        pep = matching_peptides[ind]
                        same_seq_indices = self.annot[
                            self.annot["query_sequence"] == pep
                        ].index

                        if pep == compare[1:]:
                            cleared.add(pep)
                            compare = pep
                            logging.info(f"1 {seq} {pep}")

                            if rag_flag and len(pep) == lpep - 1:
                                for i in same_seq_indices:
                                    self.annot.loc[
                                        i, "exopeptidase"
                                    ] = "Dipeptidase_seed_Aminopeptidase_activity"
                            else:
                                for i in same_seq_indices:
                                    self.annot.loc[
                                        i, "exopeptidase"
                                    ] = "Aminopeptidase_activity"
                                rag_flag = False

                        elif pep == compare[2:]:
                            cleared.add(pep)
                            compare = pep
                            rag_flag = True
                            lpep = len(pep)
                            logging.info(f"2 {seq} {pep}")

                            for i in same_seq_indices:
                                self.annot.loc[
                                    i, "exopeptidase"
                                ] = "Dipeptidase_activity"

    def visualize(self):
        """Calls Visualizer class and stores figure objects."""

        vis = Visualizer(self.df, self.annot, self.conditions, self.software, self.patterns)
        
        self.figures["General"] = vis.general()
        self.figures["CV"] = vis.cv_plot()
        self.figures["Piechart"] = vis.generate_pie_charts()
        self.figures["Heatmap"] = vis.heatmap()
        self.figures["Clustermap"] = vis.clustermap()

        if len(self.conditions) > 1:
            
            self.figures["Volcano"] = vis.volcano()
            self.figures["Fold"] = vis.fold_plot()
            self.figures["Fold_nterm"] = vis.fold_termini()

            print("Starting gallery generation...")
            logging.info("Starting gallery generation...")
            vis.gallery(stat=self.stat, cutoff=0.05, folder=self.outfolder)
            print("Finished gallery generation...")
            logging.info("Finished gallery generation...")

    def create_logos(self):
        """Create sequence logos.

        Based on logomaker package, bioconductor matrix calculation
        manual and algorithms in bioinformatics DTU course for PSSM
        construction
        """

        if len(self.conditions) > 1:

            conditions_iter = permutations(self.conditions.keys(), 2)
            for pair in conditions_iter:
                if self.stat:
                    
                    try:
                        column_name_test = f"Log10_ttest: {pair[0]}_{pair[1]}"
                        column_name_fold = f"Log2_fold_change: {pair[0]}/{pair[1]}"
                        column_test = self.annot[column_name_test]
                        column_fold = self.annot[column_name_fold]

                        ratio_name = column_name_test.split()[1].strip()

                        data = self.annot[(column_test < -1.5) & (column_fold > 1.5)]
                        data = data["p4_p4prime"].astype(str)
                        filtered = data[~data.str.contains("-|Not found|nan|X|Z|U|B|J|O")]
                        sequences = filtered.to_list()

                        if len(sequences) > 0:
                            figures = generate_logos(
                                sequences, ratio_name, self.pseudocounts, self.logo
                            )
                            self.figures[f"Logo_{ratio_name}_high"] = figures
                        else:
                            logging.debug(
                                f"Logo generation for {ratio_name} was skipped as \
                                sequences were 0"
                            )

                        data = self.annot[(column_test < -1.5) & (column_fold < -1.5)]
                        data = data["p4_p4prime"].astype(str)
                        filtered = data[~data.str.contains("-|Not found|nan|X|Z|U|B|J|O")]
                        sequences = filtered.to_list()

                        if len(sequences) > 0:
                            figures = generate_logos(
                                sequences, ratio_name, self.pseudocounts, self.logo
                            )
                            self.figures[f"Logo_{ratio_name}_low"] = figures
                        else:
                            logging.debug(
                                f"Logo generation for {ratio_name} was skipped as \
                                sequences were 0"
                            )              

                    except KeyError:
                        continue

                else:
                    column = f"Fold {pair[0]}/{pair[1]} significance"
                    data = self.annot[self.annot[column] == "significant high"]
                    data = data["p4_p4prime"].astype(str)
                    filtered = data[~data.str.contains("-|Not found|nan|X|Z|U|B|J|O")]
                    sequences = filtered.to_list()

                    if len(sequences) > 0:
                        figures = generate_logos(
                            sequences, pair[0], self.pseudocounts, self.logo
                        )
                        self.figures[f"Logo_{pair[0]}_high"] = figures
                    else:
                        logging.debug(
                            f"Logo generation for {pair[0]} was skipped as \
                            sequences were 0"
                        )

                    data = self.annot[self.annot[column] == "significant low"]
                    data = data["p4_p4prime"].astype(str)
                    filtered = data[~data.str.contains("-|Not found|nan|X|Z|U|B|J|O")]
                    sequences = filtered.to_list()

                    if len(sequences) > 0:
                        figures = generate_logos(
                            sequences, pair[0], self.pseudocounts, self.logo
                        )
                        self.figures[f"Logo_{pair[0]}_low"] = figures
                    else:
                        logging.debug(
                            f"Logo generation for {pair[0]} was skipped as \
                            sequences were 0"
                        )

        elif len(self.conditions) == 1:
            condition = list(self.conditions.keys())[0]
            data = self.annot["p4_p4prime"].astype(str)
            filtered = data[~data.str.contains("-|Not found|nan|X|Z|U|B|J|O")]
            sequences = filtered.to_list()

            if len(sequences) > 0:
                figures = generate_logos(
                    sequences, condition, self.pseudocounts, self.logo
                )
                self.figures[f"Logo_{condition}"] = figures
            else:
                logging.debug(
                    f"Logo generation for was skipped as \
                    sequences were 0"
                )

    def write_files(self):
        """Writes ouput files."""

        outfile = os.path.join(
            self.outfolder, self.outname
        )

        if self.separate:
            final_df = self.annot
        else:
            final_df = self.df.join(self.annot)

        if self.outfile_type == "csv":
            final_df.to_csv(outfile, sep=",", index=False)
        elif self.outfile_type == "tsv":
            final_df.to_csv(outfile, sep="\t", index=False)
        elif self.outfile_type == "xlsx":
            final_df.to_excel(outfile, engine="openpyxl", index=False)
        elif self.outfile_type == "json":
            final_df.to_json(outfile)
        elif self.outfile_type == "pkl":
            final_df.to_pickle(outfile, compression="infer")

        if len(self.figures) > 0:
            for k in self.figures:
                if self.figures[k] is not None:
                    for i in self.figures[k]:
                        if self.figures[k][i] is not None:
                            try:
                                # python 3.7 logos have fig.figure as None.
                                # 3.9 is fine, delete 'and' in if statement
                                if k != "Clustermap" and not k.startswith("Logo") and not k.startswith("General") and not k.startswith("Piechart"):
                                    self.figures[k][i].figure.savefig(
                                        os.path.join(self.outfolder, f"{k}_{i}.png"),
                                        format="png",
                                        dpi=300,
                                    )
                                    self.figures[k][i].figure.savefig(
                                        os.path.join(self.outfolder, f"{k}_{i}.svg"),
                                        format="svg",
                                    )
                                else:
                                    self.figures[k][i].savefig(
                                        os.path.join(self.outfolder, f"{k}_{i}.png"),
                                        format="png",
                                        dpi=300,
                                    )
                                    self.figures[k][i].savefig(
                                        os.path.join(self.outfolder, f"{k}_{i}.svg"),
                                        format="svg",
                                    )
                            except:
                                print(f"Skipped {k, self.figures[k]}")
                                logging.info(f"Skipped {k, self.figures[k]}")

        print(f"Finished. Wrote results to {self.outfolder}")
        logging.info(f"Finished. Wrote results to {self.outfolder}")

        try:
            log_basename = os.path.basename(self.logfile)
            shutil.copy(self.logfile, os.path.join(self.outfolder, log_basename))
        except: 
            pass

        shutil.make_archive(self.outfolder, "zip", self.outfolder)

        return None
