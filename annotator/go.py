# -*- coding: utf-8 -*-
"""Spyder Editor.

This is a temporary script file.
"""

import os

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats import multitest
from tqdm import tqdm


def read_basicgo():
    """Parses gene ontology terms from go-basic.obo.

    Downloaded from:
    http://geneontology.org/docs/download-ontology/#go_basic
    Returns list of dictionaries for each term
    """

    with open(
        r"C:\Users\konka\Documents\Python_Scripts\Annotator\annotator-4.0.0\data\go-basic.obo",
        "r",
    ) as fh:
        lines = fh.readlines()

    terms = {}
    term_dic = {}
    flag = False
    term = None

    for line in lines:
        line = line.strip()

        if line == "":
            terms[term] = term_dic
            flag = False
            term_dic = {}

        if flag:
            key, value = line.split(":", 1)
            if line.startswith("id"):
                term = value.strip()
            elif line.startswith("is_a"):
                if "parent" not in term_dic:
                    term_dic["parent"] = []
                    parent = value.strip().split(" ! ")[0].strip()
                    term_dic["parent"].append(parent)
                else:
                    parent = value.strip().split(" ! ")[0].strip()
                    term_dic["parent"].append(parent)
            elif line.startswith("relationship"):
                if "parent" not in term_dic:
                    term_dic["parent"] = []
                    parent = value.strip().split(" ")[1].strip()
                    term_dic["parent"].append(parent)
                else:
                    parent = value.strip().split(" ")[1].strip()
                    term_dic["parent"].append(parent)
            else:
                term_dic[key] = value.strip()

        if line == "[Term]":
            flag = True

    return terms


def read_taxgo(name="human", interaction=None):
    """Downloaded from:
    http://current.geneontology.org/products/pages/downloads.html Header
    information: http://geneontology.org/docs/go-annotation-file-gaf-
    format-2.1/.

    Trimmed header from file as it was spanning different mumber of
    lines for each organism
    """

    data = r"C:\Users\konka\Documents\Python_Scripts\Annotator\annotator-4.0.0\data"
    filename = f"goa_go_{name}.gaf"
    fh = os.path.join(data, filename)
    columns = ["database", "accesssion", "symbol", "interaction", "term"]
    frame = pd.read_csv(fh, sep="\t", header=None, names=columns, usecols=columns)

    if interaction is not None:
        frame = frame[frame["interaction"] == interaction]

    ref = {}
    ref["total"] = len(frame)
    for loc in frame.index:
        code = str(frame.loc[loc, "term"])
        if code != "nan":
            ref[code] = ref.get(code, 0) + 1

    return ref


def read_samplego(column):
    """Parse and return a dictionary of go terms for the sample."""

    sample = {}
    sample["total"] = column[~column.isna()].count()
    for loc in column.index:
        cell = str(column.loc[loc])
        if cell != "nan":
            codes = cell.split("|")
            for code in codes:
                sample[code] = sample.get(code, 0) + 1

    return sample


def map_go(frame):
    """Maps GO codes to GO names."""

    mapper = {}
    for loc in frame.index:
        cell = str(frame.loc[loc, "go_codes"])
        cell_names = str(frame.loc[loc, "go_names"])
        if cell != "nan":
            codes = cell.split("|")
            names = cell_names.split("|")
            for i in range(len(codes)):
                if codes[i] not in mapper:
                    mapper[codes[i]] = names[i]

    return mapper


df = pd.read_excel(
    r"C:\Users\konka\Documents\Python_Scripts\Annotator\annotator-4.0.0\results\results_230521151206\HUNTER_clean_annot.xlsx"
)
mapper = map_go(df)
sample = read_samplego(df["go_codes"])

lexicon = read_basicgo()
ref = read_taxgo()

for term in list(sample):
    if term in lexicon:
        parents = lexicon[term].pop("parent", None)
        if parents is not None:
            for parent in parents:
                if parent == "GO:0070062":
                    print("yes")
                popped = sample.pop(parent, None)

term = []
pval = []
sample_total = sample.pop("total")
ref_total = ref.pop("total")
for k in tqdm(sample):
    if k in ref and ref[k] > 10 and sample[k] > 5:
        tab = pd.DataFrame(
            data=[[ref[k], ref_total - ref[k]], [sample[k], sample_total - sample[k]]],
            index=["ref", "sample"],
            columns=["present", "absent"],
        )
        table = sm.stats.Table(tab)
        result = table.test_nominal_association()
        if result.pvalue < 1e-4:
            if (
                table.chi2_contribs.max().max()
                == table.chi2_contribs.loc["sample", "present"]
            ):
                term.append(k)
                pval.append(result.pvalue)

reject, new_pval, asid, abon = multitest.multipletests(
    pval, alpha=0.01, method="bonferroni"
)
names = [mapper[k] for k in sample]
enriched = pd.DataFrame(
    data=zip(term, names, new_pval), columns=["term", "names", "pval"]
)
