from __future__ import annotations

import logging
from numbers import Integral
from pathlib import Path

import numpy as np
import pandas as pd

from . import annutils


def load_available_models(accession_file: Path) -> list[str]:
    """Read available AlphaFold accession IDs from disk."""
    return annutils.read_alphafold_accessions(accession_file)


def structure_column_names(cleavagesitesize: int) -> tuple[str, str]:
    suffix = f"p{cleavagesitesize}_p{cleavagesitesize}prime"
    return (f"secondary_structure {suffix}", f"solvent_accessibility {suffix}")


def initialize_structure_columns(annot: pd.DataFrame, cleavagesitesize: int) -> tuple[str, str]:
    ss_col, sa_col = structure_column_names(cleavagesitesize)
    annot[ss_col] = np.nan
    annot[sa_col] = np.nan
    return ss_col, sa_col


def is_integer_like(value) -> bool:
    if pd.isna(value):
        return False
    if isinstance(value, Integral):
        return True
    if isinstance(value, (float, np.floating)):
        return float(value).is_integer()
    try:
        return float(value).is_integer()
    except (TypeError, ValueError):
        return False


def collect_cleavage_sites(frame: pd.DataFrame) -> dict[str, list[tuple[int, int]]]:
    """Group valid cleavage sites by accession.

    Returns mapping:
    accession -> [(row_index, cleavage_site), ...]
    """
    acc_cleavage_sites: dict[str, list[tuple[int, int]]] = {}

    for index, row in frame.iterrows():
        acc = row.get("query_accession")
        cleavage_site = row.get("p1_position")

        if pd.isna(acc) or not is_integer_like(cleavage_site):
            continue

        acc_cleavage_sites.setdefault(str(acc), []).append((index, int(cleavage_site)))

    return acc_cleavage_sites


def build_stat_columns(
    conditions: dict,
    pairwise: bool,
    multipletesting: bool,
) -> list[str]:
    """Build expected significance column names used for structural annotation."""
    cols: list[str] = []

    condition_names = list(conditions.keys())

    if len(condition_names) == 2 and not pairwise:
        base = f"Independent T-test p-value: {condition_names[0]} vs. {condition_names[1]}"
        cols.append(f"Corrected {base}" if multipletesting else base)
    elif pairwise:
        for cond_a in condition_names:
            for cond_b in condition_names:
                if cond_a == cond_b:
                    continue
                base = f"Independent T-test p-value: {cond_a} vs. {cond_b}"
                cols.append(f"Corrected {base}" if multipletesting else base)
    else:
        joined = " vs. ".join(condition_names)
        if multipletesting:
            cols.append(f"Corrected ANOVA p-value: {joined}")
        cols.append(f"ANOVA p-value: {joined}")

    # Preserve order while deduplicating.
    return list(dict.fromkeys(cols))


def calculate_structure_properties(
    acc_cleavage_sites: dict[str, list[tuple[int, int]]],
    temp_folder: Path,
    pymol_verbose: bool,
    available_models: list[str] | None,
):
    """Compute structural properties for grouped cleavage sites."""
    if not acc_cleavage_sites:
        return {}

    structure_tmp_filepath = temp_folder / "structure_properties.txt"
    return annutils.get_structure_properties(
        acc_cleavage_sites,
        structure_tmp_filepath,
        pymol_verbose,
        available_models,
    )


def apply_structure_properties(
    annot: pd.DataFrame,
    structure_properties: dict,
    ss_col: str,
    sa_col: str,
    only_when_missing: bool,
) -> None:
    for _, (index, ss, sa) in structure_properties.items():
        if index not in annot.index:
            continue

        if only_when_missing:
            if pd.notna(annot.loc[index, ss_col]) or pd.notna(annot.loc[index, sa_col]):
                continue

        annot.loc[index, ss_col] = ss
        annot.loc[index, sa_col] = sa


def annotate_significant_structure(
    annot: pd.DataFrame,
    columns: list[str],
    alpha: float,
    temp_folder: Path,
    pymol_verbose: bool,
    available_models: list[str] | None,
    ss_col: str,
    sa_col: str,
) -> None:
    for column_name in columns:
        if column_name not in annot.columns:
            logging.warning(f"Could not find statistics column '{column_name}' for structural annotation")
            continue

        subframe = annot[annot[column_name] <= alpha]
        acc_cleavage_sites = collect_cleavage_sites(subframe)
        structure_properties = calculate_structure_properties(
            acc_cleavage_sites,
            temp_folder,
            pymol_verbose,
            available_models,
        )
        apply_structure_properties(
            annot,
            structure_properties,
            ss_col,
            sa_col,
            only_when_missing=True,
        )
