import pandas as pd
import os
import numpy as np
import re
from alphabet import alphabet, background, blosum62

class PSSM:
    def __init__(
        self,
        sequences: list,
        pseudocounts=False,
    ):
        """Expects clean aligned sequences consisting only of amino acid strict
        alphabet."""

        self.sequences = sequences
        self.alphabet = alphabet
        self.blosum62 = blosum62
        self.background = background
        self.length = len(sequences[0])
        self.beta = 50
        self.alpha = len(sequences) - 1
        self.usepseudo = pseudocounts

        self.count_matrix = self.make_count_matrix()
        self.frequency_matrix = self.make_frequency_matrix()
        if pseudocounts:
            self.normalized_matrix = self.make_normalized_matrix()
        else:
            self.normalized_matrix = self.frequency_matrix
        self.weighted_matrix = self.make_weighted_matrix()

    def initialize_matrix(self):
        """Initialization of all matrices."""
        return [dict(zip(self.alphabet, [0.0] * len(self.alphabet))) for _ in range(self.length)]

    def make_count_matrix(self):
        """Count matrix initialization."""
        return [{seq[pos]: cnt + 1 for seq in self.sequences for cnt in [matrix[pos].get(seq[pos], 0)]} for pos, matrix in enumerate([self.initialize_matrix()] * self.length)]

    def make_frequency_matrix(self):
        """Create frequency matrix."""
        n = len(self.sequences)
        return [{letter: cnt / n for letter, cnt in zip(self.alphabet, [sum(seq[pos] == letter for seq in self.sequences) for letter in self.alphabet])} for pos in range(self.length)]

    def make_normalized_matrix(self):
        """Create pseudocount matrix based on blosum substitution matrix, then
        combine with the frequency matrix."""
        pseudo_matrix = [{letter_1: sum(self.frequency_matrix[pos][letter_2] * self.blosum62[letter_1][letter_2] for letter_2 in self.alphabet) for letter_1 in self.alphabet} for pos in range(self.length)]

        return [{a: (self.alpha * self.frequency_matrix[pos][a] + self.beta * pseudo_matrix[pos][a]) / (self.alpha + self.beta) for a in self.alphabet} for pos in range(self.length)]

    def make_weighted_matrix(self):
        """Create weighted matrix from normalized matrix."""
        return [{letter: 2 * (np.log(self.normalized_matrix[pos][letter] / self.background[letter]) / np.log(2)) if self.normalized_matrix[pos][letter] > 0 else 0 for letter in self.alphabet} for pos in range(self.length)]


def read_protease_file(filepath):
    """
    Reads a file with MEROPS identifiers (one per line) and returns a list of identifiers.

    Parameters:
    filepath (str): The path to the file containing the protease identifiers.

    Returns:
    list: A list of protease identifiers.
    """

    with open(filepath, 'r') as f:
        protease_identifiers = [line.strip() for line in f.readlines()]

    return protease_identifiers


def create_protease_pssm(protease, df_cleavage, df_substrate, peptide_length=8, pseudocounts=False):
    """
    Creates a PSSM for a given protease.

    Parameters:
    protease (str): The identifier of the protease to create a PSSM for.
    df_cleavage (pd.DataFrame): The dataframe of the MEROPS cleavage file.
    df_substrate (pd.DataFrame): The dataframe of the MEROPS substrate file.
    peptide_length (int): The length of the peptides to consider for the PSSM.
    pseudocounts (bool): If true, uses pseudocounts in PSSM creation.

    Returns:
    PSSM: A PSSM object for the given protease.
    """
    
    # Center position of the peptide
    imipep = int(peptide_length/2)

    # Filter df_cleavage
    df_cleavage = df_cleavage[(df_cleavage['cleavage_evidence'] == 'experimental') & 
                              (df_cleavage['cleavage_type'] == 'physiological')]

    subdf = df_cleavage[df_cleavage['code'] == protease]
    substrate_accs = subdf['uniprot_acc'].values
    positions = subdf['p1'].values

    peptides = []
    for acc, pos in zip(substrate_accs, positions):
        acc_df = df_substrate[df_substrate['uniprot_acc'] == acc]
        
        if len(acc_df) != 0:
            seq = acc_df.iloc[0]['sequence']
            
            # pos is regular numbering, so index of P1' instead of P1. That is why peptide is -4:4
            if pos > imipep and len(seq) >= pos + imipep:
                peptide = seq[pos-imipep : pos+imipep]
                peptides.append(peptide)

    # Filter out None values and sequences with non-standard amino acids
    peptides = [peptide for peptide in peptides if peptide and set(peptide) <= set(alphabet)]
    
    # Create PSSM
    pssm = PSSM(sequences=peptides, pseudocounts=pseudocounts)
    
    return pssm.weighted_matrix


def construct_pssms(protease_codes, df_cleavage, df_substrate):
    """
    Given a list of protease codes, constructs PSSMs for each and returns them in a dictionary.
    
    Parameters:
    protease_codes (list): List of protease codes.
    df_cleavage (pd.DataFrame): Dataframe from the cleavage.txt file.
    df_substrate (pd.DataFrame): Dataframe from the substrate.txt file.
    pssm_type (str): Type of PSSM matrix to be constructed.
    
    Returns:
    dict: Dictionary of PSSMs keyed by protease code.
    """

    pssms = {}
    for code in protease_codes:
        # Construct the PSSM
        pssm = create_protease_pssm(code, df_cleavage, df_substrate, peptide_length=8, pseudocounts=False)
        pssms[code] = pssm

    return pssms


def score_peptide(peptide, matrix):
    """
    Scores a peptide sequence against a Position-Specific Scoring Matrix (PSSM).

    Parameters:
    peptide (str): The peptide sequence to be scored.
    matrix (list of dicts): The PSSM used for scoring.

    Returns:
    float: The score of the peptide against the given PSSM.
    """

    # Initialize the score
    acum = 0
    
    # For each residue in the peptide, add its PSSM score to the total score
    for i in range(len(peptide)):
        acum += matrix[i][peptide[i]]
    
    # Return the final score
    return acum


def score_proteases(pssms, peptide_sequence):
    """
    Given a dictionary of PSSMs, scores a given peptide sequence against each PSSM.
    
    Parameters:
    pssms (dict): Dictionary of PSSMs keyed by protease code.
    peptide_sequence (str): Amino acid sequence of the peptide to score.

    Returns:
    str: String of protease code: score separated by '|'.
    """

    scores = []
    for code, matrix in pssms.items():
        # Score peptide
        score = score_peptide(peptide_sequence, matrix)
        # Append to scores list
        scores.append(f'{code}:{score}')

    # Return string of joined protease code: score separated by '|'
    return '|'.join(scores)
