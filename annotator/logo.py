import logomaker
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

alphabet = [
    "A",
    "R",
    "N",
    "D",
    "C",
    "Q",
    "E",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V",
]

background = {
    "A": 0.074,
    "R": 0.052,
    "N": 0.045,
    "D": 0.054,
    "C": 0.025,
    "Q": 0.034,
    "E": 0.054,
    "G": 0.074,
    "H": 0.026,
    "I": 0.068,
    "L": 0.099,
    "K": 0.058,
    "M": 0.025,
    "F": 0.047,
    "P": 0.039,
    "S": 0.057,
    "T": 0.051,
    "W": 0.013,
    "Y": 0.032,
    "V": 0.073,
}

blosum62 = {
    "A": {
        "A": 0.2901,
        "R": 0.0446,
        "N": 0.0427,
        "D": 0.041,
        "C": 0.065,
        "Q": 0.0559,
        "E": 0.0552,
        "G": 0.0783,
        "H": 0.042,
        "I": 0.0471,
        "L": 0.0445,
        "K": 0.057,
        "M": 0.0522,
        "F": 0.0338,
        "P": 0.0568,
        "S": 0.1099,
        "T": 0.073,
        "W": 0.0303,
        "Y": 0.0405,
        "V": 0.07,
    },
    "R": {
        "A": 0.031,
        "R": 0.345,
        "N": 0.0449,
        "D": 0.0299,
        "C": 0.0163,
        "Q": 0.0735,
        "E": 0.0497,
        "G": 0.0229,
        "H": 0.0458,
        "I": 0.0177,
        "L": 0.0243,
        "K": 0.1071,
        "M": 0.0321,
        "F": 0.019,
        "P": 0.0258,
        "S": 0.0401,
        "T": 0.0355,
        "W": 0.0227,
        "Y": 0.028,
        "V": 0.0219,
    },
    "N": {
        "A": 0.0256,
        "R": 0.0388,
        "N": 0.3169,
        "D": 0.069,
        "C": 0.0163,
        "Q": 0.0441,
        "E": 0.0405,
        "G": 0.0391,
        "H": 0.0534,
        "I": 0.0147,
        "L": 0.0142,
        "K": 0.0415,
        "M": 0.0201,
        "F": 0.0169,
        "P": 0.0233,
        "S": 0.0541,
        "T": 0.0434,
        "W": 0.0152,
        "Y": 0.0218,
        "V": 0.0165,
    },
    "D": {
        "A": 0.0297,
        "R": 0.031,
        "N": 0.0831,
        "D": 0.3974,
        "C": 0.0163,
        "Q": 0.0471,
        "E": 0.0902,
        "G": 0.0337,
        "H": 0.0382,
        "I": 0.0177,
        "L": 0.0152,
        "K": 0.0415,
        "M": 0.0201,
        "F": 0.0169,
        "P": 0.031,
        "S": 0.0489,
        "T": 0.0375,
        "W": 0.0152,
        "Y": 0.0187,
        "V": 0.0178,
    },
    "C": {
        "A": 0.0216,
        "R": 0.0078,
        "N": 0.009,
        "D": 0.0075,
        "C": 0.4837,
        "Q": 0.0088,
        "E": 0.0074,
        "G": 0.0108,
        "H": 0.0076,
        "I": 0.0162,
        "L": 0.0162,
        "K": 0.0086,
        "M": 0.0161,
        "F": 0.0106,
        "P": 0.0103,
        "S": 0.0175,
        "T": 0.0178,
        "W": 0.0076,
        "Y": 0.0093,
        "V": 0.0192,
    },
    "Q": {
        "A": 0.0256,
        "R": 0.0484,
        "N": 0.0337,
        "D": 0.0299,
        "C": 0.0122,
        "Q": 0.2147,
        "E": 0.0645,
        "G": 0.0189,
        "H": 0.0382,
        "I": 0.0133,
        "L": 0.0162,
        "K": 0.0535,
        "M": 0.0281,
        "F": 0.0106,
        "P": 0.0207,
        "S": 0.0332,
        "T": 0.0276,
        "W": 0.0152,
        "Y": 0.0218,
        "V": 0.0165,
    },
    "E": {
        "A": 0.0405,
        "R": 0.0523,
        "N": 0.0494,
        "D": 0.0914,
        "C": 0.0163,
        "Q": 0.1029,
        "E": 0.2965,
        "G": 0.0256,
        "H": 0.0534,
        "I": 0.0177,
        "L": 0.0202,
        "K": 0.0708,
        "M": 0.0281,
        "F": 0.019,
        "P": 0.0362,
        "S": 0.0524,
        "T": 0.0394,
        "W": 0.0227,
        "Y": 0.028,
        "V": 0.0233,
    },
    "G": {
        "A": 0.0783,
        "R": 0.0329,
        "N": 0.0652,
        "D": 0.0466,
        "C": 0.0325,
        "Q": 0.0412,
        "E": 0.035,
        "G": 0.5101,
        "H": 0.0382,
        "I": 0.0206,
        "L": 0.0213,
        "K": 0.0432,
        "M": 0.0281,
        "F": 0.0254,
        "P": 0.0362,
        "S": 0.0663,
        "T": 0.0434,
        "W": 0.0303,
        "Y": 0.0249,
        "V": 0.0247,
    },
    "H": {
        "A": 0.0148,
        "R": 0.0233,
        "N": 0.0315,
        "D": 0.0187,
        "C": 0.0081,
        "Q": 0.0294,
        "E": 0.0258,
        "G": 0.0135,
        "H": 0.355,
        "I": 0.0088,
        "L": 0.0101,
        "K": 0.0207,
        "M": 0.0161,
        "F": 0.0169,
        "P": 0.0129,
        "S": 0.0192,
        "T": 0.0138,
        "W": 0.0152,
        "Y": 0.0467,
        "V": 0.0082,
    },
    "I": {
        "A": 0.0432,
        "R": 0.0233,
        "N": 0.0225,
        "D": 0.0224,
        "C": 0.0447,
        "Q": 0.0265,
        "E": 0.0221,
        "G": 0.0189,
        "H": 0.0229,
        "I": 0.271,
        "L": 0.1154,
        "K": 0.0276,
        "M": 0.1004,
        "F": 0.0634,
        "P": 0.0258,
        "S": 0.0297,
        "T": 0.0533,
        "W": 0.0303,
        "Y": 0.0436,
        "V": 0.1646,
    },
    "L": {
        "A": 0.0594,
        "R": 0.0465,
        "N": 0.0315,
        "D": 0.028,
        "C": 0.065,
        "Q": 0.0471,
        "E": 0.0368,
        "G": 0.0283,
        "H": 0.0382,
        "I": 0.1679,
        "L": 0.3755,
        "K": 0.0432,
        "M": 0.1968,
        "F": 0.1142,
        "P": 0.0362,
        "S": 0.0419,
        "T": 0.0651,
        "W": 0.053,
        "Y": 0.0685,
        "V": 0.1303,
    },
    "K": {
        "A": 0.0445,
        "R": 0.1202,
        "N": 0.0539,
        "D": 0.0448,
        "C": 0.0203,
        "Q": 0.0912,
        "E": 0.0755,
        "G": 0.0337,
        "H": 0.0458,
        "I": 0.0236,
        "L": 0.0253,
        "K": 0.2781,
        "M": 0.0361,
        "F": 0.019,
        "P": 0.0413,
        "S": 0.0541,
        "T": 0.0454,
        "W": 0.0227,
        "Y": 0.0312,
        "V": 0.0261,
    },
    "M": {
        "A": 0.0175,
        "R": 0.0155,
        "N": 0.0112,
        "D": 0.0093,
        "C": 0.0163,
        "Q": 0.0206,
        "E": 0.0129,
        "G": 0.0094,
        "H": 0.0153,
        "I": 0.0368,
        "L": 0.0496,
        "K": 0.0155,
        "M": 0.1606,
        "F": 0.0254,
        "P": 0.0103,
        "S": 0.0157,
        "T": 0.0197,
        "W": 0.0152,
        "Y": 0.0187,
        "V": 0.0316,
    },
    "F": {
        "A": 0.0216,
        "R": 0.0174,
        "N": 0.018,
        "D": 0.0149,
        "C": 0.0203,
        "Q": 0.0147,
        "E": 0.0166,
        "G": 0.0162,
        "H": 0.0305,
        "I": 0.0442,
        "L": 0.0547,
        "K": 0.0155,
        "M": 0.0482,
        "F": 0.3869,
        "P": 0.0129,
        "S": 0.0209,
        "T": 0.0237,
        "W": 0.0606,
        "Y": 0.1308,
        "V": 0.0357,
    },
    "P": {
        "A": 0.0297,
        "R": 0.0194,
        "N": 0.0202,
        "D": 0.0224,
        "C": 0.0163,
        "Q": 0.0235,
        "E": 0.0258,
        "G": 0.0189,
        "H": 0.0191,
        "I": 0.0147,
        "L": 0.0142,
        "K": 0.0276,
        "M": 0.0161,
        "F": 0.0106,
        "P": 0.4935,
        "S": 0.0297,
        "T": 0.0276,
        "W": 0.0076,
        "Y": 0.0156,
        "V": 0.0165,
    },
    "S": {
        "A": 0.085,
        "R": 0.0446,
        "N": 0.0697,
        "D": 0.0522,
        "C": 0.0407,
        "Q": 0.0559,
        "E": 0.0552,
        "G": 0.0513,
        "H": 0.042,
        "I": 0.025,
        "L": 0.0243,
        "K": 0.0535,
        "M": 0.0361,
        "F": 0.0254,
        "P": 0.0439,
        "S": 0.2199,
        "T": 0.0927,
        "W": 0.0227,
        "Y": 0.0312,
        "V": 0.0329,
    },
    "T": {
        "A": 0.0499,
        "R": 0.0349,
        "N": 0.0494,
        "D": 0.0354,
        "C": 0.0366,
        "Q": 0.0412,
        "E": 0.0368,
        "G": 0.0297,
        "H": 0.0267,
        "I": 0.0398,
        "L": 0.0334,
        "K": 0.0397,
        "M": 0.0402,
        "F": 0.0254,
        "P": 0.0362,
        "S": 0.082,
        "T": 0.2465,
        "W": 0.0227,
        "Y": 0.028,
        "V": 0.0494,
    },
    "W": {
        "A": 0.0054,
        "R": 0.0058,
        "N": 0.0045,
        "D": 0.0037,
        "C": 0.0041,
        "Q": 0.0059,
        "E": 0.0055,
        "G": 0.0054,
        "H": 0.0076,
        "I": 0.0059,
        "L": 0.0071,
        "K": 0.0052,
        "M": 0.008,
        "F": 0.0169,
        "P": 0.0026,
        "S": 0.0052,
        "T": 0.0059,
        "W": 0.4924,
        "Y": 0.028,
        "V": 0.0055,
    },
    "Y": {
        "A": 0.0175,
        "R": 0.0174,
        "N": 0.0157,
        "D": 0.0112,
        "C": 0.0122,
        "Q": 0.0206,
        "E": 0.0166,
        "G": 0.0108,
        "H": 0.0573,
        "I": 0.0206,
        "L": 0.0223,
        "K": 0.0173,
        "M": 0.0241,
        "F": 0.0888,
        "P": 0.0129,
        "S": 0.0175,
        "T": 0.0178,
        "W": 0.0682,
        "Y": 0.3178,
        "V": 0.0206,
    },
    "V": {
        "A": 0.0688,
        "R": 0.031,
        "N": 0.027,
        "D": 0.0243,
        "C": 0.0569,
        "Q": 0.0353,
        "E": 0.0313,
        "G": 0.0243,
        "H": 0.0229,
        "I": 0.1767,
        "L": 0.0962,
        "K": 0.0328,
        "M": 0.0924,
        "F": 0.055,
        "P": 0.031,
        "S": 0.0419,
        "T": 0.071,
        "W": 0.0303,
        "Y": 0.0467,
        "V": 0.2689,
    },
}


class Logo:
    def __init__(
        self,
        sequences: list,
        name: str,
        pseudocounts=False,
    ):
        """Expects clean aligned sequences consisting only of amino acid strict
        alphabet."""

        self.sequences = sequences
        self.alphabet = alphabet
        self.blosum62 = blosum62
        self.background = background
        self.name = name
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

    def initialize_matrix(self):
        """Initalization of all matrices."""

        matrix = [0] * self.length
        for i in range(self.length):
            row = dict(zip(self.alphabet, [0.0] * len(self.alphabet)))
            matrix[i] = row

        return matrix

    def make_count_matrix(self):
        """Count matrix initialization."""

        matrix = self.initialize_matrix()

        for pos in range(self.length):
            for seq in self.sequences:

                matrix[pos][seq[pos]] += 1

        return matrix

    def make_frequency_matrix(self):
        """Create frequency matrix."""

        matrix = self.initialize_matrix()

        for pos in range(self.length):
            n = 0

            for peptide in self.sequences:
                matrix[pos][peptide[pos]] += 1
                n += 1

            for letter in alphabet:
                matrix[pos][letter] = matrix[pos][letter] / n

        return matrix

    def make_normalized_matrix(self):
        """Create pseudocount matrix based on blosum substitution matrix, then
        combine with the frequency matrix."""

        pseudo_matrix = self.initialize_matrix()

        for pos in range(self.length):
            for letter_1 in alphabet:
                for letter_2 in alphabet:
                    pseudo_matrix[pos][letter_1] += (
                        self.frequency_matrix[pos][letter_2]
                        * self.blosum62[letter_1][letter_2]
                    )

        norm_matrix = self.initialize_matrix()

        for pos in range(self.length):
            for a in self.alphabet:
                norm_matrix[pos][a] = (
                    self.alpha * self.frequency_matrix[pos][a]
                    + self.beta * pseudo_matrix[pos][a]
                ) / (self.alpha + self.beta)

        return norm_matrix

    def make_weighted_matrix(self):
        """Create weighted matrix from normalized matrix."""

        matrix = self.initialize_matrix()
        for pos in range(self.length):
            for letter in self.alphabet:
                if self.normalized_matrix[pos][letter] > 0:
                    matrix[pos][letter] = (
                        2
                        * np.log(
                            self.normalized_matrix[pos][letter]
                            / self.background[letter]
                        )
                        / np.log(2)
                    )
                else:
                    matrix[pos][letter] = 0

        return matrix

    def make_pssm(self):
        """Returns figure with weighted PSSM."""

        weighted_matrix = self.make_weighted_matrix()
        df = pd.DataFrame.from_records(weighted_matrix)
        self.weighted_matrix = df

        return self.make_logo(df, "weighted PSSM", "score")

    def make_probability(self):
        """Returns figure with frequency matrix."""

        df = pd.DataFrame.from_records(self.frequency_matrix)
        self.probability_matrix = df

        return self.make_logo(df, "probability PPSM", "percentage")

    def make_information(self):
        """Returns figure with information content PSSM.

        Includes pseudocounts
        """

        df = pd.DataFrame.from_records(self.normalized_matrix)

        tic = np.log2(len(alphabet))
        entropy = df.apply(lambda x: -x * np.log2(x)).sum(axis=1)
        fic = tic - entropy
        df = df.multiply(fic, axis=0)

        self.information_matrix = df

        return self.make_logo(df, "information content, Shannon", "information (bits)")

    def make_kullback(self):
        """Returns figure with Kullback-Leibler information content PSSM.

        Includes pseudocounts
        """

        df = pd.DataFrame.from_records(self.normalized_matrix)

        for i in df.index:
            for j in df.columns:
                df.loc[i][j] = df.loc[i][j] * np.log2(df.loc[i][j] / self.background[j])
                if np.isnan(df.loc[i][j]) or df.loc[i][j] < 0:
                    df.loc[i][j] = 0

        self.kullback_matrix = df

        return self.make_logo(df, "information content, Kullback", "information (bits)")

    def make_logo(self, frame, title, ylabel):
        """base class for all logo generation after dataframe caclucation."""

        fig, ax = plt.subplots(1, 1, figsize=[10, 6])
        logo = logomaker.Logo(
            frame,
            ax=ax,
            font_name="DejaVu Sans",
            color_scheme="NajafabadiEtAl2017",
            stack_order="big_on_top",
            center_values=False,
            flip_below=False,
            fade_below=0.7,
            shade_below=0.3,
            fade_probabilities=False,
            vpad=0.05,
            vsep=0.0,
            width=0.85,
            baseline_width=0.5,
        )

        # style using Logo methods
        logo.style_xticks(anchor=0, rotation=0)
        logo.style_spines(spines=["left", "right"], visible=False)

        if title != "weighted PSSM":
            for row in frame.index:
                for letter in frame.columns:
                    if frame.loc[row][letter] > 0.75 * frame.loc[row].sum():
                        logo.highlight_position(p=row, color="gold", alpha=0.35)

        # style using Axes methods
        logo.ax.set_xticklabels(
            x for x in ["P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'"]
        )
        logo.ax.axvline(3.5, color="goldenrod", linewidth=2, linestyle="-")
        logo.ax.set_xlim([-1, len(frame)])

        logo.ax.set_ylabel(ylabel)
        logo.ax.set_title(f"{self.name} {title} n={len(self.sequences)}")

        figure = logo.fig
        plt.close()

        return figure


def generate_logos(sequences, condition, pseudocounts, logo):
    """Generate different logos for the main class."""

    figures = {}
    pssm = Logo(sequences, condition, pseudocounts)

    if logo == "prob" or logo == "all":
        pm_logo = pssm.make_probability()
        figures["prob"] = pm_logo
    if logo == "pssm" or logo == "all":
        pssm_logo = pssm.make_pssm()
        figures["pssm"] = pssm_logo
    if logo == "shannon" or logo == "all":
        info_logo = pssm.make_information()
        figures["shannon"] = info_logo
    if logo == "kbl" or logo == "all":
        kb_logo = pssm.make_kullback()
        figures["kbl"] = kb_logo

    return figures
