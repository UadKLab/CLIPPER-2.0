import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import logomaker

from .alphabet import alphabet, background, blosum62

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

    def make_pssm(self, cleavagesitesize):
        """Returns figure with weighted PSSM."""

        weighted_matrix = self.make_weighted_matrix()
        df = pd.DataFrame.from_records(weighted_matrix)
        self.weighted_matrix = df

        return self.make_logo(df, "weighted PSSM", "score", cleavagesitesize)
    
    def make_probability(self, cleavagesitesize):
        """Returns figure with frequency matrix."""

        df = pd.DataFrame.from_records(self.frequency_matrix)
        self.probability_matrix = df

        return self.make_logo(df, "probability PPSM", "percentage", cleavagesitesize)
    

    def make_information(self, cleavagesitesize):
        """Returns figure with information content PSSM.

        Includes pseudocounts
        """

        df = pd.DataFrame.from_records(self.normalized_matrix)

        for i in df.index:
            for j in df.columns:
                df.loc[i][j] = df.loc[i][j] * np.log2(df.loc[i][j] / 0.05)
                if np.isnan(df.loc[i][j]) or df.loc[i][j] < 0:
                    df.loc[i][j] = 0


        self.information_matrix = df

        return self.make_logo(df, "information content, Shannon", "information (bits)", cleavagesitesize)


    def make_kullback(self, cleavagesitesize):
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

        return self.make_logo(df, "information content, Kullback", "information (bits)", cleavagesitesize)
    
    def make_logo(self, frame, title, ylabel, cleavagesitesize):
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

        xlabels = [x for x in [f"P{i}" for i in range(cleavagesitesize, 0, -1)] + [f"P{k}'" for k in range(1, cleavagesitesize + 1)]]
        logo.ax.set_xticklabels(xlabels)
        
        logo.ax.axvline(cleavagesitesize - 0.5, color="goldenrod", linewidth=2, linestyle="-")
        logo.ax.set_xlim([-1, len(frame)])

        logo.ax.set_ylabel(ylabel)
        logo.ax.set_title(f"{self.name} {title} n={len(self.sequences)}")

        figure = logo.fig
        plt.close()

        figure_heatmap = logo_heatmap(frame, xlabels, title)

        return figure, figure_heatmap

def logo_heatmap(df, xlabels, title, axs = None):
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Create a sample dataframe

    custom_x_labels = xlabels

    # Normalize the dataframe values to range [0, 1]
    df_normalized = (df - df.min().min()) / (df.max().max() - df.min().min())

    # Create a heatmap using seaborn with transposed data
    #plt.figure(figsize=(6, 12))

    # Set custom tick labels and center them
    #plt.xticks(ticks=np.arange(len(custom_x_labels)) + 0.5, labels=custom_x_labels, ha='center')
    #plt.yticks(ticks=np.arange(len(df.columns)) + 0.5, labels=df.columns, va='center')
    #plt.title(title)
    #sns.heatmap(df_normalized.T, cmap='RdBu_r', annot=False, cbar=False)
    plt.show()

    fig, ax = plt.subplots(1, 1, figsize=[6, 12])

    ax.set_xticks(ticks=np.arange(len(custom_x_labels)) + 0.5, labels=custom_x_labels, ha='center')
    ax.set_yticks(ticks=np.arange(len(df.columns)) + 0.5, labels=df.columns, va='center')

    plt.title(title)
    sns.heatmap(df_normalized.T, cmap='RdBu_r', annot=False, cbar=False)

    figure = fig.tight_layout()
    plt.close()

    return figure

def create_logo_helper(data, condition, pseudocounts, logo, cleavagesitesize):
    data = data[f"p{cleavagesitesize}_p{cleavagesitesize}prime"].astype(str)
    filtered = data[~data.str.contains("-|Not found|nan|X|Z|U|B|J|O")]
    sequences = filtered.to_list()

    if len(sequences) > 0:
        figures = generate_logos(sequences, condition, pseudocounts, logo, cleavagesitesize)
        return figures
    else:
        logging.debug(
            f"Logo generation for {condition} was skipped as sequences were 0"
        )
        return None


def generate_logos(sequences, condition, pseudocounts, logo, cleavagesitesize):
    """Generate different logos for the main class."""

    figures = {}
    pssm = Logo(sequences, condition, pseudocounts)
    if logo == "prob" or logo == "all":
        pm_logo = pssm.make_probability(cleavagesitesize)
        figures["prob"], figures["prob_heatmap"] = pm_logo
    if logo == "pssm" or logo == "all":
        pssm_logo = pssm.make_pssm(cleavagesitesize)
        figures["pssm"], figures["pssm_heatmap"] = pssm_logo
    if logo == "shannon" or logo == "all":
        info_logo = pssm.make_information(cleavagesitesize)
        figures["shannon"], figures["shannon_heatmap"] = info_logo
    if logo == "kbl" or logo == "all":
        kb_logo = pssm.make_kullback(cleavagesitesize)
        figures["kbl"], figures["kbl_heatmap"] = kb_logo

    return figures
