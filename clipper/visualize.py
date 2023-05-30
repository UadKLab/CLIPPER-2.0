import os
import tempfile
import gzip
import logging
import warnings
from pathlib import Path

warnings.filterwarnings('ignore')

from tqdm import tqdm

import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from umap import UMAP

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.collections as collections
from matplotlib.backends.backend_pdf import PdfPages

import pymol
from Bio import ExPASy
from Bio import SwissProt
from Bio import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from dna_features_viewer import GraphicFeature, GraphicRecord
from gprofiler import GProfiler
import networkx as nx

import annutils


alphafold_folder_name = r"\\ait-pdfs\services\BIO\Bio-Temp\Protease-Systems-Biology-temp\Kostas\CLIPPER\Datasets\Alphafold"


class Visualizer:

    """
    A class to visualize various aspects of proteomic data.

    This class provides methods to plot general statistics, create a volcano plot, plot CV values, 
    plot fold changes, fold termini, heatmaps, clustermaps, and pie charts of the data. It also 
    includes methods to generate a gallery of significant peptides, and to visualize the results 
    of PCA, UMAP, protein sequence and structure, functional enrichment, and pathway enrichment.

    Attributes:
        df (pandas.DataFrame): A dataframe containing the proteomic data.
        annot (str): The annotation.
        conditions (list): A list of conditions.
        software (str): The software used for the proteomic analysis.
        patterns (dict): A dictionary of patterns.
        pairwise (bool): A flag to indicate whether to perform pairwise analysis.

    Methods:
        general: Plots general statistics for the dataset.
        volcano: Creates a volcano plot for each condition pair.
        cv_plot: Creates a plot of CV values for all conditions.
        fold_plot: Creates a plot of fold changes of all peptides in the dataset for all conditions.
        fold_termini: Creates a plot of fold changes for internal and n-terminal peptides across all conditions.
        heatmap: Generates a heatmap of all quantification columns in the DataFrame.
        clustermap: Generates a clustermap of all quantification columns in the DataFrame.
        generate_pie_charts: Generates a series of pie charts for different peptide categories based on the DataFrame.
        gallery: Generates a PDF gallery of significant peptides from the DataFrame based on a given cutoff.
        get_significant_indices: Returns the indices of significant peptides in the DataFrame based on a given cutoff.
        pca_visualization: Performs Principal Component Analysis (PCA) on the quantification columns of the DataFrame and generates a visualization of the results.
        umap_visualization: Performs Uniform Manifold Approximation and Projection (UMAP) on the quantification columns of the DataFrame and generates a visualization of the results.
        plot_protein: Plots significant peptides for each condition on protein sequence and structure.
        plot_functional_enrichment: Performs a functional enrichment analysis and generates a bar plot of the results.
        plot_pathway_enrichment: Performs a pathway enrichment analysis and generates pathway plots for significant peptides.
    """
        
    def __init__(self, df, annot, conditions, software, patterns, pairwise=False):

        self.df = df
        self.annot = annot
        self.conditions = conditions
        self.software = software
        self.patterns = patterns
        self.pairwise = pairwise

    def general(self):

        """
        Plots general statistics for the dataset.

        Returns
        -------
        dict
            A dictionary with a key "general" and a value of a matplotlib Figure object of the plot.
        """

        acc_col = self.patterns['acc']
        seq_col = self.patterns['seq']
        mod_col = self.patterns['mod']
        quant_columns_pat = self.patterns['quant']
        label_pat = self.patterns['label']
        nter_pat = self.patterns['nterm']
        lab_nter_pat = self.patterns['nterm_label']

        data = {}
        data["proteins"] = self.df[acc_col].drop_duplicates().count()
        data["peptides"] = self.df[seq_col].drop_duplicates().count()

        if mod_col in self.df.columns:
            quant_columns = self.df.columns[self.df.columns.str.contains(quant_columns_pat)]

            if len(quant_columns) > 0:
                subframe = self.df[[seq_col, mod_col]]
                if self.software == "pd":
                    unique_mod = subframe
                    unique_mod.loc[:, mod_col] = unique_mod.loc[:, mod_col].astype(str)
                elif self.software == "sm":
                    unique_mod = subframe.drop_duplicates(subset=[seq_col, mod_col])
                quant_frame = self.df.dropna(subset=quant_columns, how="all")

                data["modified peptides"] = unique_mod[seq_col].count()
                data["labelled peptides"] = unique_mod[unique_mod[mod_col].str.contains(label_pat)][seq_col].count()
                data["nterm"] = unique_mod[unique_mod[mod_col].str.contains(nter_pat)][seq_col].count()
                data["quant proteins"] = quant_frame[acc_col].drop_duplicates().count()
                data["quant peptides"] = quant_frame[seq_col].drop_duplicates().count()
                data["quant mod peptides"] = quant_frame.drop_duplicates([seq_col, mod_col])[seq_col].count()

        ticks_no = list(range(len(data)))
        stat_names = [k for k in data]
        stat_vals = [data[k] for k in data]

        fig, ax = plt.subplots(1, 1, figsize=(12, 8))

        sns.barplot(x=ticks_no, y=stat_vals, palette="spring")
        for i in range(len(stat_vals)):
            x_pos = ticks_no[i]
            y_pos = stat_vals[i] / 2
            ax.text(x_pos, y_pos, stat_vals[i],
                fontfamily="sans-serif", ha="center", va="center", fontsize=12,)

        ax.set_yticks([int(i) for i in ax.get_yticks()])
        ax.set_xticklabels(stat_names, fontsize=14, fontfamily="sans-serif", rotation=90)
        ax.set_yticklabels([int(i) for i in ax.get_yticks()], fontsize=14, fontfamily="sans-serif")

        sns.despine(fig=fig, ax=ax)
        plt.subplots_adjust(bottom=0.3)

        fig = plt.gcf()
        plt.close()

        return {"general": fig}

    def volcano(self):

        """
        Creates a volcano plot for each condition pair.

        Returns
        -------
        dict
            A dictionary with each key being a condition pair and the value being a corresponding matplotlib Figure object of the plot.
        """

        columns_fold = self.annot.columns[self.annot.columns.str.contains("Log2_fold_change:")]
        columns_ttest = self.annot.columns[self.annot.columns.str.contains("Log10_ttest:")]
        figures = {}

        for test in columns_ttest:
            conditions = test.split()[1].split("_")

            for fold in columns_fold:

                if conditions[0] in fold and conditions[1] in fold:
                    frame = self.annot.loc[:, [fold, test]].dropna(how="any")
                    # log values are negative for ttest column, inverse
                    frame[test] = frame[test] * -1
                    frame["coding"] = np.where((~frame[fold].between(-1.5, 1.5)) & (~frame[test].between(-1.5, 1.5)), "1", "0",)
                    frame = frame.sort_values("coding")

                    if len(frame.coding.unique()) == 1:
                        colors = ["grey"]
                    else:
                        colors = ["grey", "crimson"]

                    g = sns.scatterplot(data=frame, x=fold, y=test, hue="coding", palette=colors)
                    ax = plt.gca()
                    ax.set_ylabel("- " + test)
                    xmin, xmax = ax.get_xlim()
                    ymin, ymax = ax.get_ylim()
                    ax.plot([xmin, xmax], [1.5, 1.5], "--b", alpha=0.2)
                    ax.plot([1.5, 1.5], [0, ymax], "--b", alpha=0.2)
                    ax.plot([-1.5, -1.5], [0, ymax], "--b", alpha=0.2)
                    ax.legend().set_visible(False)

                    name = fold.split(":")[1].strip()
                    name = "_".join(name.split("/"))
                    figures[name] = g
                    plt.close()

        return figures

    def cv_plot(self):

        """
        Creates a plot of CV values for all conditions.

        Returns
        -------
        dict
            A dictionary with a key "cv_plot" and a value of a matplotlib Figure object of the CV plot.
        """

        colors = sns.color_palette("husl", len(self.conditions))

        columns = self.annot.columns[self.annot.columns.str.contains("_CV")]

        if len(columns) > 0:
            fig, ax = plt.subplots()

            for i in range(len(columns)):
                g = sns.kdeplot(self.annot[columns[i]], color=colors[i], ax=ax, label=columns[i][:-3],)

            plt.legend()
            plt.xlabel("")
            plt.close()

            return {"cv_plot": g}

    def fold_plot(self):

        """
        Creates a plot of fold changes of all peptides in the dataset for all conditions.

        Returns
        -------
        dict
            A dictionary with each key being a condition and the value being a corresponding matplotlib Figure object of the plot.
        """

        columns = self.annot.columns[self.annot.columns.str.contains("Log2_fold_change:")]
        figures = {}
        colors = sns.color_palette("husl", (len(self.conditions) - 1)*len(self.conditions))

        for i in range(len(columns)):
            fig = sns.histplot(self.annot[columns[i]], kde=True, color=colors[i])
            name = columns[i].split(":")[1].strip()
            name = "_".join(name.split("/"))
            figures[name] = fig
            plt.close()

        return figures

    def fold_termini(self):

        """
        Creates a plot of fold changes for internal and n-terminal peptides across all conditions.

        Returns
        -------
        dict
            A dictionary with each key being a condition appended with 'internal_nterm' or 'natural_nterm' and the value being a corresponding matplotlib Figure object of the plot.
        """

        columns = self.annot.columns[self.annot.columns.str.contains("Log2_fold_change:")]
        figures = {}
        colors = sns.color_palette("husl", (len(self.conditions) - 1)*len(self.conditions)*2)
        
        subframe_internal = self.annot[self.annot["nterm_annot"] == "Internal"]
        subframe_natural = self.annot[self.annot["nterm_annot"] != "Internal"]
        
        ind = 0
        for i in range(len(columns)):
            
            fig = sns.histplot(subframe_internal[columns[i]], kde=True, color=colors[ind])
            ind += 1
            name = columns[i].split(":")[1].strip()
            name = "_".join(name.split("/")) + '_internal_nterm'
            figures[name] = fig
            plt.close()

            fig = sns.histplot(subframe_natural[columns[i]], kde=True, color=colors[ind])
            ind += 1
            name = columns[i].split(":")[1].strip()
            name = "_".join(name.split("/")) + '_natural_nterm'
            figures[name] = fig
            plt.close()

        return figures

    def heatmap(self):

        """
        Generate a heatmap of all quantification columns in the DataFrame.
        
        The heatmap uses logarithmic scaling with base 10. Zero values in the data are replaced with np.nan.
        Heatmap does not display y-tick labels for clarity.
        
        Returns:
            dict: A dictionary with "heatmap" as a key and the seaborn heatmap figure as the value.
        """

        columns = self.df.columns[self.df.columns.str.contains(self.patterns['quant'])]

        if len(columns) > 0:
            heat = self.df[columns]
            heat = np.log10(heat.replace(0, np.nan).dropna())
            fig = sns.heatmap(heat, yticklabels=False)
            plt.subplots_adjust(bottom=0.4)
            plt.close()

            return {"heatmap": fig}

    def clustermap(self):

        """
        Generate a clustermap of all quantification columns in the DataFrame.

        The clustermap uses logarithmic scaling with base 10. Zero values in the data are replaced with np.nan.
        Clustermap does not display y-tick labels for clarity.
        
        Returns:
            dict: A dictionary with "clustermap" as a key and the seaborn clustermap figure as the value.
        """

        columns = self.df.columns[self.df.columns.str.contains(self.patterns['quant'])]

        if len(columns) > 0:
            cluster = self.df[columns]
            cluster = np.log10(cluster.replace(0, np.nan).dropna())
            fig = sns.clustermap(cluster, yticklabels=False)
            plt.close()

            return {"clustermap": fig.fig}
    
    def generate_pie_charts(self):
        """
        Generate a series of pie charts for different peptide categories based on the DataFrame. 
        
        Pie charts represent characteristics such as the proportion of lysines that are labelled, 
        the proportion of N-termini that are acetylated, etc. The specific pie charts created 
        depend on the peptide modifications and annotations present in the data. 

        Returns:
            dict: A dictionary where keys are the names of the pie charts and values are the generated figures.
        """

        figures = {}  
        
        if self.patterns['mod'] in self.df.columns:
            # Lysine pie chart
            lysines = len(self.df[self.df[self.patterns['mod']].str.contains(r"K", na=False)])
            labelled = len(self.df[self.df[self.patterns['mod']].str.contains(self.patterns['lysine_label'], na=False)])
            if lysines - labelled > 0:
                categories = {"Free lysines": lysines - labelled,"Labelled lysines": labelled,}
                y = np.array([categories[k] for k in categories])
                if y.sum() > 0:
                    x = [k for k in categories]
                    colors = ["violet", "olivedrab"]
                    explode = [0, 0.1]
                    figures["lysine"] = create_pie_chart(y, x, colors, explode)

            # N-term pie chart
            total = len(self.df[self.patterns['mod']])
            nterm = len(self.df[self.df[self.patterns['mod']].str.contains(self.patterns['nterm'], na=False)])
            categories = {"Tryptic": total - nterm, "N-termini": nterm}
            y = np.array([categories[k] for k in categories])
            if y.sum() > 0:
                x = [k for k in categories]
                colors = ["salmon", "darkcyan"]
                explode = [0, 0.1]
                figures["nterm"] = create_pie_chart(y, x, colors, explode)

            # N-term acetylation vs labeling pie chart
            if self.software == "sm":
                term = 'Labelled'
                labelled = len(self.df[self.df[self.patterns['mod']].str.contains(self.patterns['nterm_label'], na=False)])
                acetyl = len(self.df[self.df[self.patterns['mod']].str.contains(r"Acetyl \(Protein N-term\)", na=False)])
            elif self.software == "pd":
                acetyl = len(self.df[self.df[self.patterns['mod']].str.contains(r"Acetyl \[N-Term\]", na=False)])
                if self.patterns['quant']:
                    term = 'Labelled'
                    labelled = len(self.df[self.df[self.patterns['mod']].str.contains(self.patterns['nterm_label'], na=False)])
                else:
                    term = 'Non acetylated'
                    total = len(self.df[self.df[self.patterns['mod']].str.contains(self.patterns['nterm'], na=False)])
                    labelled = total - acetyl
            
            categories = {"Acetylated": acetyl, term: labelled}
            y = np.array([categories[k] for k in categories])
            if y.sum() > 0:
                x = [k for k in categories]
                colors = ["teal", "saddlebrown"]
                explode = [0, 0.1]
                figures["modification"] = create_pie_chart(y, x, colors, explode)

            # Cysteine pie chart
            cysteines = len(self.df[self.df[self.patterns['mod']].str.contains(r"C", na=False)])
            alkylated = len(self.df[self.df[self.patterns['mod']].str.contains(r"Carbamidomethyl", na=False)])
            if cysteines - alkylated > 0: 
                categories = {"Free cysteines": cysteines - alkylated,"Carbamidomethyl": alkylated,}
                y = np.array([categories[k] for k in categories])
                if y.sum() > 0:
                    x = [k for k in categories]
                    colors = ["indianred", "cadetblue"]
                    explode = [0, 0.1]
                    figures["cysteine"] = create_pie_chart(y, x, colors, explode)


        # Exopeptidase pie chart
        if "exopeptidase" in self.annot.columns:
            exopep = (~self.annot["exopeptidase"].isna()).sum()
            noexo = self.annot["exopeptidase"].isna().sum()

            categories = {"Exopeptidase": exopep, "No activity": noexo}
            y = np.array([categories[k] for k in categories])
            if y.sum() > 0:
                x = [k for k in categories]
                colors = ["darkolivegreen", "mediumslateblue"]
                explode = [0.1, 0]
                figures["exopep"] = create_pie_chart(y, x, colors, explode)

        # N-term annotation pie charts
        column = self.annot["nterm_annot"].dropna()

        # all categories for the pie charts below
        categories = {
            "internal": (column.values == "Internal").sum(),
            "met removed": (column.values == "Met removed").sum(),
            "met intact": (column.values == "Met intact").sum(),
            "signal removed": (column.values == "Signal removed").sum(),
            "propeptide removed": (column.values == "Propeptide removed").sum(),
            "within signal": (column.values == "Cleavage within signal peptide range").sum(),
            "within propeptide": (column.values == "Cleavage within propeptide range").sum(),
            "transit removed": (column.values == "Transit peptide removed").sum(),
            "within transit peptide": (column.values == "Cleavage within transit peptide range").sum(),}

        y = np.array([categories.pop("internal"), sum(categories.values())])
        if y.sum() > 0:
            x = ["Internal", "Natural"]
            explode = [0.1, 0]
            colors = ["darkgreen","peru",]
            figures["type"] = create_pie_chart(y, x, colors, explode)

        y = np.array([categories.pop("met removed"), categories.pop("met intact"), sum(categories.values()),])
        if y.sum() > 0:
            x = ["Met removal", "Met intact", "Other"]
            colors = ["crimson", "lightpink", "seagreen"]
            explode = [0, 0.1, 0.1]
            figures["natural"] = create_pie_chart(y, x, colors, explode)

        y = np.array([categories[k] for k in categories])
        if y.sum() > 0:
            x = [k for k in categories]
            colors = ["steelblue", "goldenrod", "slategrey", "forestgreen", "violet", "tomato"]
            explode = [0, 0.1, 0, 0.1, 0, 0]
            figures["other"] = create_pie_chart(y, x, colors, explode)

        return figures

    def gallery(self, stat=False, cutoff=0.05, folder=None):

        """
        Generate a PDF gallery of significant peptides from the DataFrame based on a given cutoff.
        
        Args:
            stat (bool): A boolean flag to specify if the peptides to consider are statistically significant.
            cutoff (float): A threshold for determining the significance of peptides.
            folder (str): A string representing the path to the folder where the gallery PDF will be saved.

        Note: If the folder does not exist, it will be created.
        """

        indices = self.get_significant_indices(stat, cutoff)
        accs = self.annot.loc[indices, "query_accession"].unique()

        sns.set("paper", "ticks")
        pp = PdfPages(os.path.join(folder, "Sig_gallery.pdf"))

        for acc in tqdm(accs):

            natural = self.annot[(self.annot["query_accession"] == acc) & (self.annot["nterm_annot"].isin(["Met removed",
                            "Met intact", "Signal removed", "Propeptide removed",]))]

            internal = self.annot[(self.annot["query_accession"] == acc) & (self.annot["nterm_annot"] == "Internal")]

            columns = natural.columns[natural.columns.str.contains("_mean")]
            natural = natural.drop_duplicates("query_sequence").dropna(subset=columns)
            internal = internal.drop_duplicates("query_sequence").dropna(subset=columns)

            if len(natural) > 0 and len(internal) > 0:

                fig, ax = plt.subplots(1, 2, figsize=(12, 10))

                # Generate the mean values bar plot
                mean_values_data = get_mean_values_data(columns, natural, internal)
                generate_bar_plot(ax[0], mean_values_data, "summer", "Mean values per condition")

                # Generate the peptides with quant values bar plot
                quant_values_data = get_quant_values_data(columns, natural, internal)
                generate_bar_plot(ax[1], quant_values_data, "winter", "Peptides with quant values")

                plt.suptitle(f"Protein Uniprot ID: {acc}")
                plt.tight_layout()
                plt.savefig(pp, format="pdf")
                plt.close()

        pp.close()

    def get_significant_indices(self, stat, cutoff):

        """
        Returns the indices of significant peptides in the DataFrame based on a given cutoff.
        
        Args:
            stat (bool): A boolean flag to specify if the peptides to consider are statistically significant.
            cutoff (float): A threshold for determining the significance of peptides.
        
        Returns:
            list: A list of indices of the significant peptides.
        """
            
        if stat:
            if len(self.conditions) == 2 or self.pairwise:
                cols = self.annot.columns[self.annot.columns.str.startswith("Ttest:")]
            else:
                cols = self.annot.columns[self.annot.columns.str.startswith("ANOVA:")]
            indices = []
            for col in cols:
                ind = self.annot[self.annot[col] < cutoff].index
                indices += list(ind)
            indices = list(set(indices))
        else:
            cols = self.annot.columns[self.annot.columns.str.contains("significance")]
            indices = []
            for col in cols:
                ind = self.annot[self.annot[col].str.contains("significant", na=False)].index
                indices += list(ind)
            indices = list(set(indices))

        return indices
    
    def pca_visualization(self):

        """
        Performs Principal Component Analysis (PCA) on the quantification columns of the DataFrame 
        and generates a visualization of the results.
        
        Returns:
            dict: A dictionary with the key "PCA" and the value is the PCA figure.
        """

        # Get data and perform PCA on quantification columns
        data_columns = self.df.columns[self.df.columns.str.contains(self.patterns['quant'])]
        # Get data in arrays and replace Nan values with 0
        data = self.df[data_columns].fillna(0).values.T

        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(data)

        # Prepare a dictionary to map condition names to colors
        condition_colors = {key: plt.cm.tab10(i) for i, key in enumerate(self.conditions.keys())}

        # Create a DataFrame with PCA results, conditions, and colors
        pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2'])
        pca_df['condition'] = None
        pca_df['color'] = None

        condition_df = pd.DataFrame(index=data_columns, columns=['condition', 'color', 'annotation'])

        for condition, substrings in self.conditions.items():
            for substring in substrings:
                condition_idx = data_columns.str.contains(substring)
                condition_df.loc[condition_idx, 'condition'] = condition
                condition_colors_list = [condition_colors[condition]] * sum(condition_idx)
                condition_df.loc[condition_idx, 'color'] = condition_colors_list
                condition_df.loc[condition_idx, 'annotation'] = substring

        pca_df['condition'] = condition_df['condition'].values
        pca_df['color'] = condition_df['color'].values
        pca_df['annotation'] = condition_df['annotation'].values

        # Create a Figure object and visualize the PCA results
        fig, ax = plt.subplots()
        for condition, color in condition_colors.items():
            idx = pca_df['condition'] == condition
            ax.scatter(pca_df.loc[idx, 'PC1'], pca_df.loc[idx, 'PC2'], c=[color], label=condition)
            
            # Add text annotations above the data points
            for i, txt in enumerate(pca_df.loc[idx, 'annotation']):
                ax.annotate(txt, (pca_df.loc[idx, 'PC1'].iloc[i], pca_df.loc[idx, 'PC2'].iloc[i]),
                            textcoords="offset points", xytext=(0, 5), ha='center', fontsize=8)

        ax.set_xlabel('PC1 ({:.2f}%)'.format(pca.explained_variance_ratio_[0] * 100))
        ax.set_ylabel('PC2 ({:.2f}%)'.format(pca.explained_variance_ratio_[1] * 100))
        ax.legend()
        plt.tight_layout()
        plt.close()

        return {"PCA": fig}

    def umap_visualization(self):

        """
        Performs Uniform Manifold Approximation and Projection (UMAP) on the quantification columns 
        of the DataFrame and generates a visualization of the results.
        
        Returns:
            dict: A dictionary with the key "UMAP" and the value is the UMAP figure.
        """

        # Get data and perform UMAP on quantification columns
        data_columns = self.df.columns[self.df.columns.str.contains(self.patterns['quant'])]
        # Get data in arrays and replace Nan values with 0
        data = self.df[data_columns].fillna(0).values.T

        reducer = UMAP()
        umap_result = reducer.fit_transform(data)

        # Prepare a dictionary to map condition names to colors
        condition_colors = {key: plt.cm.tab10(i) for i, key in enumerate(self.conditions.keys())}

        # Create a DataFrame with UMAP results, conditions, and colors
        umap_df = pd.DataFrame(umap_result, columns=['UMAP1', 'UMAP2'])
        umap_df['condition'] = None
        umap_df['color'] = None

        condition_df = pd.DataFrame(index=data_columns, columns=['condition', 'color', 'annotation'])

        for condition, substrings in self.conditions.items():
            for substring in substrings:
                condition_idx = data_columns.str.contains(substring)
                condition_df.loc[condition_idx, 'condition'] = condition
                condition_colors_list = [condition_colors[condition]] * sum(condition_idx)
                condition_df.loc[condition_idx, 'color'] = condition_colors_list
                condition_df.loc[condition_idx, 'annotation'] = substring

        umap_df['condition'] = condition_df['condition'].values
        umap_df['color'] = condition_df['color'].values
        umap_df['annotation'] = condition_df['annotation'].values

        # Create a Figure object and visualize the UMAP results
        fig, ax = plt.subplots()
        for condition, color in condition_colors.items():
            idx = umap_df['condition'] == condition
            ax.scatter(umap_df.loc[idx, 'UMAP1'], umap_df.loc[idx, 'UMAP2'], c=[color], label=condition)
            
            # Add text annotations above the data points
            for i, txt in enumerate(umap_df.loc[idx, 'annotation']):
                ax.annotate(txt, (umap_df.loc[idx, 'UMAP1'].iloc[i], umap_df.loc[idx, 'UMAP2'].iloc[i]),
                            textcoords="offset points", xytext=(0, 5), ha='center', fontsize=8)

        ax.set_xlabel('UMAP1')
        ax.set_ylabel('UMAP2')
        ax.legend()
        plt.tight_layout()
        plt.close()

        return {"UMAP": fig}
    
    def plot_protein(self, cutoff=0.05, folder=None, merops=None, alphafold=None, level=None):

        """
        Plots significant peptides for each condition on protein sequence and structure.
        
        Args:
            cutoff (float): A threshold for determining the significance of peptides.
            folder (str): A string representing the path to the folder where the plots will be saved.
            merops (str): A string specifying the MEROPS database identifier for the protein.
            alphafold (str): A string specifying the AlphaFold model identifier for the protein.
            level (str): A string specifying the level of detail to include in the plot.

        Note: If the folder does not exist, it will be created.
        """

        if len(self.conditions) == 2 or self.pairwise:
            cols = self.annot.columns[self.annot.columns.str.startswith("Ttest:")]
        else:
            cols = self.annot.columns[self.annot.columns.str.startswith("ANOVA:")]

        for col in cols:
            df = self.annot[self.annot[col] < cutoff]
            accs = df["query_accession"].unique()

            # Create a PDF file to save all figures
            pdf_path = os.path.join(folder, f"{col[7:]}_sequence_plots.pdf")
            with PdfPages(pdf_path) as pp:
                for acc in tqdm(accs):
                    subframe = df[df["query_accession"] == acc]
                    if len(subframe) > 0:
                        plot_protein_figure(pp, subframe, acc, col, merops, alphafold, level)

    def plot_functional_enrichment(self, cutoff=0.05):

        """
        Performs a functional enrichment analysis and generates a bar plot of the results.
        
        Args:
            cutoff (float): A threshold for determining the significance of peptides.
        
        Returns:
            dict: A dictionary where the keys are condition names and the values are the corresponding figures.
        """

        if len(self.conditions) == 2 or self.pairwise:
            cols = self.annot.columns[self.annot.columns.str.startswith("Ttest:")]
        else:
            cols = self.annot.columns[self.annot.columns.str.startswith("ANOVA:")]

        figures = {}
        for col in cols:
            df = self.annot[self.annot[col] < cutoff]
            gene_names = df["name"].unique()
            genes = [g.split('_')[0] for g in gene_names]

            col_name = col.split(":")[1].strip()
            fig = plot_enrichment_figure(genes, col_name)
            
            figures[col_name] = fig

        return figures
    
    def plot_pathway_enrichment(self, cutoff=0.05, folder=None):

        """
        Performs a pathway enrichment analysis and generates pathway plots for significant peptides.
        
        Args:
            cutoff (float): A threshold for determining the significance of peptides.
            folder (str): A string representing the path to the folder where the pathway plots will be saved.

        Note: If the folder does not exist, it will be created.
        """

        cols = self.annot.columns[self.annot.columns.str.startswith("Ttest:")]

        for col in cols:
            col_fold = col.replace('_', '/').replace('Ttest', 'Log2_fold_change')
            subframe = self.annot[self.annot[col] < cutoff]

            # Create a PDF file to save all figures
            col_name = col.split(":")[1].strip()
            pdf_path = os.path.join(folder, f"{col_name}_pathway_plots.pdf")
            with PdfPages(pdf_path) as pp:
                plot_pathway_figures(subframe, col_fold, pp, cutoff)

    
def create_pie_chart(y, x, colors, explode):

    """
    Creates a pie chart of the given data.
    
    Parameters
    ----------
    y : pandas.Series or numpy.array
        The values for each pie slice.
    x : list
        The labels for each pie slice.
    colors : list
        The colors for each pie slice.
    explode : list
        Specifies the fraction of the radius with which to offset each pie slice.
        
    Returns
    -------
    matplotlib.figure.Figure
        The created figure object.
    """

    porcent = 100.0 * y / y.sum()
    patches, texts = plt.pie(y, explode=explode, colors=colors, startangle=90, normalize=True)
    
    labels = ["{0} - {1:1.2f} %".format(i, j) for i, j in zip(x, porcent)]
    patches, labels, dummy = zip(*sorted(zip(patches, labels, y), key=lambda x: x[2], reverse=True))
    
    plt.legend(patches, labels, loc="center left", bbox_to_anchor=(-0.5, 0.5), fontsize=8)
    plt.subplots_adjust(left=0.4)
    
    figure = plt.gcf()
    plt.close()
    
    return figure


def generate_bar_plot(ax, data, palette, title, xticklabels_rotation=90):

    """
    Generate a bar plot with the given data and parameters.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The matplotlib axes to draw on.
    data : dict
        The data to plot.
    palette : list or dict
        The colors to use for the bars.
    title : str
        The title of the plot.
    xticklabels_rotation : int, optional
        The rotation angle of the x tick labels.

    Returns
    -------
    None
    """

    names = [k for k in data.keys()]
    values = [v for v in data.values()]
    sns.barplot(x=names, y=values, palette=palette, ax=ax)
    ax.set_xticklabels(names, fontsize=12, fontfamily="sans-serif", rotation=xticklabels_rotation)
    ax.set_title(title, fontsize=11)


def get_mean_values_data(columns, natural, internal):
        
    """
    Get the mean values for the given columns in the natural and internal dataframes.

    Parameters
    ----------
    columns : list
        The columns to compute the mean of.
    natural : pandas.DataFrame
        The dataframe representing the natural data.
    internal : pandas.DataFrame
        The dataframe representing the internal data.

    Returns
    -------
    dict
        A dictionary with the column names as keys and the computed means as values.
    """

    type_term = ["natural", "internal"]
    data = {}
    for col in columns:
        for t, subframe in zip(type_term, [natural, internal]):
            data[f"{col}_{t}"] = subframe[col].mean()
    return data


def get_quant_values_data(columns, natural, internal):

    """
    Get the quant values for the given columns in the natural and internal dataframes.

    Parameters
    ----------
    columns : list
        The columns to get the values of.
    natural : pandas.DataFrame
        The dataframe representing the natural data.
    internal : pandas.DataFrame
        The dataframe representing the internal data.

    Returns
    -------
    dict
        A dictionary with the column names as keys and the retrieved values as values.
    """

    type_term = ["natural", "internal"]
    data = {}
    for t, subframe in zip(type_term, [natural, internal]):
        for row in subframe.index:
            for col in columns:
                data[f"{col[:-5]}_{t}_{subframe.loc[row, 'query_sequence']}"] = subframe.loc[row, col]
    return data


def extract_protein_features(acc, record, merops, subframe):

    """
    Extracts the protein features from the UniProt record and the MEROPS database.

    Parameters
    ----------
    acc : str
        The UniProt accession number.
    record : Bio.SeqRecord.SeqRecord
        The UniProt record.
    merops : pandas.DataFrame
        The MEROPS database dataframe.
    subframe : pandas.DataFrame
        The subframe representing the relevant data.

    Returns
    -------
    list
        A list of GraphicFeature objects representing the extracted features.
    """

    features = []
    colors_native = {"SIGNAL":'palevioletred', "PROPEP":'peru', "TRANSIT":'darkcyan'}

    for feature in record.features:
        if feature.type in ["SIGNAL", "PROPEP", "TRANSIT"]:
            start = feature.location.start.position
            end = feature.location.end.position

            if isinstance(start, int) and isinstance(end, int):
                gf = GraphicFeature(
                    start=feature.location.start.position,
                    end=feature.location.end.position,
                    strand=feature.strand,
                    color=colors_native[feature.type],
                    label=feature.type
                )
                features.append(gf)
    
    if isinstance(merops, pd.DataFrame):
        cleavage_sites = merops[merops["uniprot_acc"]== acc]["p1"].values
        codes = merops[merops["uniprot_acc"]== acc]["code"].values

        peptide_starts = set()
        for _, row in subframe.iterrows():
            try:
                start = int(row['start_pep'])
                peptide_starts.add(start)
            except:
                continue

        for site, code in zip(cleavage_sites, codes):
            site = int(site) + 1 #account for MEROPS p1 indexing (shifted by 1 to convert to p1')
            label = f"{code}: {str(site-1)}-{str(site)}"
            # Check if the cleavage site matches the start of a peptide
            if site in peptide_starts:
                gf = GraphicFeature(
                    start=site,
                    end=site,
                    strand=1,
                    color="seashell",
                    label=label
                )
                features.append(gf)

        return features


def get_pymol_image(acc, positions, colormap, vmin, vmax, alphafold):

    """
    Get the image of the protein structure with the significant positions highlighted.

    Parameters
    ----------
    acc : str
        The UniProt accession number.
    positions : list
        The positions to highlight.
    colormap : matplotlib.colors.Colormap
        The colormap to use for highlighting.
    vmin : float
        The minimum value of the colormap.
    vmax : float
        The maximum value of the colormap.
    alphafold : dict
        The AlphaFold data.

    Returns
    -------
    numpy.array or None
        The image data, or None if the AlphaFold model is not available.
    """
    
    if acc not in alphafold:
        logging.warning(f"Alphafold model for {acc} not available.")
        return None
    
    model_filename = f"AF-{acc}-F1-model_v4.cif.gz"
    alphafold_path = Path(alphafold_folder_name)
    model_path = alphafold_path / model_filename

    with tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as temp_cif:
        # Open the gzipped CIF file
        with gzip.open(model_path, 'rb') as f_in:
            # Write the decompressed contents to the temporary file
            temp_cif.write(f_in.read())
        
        temp_cif.flush()  # Ensure the file is written

        pymol.cmd.delete('all')
        pymol.cmd.load(temp_cif.name, acc)

    pymol.cmd.set('ray_trace_mode', '1')
    pymol.cmd.set('ray_trace_color', 'black')
    pymol.cmd.orient(acc)
    pymol.cmd.select('het_atoms', f'{acc} and hetatm')
    pymol.cmd.extract('het_atob', 'het_atoms')
    pymol.cmd.delete('het_atob')

    # Set the entire protein to grey
    pymol.cmd.color('grey', acc)

    positions = sorted(positions, key=lambda tup: tup[2])
    max_fold_change_pos = None
    max_fold_change = vmin

    for ind, pos in enumerate(positions):
        color = colormap((pos[2] - vmin) / (vmax - vmin))
        pymol.cmd.select('sel', f'resi {pos[0]}-{pos[1]} and {acc}')
        pymol.cmd.color('0x' + colors.to_hex(color)[1:], 'sel')

        if pos[2] > max_fold_change:
            max_fold_change = pos[2]
            max_fold_change_pos = pos

    # Orient the structure based on the highest fold change peptide
    if max_fold_change_pos is not None:
        pymol.cmd.orient(f'resi {max_fold_change_pos[0]}-{max_fold_change_pos[1]} and {acc}')
        pymol.cmd.zoom(acc, complete=0.65)
        
    tmp_file = tempfile.mktemp(suffix=".png")
    pymol.cmd.ray(1280, 720)
    pymol.cmd.png(tmp_file, dpi=600)
    img = plt.imread(tmp_file)
    os.remove(tmp_file)

    return img


def plot_protein_figure(pp, subframe, acc, col, merops, alphafold, level):

    """
    Plots the protein sequence and structure with the significant peptides highlighted, and saves the figure to the PDF file.
    Also adds existing protein features to the figure.

    Args:
        pp (PdfPages object): Object to which the figure is saved.
        subframe (DataFrame): DataFrame containing peptide information.
        acc (str): UniProt Accession of the protein of interest.
        col (str): Column name in the DataFrame to be considered for the plot.
        merops (str): Identifier for the MEROPS database.
        alphafold (bool): If True, use AlphaFold to create protein structure. 
        level (str): Level of detail in the plot.

    Returns:
        None. The function saves a figure to the PDF file represented by `pp`.
    """

    # get the protein record from UniProt
    record = SwissProt.read(ExPASy.get_sprot_raw(acc))

    # get the protein description and gene name
    desc = record.description.split(';')[0].split('=')[1]
    gene = record.gene_name[0]["Name"].split()[0] if record.gene_name else "Unknown"
    title = '-'.join([acc, gene, desc])

    # extract the protein features
    features = extract_protein_features(acc, record, merops, subframe)

    # get the significant peptides
    col_fold = col.replace('_', '/').replace('Ttest', 'Log2_fold_change')

    # Collapse duplicate peptides and take the mean of fold changes
    subframe = subframe.groupby(['query_sequence', 'start_pep', 'end_pep']).agg({col_fold: 'mean'}).reset_index()

    fold_changes = subframe[col_fold].values

    # Create a colormap centered at 0
    max_fold_change = np.max(np.abs(fold_changes))
    cmap = plt.get_cmap('coolwarm')
    norm = plt.Normalize(-max_fold_change, max_fold_change)

    for index, row in subframe.iterrows():
        peptide = row['query_sequence']
        start, end = row['start_pep'], row['end_pep']

        if isinstance(start, int) and isinstance(end, int): 
            color = cmap(norm(row[col_fold]))

            f = SeqFeature(FeatureLocation(start, end), type="peptide")
            f.qualifiers['name'] = [peptide + ' ' + '-'.join([str(start), str(end)])]

            gf = GraphicFeature(start=start, end=end, strand=+1, color=color, box_color=color, label=f.qualifiers['name'][0])
            features.append(gf)

    # create the figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 12))

    gr = GraphicRecord(sequence_length=record.sequence_length, features=features, first_index=1)
    gr.plot(ax=ax1)

    # add the colorbar to the sequence plot
    cax = fig.add_axes([0.25, 0.92, 0.5, 0.02])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation='horizontal', aspect=12)
    cb.set_label(col_fold, fontsize=12)

    # get the pymol image of the protein structure
    if alphafold and level=='both':
        positions = subframe[['start_pep', 'end_pep', col_fold]].values.tolist()
        img = get_pymol_image(acc, positions, cmap, -max_fold_change, max_fold_change, alphafold)

        if img is not None:
            ax2.imshow(img)
    
    ax2.axis('off')
    plt.suptitle(title, x=0.02, fontsize=14, horizontalalignment='left')
    plt.subplots_adjust()
    plt.savefig(pp, format="pdf")
    plt.close()


def plot_enrichment_figure(genes, col_name, organism='hsapiens'):

    """
    Perform functional enrichment analysis with gprofiler and plot the results.

    Args:
        genes (list): List of gene names to be considered for enrichment analysis.
        col_name (str): Column name to be displayed in the heatmap title.
        organism (str): Organism to be considered for enrichment analysis. Default is 'hsapiens' (Homo sapiens).

    Returns:
        fig (Figure object): The generated figure, or None if no significant results were found.
    """

    # Initialize gprofiler
    gp = GProfiler(return_dataframe=True)

    # Perform enrichment analysis
    enrichment_results = gp.profile(organism=organism, query=genes)

    # Filter significant results
    significant_results = enrichment_results[enrichment_results['p_value'] < 0.05]

    if len(significant_results) == 0:
        return None

    # Pivot data for heatmap
    heatmap_data = significant_results.pivot(index='name', columns='source', values='p_value')
    heatmap_data = -np.log10(heatmap_data)

    # Plot heatmap
    plt.figure(figsize=(10, 10))
    ax = sns.heatmap(heatmap_data, cmap="Reds", square=True, cbar_kws={'label': '-log10(p-value)'})
    ax.tick_params(axis='both', which='major', labelsize=8)
    ax.set_title(col_name)
    plt.subplots_adjust(left=0.2, bottom=0.2)
    fig = ax.get_figure()
    plt.close()

    return fig


def plot_pathway_figures(subframe, col_fold, pp, cutoff=0.05):

    """
    Perform pathway enrichment analysis with Reactome and plot the results.

    Args:
        subframe (DataFrame): DataFrame containing peptide information.
        col_fold (str): Column name representing fold changes in the DataFrame.
        pp (PdfPages object): Object to which the figure is saved.
        cutoff (float): P-value cutoff for pathway enrichment. Default is 0.05.

    Returns:
        None. The function saves a figure to the PDF file represented by `pp`.
    """

    # Get the list of proteins
    accessions = subframe["query_accession"].unique()

    # Get the enriched pathways
    enriched_pathways = annutils.get_enriched_pathways(accessions, cutoff=cutoff)
    pathways = enriched_pathways["pathways"]

    # Get the significant pathway identifiers
    significant_pathways_stIDs = [p["stId"] for p in pathways]
    logging.info(f"Significant pathways for {col_fold}: {significant_pathways_stIDs}")
    
    # Plot the pathways
    for i, pathway_stId in tqdm(enumerate(significant_pathways_stIDs)):

        # Get the pathway proteins and interaction map
        pathway_proteins, interaction_map = annutils.get_proteins_and_interactors(pathway_stId)

        if len(interaction_map) > 0:
            # Get the protein edgelist, cleavage edgelist and cleavages
            protein_edgelist, cleavage_edgelist, cleavages = annutils.construct_edgelists(subframe, interaction_map)

            # Construct the network
            network = annutils.construct_network(protein_edgelist, cleavage_edgelist, pathway_proteins)

            # Create labels and colors for proteins and cleavages
            protein_labels, protein_colors, cleavage_labels, cleavage_colors, protein_cmap, peptide_cmap = create_labels_colors(pathway_proteins, 
                                                                                                                                cleavages, subframe, accessions, col_fold)

            # Plot the network
            fig = visualize_network(network, pathway_proteins, protein_edgelist, protein_colors, protein_labels, protein_cmap, cleavages, cleavage_edgelist,
                                    cleavage_colors, cleavage_labels, peptide_cmap, col_fold, subframe, pathway_stId, pathways, i)
            fig.savefig(pp, format="pdf")
            plt.close()

        else:
            logging.warning(f"No interaction map found for pathway: {pathway_stId}")

        # save the network in a cyjs format


def create_labels_colors(proteins, cleavages, subframe, accs, col_fold):

    """
    Create labels and colors for proteins and cleavages.

    Args:
        proteins (list): List of protein identifiers.
        cleavages (list): List of cleavage site identifiers.
        subframe (DataFrame): DataFrame containing peptide information.
        accs (list): List of protein accessions considered in the study.
        col_fold (str): Column name representing fold changes in the DataFrame.

    Returns:
        Tuple containing label and color information for proteins and cleavages.
    """

    # Create a colormap for proteins and peptides
    protein_cmap = plt.cm.PuOr
    peptide_cmap = plt.cm.RdYlGn_r

    # Create a color normalizer for proteins and peptides
    protein_norm = colors.Normalize(vmin=subframe.groupby('query_accession')[col_fold].sum().min(), vmax=subframe.groupby('query_accession')[col_fold].sum().max())
    peptide_norm = colors.Normalize(vmin=subframe[col_fold].min(), vmax=subframe[col_fold].max())

    protein_labels = {}
    protein_colors = []
    cleavage_labels = {}
    cleavage_colors = []

    for prot in proteins:
        # if protein in starting accessions, color it by summed normalized intensity of its peptides
        if prot in accs:
            summed_intensity = subframe.loc[subframe['query_accession'] == prot, col_fold].sum()
            protein_colors.append(protein_cmap(protein_norm(summed_intensity)))
            protein_labels[prot] = prot
        # if protein in the pathway but not detected, color it grey
        else:
            protein_colors.append(colors.to_rgba_array('#BDB9B8'))
            protein_labels[prot] = prot

    for cleavage in cleavages:
        # Split the cleavage string to get protein, position and sequence
        prot_pos, sequence = cleavage.split(';')
        protein, position = prot_pos.split(':')

        # Query the subframe for the matching sequence
        subframe_cleavage = subframe.loc[subframe['query_sequence'] == sequence]

        if len(subframe_cleavage) > 0:
            # Calculate the average fold change for this cleavage
            avg_fold_change = subframe_cleavage[col_fold].mean()

            # Add the color and label for this cleavage
            cleavage_colors.append(peptide_cmap(peptide_norm(avg_fold_change)))
            cleavage_labels[cleavage] = cleavage

    return protein_labels, protein_colors, cleavage_labels, cleavage_colors, protein_cmap, peptide_cmap


def visualize_network(network, proteins, protein_edgelist, protein_colors, protein_labels, protein_cmap, cleavages, cleavage_edgelist, cleavage_colors, cleavage_labels, peptide_cmap, col_fold, subframe, pathway_stId, pathways, i):
    
    """
    Visualize the protein-cleavage network using matplotlib and networkx. 
    The network is represented as a graph with two types of nodes (proteins and cleavages) 
    and edges between them. 

    Args:
        network (networkx.classes.graph.Graph): Graph object representing the protein-cleavage network.
        proteins (list): List of protein identifiers.
        protein_edgelist (list): List of tuples representing edges between proteins.
        protein_colors (list): List of colors for the protein nodes.
        protein_labels (dict): Dictionary mapping proteins to their labels.
        protein_cmap (matplotlib.colors.Colormap): Colormap for proteins.
        cleavages (list): List of cleavage site identifiers.
        cleavage_edgelist (list): List of tuples representing edges between cleavage sites.
        cleavage_colors (list): List of colors for the cleavage nodes.
        cleavage_labels (dict): Dictionary mapping cleavages to their labels.
        peptide_cmap (matplotlib.colors.Colormap): Colormap for cleavages.
        col_fold (str): Column name representing fold changes in the DataFrame.
        subframe (DataFrame): DataFrame containing peptide information.
        pathway_stId (str): ID of the pathway to be visualized.
        pathways (list): List of pathways.
        i (int): Index of the current pathway in the `pathways` list.

    Returns:
        fig (matplotlib.figure.Figure): The generated figure.
    """

    # Define positions for the nodes in the network
    # Try to use the Graphviz layout, otherwise default to the spring layout
    try:
        import pygraphviz
        pos = nx.nx_agraph.graphviz_layout(network, prog='neato')
    except ImportError:
        logging.warning("pygraphviz not found, defaulting to spring layout")
        pos = nx.spring_layout(network)
    
    fig, ax = plt.subplots(3, 1, gridspec_kw={'height_ratios': [25, 3, 3]}, figsize=(20,20))
    
    all_edges = protein_edgelist + cleavage_edgelist
    all_labels = protein_labels | cleavage_labels

    prot_nodes = nx.draw_networkx_nodes(network, pos, nodelist=proteins, node_color=protein_colors, node_size=600, node_shape='o', ax=ax[0])
    cleavage_nodes = nx.draw_networkx_nodes(network, pos, nodelist=cleavages, node_color=cleavage_colors, node_size=600, node_shape='D', ax=ax[0])
    edges = nx.draw_networkx_edges(network, pos, edgelist=all_edges, arrowstyle="-", arrowsize=25, width=1, ax=ax[0])
    labels = nx.draw_networkx_labels(network, pos, all_labels, font_size=7, font_color="black", ax=ax[0])

    ax[0].axis('off')

    # Protein color bar
    pc_protein = collections.PatchCollection(edges, cmap=protein_cmap)
    pc_protein.set_array(np.array(subframe.groupby('query_accession')[col_fold].sum()))
    cb_protein = plt.colorbar(pc_protein, location='bottom', ax=ax[1], orientation='horizontal')
    cb_protein.ax.xaxis.set_ticks_position('top')
    cb_protein.ax.xaxis.set_label_position('top')
    cb_protein.set_label(f'Protein {col_fold}', fontsize=14)
    ax[1].axis('off') # turn off axis

    # Cleavage color bar
    pc_cleavage = collections.PatchCollection(edges, cmap=peptide_cmap)
    pc_cleavage.set_array(np.array(subframe[col_fold]))
    cb_cleavage = plt.colorbar(pc_cleavage, location='bottom', ax=ax[2], orientation='horizontal')
    cb_cleavage.set_label(f'Peptide {col_fold}', fontsize=14)
    ax[2].axis('off') # turn off axis

    plt.suptitle(f"{col_fold}\n{pathway_stId}: {pathways[i]['name']}\np-value: {str(pathways[i]['entities']['pValue'])}", fontsize=20)

    return fig
