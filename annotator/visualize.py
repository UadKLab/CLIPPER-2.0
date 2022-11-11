import os

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from tqdm import tqdm


class Visualizer:
    def __init__(self, df, annot, conditions, software, patterns):

        self.df = df
        self.annot = annot
        self.conditions = conditions
        self.software = software
        self.patterns = patterns

    def general(self):

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
            quant_columns = self.df.columns[
                self.df.columns.str.contains(quant_columns_pat)
            ]

            if len(quant_columns) > 0:
                subframe = self.df[[seq_col, mod_col]]
                if self.software == "pd":
                    unique_mod = subframe
                    unique_mod.loc[:, mod_col] = unique_mod.loc[:, mod_col].astype(str)
                elif self.software == "sm":
                    unique_mod = subframe.drop_duplicates(subset=[seq_col, mod_col])
                quant_frame = self.df.dropna(subset=quant_columns, how="all")

                data["modified peptides"] = unique_mod[seq_col].count()
                data["labelled peptides"] = unique_mod[
                    unique_mod[mod_col].str.contains(label_pat)
                ][seq_col].count()
                data["nterm"] = unique_mod[unique_mod[mod_col].str.contains(nter_pat)][
                    seq_col
                ].count()
                data["quant proteins"] = quant_frame[acc_col].drop_duplicates().count()
                data["quant peptides"] = quant_frame[seq_col].drop_duplicates().count()
                data["quant mod peptides"] = quant_frame.drop_duplicates(
                    [seq_col, mod_col]
                )[seq_col].count()

        ticks_no = list(range(len(data)))
        stat_names = [k for k in data]
        stat_vals = [data[k] for k in data]

        fig, ax = plt.subplots(1, 1, figsize=(12, 8))

        sns.barplot(x=ticks_no, y=stat_vals, palette="spring")
        for i in range(len(stat_vals)):
            x_pos = ticks_no[i]
            y_pos = stat_vals[i] / 2
            ax.text(
                x_pos,
                y_pos,
                stat_vals[i],
                fontfamily="sans-serif",
                ha="center",
                va="center",
                fontsize=12,
            )

        ax.set_yticks([int(i) for i in ax.get_yticks()])
        ax.set_xticklabels(
            stat_names, fontsize=14, fontfamily="sans-serif", rotation=90
        )
        ax.set_yticklabels(
            [int(i) for i in ax.get_yticks()], fontsize=14, fontfamily="sans-serif"
        )

        sns.despine(fig=fig, ax=ax)
        plt.subplots_adjust(bottom=0.3)

        fig = plt.gcf()
        plt.close()

        return {"general": fig}

    def volcano(self):

        columns_fold = self.annot.columns[
            self.annot.columns.str.contains("Log2_fold_change:")
        ]
        columns_ttest = self.annot.columns[
            self.annot.columns.str.contains("Log10_ttest:")
        ]
        figures = {}

        for test in columns_ttest:
            conditions = test.split()[1].split("_")

            for fold in columns_fold:

                if conditions[0] in fold and conditions[1] in fold:
                    frame = self.annot.loc[:, [fold, test]].dropna(how="any")
                    # log values are negative for ttest column, inverse
                    frame[test] = frame[test] * -1
                    frame["coding"] = np.where(
                        (~frame[fold].between(-1.5, 1.5))
                        & (~frame[test].between(-1.5, 1.5)),
                        "1",
                        "0",
                    )
                    frame = frame.sort_values("coding")

                    if len(frame.coding.unique()) == 1:
                        colors = ["grey"]
                    else:
                        colors = ["grey", "crimson"]

                    g = sns.scatterplot(
                        data=frame, x=fold, y=test, hue="coding", palette=colors
                    )
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

        colors = sns.color_palette("husl", len(self.conditions))

        columns = self.annot.columns[self.annot.columns.str.contains("_CV")]

        if len(columns) > 0:
            fig, ax = plt.subplots()

            for i in range(len(columns)):
                g = sns.kdeplot(
                    self.annot[columns[i]],
                    color=colors[i],
                    ax=ax,
                    label=columns[i][:-3],
                )

            plt.legend()
            plt.xlabel("")
            plt.close()

            return {"cv_plot": g}

    def fold_plot(self):

        columns = self.annot.columns[
            self.annot.columns.str.contains("Log2_fold_change:")
        ]
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

        columns = self.annot.columns[
            self.annot.columns.str.contains("Log2_fold_change:")
        ]
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
        """Heatmap of all quantification columns"""

        columns = self.df.columns[self.df.columns.str.contains(self.patterns['quant'])]

        if len(columns) > 0:
            heat = self.df[columns]
            heat = np.log10(heat.replace(0, np.nan).dropna())
            fig = sns.heatmap(heat, yticklabels=False)
            plt.subplots_adjust(bottom=0.4)
            plt.close()

            return {"heatmap": fig}

    def clustermap(self):
        """Heatmap of all quantification columns"""

        columns = self.df.columns[self.df.columns.str.contains(self.patterns['quant'])]

        if len(columns) > 0:
            cluster = self.df[columns]
            cluster = np.log10(cluster.replace(0, np.nan).dropna())
            fig = sns.clustermap(cluster, yticklabels=False)
            plt.close()

            return {"clustermap": fig.fig}

    def piechart(self):

        figures = {}

        check = False
        check_col = False
        column = self.annot["nterm_annot"].dropna()

        categories = {
            "internal": (column.values == "Internal").sum(),
            "met removed": (column.values == "Met removed").sum(),
            "met intact": (column.values == "Met intact").sum(),
            "signal removed": (column.values == "Signal removed").sum(),
            "propeptide removed": (column.values == "Propeptide removed").sum(),
            "within signal": (
                column.values == "Cleavage within signal peptide range"
            ).sum(),
            "within propeptide": (
                column.values == "Cleavage within propeptide range"
            ).sum(),
            "transit removed": (
                column.values == "Transit peptide removed"
                ).sum(),
            "within transit peptide": (
                column.values == "Cleavage within transit peptide range"
                ).sum(),
        }

        y = np.array([categories.pop("internal"), sum(categories.values())])
        for val in y:
            if isinstance(val, np.int32) and val > 0:
                check = True

        if check:
            x = ["Internal", "Natural"]
            explode = [0.1, 0]
            colors = [
                "darkgreen",
                "peru",
            ]
            porcent = 100.0 * y / y.sum()
            patches, texts = plt.pie(
                y, explode=explode, colors=colors, startangle=90, normalize=True
            )
            labels = ["{0} - {1:1.2f} %".format(i, j) for i, j in zip(x, porcent)]
            patches, labels, dummy = zip(
                *sorted(zip(patches, labels, y), key=lambda x: x[2], reverse=True)
            )
            plt.legend(
                patches,
                labels,
                loc="center left",
                bbox_to_anchor=(-0.38, 0.5),
                fontsize=8,
            )
            plt.subplots_adjust(left=0.3)
            figures["type"] = plt.gcf()
            plt.close()
            check = False

        y = np.array(
            [
                categories.pop("met removed"),
                categories.pop("met intact"),
                sum(categories.values()),
            ]
        )
        for val in y:
            if isinstance(val, np.int32) and val > 0:
                check = True

        if check:
            x = ["Met removal", "Met intact", "Other"]
            colors = ["crimson", "lightpink", "seagreen"]
            explode = [0, 0, 0.1]
            porcent = 100.0 * y / y.sum()
            patches, texts = plt.pie(
                y, explode=explode, colors=colors, startangle=90, normalize=True
            )
            labels = ["{0} - {1:1.2f} %".format(i, j) for i, j in zip(x, porcent)]
            patches, labels, dummy = zip(
                *sorted(zip(patches, labels, y), key=lambda x: x[2], reverse=True)
            )
            plt.legend(
                patches,
                labels,
                loc="center left",
                bbox_to_anchor=(-0.45, 0.5),
                fontsize=8,
            )
            plt.subplots_adjust(left=0.3)
            figures["natural"] = plt.gcf()
            plt.close()
            check = False

        y = np.array([categories[k] for k in categories])
        for val in y:
            if isinstance(val, np.int32) and val > 0:
                check = True

        if check:
            x = [k for k in categories]
            colors = ["steelblue", "goldenrod", "slategrey", "forestgreen", "violet", "tomato"]
            explode = [0, 0.1, 0, 0.1, 0, 0]
            porcent = 100.0 * y / y.sum()
            patches, texts = plt.pie(
                y, explode=explode, colors=colors, startangle=90, normalize=True
            )
            labels = ["{0} - {1:1.2f} %".format(i, j) for i, j in zip(x, porcent)]
            patches, labels, dummy = zip(
                *sorted(zip(patches, labels, y), key=lambda x: x[2], reverse=True)
            )
            plt.legend(
                patches,
                labels,
                loc="center left",
                bbox_to_anchor=(-0.52, 0.5),
                fontsize=8,
            )
            plt.subplots_adjust(left=0.3)
            figures["other"] = plt.gcf()
            plt.close()
            check = False

        if self.software == "sm":
            if self.patterns['mod'] in self.df.columns:
                term = 'Labelled'
                labelled = len(
                    self.df[
                        self.df[self.patterns['mod']].str.contains(self.patterns['nterm_label'], na=False)
                    ]
                )
                acetyl = len(
                    self.df[
                        self.df[self.patterns['mod']].str.contains(
                            r"Acetyl \(Protein N-term\)", na=False
                        )
                    ]
                )
                check_col = True
        
        elif self.software == "pd":
            if self.patterns['mod'] in self.df.columns:

                acetyl = len(
                    self.df[
                        self.df[self.patterns['mod']].str.contains(
                            r"Acetyl \[N-Term\]", na=False
                        )
                    ]
                )
                
                if self.patterns['quant']:
                    term = 'Labelled'
                    labelled = len(
                        self.df[
                            self.df[self.patterns['mod']].str.contains(
                                self.patterns['nterm_label'], na=False
                            )
                        ]
                    )
                else:
                    term = 'Non acetylated'
                    total = len(
                        self.df[
                            self.df[self.patterns['mod']].str.contains(
                                self.patterns['nterm'], na=False
                            )
                        ]
                    )
                    labelled = total - acetylated

                check_col = True

        if check_col:
            categories = {"Acetylated": acetyl, term: labelled}

            y = np.array([categories[k] for k in categories])
            for val in y:
                if isinstance(val, np.int32) and val > 0:
                    check = True

            if check:
                x = [k for k in categories]
                colors = ["teal", "saddlebrown"]
                explode = [0, 0.1]
                porcent = 100.0 * y / y.sum()
                patches, texts = plt.pie(
                    y, explode=explode, colors=colors, startangle=90, normalize=True
                )
                labels = ["{0} - {1:1.2f} %".format(i, j) for i, j in zip(x, porcent)]
                patches, labels, dummy = zip(
                    *sorted(zip(patches, labels, y), key=lambda x: x[2], reverse=True)
                )
                plt.legend(
                    patches,
                    labels,
                    loc="center left",
                    bbox_to_anchor=(-0.4, 0.5),
                    fontsize=8,
                )
                plt.subplots_adjust(left=0.3)
                figures["modification"] = plt.gcf()
                plt.close()
                check = False
                check_col = False

        if self.patterns['mod'] in self.df.columns:
            cysteines = len(
                self.df[self.df[self.patterns['mod']].str.contains(r"C", na=False)]
            )
            alkylated = len(
                self.df[
                    self.df[self.patterns['mod']].str.contains(
                        r"Carbamidomethyl", na=False
                    )
                ]
            )
            check_col = True

        if check_col:
            if cysteines - alkylated > 0: 
                categories = {
                    "Free cysteines": cysteines - alkylated,
                    "Carbamidomethyl": alkylated,
                }

                y = np.array([categories[k] for k in categories])
                for val in y:
                    if isinstance(val, np.int32) and val > 0:
                        check = True

                if check:
                    x = [k for k in categories]
                    colors = ["indianred", "cadetblue"]
                    explode = [0, 0.1]
                    porcent = 100.0 * y / y.sum()
                    patches, texts = plt.pie(
                        y, explode=explode, colors=colors, startangle=90, normalize=True
                    )
                    labels = ["{0} - {1:1.2f} %".format(i, j) for i, j in zip(x, porcent)]
                    patches, labels, dummy = zip(
                        *sorted(zip(patches, labels, y), key=lambda x: x[2], reverse=True)
                    )
                    plt.legend(
                        patches,
                        labels,
                        loc="center left",
                        bbox_to_anchor=(-0.45, 0.5),
                        fontsize=8,
                    )
                    plt.subplots_adjust(left=0.3)
                    figures["cysteine"] = plt.gcf()
                    plt.close()
                    check = False
                    check_col = False

        if self.patterns['mod'] in self.df.columns:
            lysines = len(
                self.df[self.df[self.patterns['mod']].str.contains(r"K", na=False)]
            )
            labelled = len(
                self.df[
                    self.df[self.patterns['mod']].str.contains(
                        self.patterns['lysine_label'], na=False
                    )
                ]
            )
            check_col = True

        if check_col:
            if lysines - labelled > 0:
                categories = {
                    "Free lysines": lysines - labelled,
                    "Labelled lysines": labelled,
                }

                y = np.array([categories[k] for k in categories])
                for val in y:
                    if isinstance(val, np.int32) and val > 0:
                        check = True

                if check:
                    x = [k for k in categories]
                    colors = ["violet", "olivedrab"]
                    explode = [0, 0.1]
                    porcent = 100.0 * y / y.sum()
                    patches, texts = plt.pie(
                        y, explode=explode, colors=colors, startangle=90, normalize=True
                    )
                    labels = ["{0} - {1:1.2f} %".format(i, j) for i, j in zip(x, porcent)]
                    patches, labels, dummy = zip(
                        *sorted(zip(patches, labels, y), key=lambda x: x[2], reverse=True)
                    )
                    plt.legend(
                        patches,
                        labels,
                        loc="center left",
                        bbox_to_anchor=(-0.4, 0.5),
                        fontsize=8,
                    )
                    plt.subplots_adjust(left=0.3)
                    figures["lysine"] = plt.gcf()
                    plt.close()
                    check = False
                    check_col = False

        if self.patterns['mod'] in self.df.columns:
            total = len(self.df[self.patterns['mod']])
            nterm = len(
                self.df[
                    self.df[self.patterns['mod']].str.contains(self.patterns['nterm'], na=False)
                ]
            )
            check_col = True

        if check_col:
            categories = {"Tryptic": total - nterm, "N-termini": nterm}

            y = np.array([categories[k] for k in categories])
            for val in y:
                if isinstance(val, np.int32) and val > 0:
                    check = True

            if check:
                x = [k for k in categories]
                colors = ["salmon", "darkcyan"]
                explode = [0, 0.1]
                porcent = 100.0 * y / y.sum()
                patches, texts = plt.pie(
                    y, explode=explode, colors=colors, startangle=90, normalize=True
                )
                labels = ["{0} - {1:1.2f} %".format(i, j) for i, j in zip(x, porcent)]
                patches, labels, dummy = zip(
                    *sorted(zip(patches, labels, y), key=lambda x: x[2], reverse=True)
                )
                plt.legend(
                    patches,
                    labels,
                    loc="center left",
                    bbox_to_anchor=(-0.4, 0.5),
                    fontsize=8,
                )
                plt.subplots_adjust(left=0.3)
                figures["nterm"] = plt.gcf()
                plt.close()
                check = False
                check_col = False

        if "exopeptidase" in self.annot.columns:
            exopep = (~self.annot["exopeptidase"].isna()).sum()
            noexo = self.annot["exopeptidase"].isna().sum()
            check_col = True

        if check_col:
            categories = {"Exopeptidase": exopep, "No activity": noexo}
            for val in y:
                if isinstance(val, np.int32) and val > 0:
                    check = True

            if check:
                y = np.array([categories[k] for k in categories])

                x = [k for k in categories]
                colors = ["darkolivegreen", "mediumslateblue"]
                explode = [0.1, 0]
                porcent = 100.0 * y / y.sum()
                patches, texts = plt.pie(
                    y, explode=explode, colors=colors, startangle=90, normalize=True
                )
                labels = ["{0} - {1:1.2f} %".format(i, j) for i, j in zip(x, porcent)]
                patches, labels, dummy = zip(
                    *sorted(zip(patches, labels, y), key=lambda x: x[2], reverse=True)
                )
                plt.legend(
                    patches,
                    labels,
                    loc="center left",
                    bbox_to_anchor=(-0.45, 0.5),
                    fontsize=8,
                )
                plt.subplots_adjust(left=0.3)
                figures["exopep"] = plt.gcf()
                plt.close()
                check = False
                check_col = False

        return figures

    def gallery(self, stat=False, cutoff=0.05, folder=None):

        if stat:
            cols = self.annot.columns[self.annot.columns.str.contains("Ttest:")]
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

        accs = self.annot.loc[indices, "query_accession"].unique()

        sns.set("paper", "ticks")
        pp = PdfPages(os.path.join(folder, "Sig_gallery.pdf"))

        for acc in tqdm(accs):

            natural = self.annot[
                (self.annot["query_accession"] == acc)
                & (
                    self.annot["nterm_annot"].isin(
                        [
                            "Met removed",
                            "Met intact",
                            "Signal removed",
                            "Propeptide removed",
                        ]
                    )
                )
            ]

            internal = self.annot[
                (self.annot["query_accession"] == acc)
                & (self.annot["nterm_annot"] == "Internal")
            ]

            columns = natural.columns[natural.columns.str.contains("_mean")]
            natural = natural.drop_duplicates("query_sequence").dropna(subset=columns)
            internal = internal.drop_duplicates("query_sequence").dropna(subset=columns)

            if len(natural) > 0 and len(internal) > 0:

                fig, ax = plt.subplots(1, 2, figsize=(12, 10))

                type_term = ["natural", "internal"]
                data = {}
                for col in columns:
                    for t, subframe in zip(type_term, [natural, internal]):
                        data[f"{col}_{t}"] = subframe[col].mean()

                values = [v for v in data.values()]
                names = [k for k in data.keys()]

                sns.barplot(x=names, y=values, palette="summer", ax=ax[0])
                ax[0].set_xticklabels(
                    names, fontsize=12, fontfamily="sans-serif", rotation=90
                )
                ax[0].set_title("Mean values per condition", fontsize=11)

                data = {}
                for t, subframe in zip(type_term, [natural, internal]):
                    for row in subframe.index:
                        for col in columns:
                            data[
                                f"{col[:-5]}_{t}_{subframe.loc[row, 'query_sequence']}"
                            ] = subframe.loc[row, col]

                values = [v for v in data.values()]
                names = [k for k in data.keys()]

                sns.barplot(x=names, y=values, palette="winter", ax=ax[1])
                ax[1].set_xticklabels(
                    names, fontsize=12, fontfamily="sans-serif", rotation=90
                )
                ax[1].set_title("Peptides with quant values", fontsize=11)

                plt.suptitle(f"Protein Uniprot ID: {acc}")
                plt.tight_layout()
                plt.savefig(pp, format="pdf")
                plt.close()

        pp.close()
