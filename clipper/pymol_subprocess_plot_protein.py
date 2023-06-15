import os
from ast import literal_eval
import matplotlib.pyplot as plt
from visualize import get_pymol_image
from argparse import ArgumentParser, HelpFormatter

parser = ArgumentParser(
    prog="test script",
    description="this is a script",
    formatter_class=HelpFormatter,
)

parser.add_argument(
    "-acc",
    "--accession",
    default=None,
    action="store",
    dest="acc",
    type=str,
    required=True,
    help="Accession for the current protein being processed",
)

parser.add_argument(
    "-pos",
    "--position",
    default=None,
    action="store",
    dest="pos",
    type=str,
    required=True,
    help="Position of peptide in protein",
)

parser.add_argument(
    "-_mfc",
    "--_maxfoldchange",
    default=None,
    action="store",
    dest="_mfc",
    type=float,
    required=True,
    help="Negative maximum fold change of all peptides for a protein, used to calibrate colormap",
)

parser.add_argument(
    "-mfc",
    "--maxfoldchange",
    default=None,
    action="store",
    dest="mfc",
    type=float,
    required=True,
    help="Maximum fold change of all peptides for a protein, used to calibrate colormap",
)

parser.add_argument(
    "-tf",
    "--tempfile",
    default=None,
    action="store",
    dest="tf",
    type=str,
    required=True,
    help="Path to write temp file for delivering structure properties to clipper",
)

args = parser.parse_args()


if __name__ == "__main__":

    acc = args.acc
    positions = literal_eval(args.pos)
    vmin = args._mfc
    vmax = args.mfc
    tf = args.tf
    cmap = plt.get_cmap('coolwarm')

    get_pymol_image(acc, positions, cmap, vmin, vmax, peptide_protein_plot_path = tf)

