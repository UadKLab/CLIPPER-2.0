import os
from ast import literal_eval
from annutils import calculate_structure_properties
from argparse import ArgumentParser, HelpFormatter

parser = ArgumentParser(
    prog="test script",
    description="this is a script",
    formatter_class=HelpFormatter,
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

parser.add_argument(
    "-af",
    "--alphafolder",
    default=None,
    action="store",
    dest="af",
    type=str,
    required=True,
    help="Path to alphafold folder with CIF structures",
)

parser.add_argument(
    "-csi",
    "--cleavagesiteindices",
    default=None,
    action="store",
    dest="csi",
    type=str,
    required=True,
    help="Indices for cleavage sites / accession for the current protein being processed",
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

args = parser.parse_args()


if __name__ == "__main__":

    tf = args.tf
    af = args.af
    csi = literal_eval(args.csi)
    acc = args.acc

    with open(tf, 'r') as f:
        lines = f.readlines()[0]
        structure_properties = literal_eval(lines)

    structure_properties = calculate_structure_properties(acc, cleavage_sites_indices=csi, structure_properties=structure_properties, alphafold_folder=af)

    with open(tf, 'w') as f:
        f.write(str(structure_properties))


