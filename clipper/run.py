import os
import time
import logging
import pandas as pd
from datetime import timedelta

from bin.clipper import Clipper
from bin.annutils import initialize, initialize_arguments, write_terminal_headers


def main(args=None):
    """Runs the CLIPPER annotation pipeline.

    Args:
        args (dict): Dictionary containing arguments.

    Returns:
        None
    """

    start_time = time.perf_counter()

    # If arguments are not provided, use the default values from initialize()
    if args is None:
        args = initialize_arguments()

    args = initialize(args)

    # Create an Annotator object with the provided arguments
    annotator = Clipper(args)

    # Prepare the input data for annotation
    annotator.prepare()

    # Check if the input file has already been annotated
    #args = annotator.infer_infile_annotation_status(args)

    # Perform peptide annotation
    logging.debug("Starting annotation of peptides...")
    write_terminal_headers("PEPTIDE ANNOTATION")
    
    #print(args['infileannotation'])
    #if any(value == 1 for value in args['infileannotation'].values()):
    #    logging.info(f"\nprevious annotation of the input file was detected, and {', '.join(args['infileannotation'].keys())} annotation will not be performed.\n")
    
    # load annotation dataframe, and verify, if flag is given, that previous annotation is valid. If not change argument and do annotation anyway
    annotator.initialize_annotation()


    if annotator.nomerops is False:
        annotator.read_MEROPS()
    
    # do annotation if previous annotation has not been done or is not valid
    # print(annotator.preannotated)
    # if annotator.preannotated == False:
    #     print("IN PREANNOTATED == FALSE")
    if not os.path.isfile(os.path.join(annotator.outfolder, 'annot_initial.xlsx')):
        # annotate from Uniprot
        annotator.threaded_annotate(args['threadingcores'])
        # writing annotated datafile to excel

        annotator.annot.to_excel(os.path.join(annotator.outfolder, 'annot_initial.xlsx'))
        logging.info("Finished annotation of peptides.\n")
    else:
        logging.info("Previous annotation of peptides detected, and annotation will be loaded from file.")
        annotator.annot = pd.read_excel(os.path.join(annotator.outfolder, 'annot_initial.xlsx'), index_col=0)
        logging.info("Loaded annotation from file.\n")

    # Annotate Protein Atlas data
    logging.info("Starting annotation of Protein Atlas data...")
    annotator.annotate_protein_atlas()
    logging.info("Finished annotation of Protein Atlas data.\n")

    # Perform proteoform check
    logging.info("Starting proteoform check...")
    annotator.proteoform_check()
    logging.info("Finished proteoform check.")

    # Perform exopeptidase check if specified
    if not args["noexo"]:
        logging.info("Starting exopeptidase activity check...")
        annotator.exopeptidase()
        logging.info("Finished exopeptidase activity check.")

    write_terminal_headers("Statistical analysis")

    # Perform condition-specific statistics annotation if condition file is specified
    if args["conditionfile"]:
        # Read the condition file
        annotator.read_condition_file()
        logging.info("Parsed condition file")

        # Perform general statistics annotation
        logging.info("Performing general condition statistics annotation...")
        annotator.general_conditions()
        logging.info("Finished general condition statistics annotation.")

        # Perform pairwise or all-vs-all statistical testing if specified
        if args["stat"]:
            logging.info("Performing statistical testing...")
            annotator.condition_statistics()
            logging.info("Finished statistical testing.")
            
            annotator.correct_multiple_testing()

        # Perform fold distribution check if specified
        if args["significance"]:
            logging.info("Checking fold distribution...")
            annotator.percentile_fold(0.05)
            logging.info("Finished fold distribution check.\n")

    if args["calcstructure"]:
        logging.info("Computing structural properties...")
        annotator.annotate_structure()
        logging.info("Finished computing structural properties.")

    if args["proteasefile"]:
        logging.info("Predicting protease activity...")
        annotator.predict_protease_activity()
        logging.info("Finished predicting protease activity.")

    # Generate figures if specified
    if args["visualize"]:
        write_terminal_headers('Creating visualizations')
        logging.info("Generating figures...")
        annotator.visualize()
        logging.info("Finished generating figures.")

    # Generate logos if specified
    if args["logo"] is not None:
        logging.info("Generating logos...")
        annotator.create_logos()
        logging.info("Finished generating logos.")

    # Write output files
    write_terminal_headers('Finishing up')
    logging.info("Writing output files...")
    annotator.write_files()
    logging.info("Finished writing output files")

    # Calculate and log the total running time
    end_time = time.perf_counter()
    total_time = str(timedelta(seconds=end_time - start_time))
    logging.info(f"\nCLIPPER 2.0 pipeline completed in {total_time}.\n")


if __name__ == "__main__":
    main()