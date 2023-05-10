import time
import logging
from datetime import timedelta

from annotator import Annotator
from annutils import initialize


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
        args = initialize()

    # Create an Annotator object with the provided arguments
    annotator = Annotator(args)

    # Prepare the input data for annotation
    annotator.prepare()

    # Perform peptide annotation
    logging.info("Starting annotation of peptides...")
    if not args["singlecpu"]:
        annotator.threaded_annotate()
    else:
        annotator.annotate()
    logging.info("Finished annotation of peptides.\n")

    # Perform proteoform check
    logging.info("Starting proteoform check...")
    annotator.proteoform_check()
    logging.info("Finished proteoform check.")

    # Perform exopeptidase check if specified
    if not args["noexo"]:
        logging.info("Starting exopeptidase activity check...")
        annotator.exopeptidase()
        logging.info("Finished exopeptidase activity check.\n")

    # Perform general statistics annotation if condition file is not specified
    if args["conditionfile"] is None:
        logging.info("Performing general statistics annotation...")
        annotator.general_conditions()
        logging.info("Finished general statistics annotation.\n")

    # Perform condition-specific statistics annotation if condition file is specified
    else:
        # Read the condition file
        annotator.read_condition_file()
        logging.info("Parsed condition file")

        # Perform general statistics annotation
        logging.info("Performing general statistics annotation...")
        annotator.general_conditions()
        logging.info("Finished general statistics annotation.")

        # Perform pairwise or all-vs-all statistical testing if specified
        if args["stat"]:
            logging.info("Performing statistical testing...")
            annotator.condition_statistics()
            logging.info("Finished statistical testing.")

            annotator.correct_multiple_testing()
            logging.info("Finished multiple testing correction.\n")

        # Perform fold distribution check if specified
        if args["significance"]:
            logging.info("Checking fold distribution...")
            annotator.percentile_fold(0.05)
            logging.info("Finished fold distribution check.\n")

    if args["calcstructure"]:
        logging.info("Computing structural properties...")
        annotator.annotate_structure(cutoff=0.05)
        logging.info("Finished computing structural properties.")

    # Generate figures if specified
    if args["visualize"]:
        logging.info("Generating figures...")
        annotator.visualize()
        logging.info("Finished generating figures.")

    # Generate logos if specified
    if args["logo"] is not None:
        logging.info("Generating logos...")
        annotator.create_logos()
        logging.info("Finished generating logos.\n")

    # Write output files
    logging.info("Writing output files...")
    annotator.write_files()
    logging.info("Finished writing output files")

    # Calculate and log the total running time
    end_time = time.perf_counter()
    total_time = str(timedelta(seconds=end_time - start_time))
    logging.info(f"\nCLIPPER 2.0 pipeline completed in {total_time}.\n")


if __name__ == "__main__":
    main()
