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
    logging.info("Starting annotation")
    print("Annotating peptides...")
    if not args["singlecpu"]:
        annotator.threaded_annotate()
    else:
        annotator.annotate()
    logging.info("Finished annotation")
    print("Finished annotating peptides.")

    # Perform proteoform check
    annotator.proteoform_check()
    logging.info("Finished proteoform check")

    # Perform exopeptidase check if specified
    if not args["noexo"]:
        logging.info("Starting exopeptidase check")
        print("Checking exopeptidase activity...")
        annotator.exopeptidase()
        logging.info("Finished exopeptidase check")
        print("Finished exopeptidase activity check.")

    # Perform general statistics annotation if condition file is not specified
    if args["conditionfile"] is None:
        logging.info("Starting general statistics annotation")
        print("Performing general statistics annotation...")
        annotator.general_conditions()
        logging.info("Finished general statistics annotation")
        print("Finished general statistics annotation")

    # Perform condition-specific statistics annotation if condition file is specified
    else:
        # Read the condition file
        annotator.read_condition_file()
        logging.info("Parsed condition file")

        # Perform general statistics annotation
        logging.info("Starting general statistics annotation")
        print("Performing general statistics annotation...")
        annotator.general_conditions()
        logging.info("Finished general statistics annotation")
        print("Finished general statistics annotation")

        # Perform pairwise or all-vs-all statistical testing if specified
        if args["stat"]:
            print("Performing statistical testing...")
            annotator.condition_statistics(pairwise=args["stat_pairwise"])
            logging.info("Finished statistical testing")
            print("Finished statistical testing")

        # Perform fold distribution check if specified
        if args["significance"]:
            print("Checking fold distribution...")
            annotator.percentile_fold(0.05)
            logging.info("Finished fold distribution check")
            print("Finished fold distribution check")

    # Generate figures if specified
    if args["visualize"]:
        print("Generating figures...")
        annotator.visualize()
        logging.info("Finished generating figures")
        print("Finished generating figures")

    # Generate logos if specified
    if args["logo"] is not None:
        print("Generating logos...")
        annotator.create_logos()
        logging.info("Finished generating logos")
        print("Finished generating logos")

    # Write output files
    print("Writing output files...")
    annotator.write_files()
    logging.info("Finished writing output files")

    # Calculate and log the total running time
    end_time = time.perf_counter()
    total_time = str(timedelta(seconds=end_time - start_time))
    logging.info(f"Pipeline completed in {total_time}")
    print(f"\nCLIPPER pipeline completed in {total_time}")


if __name__ == "__main__":
    main()
