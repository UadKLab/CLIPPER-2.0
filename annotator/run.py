import time
import logging
from datetime import timedelta

from annotator import Annotator
from annutils import initialize


def main(args=None):
    """Main function"""

    start = time.perf_counter()

    # run with arguments without command line. Must pass a dict with arguments
    if args is None:
        args = initialize()

    ann = Annotator(args)

    ann.prepare()

    print("Annotating...")
    if args["singlecpu"] is False:
        ann.threaded_annotate()
    else:
        ann.annotate()
    logging.info("Finished annotation")
    print("Finished annotation")

    ann.proteoform_check()
    logging.info("Finished proteoform check")

    if args["noexo"] is False:
        logging.info("Starting exopeptidase check")
        print("Starting exopeptidase check...")
        ann.exopeptidase()
        logging.info("Finished exopeptidase check")
        print("Finished exopeptidase activity check.")

    if args["conditionfile"] is None:
        print("Starting general stats annotation...")
        #ann.general_conditions()
        logging.info("Finished general stats")
        print("Finished general stats")
    else:
        ann.read_condition_file()
        logging.info("Parsed condition file")
        print("Starting general stats annotation...")
        ann.general_conditions()
        logging.info("Finished general stats")
        print("Finished general stats")

        if args["stat"]:
            print("Starting statistical testing...")
            ann.condition_statistics(pairwise=args["stat_pairwise"])
            logging.info("Finished statistical testing")
            print("Finished statistical testing")
        
        if args["significance"]:
            print("Starting fold distribution check...")
            ann.percentile_fold(0.05)
            logging.info("Finished fold distribution check")
            print("Finished fold distribution check")

    if args["visualize"]:
        print("Generating figures...")
        ann.visualize()
        logging.info("Finished figures")
        print("Finished figures")

    if args["logo"] is not None:
        print("Starting logo generation...")
        ann.create_logos()
        logging.info("Finished logos")
        print("Finished logos")

    print("Writing files...")
    ann.write_files()
    logging.info("Wrote files")

    end = time.perf_counter()
    logging.info(f"Reported success in {str(timedelta(seconds=end-start))}")
    print(f"\nReported success in {str(timedelta(seconds=end-start))}")


if __name__ == "__main__":
    main()
